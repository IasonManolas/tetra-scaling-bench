#!/usr/bin/env python3
"""Name the functions that run while the other threads have nothing to do.

The stage timers say which stage of the remesher does not speed up; they
cannot say which function inside it. This reads a `perf record` of one
multi-threaded run and answers that directly.

Every sample carries the thread it came from. The run is cut into short time
slices, and in each slice a thread counts as WORKING if it has samples in code
other than the TBB library and the kernel -- idle TBB workers spin and yield
inside exactly those two, so they must not count. A slice in which at most one
thread is working is a serial slice: whatever that one thread was running
there is serial code, and its samples, converted to milliseconds, are what that
code costs in wall time.

Caveat, stated in the output too: a thread that is retrying a lock yields in
the kernel, so a slice of heavy lock contention could read as fewer threads
working than there really were. A retrying thread also runs remeshing code
between yields, so within a 10 ms slice it normally still counts.

Usage:
    serial_profile.py perf.data [--freq 999] [--slice-ms 10] [--top 30]
                      [--out summary.txt] [--json summary.json]
"""

import argparse
import collections
import json
import re
import subprocess
import sys

# "1234/1240 5678.901234:  7f12ab symbol+0x1f (/path/to/dso)" -- perf script's
# layout for -F pid,tid,time,ip,sym,dso (without ip it prints no symbol at
# all). The symbol can contain spaces and parentheses (C++), so the dso is
# taken as the LAST parenthesised group.
LINE = re.compile(r"^\s*(\d+)/(\d+)\s+([\d.]+):\s+[0-9a-f]+\s+(.*)\s+\(([^()]*)\)\s*$")

# Where an idle TBB worker spends its time. The kernel shows as [unknown] when
# kernel addresses are hidden (kptr_restrict), as [kernel.kallsyms] otherwise.
# libc is NOT here: malloc and memcpy are real work.
IDLE_DSO = ("libtbb", "[kernel", "kallsyms", "[vdso]", "[unknown]")


def is_idle(dso, sym):
    d = dso.lower()
    return any(k in d for k in IDLE_DSO)


_DEMANGLED = {}


def demangle(sym):
    """perf leaves some very long C++ names mangled; c++filt does not."""
    if not sym.startswith("_Z"):
        return sym
    if sym not in _DEMANGLED:
        try:
            _DEMANGLED[sym] = subprocess.run(
                ["c++filt", "--no-recurse-limit", sym], capture_output=True, text=True,
                timeout=10).stdout.strip() or sym
        except Exception:
            _DEMANGLED[sym] = sym
    return _DEMANGLED[sym]


_TEMPLATE_ARGS = re.compile(r"<[^<>]*>")


def short(sym, width=110):
    """The function's own name, without the template arguments.

    A CGAL name is mostly template arguments -- the start of it names the
    triangulation type, not the function -- so cutting the full name at a
    fixed width would cut off exactly the part that identifies it.
    """
    sym = demangle(re.sub(r"\+0x[0-9a-f]+$", "", sym.strip()))
    # Innermost first, with a placeholder that has no angle brackets, so each
    # pass exposes the next level out.
    prev = None
    while prev != sym:
        prev, sym = sym, _TEMPLATE_ARGS.sub("\x00", sym)
    sym = sym.replace("\x00", "<>")
    return sym if len(sym) <= width else sym[:width - 3] + "..."


def read_samples(perf_data):
    cmd = ["perf", "script", "-i", perf_data, "-F", "pid,tid,time,ip,sym,dso"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                         text=True, errors="replace")
    for line in p.stdout:
        m = LINE.match(line)
        if m:
            pid, tid, t, sym, dso = m.groups()
            yield int(pid), int(tid), float(t), sym, dso
    p.wait()


def analyse(samples, freq, slice_s):
    slices = collections.defaultdict(lambda: collections.defaultdict(list))
    main_by_sym = collections.Counter()
    workers_by_sym = collections.Counter()
    pid0 = None
    t0 = t1 = None
    tids = set()
    names = {}
    for pid, tid, t, sym, dso in samples:
        if pid0 is None:
            pid0 = pid
        if pid != pid0:
            continue
        t0 = t if t0 is None else min(t0, t)
        t1 = t if t1 is None else max(t1, t)
        tids.add(tid)
        idle = is_idle(dso, sym)
        if idle:
            name = "(idle: %s)" % dso.split("/")[-1]
        else:
            name = names.get(sym)
            if name is None:
                name = names[sym] = short(sym)
        slices[int(t / slice_s)][tid].append((name, idle))
        if tid == pid0:
            main_by_sym[name] += 1
        else:
            workers_by_sym[name] += 1
    if t0 is None:
        sys.exit("No samples read. Was perf record run on this binary, and is "
                 "`perf script` working?")

    wall = t1 - t0
    per_slice = freq * slice_s
    # A thread is working in a slice if at least a fifth of the samples it
    # would have produced running flat out are in real code.
    need = max(1, int(round(per_slice / 5)))
    hist = collections.Counter()
    serial_by_sym = collections.Counter()
    serial_slices = 0
    for k, by_tid in slices.items():
        working = [tid for tid, xs in by_tid.items()
                   if sum(1 for _, idle in xs if not idle) >= need]
        hist[len(working)] += 1
        if len(working) <= 1:
            serial_slices += 1
            for tid in working:
                for name, idle in by_tid[tid]:
                    if not idle:
                        serial_by_sym[name] += 1
    n_slices = sum(hist.values())
    return dict(pid=pid0, threads_seen=len(tids), wall_s=wall, freq=freq,
                slice_ms=slice_s * 1e3, n_slices=n_slices,
                serial_slices=serial_slices, hist=dict(sorted(hist.items())),
                serial_by_sym=serial_by_sym, main_by_sym=main_by_sym,
                workers_by_sym=workers_by_sym)


def report(r, top):
    ms = 1e3 / r["freq"]
    wall_ms = r["wall_s"] * 1e3
    out = []
    w = out.append
    w("Per-thread profile: %d threads seen, %.2f s of samples at %d Hz, "
      "%.0f ms slices" % (r["threads_seen"], r["wall_s"], r["freq"], r["slice_ms"]))
    w("The whole process is covered, so reading the input and writing the "
      "results show up too, as stream and number-parsing functions; they are "
      "not remeshing.")
    w("")
    w("Threads working per slice (a thread 'works' if it runs code outside "
      "the TBB library and the kernel):")
    for k, n in r["hist"].items():
        w("  %3d thread(s) working: %5.1f%% of the run" % (k, 100.0 * n / r["n_slices"]))
    w("")
    w("Serial slices (at most one thread working): %.1f%% of the run, %.2f s"
      % (100.0 * r["serial_slices"] / r["n_slices"],
         r["serial_slices"] * r["slice_ms"] / 1e3))
    w("")
    w("Functions running during serial slices, by wall time they cost:")
    w("  %9s  %6s  %s" % ("ms", "% run", "function"))
    for name, n in r["serial_by_sym"].most_common(top):
        w("  %9.0f  %5.1f%%  %s" % (n * ms, 100.0 * n * ms / wall_ms, name))
    w("")
    w("All samples, main thread against the workers (a function the workers "
      "never run is serial wherever it appears):")
    w("  %9s  %9s  %8s  %s" % ("main ms", "workers ms", "workers", "function"))
    for name, n in r["main_by_sym"].most_common(top):
        if name.startswith("(idle"):
            continue
        wk = r["workers_by_sym"].get(name, 0)
        w("  %9.0f  %9.0f  %7.0f%%  %s"
          % (n * ms, wk * ms, 100.0 * wk / (n + wk), name))
    w("")
    w("Caveat: a thread retrying a lock yields in the kernel. Heavy contention "
      "can make a slice read as fewer threads working than there were.")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("perf_data")
    ap.add_argument("--freq", type=int, default=999,
                    help="the -F given to perf record")
    ap.add_argument("--slice-ms", type=float, default=10.0)
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--out")
    ap.add_argument("--json")
    a = ap.parse_args()
    r = analyse(read_samples(a.perf_data), a.freq, a.slice_ms / 1e3)
    text = report(r, a.top)
    print(text)
    if a.out:
        open(a.out, "w").write(text + "\n")
    if a.json:
        j = dict((k, v) for k, v in r.items()
                 if k not in ("serial_by_sym", "main_by_sym", "workers_by_sym"))
        j["serial_by_sym_ms"] = dict((k, v * 1e3 / a.freq) for k, v in
                                     r["serial_by_sym"].most_common(200))
        j["main_ms"] = dict((k, v * 1e3 / a.freq) for k, v in
                            r["main_by_sym"].most_common(200))
        j["workers_ms"] = dict((k, v * 1e3 / a.freq) for k, v in
                               r["workers_by_sym"].most_common(200))
        open(a.json, "w").write(json.dumps(j, indent=1))


if __name__ == "__main__":
    main()
