#!/usr/bin/env python3
"""Stage 2: the timed sweep.

Three profiles:

  full        The whole MEASUREMENT, unattended, in one command: calibrate,
              size the grid from what calibration measured, then sweep it. Use
              this when you get a single run at the machine and nobody is
              available to make decisions partway through. Quality metrics and
              the report are separate steps afterwards -- they read the saved
              meshes and do not need the machine idle.

  calibrate   Under a hard wall-clock budget (30 min by default). Finds which
              (mesh, edge factor) pairs make the 24-thread run last at least 30
              seconds, measures this machine's run-to-run noise, and sketches
              the scaling curve. Its output is what sizes the real run -- we
              cannot pick those configs in advance without knowing the machine.

  overnight   The real grid, driven by a configs file that calibrate produced
              (or a hand-written one). Resumable: a parseable result JSON means
              that cell is done and is skipped, so an interrupted run resumes
              where it stopped.

What is measured, and what is deliberately not: the timed process records only
things that cannot be recomputed later -- wall time, peak RSS, exit code, the
input's average edge length -- and writes the remeshed mesh to disk. Every
quality metric is derived afterwards from that mesh by --quality-pass. See the
README.

    python3 scripts/run_bench.py --root work --profile full --budget 43200
    python3 scripts/run_bench.py --root work --profile calibrate
    python3 scripts/run_bench.py --root work --profile overnight
    python3 scripts/run_bench.py --root work --quality-pass
"""
import argparse
import csv
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

ITERS = 10                 # remeshing iterations, as in the existing campaign
SMOOTH_CONSTRAINED = 1
TARGET_SECONDS = 30.0      # the 24-thread run must last at least this long
LADDER = [0.5, 0.4, 0.3, 0.25]   # cells grow roughly as factor^-3
# Below this, a run is dominated by process start-up and mesh I/O rather than
# remeshing, and the seq-vs-parallel ratio taken from it is meaningless. Seen
# in testing: a 0.5s config gave a 1.75x "serial ratio" on 4 cores.
ANCHOR_MIN_SECONDS = 4.0
# Floor on the time the full profile holds back for the quality pass. A sweep
# whose meshes were never measured answers half the question.
# Below this there is no point starting: wave 1 alone needs room for four arms
# on at least a couple of configs.
MIN_SWEEP_SECONDS = 600.0
THREAD_LIST = [1, 2, 4, 8, 12, 16, 24]

# --- the metrics profile ---------------------------------------------------
# Fixed configurations, named in docs/METRICS_REQUEST.md. They are hard-coded
# rather than calibrated because the point of this profile is to repeat exactly
# the configurations the earlier ladder used, with cycles and instructions added,
# so the two datasets are comparable. Everything else the smoke runs measured
# answers a different question and is deliberately left out.
METRICS_LADDER_CONFIGS = [("94665_cdt", 0.3), ("67856_cdt", 0.25),
                          ("94665_mesh3", 0.25)]
METRICS_REPS = 3
METRICS_PIN_CONFIG = ("94665_cdt", 0.3)
METRICS_PIN_THREADS = 8
METRICS_LOCK_GRIDS = [16, 24, 32, 48, 64]
METRICS_LOCK_CONFIGS = [("94665_cdt", 0.3), ("94665_mesh3", 0.25)]
METRICS_DIAG_CONFIG = ("94665_cdt", 0.3)
METRICS_DIAG_THREADS = [1, 24]
LOCK_GRID_ENV = "CGAL_TETRAHEDRAL_REMESHING_LOCK_GRID"
PERF_EVENTS = "instructions,cycles,task-clock"
# NOT a comma, even though the request writes `-x,`. perf formats its numbers in
# the machine's locale before joining them with the separator, so on a machine
# whose decimal mark is a comma the task-clock row comes out as `183,46,msec,...`
# and every field after the first shifts by one. That is silent: instructions and
# cycles still parse, task-clock reads as the number 183 followed by the event
# name "msec", and the column ends up empty for no visible reason. A semicolon
# cannot be a decimal mark in any locale.
PERF_SEP = ";"

# started_at is not decoration. Cells are measured in different epochs -- a
# calibration cell reused by the sweep can be an hour older than its own
# repeat -- and machines drift (25% between waves on the test laptop). Without
# a timestamp there is no way to tell a real effect from thermal drift after
# the fact, and no way to re-measure only the suspect rows.
CSV_HEADER = ["mesh", "factor", "arm", "threads", "rep", "wall_s", "remesh_s",
              "peak_rss_kb", "rc", "timed_out", "started_at",
              # How the machine behaved while this ran. avg_parallelism is CPU
              # time over wall time -- the cores actually working -- and the
              # frequency and throttle columns say whether the clocks held up.
              # A low speedup with avg_parallelism near the thread count is the
              # algorithm; near 1 it never went parallel; falling clocks or
              # rising throttle_events mean the machine, not the code.
              "user_s", "sys_s", "cpu_pct", "avg_parallelism",
              "freq_mean_mhz", "freq_min_mhz", "freq_max_mhz", "freq_samples",
              # The frequency of the cores that were actually WORKING. The
              # freq_* columns above average every core including idle ones, so
              # they rise with the thread count and read backwards; they are kept
              # only so old and new results sets have the same shape. Quote the
              # busy_freq_* columns instead.
              "busy_freq_mean_mhz", "busy_freq_min_mhz", "busy_freq_max_mhz",
              "busy_freq_samples",
              # Thermal-throttle entries during this cell. It was dropped from
              # this list when the busy_freq_* block was added, while run_cell
              # went on recording it -- DictWriter ignores a field the header
              # does not name, so the column vanished silently and report.py's
              # throttling check read None for every row.
              "throttle_events",
              # Work, as opposed to CPU seconds. user_s moves with the clock --
              # the same binary on the same input measured +44.6% CPU seconds
              # from 1 to 4 threads but only +17.1% instructions -- so it
              # overstates the waste roughly threefold and cannot be corrected
              # afterwards. cycles over task_clock_ms is the clock the run held.
              "instructions", "cycles", "task_clock_ms",
              # What made this cell different from the plain sweep: an empty
              # variant is an ordinary timed run. `pin` is the taskset CPU list,
              # `lock_grid` the value of CGAL_TETRAHEDRAL_REMESHING_LOCK_GRID.
              "variant", "pin", "lock_grid",
              "out_mesh", "json"]


# --------------------------------------------------------------------------
# environment capture
# --------------------------------------------------------------------------

def capture(cmd):
    try:
        return subprocess.run(cmd, text=True, capture_output=True,
                              timeout=30).stdout.strip()
    except Exception:
        return ""


def _read_first(rel):
    """One small sysfs value under /sys/devices/system/cpu, or ''."""
    try:
        return (CPU_ROOT / rel).read_text().strip()
    except Exception:
        return ""


def check_toolchain_lock(root, results_dir, allow_change):
    """Refuse to extend a results set that was measured with different binaries.

    setup.py fetches the branch tip every time it runs, which is what you want
    when starting, and exactly what you do not want on a resume: an interrupted
    12-hour run restarted after a push would rebuild, and the remaining cells
    would be measured against different code and appended to the same CSV. The
    SHAs live in env.json, which is rewritten each invocation, so afterwards
    there would be nothing to reveal the mixture.

    The first sweep into a results directory records what it used; later ones
    must match.
    """
    tc = json.loads((root / "toolchain.json").read_text())
    # The ref names are recorded as well as the commits. A SHA alone cannot say
    # WHICH branch was measured, and the branch under test has moved once
    # already (to scale24, 2026-09-22); a results set that does not name it
    # leaves anyone holding the tarball guessing.
    now = {"ours_sha": tc.get("ours_sha"), "main_sha": tc.get("main_sha"),
           "ours_ref": tc.get("ours_ref"), "main_ref": tc.get("main_ref"),
           "binaries": {n: v.get("md5") for n, v in tc.get("binaries", {}).items()}}

    lock_path = results_dir / "toolchain_lock.json"
    if not lock_path.exists():
        lock_path.write_text(json.dumps(now, indent=2))
        return

    was = json.loads(lock_path.read_text())
    # Compare only what the stored lock actually recorded. A lock written before
    # the ref names were added must not read as "the branch changed" merely
    # because it does not mention one -- that would refuse a resume which is in
    # fact perfectly consistent. Anything it did record is compared strictly.
    if all(was.get(k) == now.get(k) for k in was):
        return

    diffs = []
    for k in ("ours_ref", "ours_sha", "main_ref", "main_sha"):
        if k in was and was.get(k) != now.get(k):
            diffs.append("  %-9s %s -> %s" % (k, str(was.get(k))[:40],
                                              str(now.get(k))[:40]))
    for n, md5 in sorted(now["binaries"].items()):
        if n in was.get("binaries", {}) and was["binaries"].get(n) != md5:
            diffs.append("  %-9s %s rebuilt" % ("binary", n))

    msg = ("The code changed since this results directory was started:\n%s\n"
           % "\n".join(diffs))
    if not allow_change:
        sys.exit(msg +
                 "\nContinuing would append cells measured with different binaries to\n"
                 "the same results.csv, and nothing in the output would show it.\n"
                 "Either:\n"
                 "  - start a fresh results directory (--results-dir), or\n"
                 "  - rebuild the original commit and re-run, or\n"
                 "  - pass --allow-toolchain-change if you really mean to mix them.")
    print(msg + "Continuing because --allow-toolchain-change was given; the\n"
          "results set now mixes binaries.", file=sys.stderr)
    lock_path.write_text(json.dumps(now, indent=2))


def write_env(root, results_dir):
    tc = json.loads((root / "toolchain.json").read_text())
    governor = ""
    gp = Path("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor")
    if gp.exists():
        try:
            governor = gp.read_text().strip()
        except Exception:
            pass

    env = {
        "captured_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "nproc": os.cpu_count(),
        "cpu_governor": governor,
        # The driver matters more than the governor name. Under intel_pstate,
        # 'powersave' is the normal default and still boosts to full turbo -- it
        # is not the old ondemand-style throttling, and warning about it misleads.
        "cpufreq_driver": _read_first("cpu0/cpufreq/scaling_driver"),
        "cpufreq_max_khz": _read_first("cpu0/cpufreq/cpuinfo_max_freq"),
        "cpufreq_min_khz": _read_first("cpu0/cpufreq/cpuinfo_min_freq"),
        # 1 here means turbo is disabled outright, which caps every thread count
        # at base clock and would look exactly like poor scaling.
        "no_turbo": _read_first("intel_pstate/no_turbo"),
        "throttle_events_at_start": _read_throttle_counts(),
        "lscpu": capture(["lscpu"]),
        # A hybrid CPU's cores are not interchangeable, and which CPU number is
        # which kind is not fixed across machines. These three are what let the
        # performance and efficiency sets be identified after the fact, and are
        # what the core-pinning comparison pinned to.
        "lscpu_extended": capture(["lscpu", "-e"]),
        "cpuinfo_max_freq_khz": cpu_max_freqs_khz(),
        "core_sets": dict(zip(("performance", "efficiency"), core_sets())),
        "numactl": capture(["numactl", "--hardware"]),
        "meminfo": capture(["head", "-3", "/proc/meminfo"]),
        "compiler": capture(["c++", "--version"]).splitlines()[:1],
        "toolchain": tc,
        "jemalloc_preloaded": bool(os.environ.get("LD_PRELOAD", "")),
        "env_LD_PRELOAD": os.environ.get("LD_PRELOAD", ""),
        "env_LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
    }
    (results_dir / "env.json").write_text(json.dumps(env, indent=2))
    return env


# --------------------------------------------------------------------------
# perf, and the core sets
# --------------------------------------------------------------------------

_PERF = {}


def perf_usable():
    """Whether `perf stat` can actually count on this machine.

    Probed once, by counting a trivial command, because the two ways it fails
    are different: perf may not be installed at all, or it may be installed and
    refused by kernel.perf_event_paranoid. Either way the run must still happen
    -- a missing counter is a missing column, not a reason to lose 45 minutes of
    machine time -- so this returns a verdict rather than exiting.
    """
    if "ok" in _PERF:
        return _PERF["ok"]
    _PERF["ok"] = False
    if shutil.which("perf") is None:
        print("[perf] not on PATH; instructions/cycles/task_clock_ms will be "
              "empty. Install linux-tools to get them.", file=sys.stderr)
        return False
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".perf") as tf:
        try:
            subprocess.run(["perf", "stat", "-e", PERF_EVENTS, "-x", PERF_SEP,
                            "-o", tf.name, "true"],
                           capture_output=True, timeout=60)
            counted = read_perf(Path(tf.name))
        except Exception:
            counted = {}
    if counted.get("instructions") is None:
        para = ""
        try:
            para = Path("/proc/sys/kernel/perf_event_paranoid").read_text().strip()
        except Exception:
            pass
        print("[perf] present but cannot count (perf_event_paranoid=%s). The "
              "instructions/cycles/task_clock_ms columns will be empty.\n"
              "       To enable it:  sudo sysctl -w kernel.perf_event_paranoid=1"
              % (para or "unknown"), file=sys.stderr)
        return False
    _PERF["ok"] = True
    print("[perf] counting %s" % PERF_EVENTS)
    return True


def read_perf(path):
    """Parse one `perf stat -x,` file into the three columns we record.

    A counter that could not be read prints `<not counted>` or `<not supported>`
    where the number goes, so a value that does not parse is simply absent data.
    """
    out = {}
    try:
        text = Path(path).read_text()
    except Exception:
        return out
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(PERF_SEP)
        if len(parts) < 3:
            continue
        try:
            # A comma here is a decimal mark, not a separator: see PERF_SEP.
            value = float(parts[0].replace(",", "."))
        except ValueError:
            continue        # <not counted> / <not supported>
        event = parts[2].strip()
        if event == "instructions":
            out["instructions"] = int(value)
        elif event == "cycles":
            out["cycles"] = int(value)
        elif event == "task-clock":
            # perf reports task-clock in milliseconds.
            out["task_clock_ms"] = round(value, 3)
    return out


def cpu_max_freqs_khz():
    """Every logical CPU's maximum clock, keyed by CPU number.

    This is what separates performance cores from efficiency cores on a hybrid
    part, and it is read at run time rather than hard-coded because the CPU
    numbering of the two kinds is not fixed across machines or kernels.
    """
    out = {}
    for p in CPU_ROOT.glob("cpu[0-9]*/cpufreq/cpuinfo_max_freq"):
        try:
            out[int(p.parent.parent.name[3:])] = int(p.read_text().strip())
        except Exception:
            pass
    return out


def core_sets(n=None):
    """The fastest and the slowest CPUs, as two lists of CPU numbers.

    CPUs are grouped by their maximum clock; the highest group is the
    performance cores and the lowest the efficiency cores. On a machine whose
    cores are all the same the two groups coincide, and the caller is told so by
    getting the same list twice.
    """
    freqs = cpu_max_freqs_khz()
    if not freqs:
        return [], []
    groups = {}
    for cpu, f in freqs.items():
        groups.setdefault(f, []).append(cpu)
    fast = sorted(groups[max(groups)])
    slow = sorted(groups[min(groups)])
    if n:
        fast, slow = fast[:n], slow[:n]
    return fast, slow


# --------------------------------------------------------------------------
# one run
# --------------------------------------------------------------------------

def _kill_group(proc):
    """SIGTERM then SIGKILL the process group, and wait for it to actually die.

    Waiting matters as much as killing: returning while the group is still
    winding down hands the next measurement a machine that is not idle.
    """
    try:
        pgid = os.getpgid(proc.pid)
    except Exception:
        pgid = None

    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            if pgid is not None:
                os.killpg(pgid, sig)
            else:
                proc.send_signal(sig)
        except ProcessLookupError:
            break
        except Exception:
            break
        try:
            proc.wait(timeout=15)
            break
        except subprocess.TimeoutExpired:
            continue

    # The group's other members are not proc, so proc.wait() says nothing about
    # them. Poll until none are left, rather than trusting the signal landed.
    if pgid is not None:
        deadline = time.time() + 30
        while time.time() < deadline:
            try:
                os.killpg(pgid, 0)
            except ProcessLookupError:
                return
            except Exception:
                return
            time.sleep(0.25)
        print("Warning: process group %d did not die; later timings may be\n"
              "         contaminated. Check for stray bench_remesh processes."
              % pgid, file=sys.stderr)


CPU_ROOT = Path("/sys/devices/system/cpu")


def _read_throttle_counts():
    """Total thermal-throttle events across cores, or None if unavailable.

    A rising count during a run is the machine telling us it cut clocks to stay
    within thermal or power limits. That is the difference between 'this code
    does not scale' and 'this machine could not sustain the clocks', and the two
    have completely different answers.
    """
    # core_throttle_count only. Adding package_throttle_count to it sums two
    # different things and inflates the magnitude into something that looks
    # alarming but means nothing. Even here the magnitude is only a count of
    # throttle ENTRIES, not of severity -- what matters is whether the delta
    # across a run is zero or not.
    total = 0
    found = False
    for p in CPU_ROOT.glob("cpu*/thermal_throttle/core_throttle_count"):
        try:
            total += int(p.read_text().strip())
            found = True
        except Exception:
            pass
    return total if found else None


def _read_freqs_khz_by_cpu():
    """Current clock of every logical CPU, in kHz, keyed by CPU number."""
    out = {}
    for p in CPU_ROOT.glob("cpu[0-9]*/cpufreq/scaling_cur_freq"):
        try:
            out[int(p.parent.parent.name[3:])] = int(p.read_text().strip())
        except Exception:
            pass
    return out


def _read_busy_jiffies():
    """Per-CPU busy time from /proc/stat, keyed by CPU number.

    Busy is everything except idle and iowait. Only the DIFFERENCE between two
    reads means anything, which is why the sampler keeps the previous one.
    """
    out = {}
    try:
        for line in Path("/proc/stat").read_text().splitlines():
            if not line.startswith("cpu") or line.startswith("cpu "):
                continue
            parts = line.split()
            try:
                cpu = int(parts[0][3:])
                v = [int(x) for x in parts[1:]]
            except ValueError:
                continue
            idle = v[3] + (v[4] if len(v) > 4 else 0)
            out[cpu] = sum(v) - idle
    except Exception:
        pass
    return out


class _FreqSampler(threading.Thread):
    """Samples CPU clock while a measurement runs.

    Two numbers come out of it, and only the second is readable.

    `freq_*` is the mean over every logical CPU, idle ones included. On a 24-core
    machine running four threads, twenty idle cores sit at their minimum clock
    and drag that mean down, so it RISES with the thread count -- the opposite of
    what actually happens to the cores doing the work, and the reason the
    returned smoke data could not be used to price the clock drop.

    `busy_freq_*` is the mean over the `n_busy` CPUs that burned the most time
    between this sample and the last one -- the cores the run was actually on.
    Picking them by busy time rather than by clock matters: sorting the cores by
    their own frequency and taking the top N would select for high clocks and
    report a number that is high by construction.

    Deliberately light: a handful of small sysfs reads on a 2s period, on one
    thread, against a job that is saturating cores for tens of seconds.
    """

    def __init__(self, n_busy=None, period=2.0):
        super().__init__(daemon=True)
        self.period = period
        self.n_busy = n_busy
        # NOT self._stop: Thread._stop() is an internal method that join()
        # calls, and shadowing it with an Event makes join() raise
        # "'Event' object is not callable".
        self._stop_evt = threading.Event()
        self.samples = []       # mean MHz across every CPU, per sample
        self.busy_samples = []  # mean MHz across the busiest n_busy CPUs
        self._prev_busy = _read_busy_jiffies()

    def run(self):
        while not self._stop_evt.is_set():
            freqs = _read_freqs_khz_by_cpu()
            busy = _read_busy_jiffies()
            if freqs:
                self.samples.append(sum(freqs.values()) / len(freqs) / 1000.0)
                picked = self._busiest(busy, freqs)
                if picked:
                    self.busy_samples.append(
                        sum(freqs[c] for c in picked) / len(picked) / 1000.0)
            self._prev_busy = busy or self._prev_busy
            self._stop_evt.wait(self.period)

    def _busiest(self, busy, freqs):
        n = self.n_busy or len(freqs)
        n = max(1, min(n, len(freqs)))
        deltas = [(busy.get(c, 0) - self._prev_busy.get(c, 0), c)
                  for c in freqs if c in busy and c in self._prev_busy]
        if not deltas:
            return []
        deltas.sort(reverse=True)
        return [c for _, c in deltas[:n]]

    @staticmethod
    def _stats(samples, prefix):
        if not samples:
            return {}
        s = sorted(samples)
        return {prefix + "mean_mhz": round(sum(s) / len(s), 1),
                prefix + "min_mhz": round(s[0], 1),
                prefix + "max_mhz": round(s[-1], 1),
                prefix + "samples": len(s)}

    def stop(self):
        self._stop_evt.set()
        self.join(timeout=5)
        out = self._stats(self.samples, "freq_")
        out.update(self._stats(self.busy_samples, "busy_freq_"))
        return out


# The status file is written by a tiny shell wrapper around the binary, and is
# the only honest record of how the binary exited. Both wrappers above it lie:
# `perf stat` exits 0 for a command killed by a signal, and /usr/bin/time does
# not re-raise it either, so without this a SIGSEGV reads as a clean pass.
# $0 is the path to write, $@ is the command.
_STATUS_SH = '"$@"; s=$?; printf %s "$s" > "$0"; exit $s'


def cell_stem(mesh_id, factor, arm, threads, rep, variant=""):
    """The name every artefact of one cell shares, and its resume key.

    A variant is what made this cell different from a plain timed run -- a
    taskset CPU set, a lock-grid value, a diagnostic binary. An empty variant
    reproduces the original spelling exactly, so results directories written
    before this existed still resume.
    """
    stem = "%s_f%s_%s_t%d_r%d" % (mesh_id, factor, arm, threads, rep)
    return stem + ("_" + variant if variant else "")


def run_cell(bins, mesh_path, factor, arm, threads, rep, results_dir,
             mesh_out_dir, timeout, write_mesh, settle, done=None,
             variant="", pin_cpus=None, extra_env=None, exe=None):
    """Execute one (mesh, factor, arm, threads, rep) cell.

    Always returns a dict row. A cell that was already completed by an earlier
    invocation is returned from `done` with cached=1 and is NOT re-run and NOT
    re-appended to the CSV. Returning the cached row rather than None matters:
    calibrate reads wall times back to choose its pilot config, so a resumed
    run must still reach the same conclusions without redoing the work.
    """
    mesh_id = Path(mesh_path).stem
    stem = cell_stem(mesh_id, factor, arm, threads, rep, variant)
    json_path = results_dir / "json" / (stem + ".json")
    log_path = results_dir / "logs" / (stem + ".log")

    # Resumability, the pattern this project's thingi_sweep.sh uses: a JSON that
    # parses means the cell completed. A crashed run leaves no JSON or a broken
    # one, and gets redone.
    if json_path.exists():
        try:
            json.loads(json_path.read_text())
            if done is not None and stem in done:
                row = dict(done[stem])
                row["cached"] = 1
                return row
            # The JSON is there but the CSV row is not (CSV deleted, or the
            # process died between writing the JSON and flushing the row).
            # Fall back to the remesher's own timer, which is inside the JSON:
            # a little smaller than wall clock since it excludes I/O, but a far
            # better answer than 0.0, which would make calibrate mistake this
            # for the fastest config there is.
            d = json.loads(json_path.read_text())
            t = d.get("metrics", {}).get("Performance", {}).get("Total_Time", {}).get("Value")
            return {"mesh": mesh_id, "factor": factor, "arm": arm,
                    "threads": threads, "rep": rep,
                    "wall_s": float(t) if t is not None else 0.0,
                    "remesh_s": t, "peak_rss_kb": "", "rc": 0,
                    "timed_out": 0, "out_mesh": "", "json": str(json_path),
                    "variant": variant, "pin": ",".join(str(c) for c in pin_cpus)
                    if pin_cpus else "",
                    "lock_grid": (extra_env or {}).get(LOCK_GRID_ENV, ""),
                    "cached": 1, "recovered": 1}
        except Exception:
            json_path.unlink(missing_ok=True)

    if exe is None:
        exe = bins["bench_remesh_main"] if arm == "main" else bins["bench_remesh"]
    tag = "seq" if arm in ("seq", "main") else "par"

    cmd = [str(exe), str(mesh_path), str(ITERS), str(factor),
           str(SMOOTH_CONSTRAINED), str(threads), str(json_path), "--tag", tag]

    out_mesh = ""
    if write_mesh:
        # One output mesh per (mesh, factor, arm, threads) -- not per rep. Reps
        # exist to average timing noise; the deterministic arms produce the same
        # mesh every time, and for the parallel arm one sample per thread count
        # is what the determinism comparison needs.
        out_mesh = str(mesh_out_dir / ("%s_f%s_%s_t%d.mesh" % (mesh_id, factor, arm, threads)))
        cmd += ["--out-mesh", out_mesh]

    if settle:
        time.sleep(settle)

    rss_file = results_dir / "logs" / (stem + ".time")
    perf_file = results_dir / "perf" / (stem + ".perf")
    status_file = results_dir / "logs" / (stem + ".status")
    for f in (rss_file, perf_file, status_file):
        f.unlink(missing_ok=True)

    if pin_cpus:
        # taskset goes INSIDE the counters, so perf counts the remesher and not
        # itself, and the pinning applies to the measured process only.
        cmd = ["taskset", "-c", ",".join(str(c) for c in pin_cpus)] + cmd

    # Innermost: record the remesher's own exit status where nothing above can
    # hide it. See _STATUS_SH.
    wrapped = ["sh", "-c", _STATUS_SH, str(status_file)] + cmd

    if shutil.which("/usr/bin/time"):
        # %U and %S give CPU time, and CPU time over wall time is the average
        # number of cores actually working. Without it, a low speedup cannot be
        # told apart from the job simply not using the cores it was given.
        wrapped = ["/usr/bin/time", "-f", "%e %M %U %S %P",
                   "-o", str(rss_file)] + wrapped

    counted = perf_usable()
    if counted:
        # instructions and cycles are the only honest measure of work here:
        # CPU seconds move with the clock, which falls as cores light up.
        perf_file.parent.mkdir(parents=True, exist_ok=True)
        wrapped = ["perf", "stat", "-e", PERF_EVENTS, "-x", PERF_SEP,
                   "-o", str(perf_file)] + wrapped

    run_env = dict(os.environ)
    if extra_env:
        run_env.update((k, str(v)) for k, v in extra_env.items())

    # The whole run goes in its own process group, and a timeout kills the
    # GROUP. subprocess.run(timeout=...) only kills its direct child, which
    # here is /usr/bin/time -- the bench_remesh underneath it survives, keeps
    # every core busy, and corrupts whatever is measured next. That is not
    # hypothetical: a 24-core run recorded three cells as 240s timeouts whose
    # JSONs said "success" at 259s, and the two measurements that followed came
    # out 8% slow because an orphan was still running beside them.
    try:
        throttle_before = _read_throttle_counts()
        sampler = _FreqSampler(n_busy=threads)
        sampler.start()
    except Exception:
        throttle_before, sampler = None, _FreqSampler(n_busy=threads)

    t0 = time.time()
    started_at = time.strftime("%Y-%m-%dT%H:%M:%S")
    timed_out = False
    with open(log_path, "w") as log:
        proc = subprocess.Popen(wrapped, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True, env=run_env)
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            _kill_group(proc)
            rc = proc.returncode if proc.returncode is not None else -9
    wall = time.time() - t0

    # The wrappers' exit status is not the remesher's. Prefer what the shell
    # recorded; fall back to the outermost wrapper only if the whole group was
    # killed before it could write. A status of 128+n is a death by signal n,
    # and is normalised to a negative rc so the sweep's crash check still sees
    # it as one.
    inner_rc = None
    try:
        raw = status_file.read_text().strip()
        if raw:
            inner_rc = int(raw)
    except Exception:
        pass
    if inner_rc is None:
        inner_rc = rc
    elif inner_rc >= 128:
        inner_rc = -(inner_rc - 128)
    rc = inner_rc

    perf_cols = read_perf(perf_file) if counted else {}

    # Instrumentation must never cost a measurement. A 12-hour unattended run
    # losing a cell because a sysfs read misbehaved would be a bad trade for
    # diagnostics, so anything that goes wrong here is simply absent data.
    try:
        freq = sampler.stop()
    except Exception:
        freq = {}
    try:
        throttle_after = _read_throttle_counts()
        throttled = ("" if throttle_before is None or throttle_after is None
                     else throttle_after - throttle_before)
    except Exception:
        throttled = ""

    peak_rss, user_s, sys_s, cpu_pct = "", "", "", ""
    if rss_file.exists():
        try:
            parts = rss_file.read_text().split()
            # "%e %M %U %S %P" -> elapsed, maxrss, user, sys, cpu%
            if len(parts) >= 5:
                peak_rss, user_s, sys_s = parts[1], parts[2], parts[3]
                cpu_pct = parts[4].rstrip("%")
            elif parts:
                peak_rss = parts[-1]
        except Exception:
            pass

    # CPU time over wall time: the average number of cores actually working.
    # This is the number that separates "does not scale" from "never used the
    # cores", and it needs no trust in the scheduler or the governor.
    avg_par = ""
    try:
        if user_s and sys_s and wall > 0:
            avg_par = round((float(user_s) + float(sys_s)) / wall, 2)
    except Exception:
        pass

    if timed_out:
        # Leave nothing behind that a later run would mistake for a finished
        # cell. The binary may have written its JSON in the moments between the
        # deadline and the kill; that JSON describes a run we cut short and
        # whose timing is meaningless, and keeping it would make the cell
        # permanently skipped on resume.
        json_path.unlink(missing_ok=True)
        if out_mesh:
            Path(out_mesh).unlink(missing_ok=True)
            Path(out_mesh + ".partial").unlink(missing_ok=True)

    remesh_s = ""
    if json_path.exists():
        try:
            d = json.loads(json_path.read_text())
            remesh_s = d["metrics"]["Performance"]["Total_Time"]["Value"]
        except Exception:
            pass

    row = {"mesh": mesh_id, "factor": factor, "arm": arm, "threads": threads,
           "rep": rep, "wall_s": round(wall, 3), "remesh_s": remesh_s,
           "peak_rss_kb": peak_rss, "rc": rc, "timed_out": int(timed_out),
           "started_at": started_at, "user_s": user_s, "sys_s": sys_s,
           "cpu_pct": cpu_pct, "avg_parallelism": avg_par,
           "throttle_events": throttled,
           "instructions": "", "cycles": "", "task_clock_ms": "",
           "variant": variant,
           "pin": ",".join(str(c) for c in pin_cpus) if pin_cpus else "",
           "lock_grid": (extra_env or {}).get(LOCK_GRID_ENV, ""),
           "out_mesh": "" if timed_out else out_mesh, "json": str(json_path)}
    row.update(freq)
    row.update(perf_cols)
    return row


class Sweep:
    def __init__(self, root, results_dir, budget, args):
        self.root = root
        self.results_dir = results_dir
        self.budget = budget
        self.args = args
        self.t0 = time.time()
        self.csv_path = results_dir / "results.csv"
        new = not self.csv_path.exists()
        self.done = self._load_done()
        # Append under the header the file already has, not the one this version
        # of the script would write. A results directory started before the
        # perf columns existed has a shorter header, and writing wider rows into
        # it would silently shift every column from there on.
        fields = CSV_HEADER
        if not new:
            try:
                with open(self.csv_path, newline="") as fh:
                    old = next(csv.reader(fh))
                if old and set(old) != set(CSV_HEADER):
                    fields = old
                    missing = [c for c in CSV_HEADER if c not in old]
                    if missing:
                        print("Note: %s has an older header, so these columns "
                              "cannot be recorded here: %s.\n"
                              "      Use a fresh --results-dir to get them."
                              % (self.csv_path.name, ", ".join(missing)),
                              file=sys.stderr)
            except Exception:
                pass
        self.csv_fh = open(self.csv_path, "a", newline="")
        self.csv = csv.DictWriter(self.csv_fh, fieldnames=fields,
                                  extrasaction="ignore")
        if new:
            self.csv.writeheader()
            self.csv_fh.flush()

    def _load_done(self):
        """Prior rows, keyed the same way run_cell names its files, so a resumed
        run can read back what earlier invocations measured."""
        done = {}
        if not self.csv_path.exists():
            return done
        with open(self.csv_path, newline="") as fh:
            for r in csv.DictReader(fh):
                try:
                    stem = cell_stem(r["mesh"], r["factor"], r["arm"],
                                     int(r["threads"]), int(r["rep"]),
                                     r.get("variant") or "")
                    r["wall_s"] = float(r["wall_s"])
                    r["threads"] = int(r["threads"])
                    r["rep"] = int(r["rep"])
                    r["rc"] = int(r["rc"])
                    r["timed_out"] = int(r["timed_out"])
                    done[stem] = r
                except Exception:
                    continue
        if done:
            print("Resuming: %d cells already recorded." % len(done))
        return done

    def left(self):
        return self.budget - (time.time() - self.t0) if self.budget else float("inf")

    def is_cached(self, mesh_path, factor, arm, threads, rep, variant=""):
        """Whether this cell is already recorded. The budget guards consult this
        first: a cached cell costs no wall time, so a resumed run must not be
        stopped by the budget before it has read its earlier results back."""
        stem = cell_stem(Path(mesh_path).stem, factor, arm, threads, rep, variant)
        if stem not in self.done:
            return False
        # A CSV row is not enough. A timed-out cell has a row but its JSON was
        # deleted, so run_cell will really execute it -- and charging the budget
        # nothing for a run that is about to take the full timeout is how a
        # 30-minute calibration turns into an hour.
        return (self.results_dir / "json" / (stem + ".json")).exists()

    def affordable(self, seconds, mesh_path, factor, arm, threads, rep,
                   variant=""):
        if self.is_cached(mesh_path, factor, arm, threads, rep, variant):
            return True
        return self.left() >= seconds

    def emit(self, row):
        if row is None or row.get("cached"):
            return          # already in the CSV; do not duplicate the row
        self.csv.writerow(row)
        self.csv_fh.flush()

        # The goal requires zero crashes. A signal death means something is
        # actually wrong, so stop the whole sweep and let it be fixed rather
        # than collecting a run full of holes. A timeout is not a crash.
        if row["rc"] < 0 and not row["timed_out"]:
            print("\nCRASH: %s (rc=%d). Stopping the sweep." %
                  (row["json"], row["rc"]), file=sys.stderr)
            sys.exit(42)

    def close(self):
        self.csv_fh.close()


# --------------------------------------------------------------------------
# profiles
# --------------------------------------------------------------------------

def load_manifest(root, mesh_dir, pipelines=None):
    mp = mesh_dir / "manifest.json"
    if not mp.exists():
        sys.exit("No %s -- run scripts/prepare_meshes.py first." % mp)
    meshes = json.loads(mp.read_text())["meshes"]
    ok = [m for m in meshes.values() if m.get("cells") and Path(m["path"]).exists()]
    if not ok:
        sys.exit("manifest.json has no usable entries. If it predates the second\n"
                 "input pipeline, re-run prepare_meshes.py to regenerate it.")
    if pipelines:
        ok = [m for m in ok if m.get("pipeline") in pipelines]
    # Each entry's "key" is <id>_<pipeline>, which is also its .mesh stem, so the
    # pipeline travels with every result row for free.
    for m in ok:
        m.setdefault("key", Path(m["path"]).stem)
        m.setdefault("pipeline", "cdt")
    return sorted(ok, key=lambda m: m["cells"], reverse=True)


def interleave_pipelines(meshes, n):
    """Pick n ladder inputs covering every pipeline and spanning the size range.

    Two things this has to get right.

    First, both pipelines must appear. They are not interchangeable -- earlier
    measurements in this project had CDT and Mesh_3 inputs disagreeing on the
    sign of the sequential comparison -- so calibrating on one and running on
    both would size the grid from the wrong generator.

    Second, do NOT simply take the largest by cell count. Remeshing time does
    not track input size: measured here, a 65k-cell input ran 11.9s at its
    finest ladder rung while a 32k-cell one ran 45.5s at the same rung, because
    the target edge length is derived from the input's own average edge length
    and so depends on its geometry. Handing the ladder the n biggest inputs can
    therefore miss the bar entirely while a smaller input would have cleared it.
    Sampling across the upper size range hedges against that.
    """
    by_pipe = {}
    for m in meshes:
        by_pipe.setdefault(m["pipeline"], []).append(m)

    per_pipe = max(1, n // max(1, len(by_pipe)))
    spread = {}
    for p, ms in by_pipe.items():
        ms = sorted(ms, key=lambda m: m["cells"], reverse=True)
        # Sample from the larger half -- small inputs cannot produce a long run
        # at any edge factor -- but spread across it rather than clustering at
        # the top.
        pool = ms[:max(per_pipe, len(ms) // 2)]
        if len(pool) <= per_pipe:
            spread[p] = pool
        else:
            step = (len(pool) - 1) / float(per_pipe - 1) if per_pipe > 1 else 1
            spread[p] = [pool[int(round(i * step))] for i in range(per_pipe)]

    order = sorted(spread)
    picked, i = [], 0
    while len(picked) < n and any(spread[p] for p in order):
        p = order[i % len(order)]
        if spread[p]:
            picked.append(spread[p].pop(0))
        i += 1
    return picked


def profile_calibrate(sweep, bins, meshes, args):
    """Ordered by value, so a truncated run is still useful."""
    rd, md = sweep.results_dir, sweep.mesh_out_dir
    max_t = args.threads_max
    common = dict(results_dir=rd, mesh_out_dir=md, settle=args.settle, done=sweep.done)

    print("\n=== 1/4  size ladder: find (mesh, factor) with %d-thread wall >= %.0fs ==="
          % (max_t, TARGET_SECONDS))
    cleared = []          # (cells, mesh, factor, wall)
    measured = []         # every ladder cell that completed: (wall, mesh, path, factor)
    # Smallest input first. The ladder's job is to find where the 30s bar sits,
    # and starting at the largest input means the very first probe can blow the
    # per-run timeout and tell us nothing except "too big" -- which is what
    # happened on the 24-core run, where the largest mesh timed out at factor
    # 0.5 before anything useful had been measured.
    ladder_meshes = list(reversed(interleave_pipelines(meshes, args.calib_meshes)))
    print("  ladder inputs: %s" % ", ".join(m["key"] for m in ladder_meshes))
    for m in ladder_meshes:
        for factor in LADDER:
            if not sweep.affordable(args.run_timeout * 0.5, m["path"], factor,
                                    "par", max_t, 1):
                print("  budget spent, stopping the ladder")
                break
            row = run_cell(bins, m["path"], factor, "par", max_t, 1,
                           write_mesh=True, timeout=args.run_timeout, **common)
            sweep.emit(row)
            print("  %-10s f=%-4s  %6.1fs%s%s" %
                  (m["key"], factor, row["wall_s"],
                   "  TIMEOUT" if row["timed_out"] else "",
                   "  (cached)" if row.get("cached") else ""))
            if row["timed_out"]:
                # We have bracketed it: this factor is already too slow, so
                # descending further only wastes budget.
                break
            measured.append((row["wall_s"], m["key"], m["path"], factor))
            if row["wall_s"] >= TARGET_SECONDS:
                cleared.append((m["cells"], m["key"], m["path"], factor, row["wall_s"]))
                break

    fell_back = False
    if not cleared:
        print("\nNo config reached %.0fs at %d threads (longest: %.1fs)."
              % (TARGET_SECONDS, max_t, max(m[0] for m in measured) if measured else 0))
        print("These inputs are smaller than this machine needs. The proper fix is")
        print("bigger surfaces (relax the 10k-vertex Thingi filter), not a finer")
        print("edge factor -- the ladder already went to %s." % LADDER[-1])

        if not measured:
            print("\nNothing completed at all, so there is no grid to fall back to.")
            return None

        # Do NOT abort. On a single booked run, returning nothing because the
        # bar was missed is far worse than returning slightly noisier numbers:
        # the longest cells measured are still well clear of the noise floor,
        # and the alternative is an empty 12 hours. Note the shortfall loudly
        # and carry on with the best available.
        fell_back = True
        best = sorted(measured, reverse=True)[:args.max_configs]
        cleared = [(0, mid, mpath, mfac, mwall) for mwall, mid, mpath, mfac in best]
        print("\nFALLING BACK to the %d longest-running config(s) measured, from "
              "%.1fs down to %.1fs.\nTreat the per-config numbers as noisier than "
              "usual and read the CV column." % (len(cleared), cleared[0][4],
                                                 cleared[-1][4]))

    # The *cheapest* config that cleared the bar: enough signal, least budget.
    cleared.sort(key=lambda c: c[4])
    _, pid, ppath, pfactor, pwall = cleared[0]
    print("\n  pilot config: %s @ factor %s (%.1fs at %d threads)" %
          (pid, pfactor, pwall, max_t))

    print("\n=== 2/4  noise: 5 repeats of the pilot config ===")
    walls = [pwall]
    for rep in range(2, 7):
        if not sweep.affordable(pwall * 1.5, ppath, pfactor, "par", max_t, rep):
            print("  budget spent, stopping repeats")
            break
        row = run_cell(bins, ppath, pfactor, "par", max_t, rep,
                       write_mesh=False, timeout=args.run_timeout, **common)
        sweep.emit(row)
        walls.append(row["wall_s"])
        print("  rep %d: %6.1fs%s" %
              (rep, row["wall_s"], "  (cached)" if row.get("cached") else ""))

    cv = None
    if len(walls) > 1:
        mean = sum(walls) / len(walls)
        var = sum((w - mean) ** 2 for w in walls) / len(walls)
        cv = (var ** 0.5) / mean if mean else None
        print("  n=%d  mean=%.1fs  CV=%.2f%%" % (len(walls), mean, cv * 100))

    # Inside --profile full this stage is skipped: wave 2 of the sweep measures
    # the same curve on every config, and here it is the most expensive stage by
    # far (at 24 threads the 2-thread point alone costs ~12x the pilot). Spending
    # the calibration budget on it would starve stage 4, and stage 4 is what
    # produces the serial ratio the sweep needs to plan with.
    curve = {max_t: pwall}
    curve_threads = []
    if getattr(args, "skip_curve", False):
        print("\n=== 3/4  skipped (the sweep's wave 2 measures the curve) ===")
    else:
        print("\n=== 3/4  scaling curve on the pilot config ===")
        # Descending, so the CHEAPEST points run first. A 67s-at-24-threads
        # pilot costs roughly 800s at 2 threads; going up from 2 spends the
        # whole budget on the single most expensive point and returns a
        # one-point "curve". Working down from max_t fits several points in.
        curve_threads = sorted((x for x in THREAD_LIST if 1 < x < max_t),
                               reverse=True)

    for t in curve_threads:
        est = pwall * max_t / t
        if not sweep.affordable(est * 1.3, ppath, pfactor, "par", t, 1):
            print("  budget spent, skipping %d threads (%.0fs) and below" % (t, est))
            break
        row = run_cell(bins, ppath, pfactor, "par", t, 1,
                       write_mesh=True, timeout=args.run_timeout, **common)
        sweep.emit(row)
        curve[t] = row["wall_s"]
        print("  %2d threads: %6.1fs  (speedup vs %d-thread: %.2fx)%s" %
              (t, row["wall_s"], max_t, pwall / max(row["wall_s"], 1e-9),
               "  (cached)" if row.get("cached") else ""))

    # The anchors need a DELIBERATELY CHEAP config, not the pilot and not even
    # the fastest ladder cell. One thread costs roughly max_t times the parallel
    # wall, so a 67s-at-24-threads pilot is ~27 minutes single-threaded: it
    # cannot fit a 30-minute budget and would not even finish inside the
    # per-run timeout. On the 24-core run all three anchors timed out for
    # exactly this reason, and calibration came back with no serial ratio.
    #
    # So: smallest input, coarse edge factor. The ratio is an approximation
    # taken on a small workload either way -- what matters is that it is a
    # measurement rather than a timeout.
    anchors, anchor, serial_ratio = {}, None, None

    # Pick from what the ladder ALREADY measured: the LARGEST cell whose three
    # single-threaded arms still fit the remaining budget and the per-run
    # timeout. Largest-that-fits, not smallest: a sub-second run is fixed
    # overhead, not remeshing, and produced a nonsense 1.75x ratio in testing.
    # This also costs nothing extra, since the ladder cell is already cached.
    candidates = [c for c in sorted(measured, reverse=True)
                  if c[0] >= ANCHOR_MIN_SECONDS]
    anchor_pick = None
    for awall, aid, apath, afactor in candidates:
        need_all = awall * max_t * 3 * 1.3          # three arms at one thread
        if awall * max_t * 1.3 < args.run_timeout and sweep.left() > need_all:
            anchor_pick = (awall, aid, apath, afactor)
            break

    if anchor_pick is None:
        print("\n=== 4/4  skipped: no measured config is both long enough to time\n"
              "        reliably (>= %.0fs) and short enough to run single-threaded\n"
              "        inside the remaining budget ===" % ANCHOR_MIN_SECONDS)
    else:
        awall, aid, apath, afactor = anchor_pick
        anchor = {"mesh": aid, "factor": afactor, "wall_at_max_t": awall}
        print("\n=== 4/4  sequential anchors on %s f=%s (%.1fs at %d threads) ==="
              % (aid, afactor, awall, max_t))
        print("  one thread should take roughly %.0fs; per-run timeout is %ds"
              % (awall * max_t, args.run_timeout))
        need = awall * max_t * 1.5
        for arm in ("par", "seq", "main"):
            if not sweep.affordable(need, apath, afactor, arm, 1, 1):
                print("  budget spent, skipping %s (and the rest)" % arm)
                break
            row = run_cell(bins, apath, afactor, arm, 1, 1,
                           write_mesh=True, timeout=args.run_timeout, **common)
            sweep.emit(row)
            if row["timed_out"]:
                print("  %-5s @1: TIMED OUT at %ds -- not recorded as an anchor"
                      % (arm, args.run_timeout))
                continue
            anchors[arm] = row["wall_s"]
            print("  %-5s @1: %6.1fs%s" %
                  (arm, row["wall_s"], "  (cached)" if row.get("cached") else ""))

        # How much slower one thread is than max_t. The overnight ETA needs it:
        # the seq and par@1 arms dominate that run, and assuming they scale down
        # perfectly from the max-thread number underestimates it badly.
        #
        # Only ever computed from runs that COMPLETED. A timed-out run's wall
        # time is the timeout, not a measurement, and a ratio built from it is a
        # floor dressed up as a number.
        serial_ratio = None
        if "par" in anchors and awall > 0:
            serial_ratio = anchors["par"] / awall
            print("\n  par@1 / par@%d = %.2fx" % (max_t, serial_ratio))
        else:
            print("\n  serial ratio unknown (par@1 did not complete); the "
                  "overnight ETA will assume %dx." % max_t)

    configs = [{"mesh": c[1], "path": c[2], "factor": c[3],
                "cells": c[0], "wall_at_max_t": c[4]} for c in cleared]

    out = sweep.results_dir / "calibration.json"
    out.write_text(json.dumps({
        "target_seconds": TARGET_SECONDS,
        "met_target": not fell_back,
        "threads_max": max_t,
        "ladder": LADDER,
        "configs": configs,
        "pilot": {"mesh": pid, "factor": pfactor, "wall_at_max_t": pwall},
        "cv": cv,
        "scaling_curve_seconds": curve,
        "anchor_config": anchor,
        "anchor_walls": anchors,
        "serial_ratio": serial_ratio,
        "ladder_measured": [{"mesh": m[1], "factor": m[3], "wall": m[0]}
                            for m in sorted(measured)],
    }, indent=2))
    print("\nWrote %s -- send this back; it sizes the overnight grid." % out)
    return out


def profile_overnight(sweep, bins, meshes, args):
    calib = sweep.results_dir / "calibration.json"
    if args.configs:
        cfg = json.loads(Path(args.configs).read_text())
    elif calib.exists():
        cfg = json.loads(calib.read_text())
    else:
        sys.exit("No configs file and no calibration.json -- run --profile calibrate first.")

    configs = cfg["configs"][:args.max_configs]
    threads = [t for t in THREAD_LIST if t <= args.threads_max]
    rd, md = sweep.results_dir, sweep.mesh_out_dir
    common = dict(results_dir=rd, mesh_out_dir=md, settle=args.settle, done=sweep.done)

    # An honest ETA. The single-threaded arms dominate this run and do NOT scale
    # down perfectly from the max-thread number, so use the ratio calibration
    # actually measured; fall back to the thread count only if it is missing.
    ratio = cfg.get("serial_ratio") or float(args.threads_max)

    def est_seconds(w, arm, t):
        if arm == "par" and t == args.threads_max:
            return w
        # interpolate between 1x at max_t and `ratio`x at 1 thread
        return w * (ratio ** ((args.threads_max / max(t, 1) - 1)
                              / max(args.threads_max - 1, 1)))

    est = 0.0
    for c in configs:
        w = c["wall_at_max_t"]
        for t in threads:
            est += est_seconds(w, "par", t) * args.reps
        est += est_seconds(w, "seq", 1) * args.reps      # seq
        est += est_seconds(w, "main", 1) * args.reps     # main
    print("Configs: %d   threads: %s   reps: %d" % (len(configs), threads, args.reps))
    print("Serial ratio from calibration: %s" % (
        ("%.2fx" % ratio) if cfg.get("serial_ratio") else "unknown, assuming %.0fx" % ratio))
    print("Rough ETA: %.1f h" % (est / 3600))
    if args.budget:
        print("Budget:    %.1f h" % (args.budget / 3600))
        if est > args.budget:
            print("WARNING: the grid does not fit the budget. It will stop partway,\n"
                  "         which is safe -- re-running resumes -- but consider\n"
                  "         --max-configs or --reps to get a complete picture.")

    # Run in PRIORITY WAVES rather than config by config. Config-major order
    # means that if the budget runs out -- and on a single unattended run
    # nobody is there to notice -- the last configs have nothing at all, while
    # the first have every thread count and every repeat. Waves invert that:
    # every config gets the headline comparison before any config gets its
    # second repeat, so a truncated run is still a complete, reportable dataset,
    # just a coarser one.
    waves = [
        # 1: the headline. Speedup at full width, the parallel-overhead point,
        #    and both reference arms. Without these there is no report at all.
        ("headline", [("par", args.threads_max, 1), ("main", 1, 1),
                      ("seq", 1, 1), ("par", 1, 1)]),
        # 2: the shape of the scaling curve, cheapest thread counts first.
        ("curve", [("par", t, 1) for t in sorted(
            (x for x in threads if 1 < x < args.threads_max), reverse=True)]),
        # 3: repeats, so the headline numbers get error bars before the
        #    intermediate points get theirs.
        ("repeats", [("par", args.threads_max, r) for r in range(2, args.reps + 1)]
                    + [("main", 1, r) for r in range(2, args.reps + 1)]
                    + [("seq", 1, r) for r in range(2, args.reps + 1)]
                    + [("par", 1, r) for r in range(2, args.reps + 1)]),
        # 4: everything else.
        ("curve repeats", [("par", t, r)
                           for t in sorted((x for x in threads
                                            if 1 < x < args.threads_max), reverse=True)
                           for r in range(2, args.reps + 1)]),
    ]

    # Cheapest configs first inside a wave, for the same reason.
    ordered = sorted(configs, key=lambda c: c["wall_at_max_t"])

    for wave_name, cells in waves:
        print("\n=== wave: %s ===" % wave_name)
        for arm, t, rep in cells:
            for c in ordered:
                path, factor = c["path"], c["factor"]
                need = est_seconds(c["wall_at_max_t"], arm, t)
                if not sweep.affordable(need, path, factor, arm, t, rep):
                    print("\nBudget spent during wave '%s'. Stopping cleanly.\n"
                          "Everything from the earlier waves is complete; re-running "
                          "resumes here." % wave_name)
                    return None
                row = run_cell(bins, path, factor, arm, t, rep,
                               write_mesh=(rep == 1),
                               timeout=args.run_timeout, **common)
                sweep.emit(row)
                print("  %-14s %-5s t=%-3d r=%d  %8.1fs%s%s" %
                      (c["mesh"], arm, t, rep, row["wall_s"],
                       "  TIMEOUT" if row["timed_out"] else "",
                       "  (cached)" if row.get("cached") else ""))
    return None


def profile_full(sweep, bins, meshes, args, root, mesh_out_dir):
    """The whole MEASUREMENT in one unattended command: calibrate, then sweep.

    This exists because the benchmark machine is someone else's and we get one
    run on it. The two-step flow (calibrate, send the file back, we size the
    grid, they start the sweep) needs a person in the loop partway through;
    here the grid is sized from calibration in-process.

    Deliberately NOT included: the quality pass and the report. Those are
    post-processing -- they read the meshes this run leaves behind, they can be
    redone any number of times, and mixing them into the measurement would put
    work that does not need the machine idle inside the window that does.

    The budget is split rather than shared, so calibration cannot overrun into
    the sweep's time.
    """
    total = args.budget or 12 * 3600
    # The cap is a share, not a tenth. At 10% a short run gave calibration 90s
    # -- less than one ladder rung on a fast machine -- so it fell back before
    # measuring anything useful. --calib-budget is the real control; this only
    # stops calibration from eating a run whole.
    calib_budget = min(args.calib_budget, total * 0.5)
    sweep_budget = total - calib_budget

    if sweep_budget < MIN_SWEEP_SECONDS:
        sys.exit(
            "A %.0f s budget leaves only %.0f s for the sweep after calibration\n"
            "(%.0f s). That is not enough to measure anything. Give --budget at\n"
            "least %.0f s, or run the phases separately with --profile calibrate\n"
            "and --profile overnight."
            % (total, sweep_budget, calib_budget,
               MIN_SWEEP_SECONDS + calib_budget))

    print("=" * 70)
    print("Measurement plan for a %.1f h budget:" % (total / 3600))
    print("  calibration   %5.2f h   find configs that clear the %.0fs bar" %
          (calib_budget / 3600, TARGET_SECONDS))
    print("  sweep         %5.2f h   the grid, in priority waves" % (sweep_budget / 3600))
    print("Nothing needs to be sent back or decided mid-run.")
    print("Quality metrics and the report come afterwards, from the saved meshes.")
    print("=" * 70)

    # --- phase 1: calibration, on its own clock -----------------------------
    sweep.budget = calib_budget
    sweep.t0 = time.time()
    print("\n\n######## PHASE 1/2: CALIBRATION ########")
    # Calibration needs a much tighter per-run timeout than the sweep: its job
    # is to bracket where the 30s bar sits, and one probe allowed to run for
    # the sweep's 4-hour limit would swallow the whole calibration budget.
    sweep_timeout = args.run_timeout
    args.run_timeout = min(args.run_timeout, int(calib_budget / 3))
    args.skip_curve = True
    try:
        profile_calibrate(sweep, bins, meshes, args)
    finally:
        args.run_timeout = sweep_timeout
        args.skip_curve = False

    calib_path = sweep.results_dir / "calibration.json"
    if not calib_path.exists():
        print("\nCalibration produced no usable configs, so there is no grid to run.\n"
              "The ladder output above says why. Stopping rather than guessing.",
              file=sys.stderr)
        return None

    # --- phase 2: the sweep -------------------------------------------------
    spent = time.time() - sweep.t0
    sweep.budget = sweep_budget + max(0.0, calib_budget - spent)   # hand back any slack
    sweep.t0 = time.time()
    print("\n\n######## PHASE 2/2: SWEEP  (%.2f h) ########" % (sweep.budget / 3600))
    profile_overnight(sweep, bins, meshes, args)

    print("\n" + "=" * 70)
    print("MEASUREMENT DONE. The machine is free now -- the next step does not")
    print("need it idle, and can run whenever:\n")
    print("    python3 scripts/run_bench.py --root %s --quality-pass\n" % root)
    print("That derives every quality metric from the meshes this run saved and")
    print("writes one small tarball to send back.")
    print("=" * 70)
    return None


# --------------------------------------------------------------------------
# the metrics profile
# --------------------------------------------------------------------------

def _metrics_lookup(meshes, wanted):
    """Resolve hard-coded (mesh key, factor) pairs against the prepared inputs.

    A pair whose mesh was never prepared is reported and dropped rather than
    made fatal: three ladders minus one is still most of the answer, and a run
    that refuses to start because one input is missing wastes the machine.
    """
    by_key = dict((m["key"], m) for m in meshes)
    found, missing = [], []
    for key, factor in wanted:
        if key in by_key:
            found.append((by_key[key], factor))
        else:
            missing.append(key)
    if missing:
        print("  NOT PREPARED, skipping: %s" % ", ".join(sorted(set(missing))),
              file=sys.stderr)
    return found


def _metrics_run(sweep, bins, mesh, factor, threads, rep, label, common,
                 args, write_mesh=False, **kw):
    """One metrics cell, with the budget guard and a one-line report."""
    variant = kw.get("variant", "")
    if not sweep.affordable(args.run_timeout, mesh["path"], factor, "par",
                            threads, rep, variant):
        print("  budget spent, stopping")
        return None
    row = run_cell(bins, mesh["path"], factor, "par", threads, rep,
                   write_mesh=write_mesh, timeout=args.run_timeout,
                   **dict(common, **kw))
    sweep.emit(row)
    # A cached row comes back from the CSV, where every field is a string.
    try:
        instr = "%.1f G instr" % (float(row.get("instructions")) / 1e9)
    except (TypeError, ValueError):
        instr = "no counters"
    print("  %-26s t=%-3d r=%d  %8.1fs  %s%s%s" % (
        label, threads, rep, row["wall_s"], instr,
        "  TIMEOUT" if row["timed_out"] else "",
        "  (cached)" if row.get("cached") else ""), flush=True)
    return row


def metrics_ladder(sweep, bins, meshes, args, common):
    """1, 2, 4, 8, 12, 16, 24 threads, three reps, reps INTERLEAVED.

    Interleaved means rep 1 of every point, then rep 2, then rep 3 -- not three
    reps of one point in a row. Machines drift: 25% between measurement waves
    was measured on the test laptop here. Running all three reps of the 1-thread
    point together buries that drift inside one point's error bar and moves the
    whole point relative to the others, which is exactly the comparison the
    ladder exists to make. Interleaving spreads any drift across every point.

    Within a rep the thread counts descend, so the cheapest runs go first and a
    truncated rep still covers the wide end.
    """
    configs = _metrics_lookup(meshes, METRICS_LADDER_CONFIGS)
    if not configs:
        return
    threads = [t for t in sorted(THREAD_LIST, reverse=True)
               if t <= args.threads_max]
    print("\n=== thread ladder: %d config(s) x %s threads x %d reps ==="
          % (len(configs), threads, args.reps))
    for rep in range(1, args.reps + 1):
        for t in threads:
            for mesh, factor in configs:
                # One output mesh at each end of the first ladder, and no more.
                # Parallel against sequential quality is already settled at
                # 0.25% on cell count; this is a check that nothing regressed,
                # not a quality sweep.
                write = (rep == 1 and (mesh, factor) == configs[0]
                         and t in (1, args.threads_max))
                _metrics_run(sweep, bins, mesh, factor, t, rep,
                             "%s f=%s" % (mesh["key"], factor), common, args,
                             write_mesh=write)


def metrics_pinning(sweep, bins, meshes, args, common):
    """Unpinned against the performance cores against the efficiency cores.

    A hybrid CPU's "24 cores" are not 24 of the core the 1-thread run measured,
    so part of a scaling shortfall can simply be work landing on slower cores.
    Eight threads on eight performance cores against eight threads unpinned
    separates that from lock contention: if the pinned run wins, the scheduler's
    placement is a real part of the loss.
    """
    configs = _metrics_lookup(meshes, [METRICS_PIN_CONFIG])
    if not configs:
        return
    mesh, factor = configs[0]
    t = min(METRICS_PIN_THREADS, args.threads_max)
    fast, slow = core_sets(t)
    if not fast:
        print("\n=== core pinning: skipped, this machine exposes no "
              "cpuinfo_max_freq ===", file=sys.stderr)
        return

    arms = [("unpinned", None)]
    if fast == slow:
        print("\n=== core pinning: every core has the same maximum clock, so "
              "there is no performance/efficiency split to test. Running the "
              "unpinned arm only. ===")
    else:
        arms += [("pcore%d" % t, fast), ("ecore%d" % t, slow)]
        print("\n=== core pinning on %s f=%s at %d threads ===" %
              (mesh["key"], factor, t))
        print("  performance cores: %s" % ",".join(str(c) for c in fast))
        print("  efficiency  cores: %s" % ",".join(str(c) for c in slow))

    for rep in range(1, args.reps + 1):
        for name, cpus in arms:
            _metrics_run(sweep, bins, mesh, factor, t, rep,
                         "%s %s" % (mesh["key"], name), common, args,
                         variant=name, pin_cpus=cpus)


def metrics_lock_grid(sweep, bins, meshes, args, common):
    """Sweep CGAL_TETRAHEDRAL_REMESHING_LOCK_GRID at the widest thread count.

    The built-in grid size was chosen at four threads. The environment variable
    is read once per remesher and never on a locking path, so setting it costs
    nothing and the sweep is what says whether that choice still holds at 24.
    """
    configs = _metrics_lookup(meshes, METRICS_LOCK_CONFIGS)
    if not configs:
        return
    t = args.threads_max
    print("\n=== lock grid %s at %d threads ===" %
          (METRICS_LOCK_GRIDS, t))
    for rep in range(1, args.reps + 1):
        for grid in METRICS_LOCK_GRIDS:
            for mesh, factor in configs:
                _metrics_run(sweep, bins, mesh, factor, t, rep,
                             "%s grid=%d" % (mesh["key"], grid), common, args,
                             variant="lg%d" % grid,
                             extra_env={LOCK_GRID_ENV: grid})


def metrics_diagnostics(sweep, bins, meshes, args, common, results_dir):
    """The two instrumented binaries, once each at one thread and at the widest.

    These are diagnostics, not timings: each is the same source with one macro
    added, and their wall times are not comparable with the ladder's. What comes
    back is their stdout, saved verbatim.

    CGAL_TR_LOCKCOUNT prints its one line from a static destructor at exit, so a
    run that crashed prints nothing at all. An empty capture is therefore a
    failure and not a quiet pass, which is why the exit status is recorded
    beside it.
    """
    diag = {}
    for name in ("bench_remesh_topstage", "bench_remesh_lockcount"):
        if name in bins and Path(bins[name]).exists():
            diag[name] = Path(bins[name])
    if not diag:
        print("\n=== diagnostics: skipped, the instrumented binaries were not "
              "built. Re-run scripts/setup.py --diagnostics. ===", file=sys.stderr)
        return

    configs = _metrics_lookup(meshes, [METRICS_DIAG_CONFIG])
    if not configs:
        return
    mesh, factor = configs[0]
    diag_dir = results_dir / "diagnostics"
    diag_dir.mkdir(parents=True, exist_ok=True)
    threads = [t for t in METRICS_DIAG_THREADS if t <= args.threads_max]
    print("\n=== diagnostics on %s f=%s at %s threads ==="
          % (mesh["key"], factor, threads))

    summary = {}
    spath = diag_dir / "summary.json"
    if spath.exists():
        try:
            summary = json.loads(spath.read_text())
        except Exception:
            summary = {}

    for name, exe in sorted(diag.items()):
        variant = name.replace("bench_remesh_", "")
        for t in threads:
            row = _metrics_run(sweep, bins, mesh, factor, t, 1,
                               "%s %s" % (mesh["key"], variant), common, args,
                               variant=variant, exe=exe)
            if row is None:
                continue
            stem = cell_stem(mesh["key"], factor, "par", t, 1, variant)
            src = results_dir / "logs" / (stem + ".log")
            dst = diag_dir / (stem + ".out")
            try:
                shutil.copyfile(src, dst)
            except Exception:
                pass
            n_lines = 0
            try:
                n_lines = sum(1 for line in dst.read_text().splitlines()
                              if line.strip())
            except Exception:
                pass
            summary[stem] = {"binary": name, "threads": t, "rc": row["rc"],
                             "timed_out": row["timed_out"],
                             "stdout": dst.name, "stdout_lines": n_lines}
            if row["rc"] != 0 or n_lines == 0:
                print("  WARNING: %s exited %d and printed %d line(s). An empty "
                      "capture is a failed run, not a clean one."
                      % (stem, row["rc"], n_lines), file=sys.stderr)
    spath.write_text(json.dumps(summary, indent=2))
    print("  stdout captures and exit statuses in %s" % diag_dir)


def profile_metrics(sweep, bins, meshes, args, root, mesh_out_dir):
    """The scaling questions, and nothing else (docs/METRICS_REQUEST.md).

    About 45 minutes on a 24-core machine, against the 12 hours the full sweep
    takes, because it runs only what a scaling analysis reads: three thread
    ladders, a core-pinning comparison, a lock-grid sweep and two instrumented
    runs. No edge-factor ladder, no extra meshes, and no `seq` or `main`
    reference arms -- those run at one thread by definition and so cannot move
    with anything measured here.
    """
    common = dict(results_dir=sweep.results_dir, mesh_out_dir=mesh_out_dir,
                  settle=args.settle, done=sweep.done)
    phases = {"ladder": lambda: metrics_ladder(sweep, bins, meshes, args, common),
              "pinning": lambda: metrics_pinning(sweep, bins, meshes, args, common),
              "lockgrid": lambda: metrics_lock_grid(sweep, bins, meshes, args, common),
              "diagnostics": lambda: metrics_diagnostics(sweep, bins, meshes, args,
                                                         common, sweep.results_dir)}
    wanted = [x.strip() for x in args.metrics_phases.split(",") if x.strip()]
    unknown = [x for x in wanted if x not in phases]
    if unknown:
        sys.exit("Unknown --metrics-phases: %s. Known: %s"
                 % (", ".join(unknown), ", ".join(phases)))

    print("=" * 70)
    print("METRICS PROFILE: %s" % ", ".join(wanted))
    print("Roughly 45 min at 24 threads. The machine must be idle for all of it.")
    print("=" * 70)
    for name in wanted:
        phases[name]()

    print("\n" + "=" * 70)
    print("METRICS RUN DONE. Package it with:\n")
    print("    python3 %s --root %s --quality-pass \\\n"
          "        --results-dir %s --mesh-out-dir %s"
          % (Path(__file__).name, root, sweep.results_dir, mesh_out_dir))
    print("=" * 70)
    return None


def package(results_dir, root, mesh_out_dir=None):
    """Bundle exactly what needs to come back, so nobody has to guess."""
    import tarfile
    out = root / ("results_%s_%s.tar.gz" % (platform.node(),
                                            time.strftime("%Y%m%d-%H%M")))
    wanted = ["results.csv", "env.json", "calibration.json",
              "toolchain_lock.json", "json", "quality",
              # the metrics profile's artefacts: the raw perf counter files and
              # the instrumented runs' stdout, which is the whole point of them
              "perf", "diagnostics"]
    try:
        with tarfile.open(out, "w:gz") as tf:
            for name in wanted:
                p = results_dir / name
                if p.exists():
                    tf.add(p, arcname="results/" + name)
    except Exception as e:
        print("Could not write the bundle (%s). Send these by hand: %s"
              % (e, ", ".join(wanted)), file=sys.stderr)
        return None

    size = out.stat().st_size
    human = ("%.0f KB" % (size / 1e3)) if size < 1e6 else ("%.1f MB" % (size / 1e6))
    n = sum(1 for _ in tarfile.open(out).getnames())
    print("\n" + "=" * 70)
    print("DONE. Send back this one file:\n\n    %s   (%s, %d entries)"
          % (out, human, n))
    print("\nThe output meshes stay here as the archive -- do not delete them, so\n"
          "any further metric can be computed without re-running anything:\n"
          "    %s" % (mesh_out_dir or (root / "out_meshes")))
    print("=" * 70)
    return out


# --------------------------------------------------------------------------
# quality pass
# --------------------------------------------------------------------------

def _quality_priority(path):
    """Order the quality pass so the headline arms are measured first.

    If the pass is cut short, what survives should be the meshes the report
    leads with -- the reference arms and the full-width parallel run -- rather
    than whichever intermediate thread count happened to sort first.
    """
    name = path.stem
    if "_main_" in name:
        return (0, name)
    if "_seq_" in name:
        return (1, name)
    if "_par_t1." in name + ".":
        return (2, name)
    # remaining parallel points: widest first
    try:
        t = int(name.rsplit("_t", 1)[1])
    except Exception:
        t = 0
    return (3, -t, name)


def quality_pass(root, results_dir, mesh_out_dir, bins, force, deadline=None):
    """Derive every quality metric from the saved meshes.

    Run AFTER the sweep, on the benchmark machine, with nothing else going on.
    The colleague then returns only these JSONs (a few MB) and keeps the .mesh
    archive in place, so a metric we think of later costs a re-run of this step
    and not a re-run of the benchmark.
    """
    exe = bins["mesh_quality_report"]
    qdir = results_dir / "quality"
    qdir.mkdir(parents=True, exist_ok=True)

    meshes = sorted(mesh_out_dir.glob("*.mesh"), key=_quality_priority)
    if not meshes:
        print("No output meshes in %s -- nothing to measure." % mesh_out_dir,
              file=sys.stderr)
        return
    print("%d output meshes to measure" % len(meshes))

    ok = failed = skipped = 0
    for n, mp in enumerate(meshes, 1):
        if deadline and time.time() > deadline:
            print("\nQuality budget spent after %d of %d meshes. The meshes are all\n"
                  "still here, so re-running --quality-pass finishes the rest."
                  % (n - 1, len(meshes)))
            break
        qp = qdir / (mp.stem + ".json")
        if qp.exists() and not force:
            try:
                json.loads(qp.read_text())
                skipped += 1
                continue
            except Exception:
                qp.unlink(missing_ok=True)
        r = subprocess.run([str(exe), str(mp), str(qp)],
                           text=True, capture_output=True)
        if r.returncode == 0 and qp.exists():
            ok += 1
        else:
            failed += 1
            print("  FAILED %s: %s" % (mp.name, (r.stderr or "")[-200:]))
        if n % 10 == 0 or n == len(meshes):
            print("  [%d/%d] ok=%d skipped=%d failed=%d" %
                  (n, len(meshes), ok, skipped, failed), flush=True)

    print("\nQuality JSONs in %s" % qdir)
    package(results_dir, root, mesh_out_dir)


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(REPO_ROOT / "work"))
    ap.add_argument("--profile",
                    choices=["full", "calibrate", "overnight", "metrics"],
                    help="'full' does calibrate + sweep + quality + packaging in "
                         "one unattended command (use this if you get one run at "
                         "the machine); the others are the same phases separately. "
                         "'metrics' is the short scaling-only run described in "
                         "docs/METRICS_REQUEST.md -- about 45 min, and it does not "
                         "calibrate, since its configurations are fixed")
    ap.add_argument("--metrics-phases",
                    default="ladder,pinning,lockgrid,diagnostics",
                    help="which parts of --profile metrics to run, in order")
    ap.add_argument("--calib-budget", type=float, default=2400.0,
                    help="seconds of the full-profile budget given to calibration")
    ap.add_argument("--quality-pass", action="store_true",
                    help="derive quality from the saved meshes and exit")
    ap.add_argument("--force-quality", action="store_true")
    ap.add_argument("--results-dir")
    ap.add_argument("--mesh-out-dir")
    ap.add_argument("--mesh-dir")
    ap.add_argument("--configs", help="overnight: a configs JSON to use")
    ap.add_argument("--threads-max", type=int, default=24)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--max-configs", type=int, default=12)
    ap.add_argument("--calib-meshes", type=int, default=4,
                    help="how many inputs the ladder probes (split across pipelines)")
    ap.add_argument("--pipelines", default=None,
                    help="restrict to these input pipelines, e.g. 'cdt' or 'mesh3'; "
                         "default is everything the manifest has")
    ap.add_argument("--budget", type=float, default=None,
                    help="hard wall-clock budget in seconds "
                         "(default 1800 for calibrate, none for overnight)")
    ap.add_argument("--run-timeout", type=int, default=None,
                    help="per-run timeout in seconds "
                         "(default 240 for calibrate, 7200 for overnight)")
    ap.add_argument("--allow-toolchain-change", action="store_true",
                    help="continue even though the binaries differ from the ones "
                         "this results directory was started with (the results "
                         "set then mixes code versions)")
    ap.add_argument("--settle", type=float, default=2.0,
                    help="seconds to idle before each run")
    args = ap.parse_args()

    if not args.profile and not args.quality_pass:
        ap.error("pass --profile or --quality-pass")

    root = Path(args.root).resolve()
    tc_path = root / "toolchain.json"
    if not tc_path.exists():
        sys.exit("No %s -- run scripts/setup.py first." % tc_path)
    tc = json.loads(tc_path.read_text())
    bins = dict((n, Path(v["path"])) for n, v in tc["binaries"].items())
    # The instrumented binaries are deliberately NOT in tc["binaries"]: the
    # toolchain lock compares that dict, so adding them there would make a
    # metrics run refuse to resume a sweep, and a sweep refuse to resume after a
    # metrics run. They change no timed binary, so they do not belong in it.
    bins.update((n, Path(v["path"])) for n, v in tc.get("diag_binaries", {}).items())

    results_dir = Path(args.results_dir).resolve() if args.results_dir else root / "results"
    mesh_out_dir = Path(args.mesh_out_dir).resolve() if args.mesh_out_dir else root / "out_meshes"
    mesh_dir = Path(args.mesh_dir).resolve() if args.mesh_dir else root / "meshes"
    for d in (results_dir / "json", results_dir / "logs", results_dir / "perf",
              mesh_out_dir):
        d.mkdir(parents=True, exist_ok=True)

    if args.quality_pass:
        quality_pass(root, results_dir, mesh_out_dir, bins, args.force_quality)
        return

    if args.budget is None:
        args.budget = {"calibrate": 1800.0, "full": 12 * 3600.0,
                       # Generous against the ~45 min the work takes, so a slow
                       # machine finishes rather than being cut off mid-ladder.
                       "metrics": 3 * 3600.0}.get(args.profile, 0.0)
    if args.run_timeout is None:
        # 240s was too tight: on a 24-core run the very first ladder probe hit
        # it, and every single-threaded anchor did too, so calibration learned
        # nothing about the serial cost. A single-threaded run of a config that
        # takes ~30s on 24 threads needs room for ~15 minutes.
        args.run_timeout = 1200 if args.profile in ("calibrate", "metrics") else 14400

    check_toolchain_lock(root, results_dir, args.allow_toolchain_change)
    env = write_env(root, results_dir)
    nproc = env["nproc"] or 0
    if args.threads_max > nproc:
        print("Warning: --threads-max %d exceeds nproc %d; oversubscribing."
              % (args.threads_max, nproc), file=sys.stderr)
    # Only warn when it actually means something. Under intel_pstate,
    # 'powersave' is the stock default and still reaches full turbo; warning
    # about it sends people chasing a non-problem.
    gov, drv = env.get("cpu_governor", ""), env.get("cpufreq_driver", "")
    if gov and gov != "performance" and drv != "intel_pstate":
        print("Warning: CPU governor is '%s' with driver '%s'. Timings will be\n"
              "         noisier and low thread counts will look worse than they are."
              % (gov, drv or "unknown"), file=sys.stderr)
    if env.get("no_turbo") == "1":
        print("Warning: turbo is DISABLED (intel_pstate/no_turbo=1). Every thread\n"
              "         count is capped at base clock.", file=sys.stderr)

    pipes = [p.strip() for p in args.pipelines.split(",") if p.strip()] \
        if args.pipelines else None
    meshes = load_manifest(root, mesh_dir, pipes)
    counts = {}
    for m in meshes:
        counts[m["pipeline"]] = counts.get(m["pipeline"], 0) + 1
    print("%d prepared inputs (%s), largest %s with %d cells"
          % (len(meshes), ", ".join("%s: %d" % kv for kv in sorted(counts.items())),
             meshes[0]["key"], meshes[0]["cells"]))
    if len(counts) == 1 and not args.pipelines:
        print("Note: only the '%s' input pipeline is present. CDT and Mesh_3 inputs\n"
              "      have disagreed on the sign of the sequential comparison before,\n"
              "      so a run on one of them alone answers a narrower question than\n"
              "      it looks like. Re-run prepare_meshes.py to build both."
              % list(counts)[0], file=sys.stderr)

    sweep = Sweep(root, results_dir, args.budget, args)
    sweep.mesh_out_dir = mesh_out_dir
    started = time.time()
    try:
        if args.profile == "full":
            profile_full(sweep, bins, meshes, args, root, mesh_out_dir)
        elif args.profile == "metrics":
            profile_metrics(sweep, bins, meshes, args, root, mesh_out_dir)
        elif args.profile == "calibrate":
            profile_calibrate(sweep, bins, meshes, args)
        else:
            profile_overnight(sweep, bins, meshes, args)
    finally:
        sweep.close()

    print("\nElapsed %.0f s (%.2f h). Rows in %s"
          % (time.time() - started, (time.time() - started) / 3600, sweep.csv_path))
    if args.profile != "full":
        print("Then:  python3 scripts/run_bench.py --root %s --quality-pass" % root)


if __name__ == "__main__":
    main()
