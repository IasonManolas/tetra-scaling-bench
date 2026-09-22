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

# started_at is not decoration. Cells are measured in different epochs -- a
# calibration cell reused by the sweep can be an hour older than its own
# repeat -- and machines drift (25% between waves on the test laptop). Without
# a timestamp there is no way to tell a real effect from thermal drift after
# the fact, and no way to re-measure only the suspect rows.
CSV_HEADER = ["mesh", "factor", "arm", "threads", "rep", "wall_s", "remesh_s",
              "peak_rss_kb", "rc", "timed_out", "started_at", "out_mesh", "json"]


# --------------------------------------------------------------------------
# environment capture
# --------------------------------------------------------------------------

def capture(cmd):
    try:
        return subprocess.run(cmd, text=True, capture_output=True,
                              timeout=30).stdout.strip()
    except Exception:
        return ""


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
        "lscpu": capture(["lscpu"]),
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


def run_cell(bins, mesh_path, factor, arm, threads, rep, results_dir,
             mesh_out_dir, timeout, write_mesh, settle, done=None):
    """Execute one (mesh, factor, arm, threads, rep) cell.

    Always returns a dict row. A cell that was already completed by an earlier
    invocation is returned from `done` with cached=1 and is NOT re-run and NOT
    re-appended to the CSV. Returning the cached row rather than None matters:
    calibrate reads wall times back to choose its pilot config, so a resumed
    run must still reach the same conclusions without redoing the work.
    """
    mesh_id = Path(mesh_path).stem
    stem = "%s_f%s_%s_t%d_r%d" % (mesh_id, factor, arm, threads, rep)
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
                    "cached": 1, "recovered": 1}
        except Exception:
            json_path.unlink(missing_ok=True)

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
    rss_file.unlink(missing_ok=True)
    wrapped = cmd
    if shutil.which("/usr/bin/time"):
        wrapped = ["/usr/bin/time", "-f", "%e %M", "-o", str(rss_file)] + cmd

    # The whole run goes in its own process group, and a timeout kills the
    # GROUP. subprocess.run(timeout=...) only kills its direct child, which
    # here is /usr/bin/time -- the bench_remesh underneath it survives, keeps
    # every core busy, and corrupts whatever is measured next. That is not
    # hypothetical: a 24-core run recorded three cells as 240s timeouts whose
    # JSONs said "success" at 259s, and the two measurements that followed came
    # out 8% slow because an orphan was still running beside them.
    t0 = time.time()
    started_at = time.strftime("%Y-%m-%dT%H:%M:%S")
    timed_out = False
    with open(log_path, "w") as log:
        proc = subprocess.Popen(wrapped, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            _kill_group(proc)
            rc = proc.returncode if proc.returncode is not None else -9
    wall = time.time() - t0

    peak_rss = ""
    if rss_file.exists():
        try:
            peak_rss = rss_file.read_text().split()[-1]
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

    return {"mesh": mesh_id, "factor": factor, "arm": arm, "threads": threads,
            "rep": rep, "wall_s": round(wall, 3), "remesh_s": remesh_s,
            "peak_rss_kb": peak_rss, "rc": rc, "timed_out": int(timed_out),
            "started_at": started_at,
            "out_mesh": "" if timed_out else out_mesh, "json": str(json_path)}


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
        self.csv_fh = open(self.csv_path, "a", newline="")
        self.csv = csv.DictWriter(self.csv_fh, fieldnames=CSV_HEADER,
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
                    stem = "%s_f%s_%s_t%d_r%d" % (
                        r["mesh"], r["factor"], r["arm"],
                        int(r["threads"]), int(r["rep"]))
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

    def is_cached(self, mesh_path, factor, arm, threads, rep):
        """Whether this cell is already recorded. The budget guards consult this
        first: a cached cell costs no wall time, so a resumed run must not be
        stopped by the budget before it has read its earlier results back."""
        stem = "%s_f%s_%s_t%d_r%d" % (Path(mesh_path).stem, factor, arm, threads, rep)
        if stem not in self.done:
            return False
        # A CSV row is not enough. A timed-out cell has a row but its JSON was
        # deleted, so run_cell will really execute it -- and charging the budget
        # nothing for a run that is about to take the full timeout is how a
        # 30-minute calibration turns into an hour.
        return (self.results_dir / "json" / (stem + ".json")).exists()

    def affordable(self, seconds, mesh_path, factor, arm, threads, rep):
        if self.is_cached(mesh_path, factor, arm, threads, rep):
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


def package(results_dir, root, mesh_out_dir=None):
    """Bundle exactly what needs to come back, so nobody has to guess."""
    import tarfile
    out = root / ("results_%s_%s.tar.gz" % (platform.node(),
                                            time.strftime("%Y%m%d-%H%M")))
    wanted = ["results.csv", "env.json", "calibration.json", "json", "quality"]
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
    ap.add_argument("--profile", choices=["full", "calibrate", "overnight"],
                    help="'full' does calibrate + sweep + quality + packaging in "
                         "one unattended command (use this if you get one run at "
                         "the machine); the others are the same phases separately")
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
    ap.add_argument("--settle", type=float, default=2.0,
                    help="seconds to idle before each run")
    args = ap.parse_args()

    if not args.profile and not args.quality_pass:
        ap.error("pass --profile or --quality-pass")

    root = Path(args.root).resolve()
    tc_path = root / "toolchain.json"
    if not tc_path.exists():
        sys.exit("No %s -- run scripts/setup.py first." % tc_path)
    bins = dict((n, Path(v["path"]))
                for n, v in json.loads(tc_path.read_text())["binaries"].items())

    results_dir = Path(args.results_dir).resolve() if args.results_dir else root / "results"
    mesh_out_dir = Path(args.mesh_out_dir).resolve() if args.mesh_out_dir else root / "out_meshes"
    mesh_dir = Path(args.mesh_dir).resolve() if args.mesh_dir else root / "meshes"
    for d in (results_dir / "json", results_dir / "logs", mesh_out_dir):
        d.mkdir(parents=True, exist_ok=True)

    if args.quality_pass:
        quality_pass(root, results_dir, mesh_out_dir, bins, args.force_quality)
        return

    if args.budget is None:
        args.budget = {"calibrate": 1800.0, "full": 12 * 3600.0}.get(args.profile, 0.0)
    if args.run_timeout is None:
        # 240s was too tight: on a 24-core run the very first ladder probe hit
        # it, and every single-threaded anchor did too, so calibration learned
        # nothing about the serial cost. A single-threaded run of a config that
        # takes ~30s on 24 threads needs room for ~15 minutes.
        args.run_timeout = 1200 if args.profile == "calibrate" else 14400

    env = write_env(root, results_dir)
    nproc = env["nproc"] or 0
    if args.threads_max > nproc:
        print("Warning: --threads-max %d exceeds nproc %d; oversubscribing."
              % (args.threads_max, nproc), file=sys.stderr)
    if env["cpu_governor"] and env["cpu_governor"] != "performance":
        print("Warning: CPU governor is '%s', not 'performance'. Timings will be\n"
              "         noisier and low thread counts will look worse than they are."
              % env["cpu_governor"], file=sys.stderr)

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
