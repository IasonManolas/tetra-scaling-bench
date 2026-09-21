#!/usr/bin/env python3
"""Stage 2: the timed sweep.

Two profiles:

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
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

ITERS = 10                 # remeshing iterations, as in the existing campaign
SMOOTH_CONSTRAINED = 1
TARGET_SECONDS = 30.0      # the 24-thread run must last at least this long
LADDER = [0.5, 0.4, 0.3, 0.25]   # cells grow roughly as factor^-3
THREAD_LIST = [1, 2, 4, 8, 12, 16, 24]

CSV_HEADER = ["mesh", "factor", "arm", "threads", "rep", "wall_s", "remesh_s",
              "peak_rss_kb", "rc", "timed_out", "out_mesh", "json"]


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
    wrapped = cmd
    if shutil.which("/usr/bin/time"):
        wrapped = ["/usr/bin/time", "-f", "%e %M", "-o", str(rss_file)] + cmd

    t0 = time.time()
    timed_out = False
    try:
        with open(log_path, "w") as log:
            rc = subprocess.run(wrapped, stdout=log, stderr=subprocess.STDOUT,
                                timeout=timeout).returncode
    except subprocess.TimeoutExpired:
        rc, timed_out = -9, True
    wall = time.time() - t0

    peak_rss = ""
    if rss_file.exists():
        try:
            peak_rss = rss_file.read_text().split()[-1]
        except Exception:
            pass

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
            "out_mesh": out_mesh, "json": str(json_path)}


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
        return stem in self.done

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

def load_manifest(root, mesh_dir):
    mp = mesh_dir / "manifest.json"
    if not mp.exists():
        sys.exit("No %s -- run scripts/prepare_meshes.py first." % mp)
    meshes = json.loads(mp.read_text())["meshes"]
    ok = [m for m in meshes.values() if m.get("cdt_cells") and Path(m["path"]).exists()]
    return sorted(ok, key=lambda m: m["cdt_cells"], reverse=True)


def profile_calibrate(sweep, bins, meshes, args):
    """Ordered by value, so a truncated run is still useful."""
    rd, md = sweep.results_dir, sweep.mesh_out_dir
    max_t = args.threads_max
    common = dict(results_dir=rd, mesh_out_dir=md, settle=args.settle, done=sweep.done)

    print("\n=== 1/4  size ladder: find (mesh, factor) with %d-thread wall >= %.0fs ==="
          % (max_t, TARGET_SECONDS))
    cleared = []          # (cells, mesh, factor, wall)
    for m in meshes[:args.calib_meshes]:
        for factor in LADDER:
            if not sweep.affordable(args.run_timeout * 0.5, m["path"], factor,
                                    "par", max_t, 1):
                print("  budget spent, stopping the ladder")
                break
            row = run_cell(bins, m["path"], factor, "par", max_t, 1,
                           write_mesh=True, timeout=args.run_timeout, **common)
            sweep.emit(row)
            print("  %-10s f=%-4s  %6.1fs%s%s" %
                  (m["id"], factor, row["wall_s"],
                   "  TIMEOUT" if row["timed_out"] else "",
                   "  (cached)" if row.get("cached") else ""))
            if row["timed_out"]:
                # We have bracketed it: this factor is already too slow, so
                # descending further only wastes budget.
                break
            if row["wall_s"] >= TARGET_SECONDS:
                cleared.append((m["cdt_cells"], m["id"], m["path"], factor, row["wall_s"]))
                break

    if not cleared:
        print("\nNo config reached %.0fs at %d threads." % (TARGET_SECONDS, max_t))
        print("That is the finding: these inputs are too small for this machine.")
        print("Next step is bigger surfaces (relax the 10k-vertex Thingi filter),")
        print("not a finer edge factor -- the ladder already went to %s." % LADDER[-1])
        return None

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

    print("\n=== 3/4  scaling curve on the pilot config ===")
    curve = {max_t: pwall}
    for t in [x for x in THREAD_LIST if 1 < x < max_t]:
        est = pwall * max_t / t
        if not sweep.affordable(est * 1.3, ppath, pfactor, "par", t, 1):
            print("  budget spent, skipping %d threads and below" % t)
            break
        row = run_cell(bins, ppath, pfactor, "par", t, 1,
                       write_mesh=True, timeout=args.run_timeout, **common)
        sweep.emit(row)
        curve[t] = row["wall_s"]
        print("  %2d threads: %6.1fs  (speedup vs %d-thread: %.2fx)%s" %
              (t, row["wall_s"], max_t, pwall / max(row["wall_s"], 1e-9),
               "  (cached)" if row.get("cached") else ""))

    print("\n=== 4/4  sequential anchors, on the SMALLEST ladder config ===")
    # Deliberately not the pilot config: a 30s-at-24-threads workload can take
    # ten minutes single-threaded and would eat the whole budget on its own.
    # The overnight ETA extrapolates from this small config's ratio instead.
    small = meshes[min(args.calib_meshes, len(meshes)) - 1]
    anchors = {}
    for arm in ("par", "seq", "main"):
        if not sweep.affordable(120, small["path"], LADDER[0], arm, 1, 1):
            print("  budget spent, skipping %s" % arm)
            break
        row = run_cell(bins, small["path"], LADDER[0], arm, 1, 1,
                       write_mesh=True, timeout=args.run_timeout, **common)
        sweep.emit(row)
        anchors[arm] = row["wall_s"]
        print("  %-5s @1 on %s f=%s: %6.1fs%s" %
              (arm, small["id"], LADDER[0], row["wall_s"],
               "  (cached)" if row.get("cached") else ""))

    configs = [{"mesh": c[1], "path": c[2], "factor": c[3],
                "cells": c[0], "wall_at_max_t": c[4]} for c in cleared]

    # How much slower one thread is than max_t on the SMALL config. The
    # overnight ETA needs this: the seq and par@1 arms dominate that run, and
    # assuming they scale down perfectly from the 24-thread number would
    # underestimate it badly.
    serial_ratio = None
    if "par" in anchors and curve.get(max_t):
        small_par_max = None
        srow = run_cell(bins, small["path"], LADDER[0], "par", max_t, 1,
                        write_mesh=False, timeout=args.run_timeout, **common) \
               if sweep.affordable(60, small["path"], LADDER[0], "par", max_t, 1) else None
        if srow:
            sweep.emit(srow)
            small_par_max = srow["wall_s"]
        if small_par_max:
            serial_ratio = anchors["par"] / max(small_par_max, 1e-9)
            print("\n  par@1 / par@%d on the small config: %.2fx" % (max_t, serial_ratio))

    out = sweep.results_dir / "calibration.json"
    out.write_text(json.dumps({
        "target_seconds": TARGET_SECONDS,
        "threads_max": max_t,
        "ladder": LADDER,
        "configs": configs,
        "pilot": {"mesh": pid, "factor": pfactor, "wall_at_max_t": pwall},
        "cv": cv,
        "scaling_curve_seconds": curve,
        "anchors_small_config": {"mesh": small["id"], "factor": LADDER[0],
                                 "walls": anchors},
        "serial_ratio": serial_ratio,
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

    for c in configs:
        path, factor = c["path"], c["factor"]
        print("\n--- %s  f=%s  (%s cells in)" % (c["mesh"], factor, c.get("cells")))
        for arm, tlist in (("par", threads), ("seq", [1]), ("main", [1])):
            for t in tlist:
                for rep in range(1, args.reps + 1):
                    if not sweep.affordable(1, path, factor, arm, t, rep):
                        print("Budget spent; stopping. Re-run to resume.")
                        return None
                    row = run_cell(bins, path, factor, arm, t, rep,
                                   write_mesh=(rep == 1),
                                   timeout=args.run_timeout, **common)
                    sweep.emit(row)
                    print("  %-5s t=%-3d r=%d  %8.1fs%s" %
                          (arm, t, rep, row["wall_s"],
                           "  (cached)" if row.get("cached") else ""))
    return None


# --------------------------------------------------------------------------
# quality pass
# --------------------------------------------------------------------------

def quality_pass(root, results_dir, mesh_out_dir, bins, force):
    """Derive every quality metric from the saved meshes.

    Run AFTER the sweep, on the benchmark machine, with nothing else going on.
    The colleague then returns only these JSONs (a few MB) and keeps the .mesh
    archive in place, so a metric we think of later costs a re-run of this step
    and not a re-run of the benchmark.
    """
    exe = bins["mesh_quality_report"]
    qdir = results_dir / "quality"
    qdir.mkdir(parents=True, exist_ok=True)

    meshes = sorted(mesh_out_dir.glob("*.mesh"))
    if not meshes:
        sys.exit("No output meshes in %s" % mesh_out_dir)
    print("%d output meshes to measure" % len(meshes))

    ok = failed = skipped = 0
    for n, mp in enumerate(meshes, 1):
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
    print("Send back: results.csv, json/, quality/, env.json, calibration.json")
    print("Keep in place: %s (the mesh archive)" % mesh_out_dir)


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(REPO_ROOT / "work"))
    ap.add_argument("--profile", choices=["calibrate", "overnight"])
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
    ap.add_argument("--calib-meshes", type=int, default=3)
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
        args.budget = 1800.0 if args.profile == "calibrate" else 0.0
    if args.run_timeout is None:
        args.run_timeout = 240 if args.profile == "calibrate" else 7200

    env = write_env(root, results_dir)
    nproc = env["nproc"] or 0
    if args.threads_max > nproc:
        print("Warning: --threads-max %d exceeds nproc %d; oversubscribing."
              % (args.threads_max, nproc), file=sys.stderr)
    if env["cpu_governor"] and env["cpu_governor"] != "performance":
        print("Warning: CPU governor is '%s', not 'performance'. Timings will be\n"
              "         noisier and low thread counts will look worse than they are."
              % env["cpu_governor"], file=sys.stderr)

    meshes = load_manifest(root, mesh_dir)
    print("%d prepared meshes, largest %s with %d cells"
          % (len(meshes), meshes[0]["id"], meshes[0]["cdt_cells"]))

    sweep = Sweep(root, results_dir, args.budget, args)
    sweep.mesh_out_dir = mesh_out_dir
    try:
        if args.profile == "calibrate":
            profile_calibrate(sweep, bins, meshes, args)
        else:
            profile_overnight(sweep, bins, meshes, args)
    finally:
        sweep.close()

    print("\nElapsed %.0f s. Rows in %s" % (time.time() - sweep.t0, sweep.csv_path))
    print("Then:  python3 scripts/run_bench.py --root %s --quality-pass" % root)


if __name__ == "__main__":
    main()
