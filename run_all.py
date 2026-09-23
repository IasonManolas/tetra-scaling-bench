#!/usr/bin/env python3
"""Everything, in one command.

    python3 run_all.py --surfaces-dir /path/to/thingi10k

Runs: build -> prepare inputs -> measure (~12 h) -> derive quality -> package.
Then send back the one tarball it names at the end.

Every stage is idempotent, so if this is interrupted at any point, re-running
the same command picks up where it stopped: finished builds, prepared meshes and
completed measurement cells are all skipped.

The only stage that needs the machine otherwise idle is the measurement. It is
also the only one with a deadline, and it will not overrun it.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE / "scripts"


def hrs(seconds):
    return "%dh%02dm" % (int(seconds // 3600), int((seconds % 3600) // 60))


def stage(name, cmd, log_dir, dry):
    print("\n" + "=" * 72)
    print("STAGE: %s" % name)
    print("  " + " ".join(str(c) for c in cmd))
    print("=" * 72, flush=True)
    if dry:
        return 0.0

    t0 = time.time()
    log_path = log_dir / ("%s.log" % name.replace(" ", "_"))
    # Tee: the operator sees progress live, and the full log survives for us.
    with open(log_path, "w") as log:
        proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
        rc = proc.wait()

    took = time.time() - t0
    if rc != 0:
        sys.exit("\nStage '%s' failed (exit %d). Full log: %s\n"
                 "Fix the cause and re-run this same command -- finished stages "
                 "are skipped." % (name, rc, log_path))
    print("\n[%s took %s]" % (name, hrs(took)), flush=True)
    return took


def try_set_governor(want):
    """Best effort, never fatal.

    The CPU frequency governor decides whether the kernel scales clock speed
    with load. The usual default, 'powersave', does -- which means the same work
    takes different time depending on when it runs. Measured on a laptop here:
    25% drift between measurement waves, enough to make a 2-thread run look
    faster than a 4-thread one. 'performance' pins the clocks so runs are
    comparable. It needs root, so if it does not work we say so and continue --
    the numbers are then noisier, not wrong.
    """
    cur = []
    for p in Path("/sys/devices/system/cpu").glob("cpu*/cpufreq/scaling_governor"):
        try:
            cur.append(p.read_text().strip())
        except Exception:
            pass
    if not cur:
        print("[governor] no cpufreq on this machine; nothing to set.")
        return
    if all(g == want for g in cur):
        print("[governor] already '%s'." % want)
        return

    print("[governor] currently %s; trying to set '%s'."
          % (sorted(set(cur)), want))
    for cmd in (["cpupower", "frequency-set", "-g", want],
                ["sudo", "-n", "cpupower", "frequency-set", "-g", want]):
        if shutil.which(cmd[0]) is None:
            continue
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if r.returncode == 0:
                print("[governor] set to '%s'." % want)
                return
        except Exception:
            pass

    print("[governor] could NOT set it (needs root). Continuing anyway.\n"
          "           Timings will drift more; the report's CV column and the\n"
          "           set totals are the numbers to trust. To fix it, run:\n"
          "             sudo cpupower frequency-set -g performance", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # --off-dir is the old spelling, kept working: Thingi10K ships .stl and
    # that is the default, so the original name misdescribes the argument.
    ap.add_argument("--surfaces-dir", "--off-dir", dest="surfaces_dir",
                    required=True, metavar="DIR",
                    help="directory of Thingi10K surfaces (.stl, .off or .ply)")
    ap.add_argument("--root", default=str(HERE / "work"))
    ap.add_argument("--budget", type=float, default=12 * 3600,
                    help="seconds for the MEASUREMENT stage (default 12 h)")
    ap.add_argument("--threads-max", type=int, default=(os.cpu_count() or 24))
    ap.add_argument("--jobs", type=int, default=(os.cpu_count() or 8))
    ap.add_argument("--tbb-dir", help="directory containing TBBConfig.cmake")
    # All surfaces get a CDT (cheap, and the only honest measure of how big a
    # workload each one is). Mesh_3, the expensive pipeline, goes to the biggest
    # CDTs only. Capping by surface file size is what made a 24-core run top out
    # at 26s when a surface outside the cap would have given 67s.
    ap.add_argument("--limit", type=int, default=0,
                    help="consider only the N largest surface FILES (0 = all)")
    ap.add_argument("--mesh3-limit", type=int, default=20,
                    help="Mesh_3 partners for the N largest CDTs (0 = all)")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--governor", default="performance",
                    help="CPU governor to request; 'none' to skip")
    ap.add_argument("--smoke", action="store_true",
                    help="short validation run (~30-45 min) before committing the "
                         "machine to the real one: a few surfaces, a 15-minute "
                         "measurement, written to separate smoke_* directories so "
                         "it cannot contaminate the real results. The build and the "
                         "prepared meshes it produces are reused by the real run.")
    ap.add_argument("--metrics", action="store_true",
                    help="the scaling-only run (~45 min) instead of the 12-hour "
                         "sweep: two thread ladders with cycles and instructions, "
                         "two instrumented runs, a per-thread profile, a spatial-sort interval sweep "
                         "and a lock-grid sweep. Written to separate metrics_* "
                         "directories, like --smoke, so it cannot contaminate the "
                         "real results. It also builds the two diagnostic binaries.")
    ap.add_argument("--skip-quality", action="store_true",
                    help="stop after the measurement (quality can be run later)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the stages and exit")
    args = ap.parse_args()

    surfaces_dir = Path(args.surfaces_dir).expanduser().resolve()
    root = Path(args.root).expanduser().resolve()

    # The smoke run shares the build and the prepared meshes with the real run
    # -- those are expensive and identical -- but keeps its own results and
    # output meshes. Sharing those would mean the real run resumes on top of
    # smoke-run cells, mixing a throwaway measurement into the real dataset.
    if args.smoke and args.metrics:
        ap.error("--smoke and --metrics are different runs; pass one.")

    # The metrics run is the same build and the same prepared meshes as the real
    # one -- it measures the very configurations the real sweep measures -- but
    # it answers only the scaling questions, so it writes its own results
    # directory for the same reason --smoke does: resuming the real sweep on top
    # of these cells would mix two different measurements into one dataset.
    if args.metrics:
        if args.budget == 12 * 3600:
            args.budget = 90 * 60     # a cap; the work is ~45 min
        args.reps = 3
        results_dir = root / "metrics_results"
        mesh_out_dir = root / "metrics_out_meshes"
        print("METRICS RUN: the scaling measurements only, roughly 45 min at 24\n"
              "threads, results in %s\n" % results_dir.name)
    elif args.smoke:
        if args.limit == 0:             # untouched default
            args.limit = 6
        if args.mesh3_limit == 20:
            args.mesh3_limit = 3
        if args.budget == 12 * 3600:
            # Room for a real calibration (a single ladder rung on a fast
            # machine is ~70s) plus a sweep that actually measures something.
            args.budget = 1800
        args.reps = 1
        results_dir = root / "smoke_results"
        mesh_out_dir = root / "smoke_out_meshes"
        print("SMOKE RUN: %d surfaces, %s measurement, results in %s"
              % (args.limit, hrs(args.budget), results_dir.name))
        print("The real run reuses this build and these meshes, and writes "
              "elsewhere.\n")
    else:
        results_dir = root / "results"
        mesh_out_dir = root / "out_meshes"

    # Fail before the long stages, not four hours in.
    problems = []
    if not surfaces_dir.is_dir():
        problems.append("--surfaces-dir %s is not a directory" % surfaces_dir)
    else:
        n = sum(1 for p in surfaces_dir.iterdir()
                if p.suffix.lower() in (".stl", ".off", ".ply"))
        if n == 0:
            problems.append("no .stl/.off/.ply files in %s" % surfaces_dir)
        else:
            print("Found %d surface files in %s" % (n, surfaces_dir))
    for tool in ("cmake", "git"):
        if shutil.which(tool) is None:
            problems.append("'%s' is not on PATH" % tool)
    # Scale the disk estimate with how much is actually being built. A flat
    # figure blocks small runs for no reason: ~3 GB goes on the two CGAL clones
    # and their builds no matter what, and the rest tracks the number of
    # surfaces (roughly 0.25 GB of input mesh and 0.35 GB of output mesh each,
    # measured at ~40 bytes per cell on inputs of this size).
    n_surfaces = args.limit if args.limit else 100
    need_gb = 3.0 + n_surfaces * 0.6
    try:
        probe = root
        while not probe.exists() and probe.parent != probe:
            probe = probe.parent
        free_gb = shutil.disk_usage(probe).free / 1e9
        print("Free disk at %s: %.0f GB (this run needs roughly %.0f GB for %d "
              "surfaces)" % (probe, free_gb, need_gb, n_surfaces))
        if free_gb < need_gb:
            problems.append(
                "only %.0f GB free but roughly %.0f GB is needed for %d surfaces.\n"
                "        Free some space, point --root at a bigger disk, or lower "
                "--limit." % (free_gb, need_gb, n_surfaces))
    except Exception:
        pass

    if problems:
        # A dry run is for checking the plan, so report the problems and still
        # show it rather than exiting on the first one.
        print("\n%s\n" % ("Problems (dry run, continuing):" if args.dry_run
                          else "Cannot start:"), file=sys.stderr)
        for p in problems:
            print("  - %s" % p, file=sys.stderr)
        if not args.dry_run:
            sys.exit(1)

    root.mkdir(parents=True, exist_ok=True)
    log_dir = root / "logs_run_all"
    log_dir.mkdir(parents=True, exist_ok=True)

    n_prep = args.limit if args.limit else 100
    print("\nPlan:")
    print("  1 build          ~15 min")
    print("  2 prepare inputs ~%s   (%d largest surfaces, CDT + Mesh_3)"
          % ("10 min" if n_prep <= 5 else "1-3 h", n_prep))
    print("  3 measure        %s      MACHINE MUST BE IDLE" % hrs(args.budget))
    print("  4 quality        ~%s     (machine may be busy)"
          % ("5 min" if n_prep <= 5 else "1 h"))
    print("  5 package")
    print("\nInterrupting is safe: re-run this exact command to resume.")

    py = sys.executable or "python3"
    total0 = time.time()

    setup = [py, SCRIPTS / "setup.py", "--root", root, "--jobs", args.jobs]
    if args.tbb_dir:
        setup += ["--tbb-dir", args.tbb_dir]
    if args.metrics:
        setup += ["--diagnostics"]
    stage("1-build", setup, log_dir, args.dry_run)

    prep = [py, SCRIPTS / "prepare_meshes.py", "--root", root,
            "--surfaces-dir", surfaces_dir, "--jobs", args.jobs,
            "--limit", args.limit, "--mesh3-limit", args.mesh3_limit]
    if args.smoke:
        # Mesh_3 sizing iterates; on a big surface seven rounds is minutes each.
        # A smoke run needs the pipeline exercised, not the sizing perfected --
        # the real run re-does any pair that is still off target. The tighter
        # Mesh_3 timeout matters more than the round cap: some Thingi surfaces
        # make Mesh_3 crawl, and a smoke run must not stall on one of them.
        prep += ["--max-rounds", 3, "--mesh3-timeout", 180]
    stage("2-prepare", prep, log_dir, args.dry_run)

    if args.governor and args.governor != "none" and not args.dry_run:
        try_set_governor(args.governor)

    where = ["--results-dir", results_dir, "--mesh-out-dir", mesh_out_dir]

    measure = [py, SCRIPTS / "run_bench.py", "--root", root,
               "--profile", "metrics" if args.metrics else "full",
               "--budget", args.budget, "--threads-max", args.threads_max,
               "--reps", args.reps] + where
    if args.smoke:
        measure += ["--calib-budget", 360, "--max-configs", 2, "--run-timeout", 600]
    stage("3-measure", measure, log_dir, args.dry_run)

    if args.skip_quality:
        print("\nStopping before the quality pass, as asked. Run it with:\n"
              "    %s %s --root %s --quality-pass %s"
              % (py, SCRIPTS / "run_bench.py", root, " ".join(str(w) for w in where)))
        return

    stage("4-quality", [py, SCRIPTS / "run_bench.py", "--root", root,
                        "--quality-pass"] + where, log_dir, args.dry_run)

    if args.dry_run:
        print("\n(dry run: nothing executed)")
        return

    bundles = sorted(root.glob("results_*.tar.gz"),
                     key=lambda p: p.stat().st_mtime, reverse=True)
    print("\n" + "=" * 72)
    print("%s DONE in %s." % ("METRICS RUN" if args.metrics else
                              "SMOKE RUN" if args.smoke else "ALL",
                              hrs(time.time() - total0)))
    if bundles:
        size = bundles[0].stat().st_size
        human = ("%.0f KB" % (size / 1e3)) if size < 1e6 else ("%.1f MB" % (size / 1e6))
        print("\nSend back this one file:\n\n    %s   (%s)" % (bundles[0], human))
    if args.metrics:
        print("\nThat is the metrics run. Report it with:\n\n"
              "    python3 %s --results %s\n"
              % (SCRIPTS / "report.py", results_dir))
    elif args.smoke:
        print("\nThat is the smoke run. Once it has been checked, start the real\n"
              "one with the same command minus --smoke -- the build and the\n"
              "prepared meshes are reused, and it writes to results/ rather than\n"
              "smoke_results/, so nothing here is disturbed:\n")
        real = [a for a in sys.argv[1:] if a != "--smoke"]
        print("    nohup %s %s %s > run.log 2>&1 &"
              % (py, Path(__file__).name, " ".join(real)))
    else:
        print("\nKeep %s -- it is the archive that lets any further metric be\n"
              "computed without re-running the benchmark." % mesh_out_dir)
    print("Stage logs: %s" % log_dir)
    print("=" * 72)


if __name__ == "__main__":
    main()
