#!/usr/bin/env python3
"""Clone the two CGAL trees this benchmark compares and build all four binaries.

    ours -> github.com/IasonManolas/cgal @ scale24
    main -> github.com/CGAL/cgal @ main (tip)

CGAL is header-only for our purposes, so "building CGAL" means pointing CGAL_DIR
at a source checkout -- there is no CGAL install step. The two arms need two
separate configure directories, because CGAL_DIR is a configure-time global.

Run from anywhere:  python3 scripts/setup.py --root <workdir>
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

OURS_URL = "https://github.com/IasonManolas/cgal.git"
# The branch under test. It moved here from
# gsoc2026-Tetra_remeshing_parallel-imanolas on 2026-09-22: scale24 is that
# branch plus the three commits the scaling work needs -- CGAL_TR_TOPSTAGE,
# CGAL_TR_LOCKCOUNT and the CGAL_TETRAHEDRAL_REMESHING_LOCK_GRID override --
# and its tip was 4cc06258243 when this was written. Both the `par` and `seq`
# arms build from it; they differ by concurrency tag, not by branch.
OURS_REF = "scale24"
MAIN_URL = "https://github.com/CGAL/cgal.git"
MAIN_REF = "main"


def run(cmd, **kw):
    print("  $ " + " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def capture(cmd, cwd=None):
    return subprocess.run(cmd, cwd=cwd, check=True, text=True,
                          capture_output=True).stdout.strip()


def clone_or_update(url, ref, dest, offline):
    if dest.exists():
        if offline:
            print("[trees] %s: present, --offline so not fetching" % dest.name)
        else:
            print("[trees] %s: fetching %s" % (dest.name, ref))
            run(["git", "-C", str(dest), "fetch", "--depth", "1", "origin", ref])
            run(["git", "-C", str(dest), "checkout", "--detach", "FETCH_HEAD"])
    else:
        print("[trees] %s: cloning %s @ %s" % (dest.name, url, ref))
        # A shallow single-branch clone of CGAL is a few hundred MB rather than
        # several GB of history nobody here needs.
        run(["git", "clone", "--depth", "1", "--branch", ref, "--single-branch",
             url, str(dest)])
    sha = capture(["git", "-C", str(dest), "rev-parse", "HEAD"])
    print("[trees] %s: %s" % (dest.name, sha))
    return sha


def check_prereqs(tbb_dir):
    """Fail early and loudly, rather than 20 minutes into a build."""
    problems = []

    for tool in ("cmake", "git"):
        if shutil.which(tool) is None:
            problems.append("'%s' not found on PATH" % tool)

    if shutil.which("cmake"):
        try:
            ver = capture(["cmake", "--version"]).splitlines()[0].split()[-1]
            major, minor = (int(x) for x in ver.split(".")[:2])
            if (major, minor) < (3, 12):
                problems.append("cmake %s is too old; need >= 3.12" % ver)
        except Exception:
            pass

    # TBB is the ONE dependency this branch adds beyond what a plain CGAL build
    # needs: cgal/main's Tetrahedral_remeshing has no TBB includes at all, while
    # this branch's Parallel_tag path has nine. Everything else it newly
    # includes is Boost core or the stdlib, and boost::unordered_flat_map is
    # already used by Constrained_triangulation_3 on main.
    if tbb_dir is None:
        hits = []
        for base in ("/usr/lib", "/usr/local/lib", "/opt"):
            p = Path(base)
            if p.exists():
                hits += list(p.rglob("TBBConfig.cmake"))
                if hits:
                    break
        if not hits:
            problems.append(
                "TBB not found and --tbb-dir not given. TBB is REQUIRED for the\n"
                "      parallel arm: without it there is no --tag par and the\n"
                "      scaling measurement is impossible. Install libtbb-dev (or\n"
                "      oneAPI) and/or pass --tbb-dir <dir with TBBConfig.cmake>.")

    if problems:
        print("\nPrerequisites not met:\n", file=sys.stderr)
        for p in problems:
            print("  - %s" % p, file=sys.stderr)
        sys.exit(1)


def cmake_build(arm, cgal_dir, build_dir, tbb_dir, jobs, fresh,
                defines=None, target=None):
    if fresh and build_dir.exists():
        shutil.rmtree(build_dir)
    cfg = ["cmake", "-S", str(REPO_ROOT / "driver"), "-B", str(build_dir),
           "-DCMAKE_BUILD_TYPE=Release", "-DARM=" + arm,
           "-DCGAL_DIR=" + str(cgal_dir)]
    if tbb_dir and arm == "ours":
        cfg.append("-DTBB_DIR=" + str(tbb_dir))
    if defines:
        cfg += ["-DBENCH_DEFINES=" + defines, "-DBENCH_REMESH_ONLY=ON"]
    run(cfg)
    build = ["cmake", "--build", str(build_dir), "-j", str(jobs)]
    if target:
        build += ["--target", target]
    run(build)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(REPO_ROOT / "work"),
                    help="where to put the CGAL checkouts and build dirs")
    ap.add_argument("--ours-dir", help="use an existing checkout instead of cloning ours")
    ap.add_argument("--main-dir", help="use an existing checkout instead of cloning main")
    ap.add_argument("--tbb-dir", help="directory containing TBBConfig.cmake")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--skip-main", action="store_true",
                    help="do not clone or build the upstream reference arm. The "
                         "'main' arm is a one-thread reference and no scaling "
                         "measurement uses it, so this is how you build against a "
                         "local branch -- with --ours-dir -- without also needing a "
                         "CGAL main checkout. --profile metrics never runs it; the "
                         "full sweep does, so do not use this for a real sweep.")
    ap.add_argument("--offline", action="store_true",
                    help="do not fetch; use the checkouts as they are")
    ap.add_argument("--fresh", action="store_true",
                    help="delete the build dirs first (proves the kit is self-contained)")
    ap.add_argument("--diagnostics", action="store_true",
                    help="also build the two instrumented bench_remesh variants "
                         "(-DCGAL_TR_TOPSTAGE and -DCGAL_TR_LOCKCOUNT) that "
                         "--profile metrics runs. Adds ~5-10 min of compiling and "
                         "does not touch the timed binaries.")
    args = ap.parse_args()

    check_prereqs(args.tbb_dir)

    root = Path(args.root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    if args.ours_dir:
        ours = Path(args.ours_dir).resolve()
        ours_sha = capture(["git", "-C", str(ours), "rev-parse", "HEAD"])
        print("[trees] ours: using %s @ %s" % (ours, ours_sha))
    else:
        ours = root / "cgal-ours"
        ours_sha = clone_or_update(OURS_URL, OURS_REF, ours, args.offline)

    if args.skip_main:
        main_tree, main_sha = None, None
        print("[trees] main: SKIPPED. The reference arm will not exist, so the "
              "full sweep cannot run here -- only --profile metrics, which uses "
              "no reference arm.")
    elif args.main_dir:
        main_tree = Path(args.main_dir).resolve()
        main_sha = capture(["git", "-C", str(main_tree), "rev-parse", "HEAD"])
        print("[trees] main: using %s @ %s" % (main_tree, main_sha))
    else:
        main_tree = root / "cgal-main"
        main_sha = clone_or_update(MAIN_URL, MAIN_REF, main_tree, args.offline)

    print("\n[build] ours arm (bench_remesh, preprocess_cdt, preprocess_mesh3, "
          "mesh_quality_report)")
    cmake_build("ours", ours, root / "build_ours", args.tbb_dir, args.jobs, args.fresh)

    if not args.skip_main:
        print("\n[build] main arm (bench_remesh_main)")
        cmake_build("main", main_tree, root / "build_main", None, args.jobs,
                    args.fresh)

    # The two instrumented binaries. They are diagnostics, not timed arms: same
    # source, same Release flags, one macro added each, and each in its own
    # build directory so asking for them never rebuilds the binary the sweep is
    # measuring with.
    diag_bins = {}
    if args.diagnostics:
        for macro, name in (("CGAL_TR_TOPSTAGE", "bench_remesh_topstage"),
                            ("CGAL_TR_LOCKCOUNT", "bench_remesh_lockcount")):
            bdir = root / ("build_ours_" + macro.rsplit("_", 1)[-1].lower())
            print("\n[build] diagnostic arm (-D%s)" % macro)
            cmake_build("ours", ours, bdir, args.tbb_dir, args.jobs, args.fresh,
                        defines=macro, target="bench_remesh")
            built = bdir / "bench_remesh"
            if not built.exists():
                sys.exit("Diagnostic build %s produced no bench_remesh." % macro)
            diag_bins[name] = {"path": str(built), "md5": md5(built),
                               "define": macro}

    bins = {
        "bench_remesh":        root / "build_ours" / "bench_remesh",
        "preprocess_cdt":      root / "build_ours" / "preprocess_cdt",
        "preprocess_mesh3":    root / "build_ours" / "preprocess_mesh3",
        "mesh_quality_report": root / "build_ours" / "mesh_quality_report",
    }
    if not args.skip_main:
        bins["bench_remesh_main"] = root / "build_main" / "bench_remesh_main"
    missing = [n for n, p in bins.items() if not p.exists()]
    if missing:
        sys.exit("Build finished but these binaries are missing: %s" % missing)

    # A bench_remesh that silently lost TBB would still accept --tag par and
    # would still produce plausible-looking numbers. Prove it actually linked.
    if sys.platform.startswith("linux") and shutil.which("ldd"):
        ldd = subprocess.run(["ldd", str(bins["bench_remesh"])],
                             text=True, capture_output=True).stdout
        if "tbb" not in ldd:
            sys.exit("bench_remesh did not link TBB -- the parallel arm would be a\n"
                     "fiction. Re-run with --tbb-dir pointing at a real TBB install.")
        print("\n[check] bench_remesh links TBB: ok")

    toolchain = {
        "ours_sha": ours_sha, "ours_dir": str(ours),
        "ours_ref": args.ours_dir or OURS_REF,
        "main_sha": main_sha,
        "main_dir": str(main_tree) if main_tree else None,
        "main_ref": None if args.skip_main else MAIN_REF,
        "binaries": dict((n, {"path": str(p), "md5": md5(p)}) for n, p in bins.items()),
        # Kept OUT of "binaries" on purpose: run_bench.py's toolchain lock
        # compares that dict, so listing the diagnostic binaries there would make
        # a results directory refuse to resume merely because they were built.
        # They are never used for a timed comparison.
        "diag_binaries": diag_bins,
        "cmake": capture(["cmake", "--version"]).splitlines()[0],
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    out = root / "toolchain.json"
    out.write_text(json.dumps(toolchain, indent=2))

    print("\nDone in %.0fs. Wrote %s" % (time.time() - t0, out))
    print("Next:  python3 scripts/prepare_meshes.py --root %s --off-dir <thingi .off dir>" % root)


if __name__ == "__main__":
    main()
