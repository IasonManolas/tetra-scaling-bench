#!/usr/bin/env python3
"""Stage 1: turn Thingi10K surfaces into the .mesh inputs the sweep reads.

TWO pipelines, from the SAME surface:

  cdt    make_conforming_constrained_Delaunay_triangulation_3 -- the surface's
         own triangles are preserved as constraints.
  mesh3  make_mesh_3 -- a fresh volume mesh, sized to hit the same cell count as
         that surface's CDT.

Both are kept because they are not interchangeable: on this project's earlier
measurements the CDT set and the Mesh_3 set disagreed on the *sign* of the
sequential ours-vs-main comparison (0.95 vs 1.20). A benchmark built on one of
them alone reports a property of the input generator as if it were a property
of the remesher. Sizing Mesh_3 to the CDT's cell count is what makes the pair a
controlled comparison: the two inputs differ in element quality, not in size.

Input surfaces: Thingi10K ships .stl, which is the default here. .off and .ply
are also accepted, so a repaired/autorefined set works too -- see
--prefer-format and the README.

This step is UNTIMED and one-time. It writes manifest.json, which records each
mesh's cell count; run_bench.py picks its subset from there.

    python3 scripts/prepare_meshes.py --root work --surfaces-dir /path/to/thingi10k
"""
import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
IDS_FILE = REPO_ROOT / "meshes" / "thingi_ids.txt"

SURFACE_EXTS = (".stl", ".off", ".ply", ".STL", ".OFF", ".PLY")

# preprocess_cdt prints:  CDT construction: 12345 vertices, 67890 cells, 1.23s
CDT_LINE = re.compile(r"CDT construction:\s*(\d+)\s*vertices,\s*(\d+)\s*cells")
# preprocess_mesh3 prints: Mesh_3: 12345 vertices, 67890 cells in complex, 1.2s
MESH3_LINE = re.compile(r"Mesh_3:\s*(\d+)\s*vertices,\s*(\d+)\s*cells")

# Mesh_3 sizing is relative to the bounding-box diagonal. Cell count goes
# roughly as (1/h)^3, so one probe is enough to solve for the size that hits a
# target count:  h1 = h0 * (N_probe / N_target)^(1/3)
PROBE_REL = 0.06
# The lower clamp is a real constraint, not a formality: at 0.008 one surface
# saturated at 34k cells against a 65k target and every further round asked for
# a size the clamp refused. 0.002 is 64x the cell budget of 0.008, which puts
# the floor well clear of the targets this set needs.
REL_MIN, REL_MAX = 0.002, 0.20
MATCH_TOL = 0.15     # accept a Mesh_3 input within 15% of its CDT's cell count
MAX_ROUNDS = 7       # sizing attempts before keeping the closest one


def read_ids():
    ids = []
    for line in IDS_FILE.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            ids.append(line)
    return ids


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_surface(surface_dir, mesh_id, prefer="stl"):
    """Locate a surface for one Thingi10K id.

    Thingi10K's own release is .stl, so that is what a fresh download gives you
    and that is the default. A repaired/autorefined .off set is also usable and
    is more likely to survive CDT construction on the awkward models; pass
    --prefer-format off to take those first where both are present.
    """
    order = {
        "stl": (".stl", ".STL", ".off", ".OFF", ".ply", ".PLY"),
        "off": (".off", ".OFF", ".stl", ".STL", ".ply", ".PLY"),
    }[prefer]
    for ext in order:
        p = surface_dir / (mesh_id + ext)
        if p.exists():
            return p
    return None


def _run(cmd, timeout):
    try:
        return subprocess.run(cmd, text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None


def preprocess_cdt(exe, src, mesh_path, timeout):
    """Writes to a temp path and renames, so an interrupted run never leaves a
    truncated .mesh that a later run would accept as finished."""
    tmp = Path(str(mesh_path) + ".partial")
    r = _run([str(exe), str(src), str(tmp)], timeout)
    if r is None:
        tmp.unlink(missing_ok=True)
        return {"status": "timeout"}
    if r.returncode != 0 or not tmp.exists():
        tmp.unlink(missing_ok=True)
        return {"status": "failed", "rc": r.returncode, "stderr": (r.stderr or "")[-400:]}

    m = CDT_LINE.search(r.stdout or "")
    os.replace(tmp, mesh_path)
    return {"status": "ok", "pipeline": "cdt",
            "vertices": int(m.group(1)) if m else None,
            "cells": int(m.group(2)) if m else None}


def _mesh3_build(exe, src, out_path, rel, timeout):
    """One Mesh_3 build. Returns (cells, vertices) or None."""
    r = _run([str(exe), str(src), str(out_path), str(rel)], timeout)
    if r is None or r.returncode != 0 or not out_path.exists():
        out_path.unlink(missing_ok=True)
        return None, (r.stderr[-400:] if r is not None and r.stderr else ""), \
               ("timeout" if r is None else "failed")
    m = MESH3_LINE.search(r.stdout or "")
    if not m or int(m.group(2)) == 0:
        out_path.unlink(missing_ok=True)
        return None, (r.stdout or "")[-400:], "failed"
    return (int(m.group(2)), int(m.group(1))), "", "ok"


def preprocess_mesh3(exe, src, mesh_path, target_cells, timeout,
                     tol=MATCH_TOL, max_rounds=MAX_ROUNDS):
    """Build a Mesh_3 input sized to hit target_cells.

    Cell count against sizing is only roughly a cube law -- feature detection,
    the facet criteria and boundary effects all bend it -- so a single probe
    plus a fixed 1/3 exponent undershoots badly: measured 52% and 74% below
    target on the first two surfaces tried. Instead, fit the exponent from the
    two most recent builds and iterate until the count is within `tol`, keeping
    whichever attempt landed closest.

    This is untimed preprocessing, so the extra rounds cost only wall clock,
    never measurement validity.
    """
    tmp = Path(str(mesh_path) + ".partial")
    attempts = []          # (rel, cells)
    best = None            # (relative error, rel, cells, vertices)

    rel = PROBE_REL
    stopped_early = None
    for rnd in range(max_rounds):
        res, err, status = _mesh3_build(exe, src, tmp, rel, timeout)
        if res is None:
            if attempts:
                # Mesh_3 can fail outright at a finer sizing on an awkward
                # model even though a coarser one worked. Keep the best attempt
                # so far rather than losing the input, but record that the
                # search was cut short -- this is exactly the case where the
                # pair ends up unmatched and must not be read as controlled.
                stopped_early = "%s at rel=%.5f round %d" % (status, rel, rnd + 1)
                break
            return {"status": status, "stage": "probe", "stderr": err,
                    "facet_size_rel": rel}

        cells, verts = res
        attempts.append((rel, cells))

        if not target_cells:
            best = (0.0, rel, cells, verts)
            os.replace(tmp, mesh_path)
            break

        relerr = abs(cells - target_cells) / float(target_cells)
        if best is None or relerr < best[0]:
            best = (relerr, rel, cells, verts)
            os.replace(tmp, mesh_path)          # keep the closest attempt
        else:
            tmp.unlink(missing_ok=True)

        if relerr <= tol:
            break

        # cells ~ rel^-k. With two points we can measure k instead of assuming
        # it; with one, fall back to the ideal 3.
        k = 3.0
        if len(attempts) >= 2:
            (r0, c0), (r1, c1) = attempts[-2], attempts[-1]
            if r0 != r1 and c0 > 0 and c1 > 0:
                k = math.log(c1 / float(c0)) / math.log(r0 / float(r1))
                if not (0.5 < k < 8.0):         # nonsense fit; ignore it
                    k = 3.0
        want = rel * (cells / float(target_cells)) ** (1.0 / k)
        nxt = max(REL_MIN, min(REL_MAX, want))
        if nxt != want:
            stopped_early = ("sizing clamped at rel=%.5f (wanted %.5f); the "
                             "target needs a finer floor" % (nxt, want))
        if abs(nxt - rel) / rel < 0.01:
            break                                # converged, or stuck on a clamp
        rel = nxt

    tmp.unlink(missing_ok=True)
    if best is None or not mesh_path.exists():
        return {"status": "failed", "stage": "final"}

    relerr, rel, cells, verts = best
    return {"status": "ok", "pipeline": "mesh3", "vertices": verts, "cells": cells,
            "facet_size_rel": rel, "target_cells": target_cells,
            "size_match_error": relerr if target_cells else None,
            "rounds": len(attempts),
            "sizing_stopped_early": stopped_early,
            "matched": bool(target_cells and relerr <= tol)}


def build_pair(bins, src, mesh_dir, mesh_id, pipelines, timeout, force,
               tol=MATCH_TOL, max_rounds=MAX_ROUNDS, mesh3_timeout=None):
    """One surface -> up to two .mesh inputs. Returns {key: record}.

    The mesh3 build needs the cdt cell count as its target, so when both are
    requested they run in sequence on the same surface rather than as two
    independent jobs.
    """
    out = {}
    cdt_path = mesh_dir / ("%s_cdt.mesh" % mesh_id)
    m3_path = mesh_dir / ("%s_mesh3.mesh" % mesh_id)
    cdt_cells = None

    if "cdt" in pipelines:
        if cdt_path.exists() and not force:
            out["%s_cdt" % mesh_id] = {"status": "present"}
        else:
            res = preprocess_cdt(bins["preprocess_cdt"], src, cdt_path, timeout)
            out["%s_cdt" % mesh_id] = res
            if res["status"] == "ok":
                cdt_cells = res["cells"]

    if "mesh3" in pipelines:
        if m3_path.exists() and not force:
            out["%s_mesh3" % mesh_id] = {"status": "present"}
        else:
            # Target the CDT's count so the pair is matched on size. If the CDT
            # is not available (not requested, or it failed), fall back to the
            # probe size and record that the pair is unmatched.
            if cdt_cells is None and cdt_path.exists():
                cdt_cells = _cells_from_medit(cdt_path)
            out["%s_mesh3" % mesh_id] = preprocess_mesh3(
                bins["preprocess_mesh3"], src, m3_path, cdt_cells,
                mesh3_timeout or timeout, tol, max_rounds)
    return out


def _cells_from_medit(path):
    """Read the Tetrahedra count out of a MEDIT file without parsing the rest."""
    try:
        with open(path, "r", errors="ignore") as f:
            prev = ""
            for line in f:
                if prev.strip() == "Tetrahedra":
                    return int(line.strip())
                prev = line
    except Exception:
        pass
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(REPO_ROOT / "work"))
    # Named --off-dir when this only took .off. Thingi10K ships .stl and that is
    # now the default, so the old name misdescribes the argument; it stays as an
    # alias so existing command lines keep working.
    ap.add_argument("--surfaces-dir", "--off-dir", dest="surfaces_dir",
                    required=True, metavar="DIR",
                    help="directory of Thingi10K surfaces (.stl, .off or .ply)")
    ap.add_argument("--prefer-format", choices=("stl", "off"), default="stl",
                    help="which to take when a model exists in both formats "
                         "(default stl, which is what Thingi10K ships)")
    ap.add_argument("--match-tol", type=float, default=MATCH_TOL,
                    help="accept a Mesh_3 input within this fraction of its CDT's "
                         "cell count (default %.2f)" % MATCH_TOL)
    ap.add_argument("--max-rounds", type=int, default=MAX_ROUNDS,
                    help="Mesh_3 sizing attempts before keeping the closest")
    ap.add_argument("--mesh-dir", help="where to write .mesh (default <root>/meshes)")
    ap.add_argument("--pipelines", default="cdt,mesh3",
                    help="which input generators to build: cdt, mesh3, or both")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4) // 2))
    ap.add_argument("--timeout", type=int, default=3600,
                    help="per-build timeout in seconds, for the CDT stage")
    # Mesh_3 gets its own, much shorter limit. It is run up to --max-rounds
    # times per surface, so a shared 1-hour timeout means one pathological
    # surface can stall preprocessing for hours -- observed: a coarse probe
    # that should take seconds hitting the full 3600s. A Mesh_3 build that has
    # not finished in a few minutes at probe size is not going to produce a
    # usable matched pair, and the CDT input for that surface is unaffected.
    ap.add_argument("--mesh3-timeout", type=int, default=600,
                    help="per-build timeout for each Mesh_3 sizing round")
    ap.add_argument("--only", nargs="*", help="only these ids (for spot checks)")
    ap.add_argument("--limit", type=int, default=30,
                    help="build only the N largest surfaces (0 = all 100). The "
                         "sweep uses roughly a dozen configs and picks the "
                         "biggest, so building all of them mostly buys "
                         "preprocessing time. Default 30.")
    ap.add_argument("--budget", type=float, default=0.0,
                    help="stop starting new surfaces after this many seconds "
                         "(0 = no limit)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    pipelines = [p.strip() for p in args.pipelines.split(",") if p.strip()]
    bad = [p for p in pipelines if p not in ("cdt", "mesh3")]
    if bad:
        sys.exit("Unknown pipeline(s): %s" % bad)

    root = Path(args.root).resolve()
    surfaces_dir = Path(args.surfaces_dir).resolve()
    mesh_dir = Path(args.mesh_dir).resolve() if args.mesh_dir else root / "meshes"
    mesh_dir.mkdir(parents=True, exist_ok=True)

    tc_path = root / "toolchain.json"
    if not tc_path.exists():
        sys.exit("No %s -- run scripts/setup.py first." % tc_path)
    tc = json.loads(tc_path.read_text())
    bins = {}
    for name in ("preprocess_cdt", "preprocess_mesh3"):
        if name not in tc["binaries"]:
            sys.exit("%s is missing from toolchain.json -- re-run scripts/setup.py "
                     "(this kit gained a second input pipeline)." % name)
        bins[name] = Path(tc["binaries"][name]["path"])

    ids = args.only if args.only else read_ids()

    found = {i: find_surface(surfaces_dir, i, args.prefer_format) for i in ids}
    missing = [i for i, p in found.items() if p is None]
    if missing:
        print("%d of %d surfaces are missing from %s" % (len(missing), len(ids), surfaces_dir),
              file=sys.stderr)
        print("first few: %s" % missing[:10], file=sys.stderr)
        if len(missing) == len(ids):
            sys.exit("None of the inputs are there. Expected <id>.stl / <id>.off / "
                     "<id>.ply -- is --off-dir right?")
        print("continuing with the %d that are present\n" % (len(ids) - len(missing)),
              file=sys.stderr)
        ids = [i for i in ids if found[i] is not None]

    # Biggest surfaces first: they are the ones that produce the long-running
    # configs the scaling question needs, and if the budget runs out it should
    # run out on the small ones. File size is a fine proxy for triangle count.
    ids.sort(key=lambda i: found[i].stat().st_size, reverse=True)
    if args.limit and not args.only and len(ids) > args.limit:
        print("Building the %d largest of %d surfaces (--limit 0 for all)."
              % (args.limit, len(ids)))
        ids = ids[:args.limit]

    exts = sorted({found[i].suffix.lower() for i in ids})
    print("%d surfaces, formats: %s" % (len(ids), ", ".join(exts)))
    print("pipelines: %s   jobs: %d" % (", ".join(pipelines), args.jobs))
    if ".stl" in exts and ".off" not in exts:
        print("\nNote: these are raw Thingi10K .stl files. Many Thingi10K models are\n"
              "      not solid (self-intersecting or open), and CDT construction will\n"
              "      fail on those. The id list here is already filtered to the solid\n"
              "      subset, but if you see widespread failures, use a repaired .off\n"
              "      set instead -- see the README.\n")

    manifest_path = mesh_dir / "manifest.json"
    manifest = {}
    if manifest_path.exists() and not args.force:
        manifest = json.loads(manifest_path.read_text()).get("meshes", {})

    t0 = time.time()
    done = 0
    deadline = (t0 + args.budget) if args.budget else None
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {}
        for i in ids:
            if deadline and time.time() > deadline:
                print("Prep budget spent; %d surface(s) not started. Re-run to "
                      "continue -- finished ones are skipped."
                      % (len(ids) - len(futs)), file=sys.stderr)
                break
            futs[pool.submit(build_pair, bins, found[i], mesh_dir, i,
                             pipelines, args.timeout, args.force,
                             args.match_tol, args.max_rounds,
                             args.mesh3_timeout)] = i
        for fut in concurrent.futures.as_completed(futs):
            i = futs[fut]
            done += 1
            try:
                results = fut.result()
            except Exception as e:
                print("[%3d/%3d] %-10s EXCEPTION %s" % (done, len(futs), i, e), flush=True)
                continue
            for key, res in results.items():
                if res["status"] == "ok":
                    mp = mesh_dir / (key + ".mesh")
                    rec = {"id": i, "key": key, "pipeline": res["pipeline"],
                           "path": str(mp), "bytes": mp.stat().st_size,
                           "sha256": sha256(mp), "vertices": res.get("vertices"),
                           "cells": res.get("cells"), "source": str(found[i])}
                    for extra in ("facet_size_rel", "probe_cells", "target_cells",
                                  "size_match_error", "rounds", "matched",
                                  "sizing_stopped_early"):
                        if extra in res:
                            rec[extra] = res[extra]
                    manifest[key] = rec
                    print("[%3d/%3d] %-16s ok  cells=%s%s" % (
                        done, len(futs), key, res.get("cells"),
                        ("  (target %s)" % res["target_cells"])
                        if res.get("target_cells") else ""), flush=True)
                elif res["status"] == "present":
                    pass
                else:
                    manifest.pop(key, None)
                    print("[%3d/%3d] %-16s %s%s %s" % (
                        done, len(futs), key, res["status"],
                        " (%s)" % res["stage"] if res.get("stage") else "",
                        res.get("stderr", "")), flush=True)

    ok = {k: v for k, v in manifest.items()
          if v.get("cells") and Path(v["path"]).exists()}
    manifest_path.write_text(json.dumps({
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "surfaces_dir": str(surfaces_dir),
        "mesh_dir": str(mesh_dir),
        "pipelines": pipelines,
        "preprocess_cdt_md5": tc["binaries"]["preprocess_cdt"]["md5"],
        "preprocess_mesh3_md5": tc["binaries"]["preprocess_mesh3"]["md5"],
        "meshes": manifest,
    }, indent=2))

    by_pipe = {}
    for v in ok.values():
        by_pipe.setdefault(v["pipeline"], []).append(v)
    total_bytes = sum(v["bytes"] for v in ok.values())
    print("\n%d inputs ready (%s), %.1f GB, manifest at %s" % (
        len(ok), ", ".join("%s: %d" % (p, len(v)) for p, v in sorted(by_pipe.items())),
        total_bytes / 1e9, manifest_path))
    print("Elapsed %.0fs" % (time.time() - t0))

    # How well the Mesh_3 sizing hit its CDT target. A pair that is badly off
    # is not a controlled comparison any more, so it should be visible here
    # rather than discovered in the report.
    pairs = [v for v in ok.values()
             if v["pipeline"] == "mesh3" and v.get("target_cells")]
    if pairs:
        offs = sorted(abs(v["cells"] - v["target_cells"]) / v["target_cells"]
                      for v in pairs)
        matched = [v for v in pairs if v.get("matched")]
        print("\nMesh_3/CDT size match over %d pairs: %d within %.0f%%, median %.1f%% off"
              % (len(pairs), len(matched), 100 * args.match_tol,
                 100 * offs[len(offs) // 2]))

        bad = sorted((v for v in pairs if not v.get("matched")),
                     key=lambda v: -abs(v["cells"] - v["target_cells"]) / v["target_cells"])
        if bad:
            print("\n%d pair(s) did NOT reach the size target. These are still usable\n"
                  "inputs, but they are not a controlled CDT-vs-Mesh_3 comparison --\n"
                  "they differ in size as well as in element quality, and report.py\n"
                  "leaves them out of that comparison:" % len(bad))
            for v in bad[:10]:
                off = 100 * abs(v["cells"] - v["target_cells"]) / v["target_cells"]
                print("  %-16s %7d cells vs target %7d  (%.0f%% off, %s rounds)%s"
                      % (v["key"], v["cells"], v["target_cells"], off,
                         v.get("rounds"),
                         "  [sizing stopped: %s]" % v["sizing_stopped_early"]
                         if v.get("sizing_stopped_early") else ""))

    print("\nLargest inputs by cell count:")
    for v in sorted(ok.values(), key=lambda m: m["cells"], reverse=True)[:10]:
        print("  %-16s %9d cells  %6.1f MB" % (v["key"], v["cells"], v["bytes"] / 1e6))
    print("\nNext:  python3 scripts/run_bench.py --root %s --profile calibrate" % root)


if __name__ == "__main__":
    main()
