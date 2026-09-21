#!/usr/bin/env python3
"""Stage 1: turn Thingi10K surfaces (.off) into the .mesh inputs the sweep reads.

Only the ~18 MB of .off inputs need to travel; the ~720 MB of .mesh files are
derived, which is why this regenerates them locally instead of shipping them.

This step is UNTIMED and one-time. It writes manifest.json, which records each
mesh's cell count -- run_bench.py picks its subset by cell count from there, so
"the biggest N inputs" is data about this machine's actual meshes rather than a
list hardcoded from another machine.

    python3 scripts/prepare_meshes.py --root work --off-dir /path/to/thingi_off
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
IDS_FILE = REPO_ROOT / "meshes" / "thingi_ids.txt"

# preprocess_cdt prints e.g.
#   CDT construction: 12345 vertices, 67890 cells, 1.23s
CDT_LINE = re.compile(r"CDT construction:\s*(\d+)\s*vertices,\s*(\d+)\s*cells")


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


def find_off(off_dir, mesh_id):
    for name in ("%s.off" % mesh_id, "%s.OFF" % mesh_id):
        p = off_dir / name
        if p.exists():
            return p
    return None


def preprocess_one(exe, off_path, mesh_path, timeout):
    """Build one .mesh. Writes to a temp path and renames, so an interrupted run
    never leaves a truncated .mesh that a later run would accept as done."""
    tmp = mesh_path.with_suffix(".mesh.partial")
    try:
        r = subprocess.run([str(exe), str(off_path), str(tmp)],
                           text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        tmp.unlink(missing_ok=True)
        return {"status": "timeout"}

    if r.returncode != 0 or not tmp.exists():
        tmp.unlink(missing_ok=True)
        return {"status": "failed", "rc": r.returncode,
                "stderr": (r.stderr or "")[-400:]}

    verts = cells = None
    m = CDT_LINE.search(r.stdout or "")
    if m:
        verts, cells = int(m.group(1)), int(m.group(2))

    os.replace(tmp, mesh_path)
    return {"status": "ok", "vertices": verts, "cells": cells}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(REPO_ROOT / "work"))
    ap.add_argument("--off-dir", required=True,
                    help="directory holding the Thingi10K 'fixed' .off files")
    ap.add_argument("--mesh-dir", help="where to write .mesh (default <root>/meshes)")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4) // 2),
                    help="parallel preprocessing jobs (each is single-threaded)")
    ap.add_argument("--timeout", type=int, default=3600,
                    help="per-mesh CDT timeout in seconds")
    ap.add_argument("--only", nargs="*", help="only these ids (for spot checks)")
    ap.add_argument("--force", action="store_true",
                    help="rebuild even if the .mesh already exists")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    off_dir = Path(args.off_dir).resolve()
    mesh_dir = Path(args.mesh_dir).resolve() if args.mesh_dir else root / "meshes"
    mesh_dir.mkdir(parents=True, exist_ok=True)

    tc_path = root / "toolchain.json"
    if not tc_path.exists():
        sys.exit("No %s -- run scripts/setup.py first." % tc_path)
    exe = Path(json.loads(tc_path.read_text())["binaries"]["preprocess_cdt"]["path"])

    ids = args.only if args.only else read_ids()

    missing_off = [i for i in ids if find_off(off_dir, i) is None]
    if missing_off:
        print("%d of %d .off inputs are missing from %s" %
              (len(missing_off), len(ids), off_dir), file=sys.stderr)
        print("first few: %s" % missing_off[:10], file=sys.stderr)
        if len(missing_off) == len(ids):
            sys.exit("None of the inputs are there; is --off-dir right?")
        print("continuing with the %d that are present\n" % (len(ids) - len(missing_off)),
              file=sys.stderr)
        ids = [i for i in ids if i not in set(missing_off)]

    manifest_path = mesh_dir / "manifest.json"
    manifest = {}
    if manifest_path.exists() and not args.force:
        manifest = json.loads(manifest_path.read_text()).get("meshes", {})

    todo = []
    for i in ids:
        mp = mesh_dir / ("%s.mesh" % i)
        if not args.force and mp.exists() and i in manifest:
            continue
        todo.append(i)

    print("%d ids, %d already prepared, %d to build, %d jobs" %
          (len(ids), len(ids) - len(todo), len(todo), args.jobs))

    t0 = time.time()
    done = 0
    if todo:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futs = {}
            for i in todo:
                futs[pool.submit(preprocess_one, exe, find_off(off_dir, i),
                                 mesh_dir / ("%s.mesh" % i), args.timeout)] = i
            for fut in concurrent.futures.as_completed(futs):
                i = futs[fut]
                res = fut.result()
                done += 1
                if res["status"] == "ok":
                    mp = mesh_dir / ("%s.mesh" % i)
                    manifest[i] = {
                        "id": i,
                        "path": str(mp),
                        "bytes": mp.stat().st_size,
                        "sha256": sha256(mp),
                        "cdt_vertices": res.get("vertices"),
                        "cdt_cells": res.get("cells"),
                    }
                    print("[%3d/%3d] %-10s ok  cells=%s" %
                          (done, len(todo), i, res.get("cells")), flush=True)
                else:
                    manifest.pop(i, None)
                    print("[%3d/%3d] %-10s %s %s" %
                          (done, len(todo), i, res["status"],
                           res.get("stderr", "")), flush=True)

    ok = {k: v for k, v in manifest.items() if v.get("cdt_cells")}
    by_size = sorted(ok.values(), key=lambda m: m["cdt_cells"], reverse=True)

    manifest_path.write_text(json.dumps({
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "off_dir": str(off_dir),
        "mesh_dir": str(mesh_dir),
        "preprocess_cdt_md5": json.loads(tc_path.read_text())["binaries"]["preprocess_cdt"]["md5"],
        "meshes": manifest,
    }, indent=2))

    total_bytes = sum(m["bytes"] for m in ok.values())
    print("\n%d meshes ready, %.1f GB, manifest at %s" %
          (len(ok), total_bytes / 1e9, manifest_path))
    print("Elapsed %.0fs" % (time.time() - t0))
    print("\nLargest inputs by CDT cell count:")
    for m in by_size[:10]:
        print("  %-10s %9d cells  %6.1f MB" %
              (m["id"], m["cdt_cells"], m["bytes"] / 1e6))
    print("\nNext:  python3 scripts/run_bench.py --root %s --profile calibrate" % root)


if __name__ == "__main__":
    main()
