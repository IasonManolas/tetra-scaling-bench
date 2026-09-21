# tetra-scaling-bench

Measures how CGAL tetrahedral remeshing scales with core count, and what it does
to mesh quality, comparing three things:

| arm | tree | concurrency | threads |
|---|---|---|---|
| `main` | `CGAL/cgal` @ `main` | sequential (upstream has no concurrency tag) | 1 |
| `seq` | `IasonManolas/cgal` @ `gsoc2026-Tetra_remeshing_parallel-imanolas` | `CGAL::Sequential_tag` | 1 |
| `par` | same branch | `CGAL::Parallel_tag` | 1, 2, 4, 8, 12, 16, 24 |

`par@1` is measured deliberately: it is the only way to tell "the parallel
algorithm costs something even on one thread" apart from "it does not scale".

Everything so far has been measured on a **4-core** machine, so every scaling
claim in this project rests on 1→4 threads. This kit exists to get numbers on a
24-core machine.

---

## What you need

- Linux, a C++17 compiler, CMake ≥ 3.12, Python 3.8+, git
- **Boost** (≥ 1.81 in practice — CGAL's declared 1.74 minimum is stale, since
  `Constrained_triangulation_3` on main already uses `boost::unordered_flat_map`)
- **TBB** — `libtbb-dev`, or oneAPI. This is the *only* dependency the branch
  adds beyond a plain CGAL build: `Tetrahedral_remeshing` on `cgal/main` has no
  TBB includes at all, while this branch's `Parallel_tag` path has nine. Without
  it there is no parallel arm at all, so `setup.py` refuses to continue rather
  than quietly building a `bench_remesh` whose `--tag par` is a fiction.
- The **Thingi10K "fixed"/autorefined** dataset as `.off` files (~640 MB for the
  full set; the 100 models this kit uses are ~18 MB of it)
- ~60 GB free disk for the derived meshes and the saved outputs
- `matplotlib`, optionally, for plots in the final report

The machine should be **otherwise idle** while a measurement runs, and ideally on
the `performance` CPU governor — the scripts warn if it is not.

---

## Running it

### 1. Build (~15 min, mostly the CGAL clones)

```bash
python3 scripts/setup.py --root work
# or, if TBB is somewhere CMake cannot find:
python3 scripts/setup.py --root work --tbb-dir /opt/intel/oneapi/2025.2/lib/cmake/tbb
```

Clones both CGAL trees shallowly, builds four binaries, verifies `bench_remesh`
actually linked TBB, and records every SHA and binary md5 in
`work/toolchain.json`.

### 2. Prepare the meshes (~1–2 h, one time, untimed)

```bash
python3 scripts/prepare_meshes.py --root work --off-dir /path/to/thingi10k_off
```

Turns each `.off` into the constrained Delaunay `.mesh` the benchmark reads, and
writes `work/meshes/manifest.json` with each mesh's cell count. Parallel across
cores, restartable, and it never leaves a half-written `.mesh` behind.

Only the `.off` inputs (~18 MB) are distributed; the `.mesh` files (~720 MB) are
derived, which is why this regenerates them instead of shipping them.

### 3. Calibrate — **run this first and send the result back** (30 min)

```bash
python3 scripts/run_bench.py --root work --profile calibrate
```

The overnight grid cannot be chosen without knowing your machine. This run finds
which `(mesh, edge factor)` pairs make a **24-thread run last at least 30 s** —
far enough above the noise floor to resolve anything — by walking edge factors
down (0.5 → 0.25; output cells grow roughly as `factor^-3`, so that is ~8× the
work). It then measures run-to-run CV and sketches the scaling curve.

It runs under a **hard 30-minute budget** and does the most valuable work first,
so a truncated run is still useful. Stage 4 deliberately takes its
single-threaded anchors on the *smallest* config: a 30 s-at-24-threads workload
can take ten minutes on one thread and would eat the whole budget alone.

**Send back `work/results/calibration.json`.** That sizes the real run.

### 4. Overnight (~12 h, resumable)

```bash
nohup python3 scripts/run_bench.py --root work --profile overnight > overnight.log 2>&1 &
```

Prints an ETA (using the serial ratio calibration actually measured, not an
assumption of perfect scaling) before starting. If it is interrupted, just run it
again — it resumes. If anything dies on a signal it stops the whole sweep, so a
crash gets fixed rather than papered over.

### 5. Quality pass (after the sweep, on the same machine)

```bash
python3 scripts/run_bench.py --root work --quality-pass
```

### 6. Send back

```
work/results/results.csv
work/results/json/          per-run JSON
work/results/quality/       per-output-mesh quality JSON
work/results/env.json
work/results/calibration.json
```

A few MB. **Keep `work/out_meshes/` where it is** — that is the archive.

---

## Why quality is computed afterwards

The timed run records only what cannot be recomputed later — wall time, peak
RSS, exit code, and the input's average edge length (the target is derived from
it, so conformance is not recoverable from the output alone). It writes the
remeshed `.mesh` to disk and computes **no** quality metrics.

Every quality number is then derived from that saved mesh by
`mesh_quality_report`. So the measured process stays minimal, and a metric
nobody thought of today can be computed in three months without re-running the
benchmark — just re-run the quality pass over the archive.

One metric does **not** survive this split, and it is measured in the driver
instead: `Quality.Edge_Length` covers the c3t3's *complex/feature* edge set,
which a MEDIT file does not carry (`write_MEDIT` has named parameters for cells
and vertices, none for complex edges). After a round-trip every field of it
comes back as `-1`. Verified, not assumed. The sizing-conformance question is
answered by `Quality.Cell_Edge_Length` — every edge of every in-complex
tetrahedron, each counted once — which *is* computed offline. See
`driver/edge_metrics.h`.

---

## Reporting

```bash
python3 scripts/report.py --results work/results
```

Writes `report.md`, `all_runs.csv` and plots. Set totals are reported only over
the configs **every** arm measured; summing each arm over whatever it happens to
cover and printing the results side by side invites exactly the wrong
comparison.

---

## Layout

```
driver/           C++ sources; CMakeLists is configured twice, once per CGAL tree
  bench_remesh.cpp        timed driver, both tags in one binary, --tag picks
  bench_remesh_main.cpp   same instrument against upstream main
  mesh_quality_report.cpp offline quality from a saved .mesh
  mesh_quality.h          quality metrics, from the project's perf branch
  edge_metrics.h          the two edge-length metrics, and why they are split
  preprocess_cdt.cpp      .off -> .mesh
meshes/thingi_ids.txt     the 100 curated Thingi10K model ids
scripts/                  setup / prepare_meshes / run_bench / report
work/                     everything generated (gitignored)
```

### Driver CLI

```
bench_remesh <input.mesh> <iters> <edge_factor> <smooth_constrained> \
             <threads> <out.json> [--tag seq|par] [--out-mesh <path.mesh>]
```

Arguments 1–6 match the project's existing `benchmark_tetrahedral_remeshing`
exactly, so the parsers already in use on the measurement machines keep working
against this driver's output.

---

## Verifying the kit itself

On a machine that already has reference data:

```bash
python3 scripts/setup.py --root work --fresh          # builds from nothing
python3 scripts/prepare_meshes.py --root work --off-dir <off> --only 124534
md5sum work/meshes/124534.mesh <known-good>/124534.mesh   # must match exactly
python3 scripts/run_bench.py --root work --profile calibrate --threads-max 4 --budget 420
# interrupt it, re-run: it must resume and still reach the same conclusions
python3 scripts/run_bench.py --root work --quality-pass
python3 scripts/report.py --results work/results
```

The mesh checksum comparison is the important one: if preprocessing drifts,
every number afterwards is measured on a different input than the reference.
