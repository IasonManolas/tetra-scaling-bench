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

## TL;DR — the whole run

```bash
git clone https://github.com/IasonManolas/tetra-scaling-bench.git
cd tetra-scaling-bench

# 1) smoke run first, ~30-45 min. Send back the tarball it names.
nohup python3 run_all.py --surfaces-dir /path/to/thingi10k --smoke > smoke.log 2>&1 &

# 2) once that has been checked, the real run: same command, no --smoke
nohup python3 run_all.py --surfaces-dir /path/to/thingi10k > run.log 2>&1 &

# or, instead of (2), the scaling measurements only, ~45 min:
nohup python3 run_all.py --surfaces-dir /path/to/thingi10k --metrics > metrics.log 2>&1 &
```

**[RUNNING.md](RUNNING.md) is the operator's guide** — prerequisites and how to
check each one (`perf` is the likely blocker), the exact command for every mode,
how long each takes, what must be true while it measures, and a checklist for
telling a good run from a bad one before sending it back.

`--surfaces-dir` points at the Thingi10K `.stl` files (`.off` and `.ply` work too).

**Do the smoke run first.** It exercises the whole path — build, both input
pipelines, a real measurement, the quality pass, the tarball — on a handful of
surfaces, and above all it confirms that calibration can find configurations
that run ≥30 s at full width on *your* machine. If it cannot, we want to know
that in 30 minutes, not 12 hours in. The smoke run writes to `smoke_results/`
and `smoke_out_meshes/`, so it cannot contaminate the real dataset, and the
build and prepared meshes it produces are reused by the real run.

The real run is the whole thing: build, prepare inputs, measure (~12 h), derive
quality, package. It checks its prerequisites and free disk **before** starting
anything long, so it fails in seconds rather than four hours in. At the end it
names one `work/results_*.tar.gz` — send that back.

The machine should be otherwise idle for the measurement stage. `run_all.py`
tries to set the CPU governor to `performance` first and tells you if it could
not (that needs root; without it the timings drift more but are still usable).

**Interrupting is safe.** Every stage is idempotent — re-running the exact same
command skips finished builds, prepared meshes and completed measurement cells,
and carries on. `--dry-run` prints the plan without doing anything.

Useful knobs: `--budget` (measurement seconds, default 12 h), `--threads-max`,
`--limit` (how many surfaces to prepare), `--tbb-dir`, `--skip-quality`.

<details>
<summary>Running the stages by hand instead</summary>

```bash
python3 scripts/setup.py --root work                     # ~15 min
python3 scripts/prepare_meshes.py --root work \
        --surfaces-dir /path/to/thingi10k --jobs 24           # ~1-3 h, one time

nohup python3 scripts/run_bench.py --root work \
        --profile full --budget 43200 > run.log 2>&1 &   # ~12 h, unattended

python3 scripts/run_bench.py --root work --quality-pass   # ~1 h, afterwards
```

The quality pass is post-processing — it reads the meshes the run saved, so it
can happen whenever and the machine can be busy.

</details>

If TBB is somewhere CMake cannot find it, add
`--tbb-dir /path/to/dir/with/TBBConfig.cmake` to `setup.py`.

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
- The **Thingi10K** dataset (https://github.com/Thingi10K/Thingi10K). It ships
  `.stl`, which is what `prepare_meshes.py` takes by default. A repaired or
  autorefined `.off` set also works and survives CDT construction on more of the
  awkward models — pass `--prefer-format off` to take those where both exist.
  Only the 100 model ids in `meshes/thingi_ids.txt` are used.
- ~20 GB free disk (see **Disk** below)
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

### 2. Prepare the meshes (~2–4 h, one time, untimed)

```bash
python3 scripts/prepare_meshes.py --root work --surfaces-dir /path/to/thingi10k
```

Each surface produces **two** tetrahedral inputs:

- `<id>_cdt.mesh` — `make_conforming_constrained_Delaunay_triangulation_3`; the
  surface's own triangles are kept as constraints.
- `<id>_mesh3.mesh` — `make_mesh_3`, sized to hit the same cell count as that
  surface's CDT.

Both, because they are not interchangeable: on this project's earlier
measurements the CDT set and the Mesh_3 set **disagreed on the sign** of the
sequential ours-vs-main comparison (0.95 vs 1.20). A benchmark built on one
alone reports a property of the input generator as if it were a property of the
remesher. Matching the cell counts is what makes the pair controlled — the two
inputs then differ in element quality, not in size.

Mesh_3 sizing is iterative: cell count against sizing is only roughly a cube
law, so the script fits the exponent from its own last two builds and retries
until it is within 15% of the target, keeping the closest attempt. It prints
which pairs failed to match, and `report.py` leaves those out of the
CDT-vs-Mesh_3 comparison. Expect a few genuine failures — some surfaces make
Mesh_3 fail outright at finer sizings.

By default this builds only the **30 largest** surfaces of the 100, because the
sweep uses about a dozen configs and picks the biggest ones — building all 100
mostly buys preprocessing time. `--limit 0` builds them all; `--budget <seconds>`
caps the stage.

Restartable, parallel across cores, and it never leaves a half-written `.mesh`.
`--pipelines cdt` or `--pipelines mesh3` restricts it if you only want one.

### 3. Run it (one command, ~12 h, unattended)

```bash
nohup python3 scripts/run_bench.py --root work --profile full --budget 43200 \
  > run.log 2>&1 &
```

This is the whole measurement: it calibrates, sizes its own grid from what it
measured, and runs the sweep. **Nothing needs to be sent back mid-run and
nothing needs deciding partway.** Quality metrics and the report are separate
steps — they read the meshes this leaves behind, so they do not belong inside
the window where the machine has to stay idle.

It splits the budget rather than sharing it, so calibration cannot overrun into
the sweep:

| phase | share of 12 h | what it does |
|---|---|---|
| calibration | ~0.7 h | finds which `(mesh, edge factor)` pairs make a 24-thread run last ≥30 s, and measures this machine's noise and serial ratio |
| sweep | ~11.3 h | the grid, in priority waves |

The sweep runs in **priority waves**, not config by config. Wave 1 gives *every*
config its headline numbers (`par@24`, `par@1`, `seq`, `main`); wave 2 fills in
intermediate thread counts; wave 3 adds repeats. So if the budget runs out —
and on an unattended run nobody is there to notice — the result is a complete
dataset at coarser resolution, rather than three perfect configs and nothing for
the rest.

If it is interrupted, re-running the same command resumes where it stopped.

Useful knobs: `--threads-max` (default 24), `--reps` (3), `--max-configs` (12).

### 3c. The metrics run (~45 min)

```bash
nohup python3 run_all.py --surfaces-dir /path/to/thingi10k --metrics > metrics.log 2>&1 &
```

The scaling questions and nothing else, specified in
[docs/METRICS_REQUEST.md](docs/METRICS_REQUEST.md): three thread ladders with
instructions and cycles recorded, a core-pinning comparison, a lock-grid sweep,
and two instrumented runs. It skips everything the full sweep measures that a
scaling analysis does not read — the edge-factor ladder, the extra meshes, and
the `seq` and `main` reference arms, which run at one thread by definition and
so cannot move with anything tested here.

> **Two of the four phases need code that is not on the target branch yet.**
> `setup.py` clones `gsoc2026-Tetra_remeshing_parallel-imanolas` from GitHub, and
> at `bd32d445baa` that branch carries `CGAL_TR_LOCKCOUNT` but neither
> `CGAL_TR_TOPSTAGE` nor the `CGAL_TETRAHEDRAL_REMESHING_LOCK_GRID` environment
> variable. Until both are pushed, the lock-grid sweep measures the same build
> five times over and the stage-timing binary prints nothing. Neither failure
> announces itself in the log — the lock-grid table simply comes back flat, and
> the empty stdout capture is what the checklist in RUNNING.md §7 catches.

Like `--smoke`, it writes to its own directories (`metrics_results/`,
`metrics_out_meshes/`) and shares the build and prepared meshes with the real
run, so it cannot contaminate a full sweep. It also builds two extra
diagnostic binaries — the same source and the same Release flags with one macro
added each — into their own build directories, so the timed binary is never
rebuilt to get them.

**Set the CPU governor to `performance` first** if you can — the script warns if
it is not. On a `powersave` laptop we measured 25% drift between waves, enough
to make a 2-thread run look faster than a 4-thread one. Every row records
`started_at` so drift is detectable afterwards, and the report flags
non-monotonic scaling rather than presenting it as a result, but neither is a
substitute for a machine that holds its clock.

<details>
<summary>Running the phases separately instead</summary>

### 3a. Calibrate only (30 min)

```bash
python3 scripts/run_bench.py --root work --profile calibrate
```

The overnight grid cannot be chosen without knowing your machine. This run finds
which `(mesh, edge factor)` pairs make a **24-thread run last at least 30 s** —
far enough above the noise floor to resolve anything — by walking edge factors
down (0.5 → 0.25; output cells grow roughly as `factor^-3`, so that is ~8× the
work). It then measures run-to-run CV and sketches the scaling curve.

It runs under a **hard 30-minute budget** and does the most valuable work first,
so a truncated run is still useful. Its single-threaded anchors go on the
largest config that still *fits* single-threaded — a config sized for the 30 s
bar takes ~25 minutes on one thread and would blow the budget, while a
sub-second one measures process start-up rather than remeshing.

### 3b. Sweep, then quality

```bash
nohup python3 scripts/run_bench.py --root work --profile overnight > overnight.log 2>&1 &
python3 scripts/run_bench.py --root work --quality-pass
```

The sweep prints an ETA from the serial ratio calibration actually measured,
rather than assuming perfect scaling. If it dies on a signal it stops
everything, so a crash gets fixed rather than papered over.

</details>

### 4. Send back

`--profile full` writes a single `work/results_<host>_<date>.tar.gz` — send that.
Running the phases by hand instead, send:

```
work/results/results.csv
work/results/json/          per-run JSON
work/results/quality/       per-output-mesh quality JSON
work/results/env.json
work/results/calibration.json
```

A few MB either way. **Keep `work/out_meshes/` where it is** — that is the
archive, and it is what lets a metric nobody thought of today be computed later
without re-running the benchmark.

---

## Disk

About **10–20 GB** end to end, measured at ~40 bytes per output cell:

| | |
|---|---|
| Thingi10K surfaces | ~0.7 GB |
| 200 input `.mesh` files (100 surfaces × 2 pipelines) | ~1.5 GB |
| 2 shallow CGAL clones + build dirs | ~1.0 GB |
| Output meshes (≈12 configs × 9 arm-points × ~50 MB) | 5–11 GB |
| JSONs, logs, quality | <0.1 GB |

Pushing to edge factor 0.25 (≈8× the output cells) reaches ~45 GB.

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
  preprocess_cdt.cpp      surface -> constrained Delaunay .mesh
  preprocess_mesh3.cpp    surface -> Mesh_3 .mesh, at a chosen relative size
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
python3 scripts/prepare_meshes.py --root work --surfaces-dir <off> --only 124534
md5sum work/meshes/124534_cdt.mesh <known-good>/124534.mesh   # must match exactly
python3 scripts/run_bench.py --root work --profile calibrate --threads-max 4 --budget 420
# interrupt it, re-run: it must resume and still reach the same conclusions
python3 scripts/run_bench.py --root work --quality-pass
python3 scripts/report.py --results work/results
```

The mesh checksum comparison is the important one: if preprocessing drifts,
every number afterwards is measured on a different input than the reference.
