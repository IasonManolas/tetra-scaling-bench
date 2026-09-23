# Running this on your machine

Everything you need, in the order you need it. You should not have to ask us
anything; if a step does not behave as described here, that is a bug in the kit
and worth telling us about.

## What is being measured

The kit builds the remesher from **`IasonManolas/cgal` @ `scale24`** — tip
`4cc06258243` when this was written — and compares it against upstream
`CGAL/cgal` @ `main`. The `seq` and `par` arms are one binary built from that
one tree; they differ by concurrency tag, not by branch.

> **Two things moved on 2026-09-22.** The branch under test changed from
> `gsoc2026-Tetra_remeshing_parallel-imanolas` to `scale24`, which is that
> branch plus the three commits the scaling measurements need. And the kit
> itself gained `instructions`, `cycles` and `task_clock_ms` columns and the
> metrics run. **A `git pull` on your `main` brings both.** So a results
> tarball you produced before that date measured a different branch with a
> narrower results table, and the two are not comparable — do not merge them.
> If you are unsure which you are holding, `toolchain_lock.json` inside the
> tarball names the branch and the commit it was built from.

There are three runs in the kit, and they are separate things:

| run | what it is for | roughly |
|---|---|---|
| **smoke** | proves the kit works on your machine before you commit it to the long one | 30–45 min |
| **full sweep** | the whole dataset: scaling and mesh quality, across many meshes and edge sizes | ~14 h |
| **metrics** | the scaling measurements only — two thread ladders, two instrumented runs, two setting sweeps | ~45 min |

The metrics run is the one we are asking for now. It writes to its own
directory, so it cannot disturb a full sweep you have already done.

---

## 1. Prerequisites, and how to check each one

Run these checks before anything else. Each block says what a failure looks
like and what to do about it.

### `perf`, and permission to use it

This is the one likely to block you. The kit wraps every timed run in
`perf stat` to record instructions and cycles. CPU seconds cannot substitute
for them: the clock falls as more cores light up, so CPU seconds rise even when
the work does not.

```bash
perf stat -e instructions,cycles,task-clock true
```

**Good:** a small table with numbers next to `instructions` and `cycles`.

**Bad, "No permission to enable ...":** the kernel is refusing. Check and fix:

```bash
cat /proc/sys/kernel/perf_event_paranoid      # needs to be 2 or lower
sudo sysctl -w kernel.perf_event_paranoid=1   # takes effect immediately
```

A value of `1` is enough and is what we suggest. `3` (the default on some
Ubuntu builds) blocks it entirely. To make it survive a reboot, put
`kernel.perf_event_paranoid = 1` in `/etc/sysctl.d/99-perf.conf`.

**Bad, `perf: command not found`:**

```bash
sudo apt install linux-tools-common linux-tools-$(uname -r)
```

**Bad, `<not counted>` or `<not supported>` where the numbers should be:** the
hardware counters are not exposed to you at all. This happens inside virtual
machines and some containers, and no setting fixes it from in here.

If `perf` cannot be made to work, **run the benchmark anyway.** The
`instructions`, `cycles` and `task_clock_ms` columns will be empty and
everything else is unaffected. The run prints one line saying so at the start,
and that line is what we need to see — not silence.

### Compiler and CMake

```bash
cmake --version        # 3.12 or newer
c++ --version          # any C++17 compiler; GCC 9+ or Clang 10+ is safe
```

`setup.py` refuses to start on a CMake older than 3.12 and says so. The exact
compiler is recorded in `env.json` for us, but note what the toolchain lock
does and does not pin: it pins the two CGAL commits and the checksums of the
built binaries, **not** your compiler. So do not change compilers or update the
system toolchain partway through a run — rebuilding with a different compiler
mid-sweep would produce a results set measured with two different binaries.
When the lock sees different code, it deletes the previous results and output
meshes and starts that run over, so after a `git pull` you simply re-run the
same command.

### TBB

TBB is the one dependency this branch adds. Without it the parallel arm does
not exist and there is nothing to measure, so the build fails outright rather
than quietly producing a sequential binary.

```bash
ls /usr/lib/x86_64-linux-gnu/cmake/TBB/TBBConfig.cmake   # or anywhere on the system
```

**If it is missing:** `sudo apt install libtbb-dev`, or install oneAPI and pass
the directory holding `TBBConfig.cmake` as `--tbb-dir`.

The build checks afterwards that `bench_remesh` really linked TBB and stops if
it did not, so a silently TBB-less binary cannot reach a measurement.

### Testing a branch that is not pushed yet

`setup.py` normally clones `scale24` from GitHub. To build against a checkout
already on your disk instead — a branch still under review, or a local fix you
want to measure before pushing it — point `--ours-dir` at it and run `setup.py`
by hand before `run_all.py`:

```bash
export LD_LIBRARY_PATH=/opt/intel/oneapi/2025.2/lib/intel64/gcc4.8:$LD_LIBRARY_PATH
python3 scripts/setup.py --root ~/tetra-bench-work \
    --ours-dir /path/to/your/cgal/checkout \
    --tbb-dir /opt/intel/oneapi/2025.2/lib/cmake/tbb \
    --diagnostics --skip-main
```

Three things worth knowing about that command:

- `--ours-dir` takes the **root of a CGAL source tree** (the directory holding
  `CGALConfig.cmake`), not a build directory. Nothing is fetched and nothing in
  it is modified; it is read.
- `--skip-main` omits the upstream reference arm, so you do not also need a
  CGAL `main` checkout. `--profile metrics` never uses that arm. The full sweep
  does, so do not use `--skip-main` for a real sweep.
- `--diagnostics` builds the two instrumented binaries. Leave it off if you are
  not running the metrics mode.

The `LD_LIBRARY_PATH` line matters if your TBB came from oneAPI rather than
from the distribution: without it the binaries link but fail to start, with a
message about `libtbb.so.12` not being found. Export it in the same shell you
then run `run_all.py` from.

`run_all.py` afterwards reuses whatever `setup.py` built; it does not re-clone.

### Memory and disk

The largest run we have measured peaked at **1.8 GB** of resident memory, at 24
threads. Any machine with 8 GB has room. The check is only that nothing else
large is running at the same time:

```bash
free -g
```

Disk depends on which run:

```bash
df -h .
```

- **metrics run, meshes already prepared:** under 1 GB.
- **metrics run from nothing:** the input meshes have to be built first — see
  the note under §2 — so budget as for the full sweep.
- **full sweep:** 10–20 GB. `run_all.py` estimates the figure for your surface
  count and refuses to start if the disk cannot hold it.

### CPU governor (optional but worth it)

```bash
cpupower frequency-set -g performance    # needs root
```

`run_all.py` attempts this itself and carries on if it cannot. Under the
`intel_pstate` driver, `powersave` is the stock default and still reaches full
turbo, so it is not the problem it sounds like. On other drivers it makes
timings drift.

---

## 2. The commands

Substitute **`/path/to/thingi10k/stl`** everywhere below: it is the directory
holding the Thingi10K surface files (`.stl`, `.off` or `.ply`). Nothing else
below needs editing.

### The metrics run (~45 min of measuring) — what we are asking for

```bash
cd tetra-scaling-bench
git pull
nohup python3 run_all.py \
    --surfaces-dir /path/to/thingi10k/stl \
    --root ~/tetra-bench-work \
    --metrics > metrics.log 2>&1 &
```

It builds, prepares inputs if they are not already there, and then measures,
in this order:

- **two thread ladders** (`94665_cdt` f0.3 and `67856_cdt` f0.25) at 1, 2, 4,
  8, 12, 16 and 24 threads, three repetitions each — how the run time falls as
  threads are added, and how much work is executed at each point.
- **two instrumented runs**, one printing a line per stage of the remesher and
  one printing lock-retry counts. These are diagnostics; their stdout is what
  we want, not their timing.
- **a per-thread profile** of both ladder meshes at 24 threads: `perf record`,
  then `scripts/serial_profile.py`, which lists the functions running while
  the other threads had nothing to do. Only the summary is kept, in
  `serialprof/`; the raw recording is deleted. Needs `perf`, and is skipped
  without it.
- **a spatial-sort interval sweep** over 1, 2, 3 and 4 at 24 threads on
  `94665_cdt` f0.3.
- **a lock-grid sweep** over 16, 24, 32, 48 and 64 at 24 threads on the same
  configuration.

Core pinning, the deferred-cutoff sweep and `94665_mesh3` were answered by the
2026-09-23 run and are no longer run by default.

> **If the input meshes are not already prepared**, the preparation stage runs
> first and takes 1–3 hours. It does not need the machine idle, and it only
> happens once — a later full sweep reuses exactly the same meshes.

To run only part of it (for instance after fixing a `perf` problem, to redo
just the ladders):

```bash
python3 scripts/run_bench.py --root ~/tetra-bench-work --profile metrics \
    --metrics-phases ladder \
    --results-dir ~/tetra-bench-work/metrics_results \
    --mesh-out-dir ~/tetra-bench-work/metrics_out_meshes
```

The phase names are `ladder`, `diagnostics`, `serialprof`, `ssort` and `lockgrid` (the default), plus `pinning` and `deferred`.

### The smoke run (~30–45 min)

```bash
nohup python3 run_all.py \
    --surfaces-dir /path/to/thingi10k/stl \
    --root ~/tetra-bench-work \
    --smoke > smoke.log 2>&1 &
```

A few surfaces and a 15-minute measurement, into `smoke_results/`. Worth doing
once on a machine the kit has never run on.

### The full sweep (~14 h)

```bash
nohup python3 run_all.py \
    --surfaces-dir /path/to/thingi10k/stl \
    --root ~/tetra-bench-work > run.log 2>&1 &
```

~15 min building, 1–3 h preparing inputs, 12 h measuring, ~1 h deriving
quality. Only the measuring needs the machine idle.

---

## 3. How long, and when the machine must be idle

| stage | time | machine must be idle |
|---|---|---|
| build | ~15 min (~25 with the metrics run's two extra binaries) | no |
| prepare inputs | 1–3 h, once ever | no |
| metrics: two thread ladders | ~32 min | **yes** |
| metrics: instrumented runs | ~4 min | **yes** |
| metrics: per-thread profile | ~2 min | **yes** |
| metrics: spatial-sort interval sweep | ~5 min | **yes** |
| metrics: lock-grid sweep | ~6 min | **yes** |
| full sweep: measurement | 12 h | **yes** |
| quality pass | 5 min–1 h | no |

So the metrics run needs about **45 minutes** of an otherwise idle machine, and
the stages either side of that do not care what else is running.

These estimates come from the 2026-09-23 basquiat run, so they should be close.
The run stops starting new measurements after about 70 minutes whatever happens.
A run that takes noticeably longer is worth mentioning to us.

---

## 4. What must be true while it measures

During the timed stages — and only during those — the machine must be doing
nothing else:

- No compiling. **A build running during a timed measurement invalidates it**,
  and there is no way to tell afterwards from the numbers alone. This includes
  building anything of your own in another terminal.
- No other benchmarks, no long-running jobs, no container workloads.
- Nothing pinned to specific cores by anything else.
- Leaving an editor or a browser open is fine; they are idle.

An interactive session on the machine is fine as long as you are not running
anything in it.

---

## 5. If it is interrupted

**Re-run the identical command.** Every stage is idempotent: finished builds,
prepared meshes and completed measurement cells are all recognised and skipped,
and the run picks up where it stopped. This is true whether it was interrupted
by Ctrl-C, a reboot, or a crash.

**Do not delete `work/meshes/`.** Preparing the input meshes is the
1–3 hour stage, and it is the one thing that is expensive to lose. A resumed
run reports them as `reused` in seconds.

Safe to delete, if you want to redo one from scratch:

- `smoke_results/` and `smoke_out_meshes/` — the smoke run's output.
- `metrics_results/` and `metrics_out_meshes/` — the metrics run's output.

Never delete `results/` or `out_meshes/` from a full sweep you have not sent
back.

---

## 6. What to send, and what to keep

At the end the run names one file. **Send that one file**, nothing else:

```
~/tetra-bench-work/results_<hostname>_<date>.tar.gz
```

A few hundred KB. It contains the results table, the machine description, the
per-run JSONs, the raw `perf` counter files and the instrumented runs' stdout.

**Keep `~/tetra-bench-work/metrics_out_meshes/` and `out_meshes/` where they
are.** Those are the remeshed output meshes, and they are what lets us compute
a metric nobody thought of today without asking you to re-measure anything.
They are large; they are not something to send.

---

## 7. Before you send it: is this a good run?

Five checks. They take a minute and catch the failures that are invisible in
the log.

Set this once so the commands below fit on a line:

```bash
R=~/tetra-bench-work/metrics_results
```

**1. Everything exited cleanly.** No output at all is the pass.

```bash
awk -F, 'NR==1 {for (i=1; i<=NF; i++) c[$i]=i; next}
         $c["rc"] != 0 || $c["timed_out"] != 0 {print "BAD:", $0}' $R/results.csv
```

Any line printed here is a run that crashed or was killed, and we need to know
about it. Note that `rc` is deliberately the remesher's own exit status and not
`perf`'s: `perf stat` returns 0 even for a program killed by a signal, so
taking its status would let a crash read as a pass.

**2. The work counters are populated.** This should print a number well above
zero, and equal to the number of runs.

```bash
awk -F, 'NR==1 {for (i=1; i<=NF; i++) c[$i]=i; next}
         $c["instructions"] != "" {n++}
         END {print n+0, "of", NR-1, "runs have instruction counts"}' $R/results.csv
```

If it prints `0 of ...`, `perf` was not working — see §1. The run is still
usable but much less informative, so it is worth fixing and re-running the
ladders.

**3. The thread ladder is complete.** Each of the three configurations should
appear at all seven thread counts, three times each — 21 rows per
configuration.

```bash
awk -F, 'NR==1 {for (i=1; i<=NF; i++) c[$i]=i; next}
         $c["variant"] == "" {print $c["mesh"], $c["factor"], $c["threads"]}' \
    $R/results.csv | sort -k1,1 -k2,2 -k3,3n | uniq -c
```

The `variant` column says what made a run special — a pinned core set, a lock
grid value, a diagnostic binary — so an empty one is an ordinary ladder run.
Expect `3` in front of every line.

**4. The instrumented runs printed something.**

```bash
wc -l $R/diagnostics/*.out
```

**Every one of these files must be non-empty.** This matters more than it
looks: the lock-counting build prints its single line from a destructor that
runs as the program exits, so a run that crashed prints *nothing at all*. An
empty file is a failed run wearing the costume of a clean one. The report also
flags this, and `$R/diagnostics/summary.json` records each run's exit status
beside its output for exactly this reason.

**5. The lock-grid sweep really swept.** Five distinct grids, not one grid five
times:

```bash
awk -F, 'NR==1 {for (i=1; i<=NF; i++) c[$i]=i; next}
         $c["lock_grid"] != "" {print $c["lock_grid"]}' $R/results.csv \
    | sort -n | uniq -c
```

Expect five lines, `16 24 32 48 64`. If they all read the same number, the
build ignored `CGAL_TETRAHEDRAL_REMESHING_LOCK_GRID` — you are on a branch that
predates it, and those runs measured one configuration repeatedly. Note that an
unset run is not a distinguishing test: with the variable unset the remesher
derives the count from the thread count, and that rule returns 16 at anything up
to 4 threads anyway. The stage-timing build prints the grid it actually used, so
this settles it directly:

```bash
grep LOCKGRID $R/diagnostics/*topstage.out
```

**6. No timing that is obviously wrong.** Skim the wall-time column for a run
that took far longer than its neighbours — that is usually something else
having been running on the machine at the time. Every row records the time it
started, so a suspect row can be identified and redone rather than the whole
run being thrown away.

If checks 1 to 4 pass, send the tarball.
