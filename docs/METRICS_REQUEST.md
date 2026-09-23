# What to record on the 24-core machine

Written 2026-09-22 against `bd32d445baa`. Everything below is a measurement
request, not a code change. The two code changes are in `~/wt_scale24`.

> **Trimmed 2026-09-23.** The basquiat run of that day settled three items, and
> the kit no longer runs them by default: `94665_mesh3` (§2, §6) scaled like
> `67856_cdt` and took 114 of 161 minutes at its intended size; core pinning
> (§3) showed the scheduler already picks the performance cores; the
> deferred-cutoff sweep moved nothing. What remains runs in about 45 minutes.

## Why these

The returned data already settles the top-level split, on `94665_cdt` at edge
factor 0.3: fitting a serial fraction to the thread ladder gives 24.7% at 2
threads, 23.9% at 4 and 24.5% at 8 — flat, and measured where the run executes
only 7% more CPU than at one thread, so contention cannot be what produces it.
A serial quarter caps 24 threads at 3.59x however perfect the rest becomes.
The observed 2.48x sits below that ceiling, so there are two levers and their
sizes differ: removing ALL contention is worth 2.48x -> 3.59x, and everything
past 3.59x has to come out of the serial quarter.

What the returned data cannot say is which of those two the 12 -> 24 regression
is, because every quantity recorded at 24 threads confounds three things.

## 1. Cycles and instructions, not CPU seconds

`user_s` cannot be read as work. On the 4-core box, the same binary on
`409635_cdt_0.5` spends **+44.6% CPU seconds** at 4 threads over 1 — but
**+17.1% instructions** and **+10.0% cycles**. The rest is the clock falling
3.74 -> 2.85 GHz as more cores light up. At 24 threads that drop will be larger
still, so `avg_parallelism` and `user_s` overstate the waste by roughly
threefold and cannot be corrected after the fact.

Please wrap each timed run in:

    perf stat -e instructions,cycles,task-clock -x, -o <run>.perf <cmd>

and put `instructions`, `cycles` and `task-clock` in `results.csv`. One caution:
`perf stat` returns 0 for a run killed by a signal, so record the *inner*
command's exit status separately or a crashed run reads as a pass.

## 2. The thread ladder, on three configurations and no more

The ladder exists for `94665_cdt` at f0.3 and f0.25 only, one rep per point.
Everything above rests on it. Please run the ladder — 1, 2, 4, 8, 12, 16, 24
threads, 3 reps per point, reps interleaved rather than all reps of a point
together — on exactly these three:

| configuration | wall at 24 threads | why this one |
|---|---:|---|
| `94665_cdt` f0.3 | 27.5 s | every conclusion so far rests on it; repeating it with cycles and instructions is what makes the rest comparable |
| `67856_cdt` f0.25 | 26.1 s | the heavy end, and a different surface, so the serial fraction is not a property of one mesh |
| `94665_mesh3` f0.25 | 8.0 s | the other input pipeline, on the same surface as the first — CDT and Mesh_3 have disagreed on the sign of a comparison before |

## 3. Pin a ladder to the performance cores

The Core Ultra 9 285 is 8 performance cores plus 16 efficiency cores (the
40 MiB L2 in 12 instances is 8 single-core plus 4 clusters of 4). Efficiency
cores are far slower per core, so "24 cores" is not 24x of the core the
1-thread run measured, and part of the 12 -> 24 regression may simply be work
landing on slower cores.

To separate that from contention, please also run, on `94665_cdt` f0.3 at 8
threads, 3 reps of each of unpinned, performance-cores-only and
efficiency-cores-only:

    taskset -c <the 8 P-core ids> ... <8 threads>
    taskset -c <8 E-core ids>     ... <8 threads>

and send `lscpu -e` plus
`cat /sys/devices/system/cpu/cpu*/cpufreq/cpuinfo_max_freq` so the two sets can
be identified. If 8 threads on 8 performance cores beats 8 threads unpinned,
the scheduler placing work on efficiency cores is a real part of the loss.

## 4. Frequency of the busy cores only

`freq_mean_mhz` averages all 24 cores including idle ones, so it rises with
thread count and reads backwards — the handover already flags it. Sampling the
busiest N cores instead makes it the number that separates the clock drop in
§1 from the rest.

## 5. Two instrumented runs at 24 threads

Both binaries below are the same source with one macro added, and both are
plain Release. They are diagnostics: run them once each, not as part of the
timed sweep.

- `-DCGAL_TR_TOPSTAGE` prints one `TOPSTAGE` line per top-level stage of
  `remesh()` with its wall and CPU time: the six elementary operations, the
  serial spatial-sort rebuild inside `split()`, the smoothing's `refresh()`,
  the edge scan of `resolution_reached()`, and `postprocess()`. The stages sum
  to within 0.2% of the total, so nothing hides between them. Run it at 1 and
  at 24 threads: any stage whose wall time does not fall is what the serial
  quarter is made of. At 4 threads on the 4-core box no stage is serial — every
  one scales 1.7x to 2.9x — so the quarter is spread through all of them, and
  whether that is still true at 24 is the open question.

- `-DCGAL_TR_LOCKCOUNT` prints retries per elementary operation. Same
  configuration, same two thread counts. At 4 threads it is about 1.4; what it
  reaches at 24 says whether the lock grid is the binding constraint.

## 6. Sweep the lock grid

`~/wt_scale24` makes the grid size settable at run time through
`CGAL_TETRAHEDRAL_REMESHING_LOCK_GRID`, read once per remesher and never on a
locking path. Please sweep 16, 24, 32, 48, 64 at 24 threads, 3 reps, on
`94665_cdt` f0.3 and `94665_mesh3` f0.25. The built-in default was chosen at
FOUR threads and the new rule extrapolates it; the sweep is what says whether
the extrapolation is right.

## 7. What NOT to run, and why

Roughly 45 minutes of machine time covers everything above. Most of what the
smoke runs measured does not earn its place in it:

- **Edge factors 0.5 and 0.4, on every mesh.** They answer how wall time grows
  with mesh size at a fixed thread count. That is a different question, it is
  already answered, and none of those runs says anything about scaling.
- **`94147_cdt`, `94157_mesh3`, `94662_mesh3`, `509312_mesh3`, in full.** Each
  would add another point to a size trend nobody is asking about. None is
  longer at 24 threads than the three chosen, and the Mesh_3 pipeline is
  already represented by `94665_mesh3`.
- **The `main` and `seq` reference arms.** They are the headline
  speedup-against-CGAL numbers, they are already measured on `67856_cdt` and
  `94665_cdt` at both factors, and they run at one thread by definition — so
  they cannot move with anything being tested here. Two of them cost 9 minutes
  for nothing. Run none.
- **Quality reports on every configuration.** Parallel against sequential is
  0.25% on cell count and the slivers agree to 0.04 percentage points; that
  question is settled. One quality report at each end of one ladder is enough
  to confirm nothing regressed.

The budget, if it helps in planning: about 27 minutes for the three ladders,
9 for the lock-grid sweep, 5 for the core pinning and 3 for the two
instrumented runs.
