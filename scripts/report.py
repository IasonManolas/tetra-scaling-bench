#!/usr/bin/env python3
"""Stage 3: turn the returned raw data into a report.

Runs on whichever machine is doing the analysis, not on the benchmark machine.
Reads results.csv, the per-run JSONs and the per-mesh quality JSONs, joins them
into one tidy table, and writes all_runs.csv, report.md and (if matplotlib is
available) a few PNGs.

    python3 scripts/report.py --results work/results --out work/report

Adding a metric later does not mean re-running the benchmark: re-run
`run_bench.py --quality-pass --force-quality` over the archived output meshes,
then re-run this.
"""
import argparse
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ARM_ORDER = ["main", "seq", "par"]


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------

def load_rows(results_dir):
    csv_path = results_dir / "results.csv"
    if not csv_path.exists():
        sys.exit("No %s" % csv_path)

    qdir = results_dir / "quality"
    rows = []
    with open(csv_path, newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                r["threads"] = int(r["threads"])
                r["rep"] = int(r["rep"])
                r["wall_s"] = float(r["wall_s"])
                r["factor"] = float(r["factor"])
                r["rc"] = int(r["rc"])
                r["timed_out"] = int(r["timed_out"])
                r["peak_rss_kb"] = int(r["peak_rss_kb"]) if r["peak_rss_kb"] else None
                r["remesh_s"] = float(r["remesh_s"]) if r["remesh_s"] else None
                for k in ("avg_parallelism", "freq_mean_mhz", "freq_min_mhz",
                          "freq_max_mhz", "user_s", "sys_s",
                          "busy_freq_mean_mhz", "busy_freq_min_mhz",
                          "busy_freq_max_mhz",
                          "instructions", "cycles", "task_clock_ms"):
                    r[k] = float(r[k]) if r.get(k) else None
                # What made this row different from a plain timed run: a taskset
                # CPU set, a lock-grid value, a diagnostic binary. Rows with a
                # variant are NOT part of the thread ladder and must not be
                # aggregated into it.
                r["variant"] = (r.get("variant") or "").strip()
                r["pin"] = (r.get("pin") or "").strip()
                r["lock_grid"] = int(r["lock_grid"]) if r.get("lock_grid") else None
                r["throttle_events"] = (int(r["throttle_events"])
                                        if r.get("throttle_events") else None)
            except Exception:
                continue

            # from the run JSON: the input's average edge length and the derived
            # target, which the output mesh alone cannot tell us
            rj = Path(r["json"])
            if rj.exists():
                try:
                    d = json.loads(rj.read_text())
                    info = d.get("run_metadata", {}).get("run_info", {})
                    r["target_edge_length"] = info.get("Edge Length")
                    r["avg_edge_length_in"] = info.get("Avg_edge_length")
                    r["cells_in"] = info.get("Mesh", {}).get("Cells")
                    r["feature_edge_mean"] = (d.get("metrics", {}).get("Quality", {})
                                               .get("Edge_Length", {}).get("Mean"))
                except Exception:
                    pass

            # from the quality JSON, keyed by the output mesh's stem
            qstem = "%s_f%s_%s_t%d" % (r["mesh"], r["factor"], r["arm"], r["threads"])
            qp = qdir / (qstem + ".json")
            if qp.exists():
                try:
                    q = json.loads(qp.read_text())["metrics"]["Quality"]
                    r["cells_out"] = q.get("Cell_Count", {}).get("Value")
                    r["min_dihedral"] = q.get("Dihedral_Angle", {}).get("Minimum")
                    r["mean_dihedral"] = q.get("Dihedral_Angle", {}).get("Mean")
                    r["median_min_dihedral"] = q.get("Min_Dihedral_Angle", {}).get("Median")
                    r["sliver_5deg"] = q.get("Sliver_Fraction", {}).get("Below_5_Deg")
                    r["sliver_10deg"] = q.get("Sliver_Fraction", {}).get("Below_10_Deg")
                    r["min_cell_volume"] = q.get("Cell_Volume", {}).get("Minimum")
                    r["radius_ratio_min"] = q.get("Radius_Radius_Ratio", {}).get("Minimum")
                    cel = q.get("Cell_Edge_Length", {})
                    r["cell_edge_mean"] = cel.get("Mean")
                    r["cell_edge_median"] = cel.get("Median")
                    r["cell_edge_std"] = cel.get("Std_Dev")
                except Exception:
                    pass

            if r.get("cell_edge_mean") and r.get("target_edge_length"):
                r["edge_conformance"] = r["cell_edge_mean"] / r["target_edge_length"]

            # prepare_meshes.py names every input <thingi id>_<pipeline>.mesh, so
            # the generator that produced the input travels with the row.
            name = r["mesh"]
            if "_" in name:
                base, _, suffix = name.rpartition("_")
                if suffix in ("cdt", "mesh3"):
                    r["thingi_id"], r["pipeline"] = base, suffix
            r.setdefault("thingi_id", name)
            r.setdefault("pipeline", "unknown")

            rows.append(r)
    return rows


def key(r):
    return (r["mesh"], r["factor"], r["arm"], r["threads"])


def plain(rows):
    """The ordinary timed runs: everything that is not a variant."""
    return [r for r in rows if not r.get("variant")]


def _med(rs, field):
    """Median of a field over the repeats that have it, or None."""
    vals = [r[field] for r in rs if r.get(field) is not None]
    return statistics.median(vals) if vals else None


def aggregate(rows):
    """Median wall time and CV per (mesh, factor, arm, threads)."""
    groups = defaultdict(list)
    for r in plain(rows):
        if r["rc"] == 0 and not r["timed_out"]:
            groups[key(r)].append(r)

    agg = {}
    for k, rs in groups.items():
        walls = [r["wall_s"] for r in rs]
        mean = statistics.fmean(walls)
        cv = (statistics.pstdev(walls) / mean) if len(walls) > 1 and mean else None
        best = min(rs, key=lambda r: r["wall_s"])
        agg[k] = {
            "n": len(walls),
            "median_wall": statistics.median(walls),
            "min_wall": min(walls),
            "cv": cv,
            "peak_rss_kb": max((r["peak_rss_kb"] or 0) for r in rs) or None,
            "quality": {f: best.get(f) for f in (
                "cells_out", "min_dihedral", "mean_dihedral",
                "median_min_dihedral", "sliver_5deg", "sliver_10deg",
                "min_cell_volume", "radius_ratio_min", "cell_edge_mean",
                "cell_edge_median", "edge_conformance")},
            "cells_in": best.get("cells_in"),
            "target_edge_length": best.get("target_edge_length"),
            "avg_parallelism": best.get("avg_parallelism"),
            "freq_mean_mhz": best.get("freq_mean_mhz"),
            "freq_min_mhz": best.get("freq_min_mhz"),
            "busy_freq_mean_mhz": _med(rs, "busy_freq_mean_mhz"),
            "busy_freq_min_mhz": _med(rs, "busy_freq_min_mhz"),
            "instructions": _med(rs, "instructions"),
            "cycles": _med(rs, "cycles"),
            "task_clock_ms": _med(rs, "task_clock_ms"),
            "throttle_events": max((r.get("throttle_events") or 0) for r in rs) or None,
        }
    return agg


# --------------------------------------------------------------------------
# report sections
# --------------------------------------------------------------------------

def fmt(v, spec="%.3f"):
    if v is None:
        return "-"
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return "-"
    try:
        return spec % v
    except Exception:
        return str(v)


def section_env(results_dir, out):
    ep = results_dir / "env.json"
    if not ep.exists():
        return
    env = json.loads(ep.read_text())
    out.append("## Machine\n")
    out.append("| | |")
    out.append("|---|---|")
    out.append("| host | `%s` |" % env.get("hostname"))
    out.append("| nproc | %s |" % env.get("nproc"))
    out.append("| governor | %s (driver %s) |" % (env.get("cpu_governor") or "unknown",
                                                  env.get("cpufreq_driver") or "unknown"))
    if env.get("cpufreq_max_khz"):
        out.append("| clock range | %.0f - %.0f MHz |"
                   % (float(env["cpufreq_min_khz"]) / 1000,
                      float(env["cpufreq_max_khz"]) / 1000))
    out.append("| platform | %s |" % env.get("platform"))
    tc = env.get("toolchain", {})
    out.append("| ours | `%s` @ `%s` |" % (tc.get("ours_ref"), (tc.get("ours_sha") or "")[:12]))
    out.append("| main | `%s` @ `%s` |" % (tc.get("main_ref"), (tc.get("main_sha") or "")[:12]))
    out.append("")
    # Only flag the governor when it can actually matter. Under intel_pstate,
    # 'powersave' is the stock default and still reaches full turbo.
    gov, drv = env.get("cpu_governor"), env.get("cpufreq_driver")
    if gov and gov != "performance" and drv != "intel_pstate":
        out.append("> The CPU governor was `%s` with driver `%s`. Low thread counts "
                   "look worse than they are under a scaling governor, so the speedup "
                   "numbers below are, if anything, optimistic.\n" % (gov, drv or "?"))
    if env.get("no_turbo") == "1":
        out.append("> **Turbo was disabled** on this machine, so every thread count "
                   "ran at base clock.\n")


def section_scaling(agg, out, threads_max):
    out.append("## Scaling\n")
    configs = sorted({(m, f) for (m, f, a, t) in agg}, key=lambda c: str(c))
    for (mesh, factor) in configs:
        par_t = sorted(t for (m, f, a, t) in agg if m == mesh and f == factor and a == "par")
        if not par_t:
            continue
        base1 = agg.get((mesh, factor, "par", 1), {}).get("median_wall")
        seq = agg.get((mesh, factor, "seq", 1), {}).get("median_wall")
        mainw = agg.get((mesh, factor, "main", 1), {}).get("median_wall")
        cells_in = agg[(mesh, factor, "par", par_t[0])].get("cells_in")

        out.append("### %s, edge factor %s (input %s cells)\n" % (mesh, factor, cells_in))
        out.append("| threads | wall (s) | n | CV | vs par@1 | vs seq | vs main | eff. | peak RSS (MB) |")
        out.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for t in par_t:
            a = agg[(mesh, factor, "par", t)]
            w = a["median_wall"]
            sp1 = base1 / w if base1 and w else None
            out.append("| %d | %s | %d | %s | %s | %s | %s | %s | %s |" % (
                t, fmt(w, "%.1f"), a["n"],
                fmt(a["cv"] and a["cv"] * 100, "%.2f%%"),
                fmt(sp1, "%.2fx"),
                fmt(seq / w if seq and w else None, "%.2fx"),
                fmt(mainw / w if mainw and w else None, "%.2fx"),
                fmt(sp1 / t if sp1 else None, "%.2f"),
                fmt(a["peak_rss_kb"] and a["peak_rss_kb"] / 1024, "%.0f")))
        # More threads taking longer is physically possible past the scaling
        # knee, but more threads being FASTER at a lower count than a higher one
        # by a wide margin usually means drift between measurement epochs, not a
        # property of the code. Say so rather than letting it read as a result.
        walls = [(t, agg[(mesh, factor, "par", t)]["median_wall"]) for t in par_t]
        inversions = [(t1, w1, t2, w2)
                      for (t1, w1), (t2, w2) in zip(walls, walls[1:])
                      if w2 > w1 * 1.05]
        if inversions:
            t1, w1, t2, w2 = inversions[0]
            out.append("")
            out.append("> **Non-monotonic**: %d threads (%.1fs) is slower than %d "
                       "(%.1fs). Past the scaling knee this can be real; a large "
                       "gap more often means the two were measured at different "
                       "times on a drifting machine. Check `started_at` in "
                       "`all_runs.csv` and the CV column before quoting these."
                       % (t2, w2, t1, w1))

        for arm, w in (("seq", seq), ("main", mainw)):
            if w:
                a = agg[(mesh, factor, arm, 1)]
                out.append("| %s@1 | %s | %d | %s | - | - | - | - | %s |" % (
                    arm, fmt(w, "%.1f"), a["n"],
                    fmt(a["cv"] and a["cv"] * 100, "%.2f%%"),
                    fmt(a["peak_rss_kb"] and a["peak_rss_kb"] / 1024, "%.0f")))
        out.append("")


def section_totals(agg, out, threads_max):
    """Total time over the whole set.

    This project's SETS.md treats the TOTAL over the set as the reading, not the
    median of per-config ratios -- on a small number of configs the median is
    weak and has misread changes before. Both are printed; the total is the one
    to quote.
    """
    out.append("## Set totals\n")
    all_configs = sorted({(m, f) for (m, f, a, t) in agg})
    arms = [("main", 1), ("seq", 1)] + [("par", t) for t in
                                        sorted({t for (m, f, a, t) in agg if a == "par"})]

    # Totals are only comparable over a COMMON set of configs. Summing each arm
    # over whatever it happens to cover and putting the results in one column
    # invites exactly the wrong comparison, so restrict to the configs every
    # listed arm measured and say which those are.
    common = [c for c in all_configs
              if all((c[0], c[1], arm, t) in agg for arm, t in arms)]
    if not common:
        out.append("No single config was measured by every arm, so there is no "
                   "comparable total. Per-config numbers above are the reading.\n")
        # Report which arm is the bottleneck, so the gap is actionable.
        for arm, t in arms:
            n = sum(1 for c in all_configs if (c[0], c[1], arm, t) in agg)
            out.append("- `%s`: %d/%d configs" %
                       (arm if arm != "par" else "par@%d" % t, n, len(all_configs)))
        out.append("")
        return

    out.append("Over the %d config(s) measured by every arm: %s\n" %
               (len(common), ", ".join("%s f=%s" % c for c in common)))
    out.append("| arm | total wall (s) | vs main | vs seq |")
    out.append("|---|---:|---:|---:|")
    totals = {}
    for arm, t in arms:
        totals[(arm, t)] = sum(agg[(c[0], c[1], arm, t)]["median_wall"] for c in common)
    base_main = totals.get(("main", 1))
    base_seq = totals.get(("seq", 1))
    for arm, t in arms:
        tot = totals[(arm, t)]
        label = arm if arm != "par" else "par@%d" % t
        out.append("| %s | %s | %s | %s |" % (
            label, fmt(tot, "%.1f"),
            fmt(base_main / tot if base_main and tot else None, "%.2fx"),
            fmt(base_seq / tot if base_seq and tot else None, "%.2fx")))
    out.append("")
    if len(common) < len(all_configs):
        out.append("_%d of %d configs are excluded from this table because not "
                   "every arm covers them; they still appear per-config above._\n"
                   % (len(all_configs) - len(common), len(all_configs)))


def section_machine_behaviour(agg, out):
    """Why the speedup is what it is: the code, or the machine?

    A speedup well below the thread count has three very different causes, and
    wall time alone cannot tell them apart:

      - avg_parallelism near the thread count, clocks steady
            -> the cores were used and held their speed. The shortfall is the
               algorithm. This is the finding the benchmark exists to produce.
      - avg_parallelism far below the thread count
            -> the work never went parallel: serial sections, lock contention,
               or threads starved. Not a clock problem.
      - clocks falling as threads rise, or throttle_events rising
            -> the machine could not sustain the frequency. All-core turbo is
               always below single-core turbo, so SOME drop is normal and
               expected; throttle events are not.
    """
    rows = [(k, a) for k, a in agg.items() if k[2] == "par"
            and a.get("avg_parallelism")]
    if not rows:
        return

    out.append("## Cores actually used, and clocks held\n")
    out.append("`used` is CPU time over wall time — the average number of cores "
               "genuinely working. `eff` is `used` over the threads requested: "
               "high `eff` with low speedup means the cores were busy but not "
               "productive, low `eff` means they were never taken. Some clock "
               "drop as threads rise is normal, since all-core turbo sits below "
               "single-core turbo. `throttle` counts throttle entries during the "
               "run — the magnitude is not a severity, but anything above 0 "
               "means the machine was cutting clocks.\n")
    out.append("`busy clock` is the mean clock of the N cores that were actually "
               "working, N being the run's thread count. `all-core clock` is the "
               "mean over every core, idle ones included — it RISES with the "
               "thread count as idle cores stop dragging it down, which reads "
               "backwards, and it is shown only so old results sets can be "
               "compared. Quote the busy column.\n")
    out.append("| config | threads | wall (s) | used | eff | busy clock (MHz) | "
               "busy min | all-core clock (MHz) | throttle |")
    out.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")

    for k, a in sorted(rows, key=lambda kv: (kv[0][0], kv[0][1], kv[0][3])):
        mesh, factor, _, t = k
        used = a["avg_parallelism"]
        out.append("| %s f=%s | %d | %s | %s | %s | %s | %s | %s | %s |" % (
            mesh, factor, t, fmt(a["median_wall"], "%.1f"), fmt(used, "%.1f"),
            fmt(100.0 * used / t if t else None, "%.0f%%") if used else "-",
            fmt(a.get("busy_freq_mean_mhz"), "%.0f"),
            fmt(a.get("busy_freq_min_mhz"), "%.0f"),
            fmt(a.get("freq_mean_mhz"), "%.0f"),
            a.get("throttle_events") if a.get("throttle_events") is not None else "-"))
    out.append("")

    # Say the conclusion rather than leaving it in the table.
    worst = None
    for k, a in rows:
        t = k[3]
        if t and t > 1:
            e = a["avg_parallelism"] / t
            if worst is None or e < worst[0]:
                worst = (e, k, a)
    throttled = [k for k, a in rows if (a.get("throttle_events") or 0) > 0]

    if throttled:
        out.append("> **The machine throttled** during %d configuration(s), so "
                   "some of the shortfall is thermal or power limiting rather "
                   "than the code. Those rows cannot be read as scaling results.\n"
                   % len(throttled))
    if worst and worst[0] < 0.6:
        e, k, a = worst
        out.append("> At %d threads on `%s f=%s`, only **%.1f of %d cores** were "
                   "working on average (%.0f%%). The cores were available and "
                   "were not taken, so this is the algorithm — serial sections, "
                   "contention or starvation — not clock behaviour.\n"
                   % (k[3], k[0], k[1], a["avg_parallelism"], k[3], 100 * e))
    elif worst and worst[0] >= 0.8 and not throttled:
        out.append("> Cores were used well (worst case %.0f%% of the threads "
                   "requested) and no throttling was recorded, so the speedup "
                   "shortfall is real algorithmic scaling, not the machine.\n"
                   % (100 * worst[0]))




def amdahl_serial_fraction(speedup, threads):
    """The serial fraction f that would produce exactly this speedup.

    Amdahl's law says S = 1 / (f + (1-f)/p). Solved for f at a measured S and p.
    Fitting it separately at every thread count is the useful part: a fraction
    that stays flat as p rises is a genuine serial section, while one that
    climbs is a cost that grows with the thread count -- contention -- wearing a
    serial section's clothes.
    """
    if not speedup or threads is None or threads <= 1 or speedup <= 0:
        return None
    f = (1.0 / speedup - 1.0 / threads) / (1.0 - 1.0 / threads)
    return f


def implied_ghz(cycles, task_clock_ms):
    """The clock the run actually held: cycles over CPU time on the cores."""
    if not cycles or not task_clock_ms:
        return None
    return cycles / (task_clock_ms / 1000.0) / 1e9


def section_thread_ladder(agg, out):
    """Speedup, the serial fraction behind it, and the work that was executed.

    CPU seconds are not reported here, and deliberately. The same binary on the
    same input measured +44.6% CPU seconds from one thread to four but only
    +17.1% instructions and +10.0% cycles -- the rest is the clock falling as
    more cores light up. Instructions and cycles are the work; `clock` is the
    part CPU seconds would have blamed on the code.
    """
    usable = [k for k in agg if k[2] == "par"]
    if not usable:
        return
    out.append("## Thread ladder: speedup, serial fraction and work\n")
    out.append("`f` is the serial fraction that fits that one point under "
               "Amdahl's law, and `ceiling` is the speedup that fraction would "
               "allow at the widest thread count measured. `instr` and `cycles` "
               "are relative to the same configuration at one thread. `clock` is "
               "cycles over task-clock — the frequency the working cores held.\n")

    configs = sorted({(m, f) for (m, f, a, t) in usable})
    for (mesh, factor) in configs:
        ts = sorted(t for (m, f, a, t) in usable if m == mesh and f == factor)
        base = agg.get((mesh, factor, "par", 1))
        if len(ts) < 2:
            continue
        widest = ts[-1]
        out.append("### %s, edge factor %s\n" % (mesh, factor))
        if not base:
            out.append("_No 1-thread run, so speedup and the work ratios cannot "
                       "be computed for this configuration._\n")
            continue
        out.append("| threads | wall (s) | n | speedup | f | ceiling at %d | "
                   "instr vs 1t | cycles vs 1t | clock (GHz) |" % widest)
        out.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for t in ts:
            a = agg[(mesh, factor, "par", t)]
            sp = base["median_wall"] / a["median_wall"] if a["median_wall"] else None
            f = amdahl_serial_fraction(sp, t)
            ceiling = (1.0 / (f + (1 - f) / widest)) if f not in (None, 1.0) else None
            ir = (a["instructions"] / base["instructions"]
                  if a.get("instructions") and base.get("instructions") else None)
            cr = (a["cycles"] / base["cycles"]
                  if a.get("cycles") and base.get("cycles") else None)
            out.append("| %d | %s | %d | %s | %s | %s | %s | %s | %s |" % (
                t, fmt(a["median_wall"], "%.1f"), a["n"],
                fmt(sp, "%.2fx"),
                fmt(f and f * 100, "%.1f%%"),
                fmt(ceiling, "%.2fx"),
                fmt(ir and (ir - 1) * 100, "%+.1f%%"),
                fmt(cr and (cr - 1) * 100, "%+.1f%%"),
                fmt(implied_ghz(a.get("cycles"), a.get("task_clock_ms")), "%.2f")))
        out.append("")

        if not base.get("instructions"):
            out.append("> No `instructions`/`cycles` for this configuration, so "
                       "the work columns are empty. Either `perf` was missing on "
                       "the benchmark machine or `kernel.perf_event_paranoid` "
                       "refused it; the run log says which.\n")
            continue
        fits = [(t, amdahl_serial_fraction(
            base["median_wall"] / agg[(mesh, factor, "par", t)]["median_wall"], t))
            for t in ts if t > 1
            and agg[(mesh, factor, "par", t)]["median_wall"]]
        fits = [(t, f) for t, f in fits if f is not None]
        if len(fits) >= 3:
            low = [f for t, f in fits if t <= 8]
            high = [f for t, f in fits if t > 8]
            if low and high and statistics.fmean(high) > statistics.fmean(low) * 1.2:
                out.append("> The fitted serial fraction is %.1f%% up to 8 threads "
                           "and %.1f%% above it. A fraction that only climbs at the "
                           "wide end is a cost that grows with the thread count, "
                           "not a serial section — compare the `cycles vs 1t` "
                           "column at the same rows.\n"
                           % (100 * statistics.fmean(low), 100 * statistics.fmean(high)))


def section_variants(rows, out):
    """Core pinning and the lock-grid sweep.

    Both are one configuration measured several ways at a fixed thread count, so
    the comparison is within the block and never against the ladder.
    """
    var = [r for r in rows if r.get("variant") and r["rc"] == 0
           and not r["timed_out"]]
    if not var:
        return

    # --- core pinning ---------------------------------------------------
    pin_rows = [r for r in var if r["variant"].startswith(("unpinned", "pcore",
                                                           "ecore"))]
    if pin_rows:
        out.append("## Core pinning\n")
        out.append("The same configuration at the same thread count, placed three "
                   "ways. On a hybrid CPU the cores are not interchangeable, so a "
                   "pinned run beating the unpinned one means the scheduler's "
                   "placement is part of the loss — a different problem from lock "
                   "contention, with a different fix. `pin` is the CPU list "
                   "`taskset` was given.\n")
        out.append("| config | threads | placement | pin | n | wall (s) | vs unpinned | "
                   "instr | cycles | clock (GHz) |")
        out.append("|---|---:|---|---|---:|---:|---:|---:|---:|---:|")
        groups = defaultdict(list)
        for r in pin_rows:
            groups[(r["mesh"], r["factor"], r["threads"], r["variant"])].append(r)
        base = {}
        for k, rs in groups.items():
            if k[3] == "unpinned":
                base[k[:3]] = statistics.median(r["wall_s"] for r in rs)
        for k in sorted(groups, key=lambda k: (k[0], k[1], k[2], k[3])):
            rs = groups[k]
            w = statistics.median(r["wall_s"] for r in rs)
            b = base.get(k[:3])
            out.append("| %s f=%s | %d | %s | `%s` | %d | %s | %s | %s | %s | %s |" % (
                k[0], k[1], k[2], k[3], rs[0].get("pin") or "-", len(rs),
                fmt(w, "%.1f"), fmt(b / w if b and w else None, "%.2fx"),
                fmt(_gmed(rs, "instructions"), "%.3g"),
                fmt(_gmed(rs, "cycles"), "%.3g"),
                fmt(implied_ghz(_gmed(rs, "cycles"), _gmed(rs, "task_clock_ms")),
                    "%.2f")))
        out.append("")

    # --- lock grid ------------------------------------------------------
    lg_rows = [r for r in var if r.get("lock_grid")]
    if lg_rows:
        out.append("## Lock grid\n")
        out.append("`CGAL_TETRAHEDRAL_REMESHING_LOCK_GRID` is the number of lock "
                   "cells per axis, read once per remesher. The built-in default "
                   "was chosen at four threads, so this says whether it still "
                   "holds at the widest.\n")
        out.append("| config | threads | grid | n | wall (s) | vs best | instr | "
                   "cycles |")
        out.append("|---|---:|---:|---:|---:|---:|---:|---:|")
        groups = defaultdict(list)
        for r in lg_rows:
            groups[(r["mesh"], r["factor"], r["threads"], r["lock_grid"])].append(r)
        walls = {}
        for k, rs in groups.items():
            walls[k] = statistics.median(r["wall_s"] for r in rs)
        best = {}
        for k, w in walls.items():
            cur = best.get(k[:3])
            if cur is None or w < cur[1]:
                best[k[:3]] = (k[3], w)
        for k in sorted(groups):
            rs, w = groups[k], walls[k]
            bgrid, bw = best[k[:3]]
            out.append("| %s f=%s | %d | %d | %d | %s | %s | %s | %s |" % (
                k[0], k[1], k[2], k[3], len(rs), fmt(w, "%.1f"),
                fmt(bw / w, "%.2fx"),
                fmt(_gmed(rs, "instructions"), "%.3g"),
                fmt(_gmed(rs, "cycles"), "%.3g")))
        out.append("")
        for cfg, (bgrid, bw) in sorted(best.items()):
            out.append("- `%s f=%s` at %d threads is fastest at grid **%d** (%.1f s)."
                       % (cfg[0], cfg[1], cfg[2], bgrid, bw))
        out.append("")


def _gmed(rs, field):
    vals = [r[field] for r in rs if r.get(field) is not None]
    return statistics.median(vals) if vals else None


def section_diagnostics(results_dir, out):
    """What the two instrumented binaries printed, and whether they printed.

    CGAL_TR_LOCKCOUNT reports from a static destructor at exit, so a run that
    died prints nothing at all. An empty capture is a failed run, which is why
    the exit status is recorded beside it and checked here.
    """
    ddir = results_dir / "diagnostics"
    spath = ddir / "summary.json"
    if not spath.exists():
        return
    try:
        summary = json.loads(spath.read_text())
    except Exception:
        return
    out.append("## Instrumented runs\n")
    out.append("| run | binary | threads | exit | lines of output |")
    out.append("|---|---|---:|---:|---:|")
    bad = []
    for stem in sorted(summary):
        d = summary[stem]
        out.append("| `%s` | `%s` | %s | %s | %s |" % (
            stem, d.get("binary"), d.get("threads"), d.get("rc"),
            d.get("stdout_lines")))
        if d.get("rc") != 0 or not d.get("stdout_lines"):
            bad.append(stem)
    out.append("")
    if bad:
        out.append("> **%d instrumented run(s) produced nothing usable**: %s. "
                   "`CGAL_TR_LOCKCOUNT` prints its line from a static destructor "
                   "at exit, so a crash prints nothing — an empty capture is a "
                   "failure, not a clean run.\n" % (len(bad), ", ".join(bad)))
    out.append("The captures themselves are in `%s`, verbatim.\n" % ddir.name)


def load_manifest(results_dir, explicit=None):
    """The input manifest, if it can be found, keyed by .mesh stem.

    It carries whether each Mesh_3 input actually hit its CDT's cell count,
    which decides whether a CDT-vs-Mesh_3 row is a controlled comparison.
    """
    cands = [Path(explicit)] if explicit else []
    cands += [results_dir.parent / "meshes" / "manifest.json",
              results_dir / "manifest.json"]
    for cand in cands:
        if cand.exists():
            try:
                return json.loads(cand.read_text()).get("meshes", {})
            except Exception:
                pass
    return {}


def section_pipelines(agg, rows, out, manifest):
    """The ours-vs-main comparison, computed separately per input generator.

    These two have disagreed on the SIGN of that comparison before (0.95 on one
    set, 1.20 on the other), which is why both are measured. Reporting only the
    pooled number would average that disagreement away and hide the one thing
    this split exists to show.
    """
    pipe_of = {}
    for r in rows:
        pipe_of[r["mesh"]] = r.get("pipeline", "unknown")
    pipes = sorted({p for p in pipe_of.values() if p != "unknown"})
    if len(pipes) < 2:
        return

    # A Mesh_3 input that missed its CDT's cell count differs in SIZE as well as
    # in element quality, so including it turns a controlled comparison into a
    # confounded one. Drop those here and say how many were dropped.
    unmatched = set()
    for stem, rec in manifest.items():
        if rec.get("pipeline") == "mesh3" and rec.get("target_cells") \
                and not rec.get("matched"):
            unmatched.add(stem)
            unmatched.add(stem.replace("_mesh3", "_cdt"))   # drop its partner too

    out.append("## CDT vs Mesh_3 inputs\n")
    out.append("The same surfaces, tetrahedralized two ways, with Mesh_3 sized to "
               "match each CDT's cell count. Each row is a total over the configs "
               "of that pipeline measured by every arm in the row.\n")
    if unmatched and manifest:
        out.append("_%d input(s) are excluded because the Mesh_3 build missed its "
                   "CDT's cell count, so the pair differs in size as well as "
                   "quality: %s._\n"
                   % (len(unmatched), ", ".join(sorted(unmatched))))
    out.append("| pipeline | configs | main (s) | seq (s) | seq vs main | par@max (s) | par vs main |")
    out.append("|---|---:|---:|---:|---:|---:|---:|")

    maxt = max((t for (m, f, a, t) in agg if a == "par"), default=1)
    disagree = []
    for p in pipes:
        cfgs = sorted({(m, f) for (m, f, a, t) in agg
                       if pipe_of.get(m) == p and m not in unmatched})
        usable = [c for c in cfgs
                  if (c[0], c[1], "main", 1) in agg
                  and (c[0], c[1], "seq", 1) in agg]
        if not usable:
            out.append("| %s | 0 | - | - | - | - | - |" % p)
            continue
        tm = sum(agg[(c[0], c[1], "main", 1)]["median_wall"] for c in usable)
        ts = sum(agg[(c[0], c[1], "seq", 1)]["median_wall"] for c in usable)
        with_par = [c for c in usable if (c[0], c[1], "par", maxt) in agg]
        tp = sum(agg[(c[0], c[1], "par", maxt)]["median_wall"] for c in with_par) \
            if len(with_par) == len(usable) else None
        ratio_seq = tm / ts if ts else None
        disagree.append((p, ratio_seq))
        out.append("| %s | %d | %s | %s | %s | %s | %s |" % (
            p, len(usable), fmt(tm, "%.1f"), fmt(ts, "%.1f"),
            fmt(ratio_seq, "%.3fx"),
            fmt(tp, "%.1f") if tp else "-",
            fmt(tm / tp if tp else None, "%.2fx")))
    out.append("")

    known = [(p, r) for p, r in disagree if r]
    if len(known) >= 2:
        lo, hi = min(r for _, r in known), max(r for _, r in known)
        if (lo - 1.0) * (hi - 1.0) < 0:
            out.append("> **The two pipelines disagree on the sign.** `seq` is "
                       "%s on one generator and %s on the other (%s). The pooled "
                       "number is not a meaningful summary here; quote the two "
                       "separately.\n" % (
                           "faster than main" if hi > 1 else "slower than main",
                           "slower" if lo < 1 else "faster",
                           ", ".join("%s %.3fx" % (p, r) for p, r in known)))
        else:
            out.append("> Both pipelines agree in direction (%s), so the pooled "
                       "comparison is safe to quote.\n"
                       % ", ".join("%s %.3fx" % (p, r) for p, r in known))


def section_quality(agg, out):
    out.append("## Quality\n")
    out.append("All of this is derived from the saved output meshes, not from the "
               "timed run. `edge conformance` is the mean tetrahedron edge length "
               "divided by the target edge length: 1.00 is on target.\n")
    out.append("| config | arm | cells out | min dihedral | median min dihedral | "
               "slivers <5deg | slivers <10deg | edge conformance |")
    out.append("|---|---|---:|---:|---:|---:|---:|---:|")
    for k in sorted(agg, key=lambda k: (k[0], k[1], ARM_ORDER.index(k[2]) if k[2] in ARM_ORDER else 9, k[3])):
        mesh, factor, arm, t = k
        q = agg[k]["quality"]
        if q.get("cells_out") is None:
            continue
        label = arm if arm != "par" else "par@%d" % t
        out.append("| %s f=%s | %s | %s | %s | %s | %s | %s | %s |" % (
            mesh, factor, label, fmt(q["cells_out"], "%.0f"),
            fmt(q["min_dihedral"], "%.4g"), fmt(q["median_min_dihedral"], "%.2f"),
            fmt(q["sliver_5deg"], "%.4f"), fmt(q["sliver_10deg"], "%.4f"),
            fmt(q["edge_conformance"], "%.3f")))
    out.append("")


def section_validity(rows, agg, out):
    out.append("## Validity and determinism\n")

    bad = [r for r in rows if r["rc"] != 0 or r["timed_out"]]
    if bad:
        out.append("**%d runs did not complete cleanly:**\n" % len(bad))
        for r in bad[:25]:
            out.append("- `%s` f=%s %s@%d rep %d: rc=%d%s" % (
                r["mesh"], r["factor"], r["arm"], r["threads"], r["rep"],
                r["rc"], " (timeout)" if r["timed_out"] else ""))
        out.append("")
    else:
        out.append("Every run exited 0 and none timed out.\n")

    degenerate = [k for k, a in agg.items()
                  if a["quality"].get("min_dihedral") is not None
                  and a["quality"]["min_dihedral"] <= 0.0]
    if degenerate:
        out.append("**Degenerate output (min dihedral <= 0):**\n")
        for k in degenerate:
            out.append("- %s f=%s %s@%d" % k)
        out.append("")

    # Cross-thread agreement. The parallel arm is not expected to be
    # bit-identical across thread counts -- scheduling changes the order of
    # equally valid operations -- so the question is whether the spread is small
    # and stable, not whether it is zero. The seq arm is the reference.
    out.append("### Cell count across thread counts\n")
    out.append("| config | seq | par min | par max | spread | vs seq |")
    out.append("|---|---:|---:|---:|---:|---:|")
    configs = sorted({(m, f) for (m, f, a, t) in agg})
    for (mesh, factor) in configs:
        pars = [agg[k]["quality"]["cells_out"] for k in agg
                if k[0] == mesh and k[1] == factor and k[2] == "par"
                and agg[k]["quality"].get("cells_out")]
        if not pars:
            continue
        seqc = agg.get((mesh, factor, "seq", 1), {}).get("quality", {}).get("cells_out")
        lo, hi = min(pars), max(pars)
        spread = (hi - lo) / lo if lo else None
        out.append("| %s f=%s | %s | %s | %s | %s | %s |" % (
            mesh, factor, fmt(seqc, "%.0f"), fmt(lo, "%.0f"), fmt(hi, "%.0f"),
            fmt(spread and spread * 100, "%.2f%%"),
            fmt(seqc and (hi - seqc) / seqc * 100, "%+.2f%%")))
    out.append("")


# --------------------------------------------------------------------------
# plots
# --------------------------------------------------------------------------

def make_plots(agg, out_dir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return []

    made = []
    configs = sorted({(m, f) for (m, f, a, t) in agg})

    # speedup vs threads, one line per config
    fig, ax = plt.subplots(figsize=(7, 5))
    plotted = 0
    for (mesh, factor) in configs:
        ts = sorted(t for (m, f, a, t) in agg if m == mesh and f == factor and a == "par")
        base = agg.get((mesh, factor, "par", 1), {}).get("median_wall")
        if not base or len(ts) < 2:
            continue
        xs = [t for t in ts]
        ys = [base / agg[(mesh, factor, "par", t)]["median_wall"] for t in ts]
        ax.plot(xs, ys, marker="o", label="%s f=%s" % (mesh, factor))
        plotted += 1
    if plotted:
        lim = max(t for (m, f, a, t) in agg if a == "par")
        ax.plot([1, lim], [1, lim], "--", color="0.6", label="ideal")
        ax.set_xlabel("threads")
        ax.set_ylabel("speedup vs par@1")
        ax.set_title("Strong scaling")
        ax.legend(fontsize=7)
        ax.grid(alpha=.3)
        p = out_dir / "speedup_vs_threads.png"
        fig.tight_layout(); fig.savefig(p, dpi=140); made.append(p)
    plt.close(fig)

    # speedup at the highest thread count vs input size
    fig, ax = plt.subplots(figsize=(7, 5))
    maxt = max((t for (m, f, a, t) in agg if a == "par"), default=None)
    pts = []
    for (mesh, factor) in configs:
        base = agg.get((mesh, factor, "par", 1), {}).get("median_wall")
        top = agg.get((mesh, factor, "par", maxt), {})
        if base and top.get("median_wall") and top.get("cells_in"):
            pts.append((top["cells_in"], base / top["median_wall"], mesh, factor))
    if pts:
        ax.scatter([p[0] for p in pts], [p[1] for p in pts])
        for c, s, mesh, factor in pts:
            ax.annotate("%s f=%s" % (mesh, factor), (c, s), fontsize=6,
                        xytext=(3, 3), textcoords="offset points")
        ax.set_xscale("log")
        ax.set_xlabel("input cells")
        ax.set_ylabel("speedup at %d threads" % maxt)
        ax.set_title("Does it scale on big inputs?")
        ax.grid(alpha=.3)
        p = out_dir / "speedup_vs_size.png"
        fig.tight_layout(); fig.savefig(p, dpi=140); made.append(p)
    plt.close(fig)

    return made


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="the results directory")
    ap.add_argument("--out", help="where to write the report (default <results>/../report)")
    ap.add_argument("--manifest",
                    help="the input manifest.json, if it is not beside the "
                         "results directory")
    args = ap.parse_args()

    results_dir = Path(args.results).resolve()
    out_dir = Path(args.out).resolve() if args.out else results_dir.parent / "report"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(results_dir)
    if not rows:
        sys.exit("No usable rows in %s" % results_dir)
    agg = aggregate(rows)

    # tidy table, one row per run
    fields = sorted({k for r in rows for k in r})
    tidy = out_dir / "all_runs.csv"
    with open(tidy, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    threads_max = max((r["threads"] for r in rows), default=1)
    out = ["# Parallel tetrahedral remeshing: scaling and quality\n",
           "%d runs, %d distinct (mesh, factor, arm, threads) cells.\n"
           % (len(rows), len(agg))]
    # If calibration never found a config that ran long enough, every timing
    # below sits closer to the noise floor than intended. That has to be said
    # once at the top, not left for the reader to infer from the CV column.
    cpath = results_dir / "calibration.json"
    if cpath.exists():
        try:
            calib = json.loads(cpath.read_text())
            if calib.get("met_target") is False:
                out.append("> **Calibration missed its target.** No configuration "
                           "reached %.0f s at %d threads, so the sweep fell back to "
                           "the longest runs available. Everything below is noisier "
                           "than intended; check the CV column before quoting any "
                           "single number, and prefer the set totals.\n"
                           % (calib.get("target_seconds", 30),
                              calib.get("threads_max", 24)))
        except Exception:
            pass

    section_env(results_dir, out)
    section_thread_ladder(agg, out)
    section_variants(rows, out)
    section_diagnostics(results_dir, out)
    section_scaling(agg, out, threads_max)
    section_totals(agg, out, threads_max)
    section_machine_behaviour(agg, out)
    section_pipelines(agg, rows, out, load_manifest(results_dir, args.manifest))
    section_quality(agg, out)
    section_validity(rows, agg, out)

    plots = make_plots(agg, out_dir)
    if plots:
        out.append("## Plots\n")
        for p in plots:
            out.append("![%s](%s)\n" % (p.stem, p.name))
    else:
        out.append("_(matplotlib not available, so no plots were drawn.)_\n")

    report = out_dir / "report.md"
    report.write_text("\n".join(out))
    print("Wrote %s" % report)
    print("Wrote %s (%d rows)" % (tidy, len(rows)))
    for p in plots:
        print("Wrote %s" % p)


if __name__ == "__main__":
    main()
