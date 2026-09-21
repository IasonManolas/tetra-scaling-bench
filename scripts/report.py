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


def aggregate(rows):
    """Median wall time and CV per (mesh, factor, arm, threads)."""
    groups = defaultdict(list)
    for r in rows:
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
    out.append("| governor | %s |" % (env.get("cpu_governor") or "unknown"))
    out.append("| platform | %s |" % env.get("platform"))
    tc = env.get("toolchain", {})
    out.append("| ours | `%s` @ `%s` |" % (tc.get("ours_ref"), (tc.get("ours_sha") or "")[:12]))
    out.append("| main | `%s` @ `%s` |" % (tc.get("main_ref"), (tc.get("main_sha") or "")[:12]))
    out.append("")
    gov = env.get("cpu_governor")
    if gov and gov != "performance":
        out.append("> The CPU governor was `%s`, not `performance`. Low thread counts "
                   "look worse than they are under a scaling governor, so the speedup "
                   "numbers below are, if anything, optimistic.\n" % gov)


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
    section_scaling(agg, out, threads_max)
    section_totals(agg, out, threads_max)
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
