#pragma once
// Edge-length metrics, split by WHERE they can be computed.
//
// `mesh_quality.h` reports `Quality.Edge_Length` over the c3t3's *complex*
// edges (the feature/constrained edge set). That set is built in memory by the
// remesher and is NOT carried by a MEDIT file: write_MEDIT has named parameters
// for cells and vertices but none for the complex edge set, so after
// read_MEDIT + rescan_after_load_of_triangulation() the complex edge set is
// empty and every `Edge_Length` field comes back as the -1 sentinel. Measured,
// not assumed: an offline pass over a saved mesh reports Edge_Length = -1 while
// an in-process run on the same config reports mean 3.13.
//
// So the two metrics below live on opposite sides of the file:
//
//   append_complex_edge_lengths()  must run IN the benchmark driver, in
//                                  process, before the mesh is written out.
//                                  It is O(#complex edges) -- a small fraction
//                                  of the cell count -- so it costs nothing
//                                  next to the remeshing it follows, and it is
//                                  called after the timer has stopped.
//
//   append_cell_edge_lengths()     runs OFFLINE in mesh_quality_report, over
//                                  every edge of the complex's tetrahedra.
//                                  This is the one the sizing-conformance
//                                  question actually needs: how close the
//                                  output's edges are to the target length.
//                                  `Quality.Edge_Length` cannot answer that --
//                                  it only ever looked at feature edges.
//
// Both emit into the same `metrics.Quality` object as mesh_quality.h, so the
// post-processing sees one flat namespace.

#include "bench_common.h"

#include <CGAL/squared_distance_3.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <nlohmann/json.hpp>

namespace bench_edges {

namespace detail {

inline void append_stats(nlohmann::json& results_json,
                         const std::string& metric_name,
                         std::vector<double>& lengths)
{
  using benchmarking::append_metric_result;

  if(lengths.empty()) {
    // -1 is mesh_quality.h's `metric_undefined`; keep the same sentinel so the
    // post-processing has one rule for "not measurable here".
    append_metric_result(results_json, "Quality", metric_name, "Count", 0);
    append_metric_result(results_json, "Quality", metric_name, "Minimum", -1.);
    append_metric_result(results_json, "Quality", metric_name, "Mean", -1.);
    append_metric_result(results_json, "Quality", metric_name, "Median", -1.);
    append_metric_result(results_json, "Quality", metric_name, "Maximum", -1.);
    append_metric_result(results_json, "Quality", metric_name, "Std_Dev", -1.);
    return;
  }

  const double total = std::accumulate(lengths.begin(), lengths.end(), 0.);
  const double mean = total / lengths.size();

  double sq = 0.;
  for(const double l : lengths) sq += (l - mean) * (l - mean);
  const double stddev = std::sqrt(sq / lengths.size());

  const std::size_t mid = lengths.size() / 2;
  std::nth_element(lengths.begin(), lengths.begin() + mid, lengths.end());
  const double median = lengths[mid];

  append_metric_result(results_json, "Quality", metric_name, "Count",
                       static_cast<double>(lengths.size()));
  append_metric_result(results_json, "Quality", metric_name, "Minimum",
                       *std::min_element(lengths.begin(), lengths.end()));
  append_metric_result(results_json, "Quality", metric_name, "Mean", mean);
  append_metric_result(results_json, "Quality", metric_name, "Median", median);
  append_metric_result(results_json, "Quality", metric_name, "Maximum",
                       *std::max_element(lengths.begin(), lengths.end()));
  append_metric_result(results_json, "Quality", metric_name, "Std_Dev", stddev);
}

} // namespace detail

/// Feature/complex edges. Driver-side only -- see the header comment.
template <typename C3t3>
void append_complex_edge_lengths(const C3t3& c3t3, nlohmann::json& results_json)
{
  std::vector<double> lengths;
  lengths.reserve(c3t3.number_of_edges_in_complex());

  for(auto eit = c3t3.edges_in_complex_begin();
      eit != c3t3.edges_in_complex_end(); ++eit)
  {
    const auto& p1 = eit->first->vertex(eit->second)->point();
    const auto& p2 = eit->first->vertex(eit->third)->point();
    lengths.push_back(std::sqrt(CGAL::squared_distance(p1, p2)));
  }

  detail::append_stats(results_json, "Edge_Length", lengths);
}

/// Every edge of every in-complex tetrahedron, each counted once.
/// This is what sizing conformance is measured against.
template <typename C3t3>
void append_cell_edge_lengths(const C3t3& c3t3, nlohmann::json& results_json)
{
  const auto& tr = c3t3.triangulation();

  std::vector<double> lengths;
  for(auto eit = tr.finite_edges_begin(); eit != tr.finite_edges_end(); ++eit)
  {
    // An edge belongs to the complex if any tetrahedron around it does. Walking
    // cells and taking their 6 edges instead would count each edge once per
    // incident cell, which silently weights interior edges by their valence and
    // would bias the mean.
    bool in_complex = false;
    auto circ = tr.incident_cells(*eit), done = circ;
    if(circ != nullptr) {
      do {
        if(c3t3.is_in_complex(circ)) { in_complex = true; break; }
      } while(++circ != done);
    }
    if(!in_complex) continue;

    const auto& p1 = eit->first->vertex(eit->second)->point();
    const auto& p2 = eit->first->vertex(eit->third)->point();
    lengths.push_back(std::sqrt(CGAL::squared_distance(p1, p2)));
  }

  detail::append_stats(results_json, "Cell_Edge_Length", lengths);
}

} // namespace bench_edges
