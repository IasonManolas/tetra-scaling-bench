// Timed driver for the 24-core scaling benchmark.
//
// Usage: bench_remesh <input.mesh> <iters> <edge_factor> <smooth_constrained>
//                     <threads> <out.json> [--tag seq|par] [--out-mesh <path>]
//
// Argv 1-6 are deliberately identical to the project's existing
// benchmark_tetrahedral_remeshing, so the parsers already in use keep working.
//
// This driver records ONLY what cannot be recomputed later: the remeshing wall
// time, the input's average edge length (and hence the target), the input
// counts, and the run's identity. It computes NO quality metrics -- it writes
// the remeshed mesh instead, and `mesh_quality_report` derives every quality
// number from that file afterwards. On a multi-million-cell mesh the metric
// pass is not free, and everything it produces is a function of the output.

#include "bench_common.h"
#include "edge_metrics.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/IO/File_medit.h>

#include <tbb/global_control.h>
#include <tbb/version.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <nlohmann/json.hpp>

using K  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
using Cb = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;

// On this branch the concurrency choice is a TEMPLATE PARAMETER of the TDS, not
// a compile-time macro, so both instantiations can live in one binary and be
// selected at run time by --tag. That is the whole reason there is no second
// build here.
template <typename ConcurrencyTag>
using Tds3 = CGAL::Triangulation_data_structure_3<Vb, Cb, ConcurrencyTag>;
template <typename ConcurrencyTag>
using T3   = CGAL::Triangulation_3<K, Tds3<ConcurrencyTag>>;
template <typename ConcurrencyTag>
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<T3<ConcurrencyTag>, int, int>;

template <typename ConcurrencyTag>
static int run(const std::string& input,
               int num_iterations,
               double edge_factor,
               bool smooth_constrained_edges,
               int num_threads,
               const std::string& out_json,
               const std::string& out_mesh,
               nlohmann::json& results_json)
{
  using namespace benchmarking;

  C3t3<ConcurrencyTag> c3t3;
  auto& tr = c3t3.triangulation();

  std::ifstream is(input, std::ios_base::in);
  // allow_non_manifold: read_MEDIT rejects a non-manifold triangulation by
  // default, and that refusal is fatal to a whole measurement cell. It costs
  // nothing on a manifold input and keeps one awkward mesh from taking an arm
  // down with it. (The same flag is required in mesh_quality_report -- see the
  // note there; remeshed Mesh_3 outputs genuinely are non-manifold.)
  if(!CGAL::IO::read_MEDIT(is, tr, CGAL::parameters::allow_non_manifold(true)))
    fatal_error(std::string("Could not read input mesh '") + input + "'");

  const std::string input_name = std::filesystem::path(input).stem().string();
  write_triangulation_info(results_json, tr, input_name);

  // Recorded because it is NOT recoverable from the output mesh: the target is
  // derived from the INPUT's average edge length, so edge-length conformance
  // cannot be judged later without it.
  const double avg_edge_length = compute_average_edge_length(tr);
  if(avg_edge_length <= 0.0)
    fatal_error("Could not compute average edge length.");
  const double target_edge_length = avg_edge_length * edge_factor;

  append_run_info(results_json, "Avg_edge_length", avg_edge_length);
  append_run_info(results_json, "Edge_factor", edge_factor);
  append_run_info(results_json, "Edge Length", target_edge_length);
  append_run_info(results_json, "Num_iterations", num_iterations);
  append_run_info(results_json, "Smooth_constrained_edges", smooth_constrained_edges);

  CGAL::Real_timer t;
  t.start();
  CGAL::tetrahedral_isotropic_remeshing(
      c3t3, target_edge_length,
      CGAL::parameters::number_of_iterations(num_iterations)
                       .smooth_constrained_edges(smooth_constrained_edges));
  t.stop();

  append_metric_result(results_json, "Performance", "Total_Time", "Value", t.time());
  append_metric_result(results_json, "Performance", "Memory", "Value",
                       CGAL::Memory_sizer().virtual_size() >> 20);

  // Feature-edge lengths must be taken here: they do not survive the MEDIT
  // round-trip, so the offline pass cannot recover them. See edge_metrics.h.
  bench_edges::append_complex_edge_lengths(c3t3, results_json);

  // Written after the timer stops, so the I/O is outside the measurement.
  if(!out_mesh.empty())
  {
    // Write to a temporary and rename, so a killed run never leaves a truncated
    // file that the quality pass would later read as a valid mesh.
    const std::filesystem::path final_path = std::filesystem::absolute(out_mesh);
    std::filesystem::create_directories(final_path.parent_path());
    const std::filesystem::path tmp_path = final_path.string() + ".partial";

    CGAL::Real_timer tw;
    tw.start();
    {
      std::ofstream os(tmp_path);
      if(!os) fatal_error("Cannot open '" + tmp_path.string() + "' for writing.");
      CGAL::IO::write_MEDIT(os, tr);
    }
    tw.stop();
    std::filesystem::rename(tmp_path, final_path);

    append_run_info(results_json, "Out_mesh", final_path.string());
    append_metric_result(results_json, "Performance", "Write_Time", "Value", tw.time());
  }

  append_execution_status(results_json, "success");
  write_results_json(results_json, out_json);
  return 0;
}

int main(int argc, char** argv)
{
  using namespace benchmarking;
  nlohmann::json results_json;
  std::cout << std::setprecision(17);
  std::cerr << std::setprecision(17);

  if(argc < 7) {
    fatal_error(std::string("Usage: ") + argv[0] +
                " <input_mesh> <num_iterations> <remeshing_target_edge_factor>"
                " <smooth_constrained_edges> <num_threads> <results_json_path>"
                " [--tag seq|par] [--out-mesh <path.mesh>]");
  }

  const std::string input   = argv[1];
  const int num_iterations  = std::stoi(argv[2]);
  const double edge_factor  = std::stod(argv[3]);
  const bool smooth_ce      = std::stoi(argv[4]) != 0;
  const int num_threads     = std::stoi(argv[5]);
  const std::string out_json = argv[6];

  std::string tag = "par";
  std::string out_mesh;
  for(int i = 7; i < argc; ++i) {
    const std::string a = argv[i];
    if(a == "--tag" && i + 1 < argc)            tag = argv[++i];
    else if(a == "--out-mesh" && i + 1 < argc)  out_mesh = argv[++i];
    else fatal_error("Unknown argument '" + a + "'");
  }
  if(tag != "seq" && tag != "par")
    fatal_error("--tag must be 'seq' or 'par', got '" + tag + "'");
  if(num_threads <= 0)
    fatal_error("num_threads must be positive.");
  if(tag == "seq" && num_threads != 1)
    std::cerr << "Warning: --tag seq ignores num_threads (" << num_threads << ")" << std::endl;

  std::filesystem::create_directories(
      std::filesystem::path(std::filesystem::absolute(out_json)).parent_path());

  append_run_info(results_json, "Arm_tag", tag);
  append_run_info(results_json, "Num_threads", tag == "seq" ? 1 : num_threads);
  append_run_info(results_json, "TBB_version", TBB_VERSION_STRING);

  std::cout << "tag=" << tag << " threads=" << (tag == "seq" ? 1 : num_threads) << std::endl;

  if(tag == "seq")
    return run<CGAL::Sequential_tag>(input, num_iterations, edge_factor, smooth_ce,
                                     1, out_json, out_mesh, results_json);

  // global_control must outlive the remeshing call.
  tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);
  return run<CGAL::Parallel_tag>(input, num_iterations, edge_factor, smooth_ce,
                                 num_threads, out_json, out_mesh, results_json);
}
