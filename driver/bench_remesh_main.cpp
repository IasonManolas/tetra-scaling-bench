// The reference arm: the same timed driver, built against upstream cgal/main.
//
// Usage: bench_remesh_main <input.mesh> <iters> <edge_factor> <smooth_constrained>
//                          <threads> <out.json> [--tag seq] [--out-mesh <path>]
//
// Upstream's Remeshing triangulation has no concurrency tag parameter, so this
// is sequential by construction. `--tag` is accepted and ignored so that
// run_bench.py can build one command line for every arm; `threads` must be 1.
//
// Kept as a separate translation unit from bench_remesh.cpp because it compiles
// against a DIFFERENT CGAL tree (a different -DCGAL_DIR), not because the body
// differs in any interesting way. Everything else -- the JSON schema, the
// timer placement, the atomic mesh write -- is identical on purpose, so the
// reference arm is measured by exactly the same instrument.

#include "bench_common.h"
#include "edge_metrics.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/IO/File_medit.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <nlohmann/json.hpp>

using K   = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb  = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
using Cb  = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
using T3  = CGAL::Triangulation_3<K, CGAL::Triangulation_data_structure_3<Vb, Cb>>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<T3, int, int>;

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
                " [--tag seq] [--out-mesh <path.mesh>]");
  }

  const std::string input    = argv[1];
  const int num_iterations   = std::stoi(argv[2]);
  const double edge_factor   = std::stod(argv[3]);
  const bool smooth_ce       = std::stoi(argv[4]) != 0;
  const int num_threads      = std::stoi(argv[5]);
  const std::string out_json = argv[6];

  std::string out_mesh;
  for(int i = 7; i < argc; ++i) {
    const std::string a = argv[i];
    if(a == "--tag" && i + 1 < argc)           ++i;   // accepted and ignored
    else if(a == "--out-mesh" && i + 1 < argc) out_mesh = argv[++i];
    else fatal_error("Unknown argument '" + a + "'");
  }
  if(num_threads != 1)
    std::cerr << "Warning: upstream main is sequential; ignoring num_threads="
              << num_threads << std::endl;

  std::filesystem::create_directories(
      std::filesystem::path(std::filesystem::absolute(out_json)).parent_path());

  append_run_info(results_json, "Arm_tag", "main");
  append_run_info(results_json, "Num_threads", 1);

  C3t3 c3t3;
  T3& tr = c3t3.triangulation();
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

  const double avg_edge_length = compute_average_edge_length(tr);
  if(avg_edge_length <= 0.0)
    fatal_error("Could not compute average edge length.");
  const double target_edge_length = avg_edge_length * edge_factor;

  append_run_info(results_json, "Avg_edge_length", avg_edge_length);
  append_run_info(results_json, "Edge_factor", edge_factor);
  append_run_info(results_json, "Edge Length", target_edge_length);
  append_run_info(results_json, "Num_iterations", num_iterations);
  append_run_info(results_json, "Smooth_constrained_edges", smooth_ce);

  CGAL::Real_timer t;
  t.start();
  CGAL::tetrahedral_isotropic_remeshing(
      c3t3, target_edge_length,
      CGAL::parameters::number_of_iterations(num_iterations)
                       .smooth_constrained_edges(smooth_ce));
  t.stop();

  append_metric_result(results_json, "Performance", "Total_Time", "Value", t.time());
  append_metric_result(results_json, "Performance", "Memory", "Value",
                       CGAL::Memory_sizer().virtual_size() >> 20);

  // Feature-edge lengths must be taken here: they do not survive the MEDIT
  // round-trip, so the offline pass cannot recover them. See edge_metrics.h.
  bench_edges::append_complex_edge_lengths(c3t3, results_json);

  if(!out_mesh.empty())
  {
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
