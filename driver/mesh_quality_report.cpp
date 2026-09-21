// Computes the mesh_quality.h metrics for an existing MEDIT .mesh file, so that
// meshes produced by other tools (ParMMG) are measured by the same code as ours.
//
// Usage: mesh_quality_report <input.mesh> <results_json_path>
//
// Only cells with a non-zero subdomain index are measured. This is CGAL's own
// definition of complex membership, and it is also the ref!=0 filter ParMMG
// meshes need: their exterior hull-fill tets are written with ref 0.

#include "bench_common.h"
#include "mesh_quality.h"
#include "edge_metrics.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/IO/File_medit.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <nlohmann/json.hpp>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
using Cb = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
using T3 = CGAL::Triangulation_3<K, CGAL::Triangulation_data_structure_3<Vb, Cb>>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<T3, int, int>;

int main(int argc, char** argv)
{
  using namespace benchmarking;
  using nlohmann::json;

  std::cout << std::setprecision(17);
  if(argc != 3 && argc != 4) {
    fatal_error(std::string("Usage: ") + argv[0] +
                " <input.mesh> <results_json_path> [subdomain_index]");
  }
  const std::string input = argv[1];
  const std::string results_json_path = argv[2];
  // measure one subdomain only; useful when a single subdomain was remeshed and
  // the others, left as the mesher produced them, would dominate the extremes
  const bool one_subdomain = (argc == 4);
  const int subdomain = one_subdomain ? std::stoi(argv[3]) : 0;

  std::filesystem::create_directories(
    std::filesystem::path(std::filesystem::absolute(results_json_path)).parent_path());

  C3t3 c3t3;
  T3& tr = c3t3.triangulation();
  std::ifstream is(input, std::ios_base::in);
  // allow_non_manifold is required, not defensive. read_MEDIT refuses a
  // non-manifold triangulation by default, and remeshing a Mesh_3-derived
  // input produces outputs that trip it: measured here, every Mesh_3 output
  // failed to read while every CDT output of the same run read fine. Refusing
  // them would silently drop one of the two input pipelines from the quality
  // half of the report -- the half the pipeline comparison depends on.
  if(!CGAL::IO::read_MEDIT(is, tr, CGAL::parameters::allow_non_manifold(true)))
    fatal_error(std::string("Could not read input mesh '") + input + "'");

  // the metrics are computed over the cells of the complex, so restricting to
  // one subdomain is a matter of putting every other cell outside of it
  if(one_subdomain)
  {
    for(auto c : tr.finite_cell_handles())
    {
      if(c->subdomain_index() != subdomain)
        c->set_subdomain_index(0);
    }
  }

  // read_MEDIT fills the triangulation and its per-cell subdomain indices, but
  // leaves the complex bookkeeping empty; this rebuilds it, counting only the
  // non-zero-subdomain cells.
  c3t3.rescan_after_load_of_triangulation();

  std::cout << "Vertices: " << tr.number_of_vertices() << "\n"
            << "Finite cells in file: " << tr.number_of_finite_cells() << "\n";
  if(one_subdomain)
    std::cout << "Cells in complex (subdomain == " << subdomain << "): ";
  else
    std::cout << "Cells in complex (subdomain != 0): ";
  std::cout << c3t3.number_of_cells_in_complex() << "\n"
            << "Facets in complex: " << c3t3.number_of_facets_in_complex() << std::endl;

  if(c3t3.number_of_cells_in_complex() == 0)
    fatal_error("No cells in complex - every metric would be undefined.");

  json results_json;
  const std::string input_name = std::filesystem::path(input).stem().string();
  write_triangulation_info(results_json, tr, input_name);
  append_run_info(results_json, "Source", input);
  append_run_info(results_json, "Measured_subdomain",
                  one_subdomain ? std::to_string(subdomain) : std::string("all (!= 0)"));
  generate_quality_metrics(c3t3, results_json);

  // Quality.Edge_Length above is the complex/feature edge set, which a MEDIT
  // file does not carry; it reads back as the -1 sentinel here. The sizing
  // conformance question is answered by the tetrahedra own edges instead.
  bench_edges::append_cell_edge_lengths(c3t3, results_json);

  std::ofstream os(results_json_path);
  os << std::setw(2) << results_json << std::endl;
  std::cout << "Wrote results to " << results_json_path << std::endl;
  return 0;
}
