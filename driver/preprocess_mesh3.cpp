// Companion to thingi10k_preprocess_cdt.cpp: converts the SAME Thingi10K .off
// surface into a Mesh_3 volume mesh instead of a conforming constrained
// Delaunay triangulation, so each surface yields a MATCHED PAIR of inputs.
//
// Why: the project's acceptance set is entirely Mesh_3 while the 100-mesh
// Thingi sweep is entirely CDT, and the two disagree on the SIGN of the
// sequential ours-vs-main comparison (0.95 vs 1.20). Matched pairs make input
// quality a measurable axis instead of a confound between datasets.
//
// Usage: thingi10k_preprocess_mesh3 <input.off> <output.mesh> [facet_size_rel]
//   facet_size_rel: facet/cell size as a fraction of the bounding-box diagonal
//                   (default 0.02). Sizing is set RELATIVE so meshes of very
//                   different scale get comparable element counts.

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

using K       = CGAL::Exact_predicates_inexact_constructions_kernel;
using Surface = CGAL::Surface_mesh<K::Point_3>;
using Domain  = CGAL::Polyhedral_mesh_domain_with_features_3<K, Surface>;

// Default Mesh_3 bases. Mesh_3 builds a REGULAR triangulation, so the
// Remeshing_* bases cannot be substituted here; the documented route is to
// mesh with the default bases and then hand the result to
// convert_to_triangulation_3(), exactly as CGAL's own
// mesh_and_remesh_with_adaptive_sizing.cpp example does.
using Tr   = CGAL::Mesh_triangulation_3<Domain>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;
using Criteria = CGAL::Mesh_criteria_3<Tr>;

// Triangulation for remeshing / MEDIT output
using T3 = CGAL::Triangulation_3<Tr::Geom_traits, Tr::Triangulation_data_structure>;

int main(int argc, char** argv)
{
  if(argc != 3 && argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input.off> <output.mesh> [facet_size_rel]\n";
    return EXIT_FAILURE;
  }
  const std::string input = argv[1], output = argv[2];
  const double rel = (argc == 4) ? std::stod(argv[3]) : 0.02;

  Surface surface;
  if(!CGAL::IO::read_polygon_mesh(input, surface) || surface.is_empty()) {
    std::cerr << "Error: cannot read " << input << "\n";
    return EXIT_FAILURE;
  }
  if(!CGAL::is_triangle_mesh(surface)) {
    std::cerr << "Error: not a triangle mesh: " << input << "\n";
    return EXIT_FAILURE;
  }

  const CGAL::Bbox_3 bb = CGAL::Polygon_mesh_processing::bbox(surface);
  const double diag = std::sqrt(CGAL::square(bb.xmax()-bb.xmin())
                              + CGAL::square(bb.ymax()-bb.ymin())
                              + CGAL::square(bb.zmax()-bb.zmin()));
  const double h = rel * diag;

  Domain domain(surface);
  domain.detect_features();

  Criteria criteria(CGAL::parameters::edge_size(h)
                                     .facet_angle(25)
                                     .facet_size(h)
                                     .facet_distance(h/10.)
                                     .cell_radius_edge_ratio(3)
                                     .cell_size(h));

  CGAL::Real_timer t; t.start();
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      CGAL::parameters::no_perturb().no_exude());
  t.stop();

  if(c3t3.triangulation().number_of_cells() == 0) {
    std::cerr << "Error: Mesh_3 produced an empty triangulation for " << input << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "Mesh_3: " << c3t3.triangulation().number_of_vertices() << " vertices, "
            << c3t3.number_of_cells_in_complex() << " cells in complex, "
            << t.time() << "s (h=" << h << ")" << std::endl;

  std::ofstream out(output);
  if(!out) { std::cerr << "Error: cannot open " << output << "\n"; return EXIT_FAILURE; }
  T3 tr = CGAL::convert_to_triangulation_3(std::move(c3t3));
  CGAL::IO::write_MEDIT(out, tr);
  return EXIT_SUCCESS;
}
