// Stage 1 of the Thingi10K benchmarking pipeline: converts a surface mesh
// (Thingi10K "fixed"/autorefined .off, faces already conforming) into an
// initial constrained tetrahedralization, written out as a MEDIT .mesh file.
//
// This is a one-time, untimed preprocessing step per input mesh. It does NOT
// run tetrahedral_isotropic_remeshing -- that happens later, timed, in the
// stage-2 benchmark executables (benchmark_original_tetrahedral_remeshing,
// benchmark_refactored_sequential_tetrahedral_remeshing,
// benchmark_refactored_parallel_tetrahedral_remeshing), which read the
// .mesh file produced here.
//
// Usage: thingi10k_preprocess_cdt <input.off> <output.mesh>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Conforming_constrained_Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_cell_base_3.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_vertex_base_3.h>

#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Real_timer.h>

#include <fstream>
#include <iostream>
#include <string>

using K    = CGAL::Exact_predicates_inexact_constructions_kernel;

using Vbb  = CGAL::Tetrahedral_remeshing::Remeshing_vertex_base_3<K>;
using Vb   = CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3<K, Vbb>;

using Cbb  = CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3<K>;
using Cb   = CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3<K, Cbb>;

using Tds  = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Tr   = CGAL::Triangulation_3<K, Tds>;
using CCDT = CGAL::Conforming_constrained_Delaunay_triangulation_3<K, Tr>;

using CCDT_Tr = CCDT::Triangulation;
using Triangulation_3 = CGAL::Triangulation_3<K, CCDT_Tr::Triangulation_data_structure>;

int main(int argc, char** argv)
{
  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input.off> <output.mesh>" << std::endl;
    return EXIT_FAILURE;
  }
  const std::string input = argv[1];
  const std::string output = argv[2];

  CGAL::Surface_mesh<K::Point_3> mesh;
  if(!CGAL::IO::read_polygon_mesh(input, mesh))
  {
    std::cerr << "Error: cannot read file " << input << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Read " << mesh.number_of_vertices() << " vertices, "
            << mesh.number_of_faces() << " faces from " << input << std::endl;

  CGAL::Real_timer t;
  t.start();

  CCDT ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3<CCDT>(mesh);

  namespace Tet_remesh = CGAL::Tetrahedral_remeshing;
  Tr tr = Tet_remesh::get_remeshing_triangulation(std::move(ccdt));

  t.stop();

  std::cout << "CDT construction: " << tr.number_of_vertices() << " vertices, "
            << tr.number_of_cells() << " cells, " << t.time() << "s" << std::endl;

  std::ofstream out(output);
  if(!out) {
    std::cerr << "Error: cannot open output file " << output << std::endl;
    return EXIT_FAILURE;
  }
  CGAL::IO::write_MEDIT(out, tr);

  return EXIT_SUCCESS;
}
