# Build file for the 24-core scaling benchmark kit.
#
# This same file is configured TWICE, into two build directories, because
# CGAL_DIR is a configure-time global and the reference arm must compile against
# a different CGAL tree:
#
#   cmake -S driver -B build/ours -DCGAL_DIR=<our branch worktree> -DARM=ours
#   cmake -S driver -B build/main -DCGAL_DIR=<cgal main worktree> -DARM=main
#
# ARM=ours builds bench_remesh (both concurrency tags, needs TBB),
# preprocess_cdt and mesh_quality_report. ARM=main builds only
# bench_remesh_main, which is sequential and needs no TBB.
#
# scripts/setup.py drives both configurations; see the README.

cmake_minimum_required(VERSION 3.12...3.31)
project(tetra_scaling_bench CXX)

set(ARM "ours" CACHE STRING "Which arm to build: ours | main")
set_property(CACHE ARM PROPERTY STRINGS ours main)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

find_package(CGAL REQUIRED)

# Vendored so the kit builds with no network access and no package manager.
add_library(bench_json INTERFACE)
target_include_directories(bench_json INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/external")

function(bench_add_target name src)
  add_executable(${name} "${src}")
  target_include_directories(${name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(${name} PRIVATE CGAL::CGAL bench_json)
  target_compile_features(${name} PRIVATE cxx_std_17)
endfunction()

if(ARM STREQUAL "ours")
  # The parallel arm needs TBB. This is REQUIRED rather than optional on
  # purpose: silently falling back to a TBB-less build would produce a
  # bench_remesh whose --tag par is a lie, and the whole scaling measurement
  # would be meaningless while still looking like it worked.
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)

  bench_add_target(bench_remesh        bench_remesh.cpp)
  bench_add_target(preprocess_cdt      preprocess_cdt.cpp)
  bench_add_target(preprocess_mesh3    preprocess_mesh3.cpp)
  bench_add_target(mesh_quality_report mesh_quality_report.cpp)

  if(NOT TARGET CGAL::TBB_support)
    message(FATAL_ERROR "TBB was found but CGAL::TBB_support was not created; "
                        "check CGAL_DIR and TBB_DIR.")
  endif()
  target_link_libraries(bench_remesh PRIVATE CGAL::TBB_support)

elseif(ARM STREQUAL "main")
  # Upstream's remeshing triangulation has no concurrency tag, so this arm is
  # sequential by construction and deliberately links no TBB.
  bench_add_target(bench_remesh_main bench_remesh_main.cpp)

else()
  message(FATAL_ERROR "ARM must be 'ours' or 'main', got '${ARM}'")
endif()
