cmake_minimum_required(VERSION 3.16)
project(simio LANGUAGES C CXX)

option(SIMIO_BUILD_EXAMPLES "Build examples" ON)
option(SIMIO_BUILD_TESTS "Build tests" ON)
option(SIMIO_DEBUG "Enable debug assertions/checks" ON)
option(SIMIO_TRACE "Enable verbose per-thread trace logs in example" OFF)

find_package(Threads REQUIRED)

add_library(simio INTERFACE)
target_include_directories(simio
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_compile_features(simio INTERFACE cxx_std_17)
target_link_libraries(simio INTERFACE Threads::Threads)
# GCC < 9 needs explicit std::filesystem linkage.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
   CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    target_link_libraries(simio INTERFACE stdc++fs)
endif()
if(SIMIO_DEBUG)
    target_compile_definitions(simio INTERFACE SIMIO_DEBUG=1)
else()
    target_compile_definitions(simio INTERFACE SIMIO_DEBUG=0)
endif()

add_library(reader_xtc STATIC
    reader/src/xtc.c
    reader/src/stub_trr.c
)
target_include_directories(reader_xtc PUBLIC reader/include)

add_library(simio_properties STATIC
    src/properties/DensityXZ.cpp
    src/properties/DipoleXZ.cpp
    src/properties/GatingCenterPlane.cpp
    src/properties/GatingPbcWrap.cpp
)
target_include_directories(simio_properties
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_compile_features(simio_properties PUBLIC cxx_std_17)
target_link_libraries(simio_properties PUBLIC simio)

if(SIMIO_BUILD_EXAMPLES)
    add_executable(simio_example examples/main.cpp)
    target_link_libraries(simio_example PRIVATE simio)
    if(SIMIO_TRACE)
        target_compile_definitions(simio_example PRIVATE SIMIO_TRACE=1)
    else()
        target_compile_definitions(simio_example PRIVATE SIMIO_TRACE=0)
    endif()

    add_executable(simio_xtc_pipeline examples/xtc_pipeline_real.cpp)
    target_link_libraries(simio_xtc_pipeline PRIVATE simio reader_xtc)

    add_executable(simio_xtc_density_x examples/xtc_density_x.cpp)
    target_link_libraries(simio_xtc_density_x PRIVATE simio reader_xtc)

    add_executable(simio_xtc_dipole_x examples/xtc_dipole_x.cpp)
    target_link_libraries(simio_xtc_dipole_x PRIVATE simio reader_xtc)

    add_executable(simio_xtc_coord_x examples/xtc_coord_x.cpp)
    target_link_libraries(simio_xtc_coord_x PRIVATE simio reader_xtc)

    add_executable(simio_xtc_gating_flux examples/xtc_gating_flux.cpp)
    target_link_libraries(simio_xtc_gating_flux PRIVATE simio simio_properties reader_xtc)

    add_executable(simio_xtc_all_properties
        examples/xtc_all_properties.cpp
        examples/all_props/channel_count_xz.cpp
        examples/all_props/density_x.cpp
        examples/all_props/density_z_in_x_channel.cpp
        examples/all_props/dipole_x.cpp
        examples/all_props/dipole_z_in_x_channel.cpp
        examples/all_props/coord_x.cpp
        examples/all_props/gating_flux.cpp
        examples/all_props/jump_msd.cpp
    )
    target_link_libraries(simio_xtc_all_properties PRIVATE simio reader_xtc)

    add_executable(simio_layered_demo examples/layered_arch_demo.cpp)
    target_link_libraries(simio_layered_demo PRIVATE simio simio_properties)
endif()

if(SIMIO_BUILD_TESTS)
    enable_testing()

    add_executable(simio_smoke_test tests/smoke_test.cpp)
    target_link_libraries(simio_smoke_test PRIVATE simio)
    add_test(NAME simio_smoke_test COMMAND simio_smoke_test)

    add_executable(simio_system_5w_3na_3cl_test tests/system_5w_3na_3cl_test.cpp)
    target_link_libraries(simio_system_5w_3na_3cl_test PRIVATE simio)
    add_test(NAME simio_system_5w_3na_3cl_test COMMAND simio_system_5w_3na_3cl_test)

    add_executable(simio_pbc_boundary_crossing_test tests/pbc_boundary_crossing_test.cpp)
    target_link_libraries(simio_pbc_boundary_crossing_test PRIVATE simio)
    add_test(NAME simio_pbc_boundary_crossing_test COMMAND simio_pbc_boundary_crossing_test)

    add_executable(simio_layered_flux_gating_test tests/layered_flux_gating_test.cpp)
    target_link_libraries(simio_layered_flux_gating_test PRIVATE simio simio_properties)
    add_test(NAME simio_layered_flux_gating_test COMMAND simio_layered_flux_gating_test)

    add_executable(simio_layered_density_roi_test tests/layered_density_roi_test.cpp)
    target_link_libraries(simio_layered_density_roi_test PRIVATE simio simio_properties)
    add_test(NAME simio_layered_density_roi_test COMMAND simio_layered_density_roi_test)
endif()
