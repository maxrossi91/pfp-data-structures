#add_executable(pfp_ds_test pfp_ds_test.cpp)
#target_link_libraries(pfp_ds_test common strdup gsacak sdsl malloc_count)

add_executable(ds_building ds_building.cpp)
target_link_libraries(ds_building common strdup pfp gsacak sdsl malloc_count)

add_executable(ds_building64 ds_building.cpp)
target_link_libraries(ds_building64 common strdup pfp gsacak64 sdsl malloc_count)
target_compile_options(ds_building64 PUBLIC -DM64)

add_executable(ds_building_wt_sdsl ds_building_wt_sdsl.cpp)
target_link_libraries(ds_building_wt_sdsl common strdup pfp gsacak sdsl malloc_count)

add_executable(ds_building_wt_sdsl64 ds_building_wt_sdsl.cpp)
target_link_libraries(ds_building_wt_sdsl64 common strdup pfp gsacak64 sdsl malloc_count)
target_compile_options(ds_building_wt_sdsl64 PUBLIC -DM64)

add_executable(ds_building_wt_sdsl_2 ds_building_wt_sdsl_2.cpp)
target_link_libraries(ds_building_wt_sdsl_2 common strdup pfp gsacak sdsl malloc_count)

add_executable(ds_building_wt_sdsl_264 ds_building_wt_sdsl_2.cpp)
target_link_libraries(ds_building_wt_sdsl_264 common strdup pfp gsacak64 sdsl malloc_count)
target_compile_options(ds_building_wt_sdsl_264 PUBLIC -DM64)

add_executable(sdsl_baseline_cst sdsl_baseline_cst.cpp)
target_link_libraries(sdsl_baseline_cst common strdup sdsl divsufsort divsufsort64 malloc_count)

add_executable(sdsl_baseline_csa sdsl_baseline_csa.cpp)
target_link_libraries(sdsl_baseline_csa common strdup sdsl divsufsort divsufsort64 malloc_count)

add_executable(sdsl_wt_test sdsl_wt_test.cpp)
target_link_libraries(sdsl_wt_test common strdup sdsl divsufsort divsufsort64 malloc_count)

add_executable(lce_test lce_test.cpp)
target_link_libraries(lce_test common strdup pfp gsacak sdsl divsufsort divsufsort64 malloc_count gtest_main pthread)

add_executable(sa_test sa_test.cpp)
target_link_libraries(sa_test common strdup pfp gsacak sdsl divsufsort divsufsort64 malloc_count gtest_main pthread)

add_executable(pfp_test pfp_test.cpp)
target_link_libraries(pfp_test common pfp gsacak sdsl divsufsort divsufsort64 malloc_count gtest_main pthread)

add_executable(lce_benchmark lce_benchmark.cpp)
target_link_libraries(lce_benchmark common strdup pfp gsacak sdsl divsufsort divsufsort64 malloc_count benchmark pthread)

add_executable(sa_benchmark sa_benchmark.cpp)
target_link_libraries(sa_benchmark common strdup pfp gsacak sdsl divsufsort divsufsort64 benchmark pthread)

add_executable(pfp_benchmarks pfp_benchmarks.cpp)
target_link_libraries(pfp_benchmarks common strdup pfp gsacak sdsl divsufsort divsufsort64 benchmark pthread)

add_executable(pfp_load_test pfp_load_test.cpp)
target_link_libraries(pfp_load_test common strdup pfp gsacak sdsl malloc_count)

add_executable(sdsl_load_test sdsl_load_test.cpp)
target_link_libraries(sdsl_load_test common strdup pfp gsacak sdsl malloc_count sdsl divsufsort divsufsort64)
