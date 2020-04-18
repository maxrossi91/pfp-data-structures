#add_executable(pfp_ds_test pfp_ds_test.cpp)
#target_link_libraries(pfp_ds_test common strdup gsacak sdsl malloc_count)

add_executable(ds_building ds_building.cpp)
target_link_libraries(ds_building common strdup pfp gsacak sdsl malloc_count)

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

add_executable(lce_benchmark lce_benchmark.cpp)
target_link_libraries(lce_benchmark common strdup pfp gsacak sdsl divsufsort divsufsort64 malloc_count benchmark pthread)