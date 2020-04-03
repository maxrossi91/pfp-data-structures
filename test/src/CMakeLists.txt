add_executable(pfp_ds_test pfp_ds_test.cpp)
target_link_libraries(pfp_ds_test common gsacak sdsl malloc_count)

add_executable(sdsl_baseline_cst sdsl_baseline_cst.cpp)
target_link_libraries(sdsl_baseline_cst common sdsl divsufsort divsufsort64 malloc_count)
