#add_executable(pfp_ds_test pfp_ds_test.cpp)
#target_link_libraries(pfp_ds_test common gsacak sdsl malloc_count)

add_executable(ds_building ds_building.cpp)
target_link_libraries(ds_building common pfp gsacak sdsl malloc_count)

add_executable(sdsl_baseline_cst sdsl_baseline_cst.cpp)
target_link_libraries(sdsl_baseline_cst common sdsl divsufsort divsufsort64 malloc_count)

add_executable(sdsl_baseline_csa sdsl_baseline_csa.cpp)
target_link_libraries(sdsl_baseline_csa common sdsl divsufsort divsufsort64 malloc_count)
