include(FetchContent)

## Add Big-BWT
FetchContent_Declare(
  bigbwt
  GIT_REPOSITORY https://github.com/alshai/Big-BWT.git
  )
  
FetchContent_GetProperties(bigbwt)
if(NOT bigbwt_POPULATED)
  FetchContent_Populate(bigbwt)
  add_subdirectory(${bigbwt_SOURCE_DIR} ${bigbwt_BINARY_DIR})

  # add_library(gsacak OBJECT ${bigbwt_SOURCE_DIR}/gsa/gsacak.c ${bigbwt_SOURCE_DIR}/gsa/gsacak.h)
  # target_include_directories(gsacak PUBLIC "${bigbwt_SOURCE_DIR}/gsa")
  
  add_library(malloc_count OBJECT ${bigbwt_SOURCE_DIR}/malloc_count.c ${bigbwt_SOURCE_DIR}/malloc_count.h)
  target_link_libraries(malloc_count dl)
  target_include_directories(malloc_count PUBLIC "${bigbwt_SOURCE_DIR}")
  
  endif()
  
## Add gsacak
FetchContent_Declare(
  gsacak
  GIT_REPOSITORY https://github.com/felipelouza/gsa-is.git
)

FetchContent_GetProperties(gsacak)
if(NOT gsacak_POPULATED)
    FetchContent_Populate(gsacak)
    add_library(gsacak OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak PUBLIC "${gsacak_SOURCE_DIR}")
    
    add_library(gsacak64 OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak64 PUBLIC "${gsacak_SOURCE_DIR}")
    target_compile_options(gsacak64 PUBLIC -DM64)
  # add_subdirectory(${gsacak_SOURCE_DIR} ${gsacak_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()


## Google benchmark
FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        v1.5.0
  )
  
FetchContent_GetProperties(benchmark)
if(NOT benchmark_POPULATED)
  FetchContent_Populate(benchmark)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable testing of the benchmark library." FORCE)
  add_subdirectory(${benchmark_SOURCE_DIR} ${benchmark_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()


## Add sdsl
FetchContent_Declare(
  sdsl
  # GIT_REPOSITORY https://github.com/elarielcl/sdsl-lite.git
  GIT_REPOSITORY https://github.com/simongog/sdsl-lite
)

FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
  FetchContent_Populate(sdsl)

  set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")
  
  add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

## Add divsuffsort
FetchContent_Declare(
  divsufsort
  GIT_REPOSITORY https://github.com/simongog/libdivsufsort.git
  GIT_TAG        2.0.1
)

FetchContent_GetProperties(divsufsort)
if(NOT divsufsort_POPULATED)
  FetchContent_Populate(divsufsort)

  set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
  set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
  set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")

  add_subdirectory(${divsufsort_SOURCE_DIR} ${divsufsort_BINARY_DIR} EXCLUDE_FROM_ALL)

  target_include_directories(divsufsort PUBLIC "${divsufsort_BINARY_DIR}/include")
  target_include_directories(divsufsort64 PUBLIC "${divsufsort_BINARY_DIR}/include")
endif()





# ## Add BT-CST
# FetchContent_Declare(
#   bt_cst
#   GIT_REPOSITORY https://github.com/maxrossi91/BT-CST.git
#   # GIT_REPOSITORY https://github.com/simongog/sdsl-lite
# )

# FetchContent_GetProperties(bt_cst)
# if(NOT bt_cst_POPULATED)
#   FetchContent_Populate(bt_cst)

  
#   add_subdirectory(${bt_cst_SOURCE_DIR} ${bt_cst_BINARY_DIR} EXCLUDE_FROM_ALL)
#   target_include_directories(blocktree PUBLIC "${bt_cst_SOURCE_DIR}/includes")
#   message("${bt_cst_SOURCE_DIR}/includes : "${bt_cst_SOURCE_DIR}/includes)
# endif()



# ## Add googletest
# FetchContent_Declare(
#   googletest
#   GIT_REPOSITORY https://github.com/google/googletest.git
#   GIT_TAG        release-1.8.0
# )

# FetchContent_GetProperties(googletest)
# if(NOT googletest_POPULATED)
#   FetchContent_Populate(googletest)
#   add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
#   target_include_directories(gtest_main PUBLIC "${googletest_BINARY_DIR}/include")
# endif()



# ## Add sdsl-lite
# set(sdsl-lite_cmake_file ${CMAKE_CURRENT_LIST_DIR}/sdsl-lite/CMakeLists.txt)
# if(NOT EXISTS ${sdsl-lite_cmake_file})
#     EXECUTE_PROCESS(COMMAND git submodule update --init thirdparty/sdsl-lite
#                     WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
#                     OUTPUT_QUIET
#             )
# endif(NOT EXISTS ${sdsl-lite_cmake_file})

# set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")

# add_subdirectory(sdsl-lite)

# target_include_directories(divsufsort PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/external/libdivsufsort/include")
# target_include_directories(divsufsort64 PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/external/libdivsufsort/include")
# target_include_directories(gtest_main PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/external/googletest/include")


# ## Add BT-CST
# set(BT-CST_cmake_file ${CMAKE_CURRENT_LIST_DIR}/BT-CST/CMakeLists.txt)
# if(NOT EXISTS ${BT-CST_cmake_file})
#     EXECUTE_PROCESS(COMMAND git submodule update --init thirdparty/BT-CST
#                     WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
#                     OUTPUT_QUIET
#             )
# endif(NOT EXISTS ${BT-CST_cmake_file})

# add_subdirectory(BT-CST)

# target_include_directories(blocktree PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/BT-CST/include")
