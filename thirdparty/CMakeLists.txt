## Add Big-BWT
set(Big-BWT_cmake_file ${CMAKE_CURRENT_LIST_DIR}/Big-BWT/CMakeLists.txt)
if(NOT EXISTS ${Big-BWT_cmake_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init thirdparty/Big-BWT
                    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
                    OUTPUT_QUIET
            )
endif(NOT EXISTS ${Big-BWT_cmake_file})

add_subdirectory(Big-BWT)

add_library(gsacak OBJECT Big-BWT/gsa/gsacak.c Big-BWT/gsa/gsacak.h)
target_include_directories(gsacak PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/Big-BWT/gsa")
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/Big-BWT/gsa/gsacak.h" "${CMAKE_CURRENT_BINARY_DIR}/Big-BWT/gsa/gsacak.h" COPYONLY )



## Add sdsl-lite
set(sdsl-lite_cmake_file ${CMAKE_CURRENT_LIST_DIR}/sdsl-lite/CMakeLists.txt)
if(NOT EXISTS ${sdsl-lite_cmake_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init thirdparty/sdsl-lite
                    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
                    OUTPUT_QUIET
            )
endif(NOT EXISTS ${sdsl-lite_cmake_file})

set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")

add_subdirectory(sdsl-lite)

target_include_directories(divsufsort PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/external/libdivsufsort/include")
