
set(COMMON_SOURCES common.hpp)

add_library(common OBJECT ${COMMON_SOURCES})
target_link_libraries(common)

set(STRDUP_SOURCES strdup.hpp)
add_library(strdup OBJECT ${STRDUP_SOURCES})
target_link_libraries(strdup)