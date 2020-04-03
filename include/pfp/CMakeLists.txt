set(PFP_SOURCES dictionary.hpp
                parse.hpp
                wt.hpp
                lce_support.hpp
                sa_support.hpp
                pfp.hpp)

add_library(pfp OBJECT ${PFP_SOURCES})
target_link_libraries(pfp common)
