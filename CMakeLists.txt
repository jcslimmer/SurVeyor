cmake_minimum_required(VERSION 2.8)
project(SurVeyor)

set(CMAKE_CXX_STANDARD 11)
if (NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    message("Building with -march=native")
endif()
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(htslib-1.18)
include_directories(htslib-1.18/include)

add_executable(survindel2_reads_categorizer SurVIndel2/reads_categorizer.cpp)
add_executable(survindel2_clip_consensus_builder libs/ssw.c libs/ssw_cpp.cpp SurVIndel2/clip_consensus_builder.cpp)
add_executable(survindel2_merge_identical_calls SurVIndel2/merge_identical_calls.cpp)
add_executable(survindel2_dp_clusterer libs/ssw.c libs/ssw_cpp.cpp libs/kdtree.c SurVIndel2/dp_clusterer.cpp)
add_executable(survindel2_normalise SurVIndel2/normalise.cpp)

add_executable(insurveyor_reads_categorizer INSurVeyor/reads_categorizer.cpp)
add_executable(insurveyor_clip_consensus_builder INSurVeyor/clip_consensus_builder.cpp)
add_executable(insurveyor_call_insertions libs/ssw.c libs/ssw_cpp.cpp INSurVeyor/call_insertions.cpp)
add_executable(insurveyor_dc_remapper libs/ssw.c libs/ssw_cpp.cpp INSurVeyor/dc_remapper.cpp)
add_executable(insurveyor_add_filtering_info INSurVeyor/add_filtering_info.cpp)
add_executable(insurveyor_filter INSurVeyor/filter.cpp)
add_executable(insurveyor_normalise INSurVeyor/normalise.cpp)

add_executable(reads_categorizer src/reads_categorizer.cpp)

add_executable(compare libs/ssw.c libs/ssw_cpp.cpp utils/compare.cpp)

find_library(HTS_LIB hts PATHS htslib-1.18/lib NO_DEFAULT_PATH)
target_link_libraries(survindel2_reads_categorizer "${HTS_LIB}")
target_link_libraries(survindel2_clip_consensus_builder "${HTS_LIB}")
target_link_libraries(survindel2_merge_identical_calls "${HTS_LIB}")
target_link_libraries(survindel2_dp_clusterer "${HTS_LIB}")
target_link_libraries(survindel2_normalise "${HTS_LIB}")

target_link_libraries(insurveyor_reads_categorizer "${HTS_LIB}")
target_link_libraries(insurveyor_clip_consensus_builder "${HTS_LIB}")
target_link_libraries(insurveyor_call_insertions "${HTS_LIB}")
target_link_libraries(insurveyor_dc_remapper "${HTS_LIB}")
target_link_libraries(insurveyor_add_filtering_info "${HTS_LIB}")
target_link_libraries(insurveyor_filter "${HTS_LIB}")
target_link_libraries(insurveyor_normalise "${HTS_LIB}")

target_link_libraries(reads_categorizer "${HTS_LIB}")

target_link_libraries(compare "${HTS_LIB}")
