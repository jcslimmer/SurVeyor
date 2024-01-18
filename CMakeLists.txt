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

set(CMAKE_OSX_ARCHITECTURES "x86_64")

include_directories(htslib-1.18)
include_directories(htslib-1.18/include)

add_executable(normalise utils/normalise.cpp)

add_executable(reads_categorizer src/reads_categorizer.cpp)
add_executable(clip_consensus_builder libs/ssw.c libs/ssw_cpp.cpp src/clip_consensus_builder.cpp)
add_executable(find_svs_from_sr_consensuses libs/ssw.c libs/ssw_cpp.cpp src/find_svs_from_sr_consensuses.cpp)
add_executable(merge_identical_calls src/merge_identical_calls.cpp)
add_executable(dp_clusterer libs/ssw.c libs/ssw_cpp.cpp src/dp_clusterer.cpp)
add_executable(insertions_assembler libs/ssw.c libs/ssw_cpp.cpp src/insertions_assembler.cpp)
add_executable(add_filtering_info src/add_filtering_info.cpp)

add_executable(compare libs/ssw.c libs/ssw_cpp.cpp utils/compare.cpp)

find_library(HTS_LIB hts PATHS htslib-1.18/lib NO_DEFAULT_PATH)

target_link_libraries(reads_categorizer "${HTS_LIB}")
target_link_libraries(clip_consensus_builder "${HTS_LIB}")
target_link_libraries(find_svs_from_sr_consensuses "${HTS_LIB}")
target_link_libraries(merge_identical_calls "${HTS_LIB}")
target_link_libraries(dp_clusterer "${HTS_LIB}")
target_link_libraries(insertions_assembler "${HTS_LIB}")
target_link_libraries(add_filtering_info "${HTS_LIB}")

target_link_libraries(compare "${HTS_LIB}")
target_link_libraries(normalise "${HTS_LIB}")
