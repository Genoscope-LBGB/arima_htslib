cmake_minimum_required(VERSION 3.14)
project(arima_htslib LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    htslib
    GIT_REPOSITORY https://github.com/samtools/htslib.git
    GIT_TAG 1.19
)
FetchContent_MakeAvailable(htslib)

add_subdirectory(filter_five_end)
add_subdirectory(get_stats)
add_subdirectory(two_read_bam_combiner)
