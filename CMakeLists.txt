cmake_minimum_required(VERSION 3.20)
project(FastANI
    VERSION 1.34
    DESCRIPTION "Fast computation of whole-genome Average Nucleotide Identity (ANI)."
    LANGUAGES CXX)
#
# Main executable
add_executable(FastANI src/cgi/core_genome_identity.cpp src/cgi/main.cpp)
target_include_directories(FastANI PUBLIC "src/")
#
# by default, static linking
option(GSL_SHARED "Build using shared libraries" OFF)
option(BUILD_TESTING "Build tests" ON) # build tests
#
# zlib dependency
find_package(ZLIB 1.1 REQUIRED)
if(${ZLIB_FOUND})
    message("Using zlib library v${ZLIB_VERSION_STRING}")
endif(${ZLIB_FOUND})
target_link_libraries(FastANI PUBLIC ZLIB::ZLIB)
#
# openmp dependency
find_package(OpenMP)
if(${OpenMP_CXX_FOUND})
    message("Using OpenMP_CXX v${OpenMP_CXX_VERSION}")
    target_link_libraries(FastANI PUBLIC ${OpenMP_CXX_LIBRARIES})
    target_include_directories(FastANI PUBLIC ${OpenMP_CXX_INCLUDE_DIRS})
    target_compile_options(FastANI PUBLIC ${OpenMP_CXX_FLAGS})
else()
    message("NOTE: OpenMP not found; Compiling as a single threaded app...")
endif(${OpenMP_CXX_FOUND})
#
# GSL or Boost math
find_package(GSL 1.6)
if(${GSL_FOUND})
    message("GSL v${GSL_VERSION} Found: ${GSL_LIBDIR}")
    if(${GSL_SHARED})
        target_link_libraries(FastANI PUBLIC ${GSL_LIBRARIES})
    else()
        # GSL_LIBDIR is not documented
        target_link_libraries(FastANI PUBLIC
                ${GSL_LIBDIR}/libgsl.a
                ${GSL_LIBDIR}/libgslcblas.a)
    endif(${GSL_SHARED})
    target_include_directories(FastANI PUBLIC ${GSL_INCLUDE_DIRS})
else()
    set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
    set(Boost_USE_DEBUG_LIBS        OFF)  # ignore debug libs and
    set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
    set(Boost_USE_MULTITHREADED      ON)
    find_package(Boost 1.45 REQUIRED COMPONENTS math_c99)
    target_link_libraries(FastANI PUBLIC ${Boost_MATCH_C99_LIBRARY})
    target_include_directories(FastANI PUBLIC ${Boost_INCLUDE_DIRS})
    target_compile_definitions(FastANI PUBLIC USE_BOOST=1)
endif(${GSL_FOUND})
#
#
if(${BUILD_TESTING})
    add_subdirectory(ext/Catch2)
    include(CTest)
    include(Catch)
    add_executable(FastANITest tests/fastani_tests.cpp src/cgi/core_genome_identity.cpp)
    target_include_directories(FastANITest PUBLIC "src/"
                            ${OpenMP_CXX_INCLUDE_DIRS})
    target_compile_options(FastANITest PUBLIC
            --coverage -g -O0 -fprofile-arcs -ftest-coverage ${OpenMP_CXX_FLAGS})
    target_link_libraries(FastANITest PRIVATE
                          Catch2::Catch2 Catch2::Catch2WithMain
                          ZLIB::ZLIB GSL::gsl GSL::gslcblas
                          ${OpenMP_CXX_LIBRARIES}
                          gcov)
    catch_discover_tests(FastANITest)
    file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    add_custom_target(lcov lcov -c -d .. -o FastANITest.out
                        COMMAND genhtml -o ../lcov/ FastANITest.out
                        DEPENDS FastANITest)
    add_custom_target(lcov2 lcov -c -d ./ -o FastANITest.info
                        COMMAND lcov -e FastANITest.info "*/FastANI/src/*/*.cpp" "*/FastANI/src/*/*.hpp"  > coverage.info
                        COMMAND lcov --list coverage.info |tee coverage_list.txt
                        COMMAND lcov --summary coverage.info |tee coverage_summary.txt
                        DEPENDS FastANITest)
endif(${BUILD_TESTING})

set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_DIRECTORY ${PROJECT_BINARY_DIR}/package)

set(CPACK_SOURCE_IGNORE_FILES .git/ .github/ 
    .vscode/ .cache/ .mypy_cache/
    .idea/ complie_comands.json
    lcov/ lcov0/ lcov2/ data/ cmake-build-debug/
    _CPack_Packages/ build.sh/ build/  build.clang/
    ${CMAKE_BINARY_DIR}/ ${PROJECT_BINARY_DIR}/)

include(CPack)
