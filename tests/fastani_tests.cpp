
#include <functional>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include "map/include/map_parameters.hpp"
#include "map/include/parseCmdArgs.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"
#include "cgi/include/computeCoreIdentity.hpp" 

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Single Threaded Pair Test", "[single threaded pair]" ) {
    skch::Parameters parameters;
    // -q [QUERY_GENOME] -r [REFERENCE_GENOME] -o [OUTPUT_FILE] 
    const char *argv[] =  {"t1", 
                    "-q", "data/Escherichia_coli_str_K12_MG1655.fna",
                    "-r", "data/Shigella_flexneri_2a_01.fna",
                    "-o", "stt-test.txt"};
    skch::parseandSave(7, (char**)argv, parameters);
    std::string fileName = parameters.outFileName;
    //
    //Build the sketch for reference
    skch::Sketch referSketch(parameters);
    //
    // Map Results
    skch::MappingResultsVector_t mapResults;
    uint64_t totalQueryFragments = 0;
    auto queryno = 0;
    auto fn = std::bind(skch::Map::insertL2ResultsToVec, 
                        std::ref(mapResults), std::placeholders::_1);
    skch::Map mapper = skch::Map(parameters, referSketch,
                                 totalQueryFragments, queryno, fn);
    //
    //
    //Final output vector of ANI computation
    std::vector<cgi::CGI_Results> finalResults;
    cgi::computeCGI(parameters, mapResults, mapper, referSketch,
                    totalQueryFragments, queryno, fileName, finalResults);

    REQUIRE(finalResults.size() == 1);
    REQUIRE(finalResults[0].identity > 97.663);  
    REQUIRE(finalResults[0].identity < 97.665);
    REQUIRE(finalResults[0].countSeq == 1322);
    REQUIRE(finalResults[0].totalQueryFragments == 1547);
    //
    // Write output
    // name of genome -> length
    // std::unordered_map <std::string, uint64_t> genomeLengths;
    // cgi::computeGenomeLengths(parameters, genomeLengths);

    //report output in file
    //cgi::outputCGI (parameters, genomeLengths, finalResults, fileName);

    // REQUIRE( Factorial(1) == 1 );
    // REQUIRE( Factorial(2) == 2 );
    // REQUIRE( Factorial(3) == 6 );
    // REQUIRE( Factorial(10) == 3628800 );
}


TEST_CASE( "Single Threaded Multi Test", "[single threaded multi]" ) {
    // -q [QUERY_GENOME] -r [REFERENCE_GENOME] -o [OUTPUT_FILE] 
    const char *argv[] =  {"t12", 
                    "--ql", "data/D4/multiq.txt",
                    "-r", "data/D4/2000031001.LargeContigs.fna",
                    "-o", "mtt-test.txt"};
    //
    skch::Parameters parameters;
    skch::parseandSave(7, (char**)argv, parameters);
    std::string fileName = parameters.outFileName;
    //
    //Build the sketch for reference
    skch::Sketch referSketch(parameters);
    //
    // Map Results
    skch::MappingResultsVector_t mapResults;
    uint64_t totalQueryFragments = 0;
    std::vector<cgi::CGI_Results> finalResults;
    for(uint64_t queryno = 0; queryno < parameters.querySequences.size(); queryno++) {
        auto fn = std::bind(skch::Map::insertL2ResultsToVec, 
                            std::ref(mapResults), std::placeholders::_1);
        skch::Map mapper = skch::Map(parameters, referSketch,
                                    totalQueryFragments, queryno, fn);

        cgi::computeCGI(parameters, mapResults, mapper, referSketch,
                        totalQueryFragments, queryno, fileName, finalResults);
    }
    //
    //
    //Final output vector of ANI computation

    // REQUIRE(finalResults.size() == 1);
    // REQUIRE(finalResults[0].identity > 97.663);  
    // REQUIRE(finalResults[0].identity < 97.665);
    // REQUIRE(finalResults[0].countSeq == 1322);
    // REQUIRE(finalResults[0].totalQueryFragments == 1547);
    //
    // Write output
    // name of genome -> length
    std::unordered_map <std::string, uint64_t> genomeLengths;
    cgi::computeGenomeLengths(parameters, genomeLengths);

    //report output in file
    //cgi::outputCGI (parameters, genomeLengths, finalResults, fileName);
    REQUIRE(finalResults.size() == 2);
    //data/2000031004.LargeContigs.fna        data/2000031001.LargeContigs.fna        99.9493 1708    1722
    REQUIRE(finalResults[0].identity > 99.9492);
    REQUIRE(finalResults[0].identity < 99.9494);
    REQUIRE(finalResults[0].countSeq == 1708);
    REQUIRE(finalResults[0].totalQueryFragments == 1722);
    //data/2000031006.LargeContigs.fna        data/2000031001.LargeContigs.fna        99.9973 1259    3517
    REQUIRE(finalResults[1].identity > 99.9972);
    REQUIRE(finalResults[1].identity < 99.9974);
    REQUIRE(finalResults[1].countSeq == 1259);
    REQUIRE(finalResults[1].totalQueryFragments == 3517);
 
    //REQUIRE( Factorial(1) == 1 );
    //REQUIRE( Factorial(2) == 2 );
    // REQUIRE( Factorial(3) == 6 );
    // REQUIRE( Factorial(10) == 3628800 );
}
