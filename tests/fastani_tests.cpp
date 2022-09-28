
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

bool compareCGI(const cgi::CGI_Results& i1, const cgi::CGI_Results& i2)
{
    return i1.qryGenomeId == i2.qryGenomeId ? 
        i1.refGenomeId < i2.refGenomeId : i1.qryGenomeId < i2.qryGenomeId;
}

TEST_CASE( "Single Threaded Pair Query Ref", "[single threaded pair]" ) {
    skch::Parameters parameters;
    // -q [QUERY_GENOME] -r [REFERENCE_GENOME] -o [OUTPUT_FILE] 
    const char *argv[] =  {"single-pair", 
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

    //
    // Write output
    // name of genome -> length
    // std::unordered_map <std::string, uint64_t> genomeLengths;
    // cgi::computeGenomeLengths(parameters, genomeLengths);

    //report output in file
    //cgi::outputCGI (parameters, genomeLengths, finalResults, fileName);

    REQUIRE(finalResults.size() == 1);
    REQUIRE(finalResults[0].identity > 97.663);  
    REQUIRE(finalResults[0].identity < 97.665);
    REQUIRE(finalResults[0].countSeq == 1322);
    REQUIRE(finalResults[0].totalQueryFragments == 1547);
}


TEST_CASE( "Single Threaded Multi Query", "[single threaded multi query]" ) {
    // -q [QUERY_GENOME] -r [REFERENCE_GENOME] -o [OUTPUT_FILE] 
    const char *argv[] =  {"single-ref-multi-query", 
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
    // Write output
    // name of genome -> length
    //std::unordered_map <std::string, uint64_t> genomeLengths;
    //cgi::computeGenomeLengths(parameters, genomeLengths);
    //
    //report output in file
    //cgi::outputCGI (parameters, genomeLengths, finalResults, fileName);
    std::sort(finalResults.begin(), finalResults.end(), compareCGI);
    for(auto fx : finalResults){
        std::cout << fx.qryGenomeId << " " << fx.refGenomeId << " " << 
          fx.identity << " " << fx.countSeq << " " << fx.totalQueryFragments <<
          std::endl;
    }
    REQUIRE(finalResults.size() == 2);
    //data/2000031004.LargeContigs.fna
    //data/2000031001.LargeContigs.fna        99.9493 1708    1722
    REQUIRE(finalResults[0].identity > 99.9492);
    REQUIRE(finalResults[0].identity < 99.9494);
    REQUIRE(finalResults[0].countSeq == 1708);
    REQUIRE(finalResults[0].totalQueryFragments == 1722);
    //data/2000031006.LargeContigs.fna
    //data/2000031001.LargeContigs.fna        99.9973 1259    3517
    REQUIRE(finalResults[1].identity > 99.9972);
    REQUIRE(finalResults[1].identity < 99.9974);
    REQUIRE(finalResults[1].countSeq == 1259);
    REQUIRE(finalResults[1].totalQueryFragments == 3517);
 
}

TEST_CASE( "Single Threaded Multi Ref.", "[single threaded multi ref.]" ) {
    // -q [QUERY_GENOME] -r [REFERENCE_GENOME] -o [OUTPUT_FILE] 
    const char *argv[] =  {"t12", 
                    "-q", "data/D4/2000031001.LargeContigs.fna",
                    "--rl", "data/D4/multiref.txt",
                    "-o", "mtt-test-rq.txt"};
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
    //Final output vector of ANI computation
    // Write output
    // name of genome -> length
    //std::unordered_map <std::string, uint64_t> genomeLengths;
    //cgi::computeGenomeLengths(parameters, genomeLengths);
    //
    //report output in file
    //cgi::outputCGI (parameters, genomeLengths, finalResults, fileName);
    std::sort(finalResults.begin(), finalResults.end(), compareCGI);
    for(auto fx : finalResults){
        std::cout << fx.qryGenomeId << " " << fx.refGenomeId << " " << 
          fx.identity << " " << fx.countSeq << " " << fx.totalQueryFragments <<
          std::endl;
    }
    REQUIRE(finalResults.size() == 2);
    // data/D4/2000031001.LargeContigs.fna 
    // data/D4/2000031004.LargeContigs.fna     99.9759 1704    1711
    REQUIRE(finalResults[0].identity > 99.9758);
    REQUIRE(finalResults[0].identity < 99.9760);
    REQUIRE(finalResults[0].countSeq == 1704);
    REQUIRE(finalResults[0].totalQueryFragments == 1711);
    //data/D4/2000031001.LargeContigs.fna
    //data/D4/2000031006.LargeContigs.fna     99.9867 1703    1711
    REQUIRE(finalResults[1].identity > 99.9866);
    REQUIRE(finalResults[1].identity < 99.9868);
    REQUIRE(finalResults[1].countSeq == 1703);
    REQUIRE(finalResults[1].totalQueryFragments == 1711);
}


TEST_CASE( "Single Threaded Multi Q. Multi Ref.", "[single threaded multi q. multi ref.]" ) {
    // -q [QUERY_GENOME] -r [REFERENCE_GENOME] -o [OUTPUT_FILE] 
    const char *argv[] =  {"multiq-multi-ref", 
                    "--ql", "data/D4/multiq2.txt",
                    "--rl", "data/D4/multiref2.txt",
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
    //Final output vector of ANI computation
    // Write output
    // name of genome -> length
    // std::unordered_map <std::string, uint64_t> genomeLengths;
    // cgi::computeGenomeLengths(parameters, genomeLengths);
    //
    //report output in file
    // cgi::outputCGI (parameters, genomeLengths, finalResults, fileName);
    std::sort(finalResults.begin(), finalResults.end(), compareCGI);
    for(auto fx : finalResults){
        std::cout << fx.qryGenomeId << " " << fx.refGenomeId << " " << 
          fx.identity << " " << fx.countSeq << " " << fx.totalQueryFragments <<
          std::endl;
    }
    REQUIRE(finalResults.size() == 6);
    // data/D4/2000031001.LargeContigs.fna
    // data/D4/2000031008.LargeContigs.fna     99.9785 1700    1711
    REQUIRE(finalResults[0].identity > 99.9784);
    REQUIRE(finalResults[0].identity < 99.9786);
    REQUIRE(finalResults[0].countSeq == 1700);
    REQUIRE(finalResults[0].totalQueryFragments == 1711);
    // data/D4/2000031001.LargeContigs.fna
    // data/D4/2000031009.LargeContigs.fna     99.9795 1704    1711
    REQUIRE(finalResults[1].identity > 99.9794);
    REQUIRE(finalResults[1].identity < 99.9796);
    REQUIRE(finalResults[1].countSeq == 1704);
    REQUIRE(finalResults[1].totalQueryFragments == 1711);
    // data/D4/2000031004.LargeContigs.fna 
    // data/D4/2000031008.LargeContigs.fna     99.9943 1223    3433
    REQUIRE(finalResults[2].identity > 99.9942);
    REQUIRE(finalResults[2].identity < 99.9944);
    REQUIRE(finalResults[2].countSeq == 1223);
    REQUIRE(finalResults[2].totalQueryFragments == 3433);
    // data/D4/2000031004.LargeContigs.fna
    // data/D4/2000031009.LargeContigs.fna     99.9958 1248    3433
    REQUIRE(finalResults[3].identity > 99.9957);
    REQUIRE(finalResults[3].identity < 99.9959);
    REQUIRE(finalResults[3].countSeq == 1248);
    REQUIRE(finalResults[3].totalQueryFragments == 3433);
    // data/D4/2000031006.LargeContigs.fna 
    // data/D4/2000031008.LargeContigs.fna     99.9967 1110    5228
    REQUIRE(finalResults[4].identity > 99.9966);
    REQUIRE(finalResults[4].identity < 99.9968);
    REQUIRE(finalResults[4].countSeq == 1110);
    REQUIRE(finalResults[4].totalQueryFragments == 5228);
    // data/D4/2000031006.LargeContigs.fna
    // data/D4/2000031009.LargeContigs.fna     99.9977 1094    5228
    REQUIRE(finalResults[5].identity > 99.9976);
    REQUIRE(finalResults[5].identity < 99.9978);
    REQUIRE(finalResults[5].countSeq == 1094);
    REQUIRE(finalResults[5].totalQueryFragments == 5228); 
}
