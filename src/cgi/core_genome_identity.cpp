/**
 * @file    core_genome_identity.cpp
 * @ingroup src
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>
#include <omp.h>

//Own includes
#include "map/include/map_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/parseCmdArgs.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"
#include "map/include/commonFunc.hpp"
#include "cgi/include/computeCoreIdentity.hpp" 

int main(int argc, char** argv)
{
  /*
   * Make sure env variable MALLOC_ARENA_MAX is unset 
   * for efficient multi-threaded execution
   */
  unsetenv((char *)"MALLOC_ARENA_MAX");
  using namespace std::placeholders;  // for _1, _2, _3...

  //Parse command line arguements   
  skch::Parameters parameters;        //sketching and mapping parameters
  skch::parseandSave(argc, argv, parameters);

  std::string fileName = parameters.outFileName;

  //To redirect Mashmap's mapping output to null fs, using file name for CGI output
  parameters.outFileName = "/dev/null";

  //Set up for parallel execution
  omp_set_num_threads( parameters.threads ); 
  std::vector <skch::Parameters> parameters_split (parameters.threads);
  cgi::splitReferenceGenomes (parameters, parameters_split);

  //Final output vector of ANI computation
  std::vector<cgi::CGI_Results> finalResults;

#pragma omp parallel for schedule(static,1)
  for (uint64_t i = 0; i < parameters.threads; i++)
  {
    if ( omp_get_thread_num() == 0)
      std::cerr << "INFO [thread 0], skch::main, Count of threads executing parallel_for : " << omp_get_num_threads() << std::endl;

    //start timer
    auto t0 = skch::Time::now();

    //Build the sketch for reference
    skch::Sketch referSketch(parameters_split[i]);

    std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;

    if ( omp_get_thread_num() == 0)
      std::cerr << "INFO [thread 0], skch::main, Time spent sketching the reference : " << timeRefSketch.count() << " sec" << std::endl;

    //Final output vector of ANI computation
    std::vector<cgi::CGI_Results> finalResults_local;

    //Loop over query genomes
    for(uint64_t queryno = 0; queryno < parameters_split[i].querySequences.size(); queryno++)
    {
      t0 = skch::Time::now();

      skch::MappingResultsVector_t mapResults;
      uint64_t totalQueryFragments = 0;

      auto fn = std::bind(skch::Map::insertL2ResultsToVec, std::ref(mapResults), _1);
      skch::Map mapper = skch::Map(parameters_split[i], referSketch, totalQueryFragments, queryno, fn);

      std::chrono::duration<double> timeMapQuery = skch::Time::now() - t0;

      if ( omp_get_thread_num() == 0)
        std::cerr << "INFO [thread 0], skch::main, Time spent mapping fragments in query #" << queryno + 1 <<  " : " << timeMapQuery.count() << " sec" << std::endl;

      t0 = skch::Time::now();

      cgi::computeCGI(parameters_split[i], mapResults, mapper, referSketch, totalQueryFragments, queryno, fileName, finalResults_local);

      std::chrono::duration<double> timeCGI = skch::Time::now() - t0;

      if ( omp_get_thread_num() == 0)
        std::cerr << "INFO [thread 0], skch::main, Time spent post mapping : " << timeCGI.count() << " sec" << std::endl;
    }

    cgi::correctRefGenomeIds (finalResults_local);

#pragma omp critical
    {
      finalResults.insert (finalResults.end(), finalResults_local.begin(), finalResults_local.end());
    }

#pragma omp critical
    {
      std::cerr << "INFO [thread " << omp_get_thread_num() << "], skch::main, ready to exit the loop" << std::endl;
    }
  }

  std::cerr << "INFO, skch::main, parallel_for execution finished" << std::endl;

  std::unordered_map <std::string, uint64_t> genomeLengths;    // name of genome -> length
  cgi::computeGenomeLengths(parameters, genomeLengths);

  //report output in file
  cgi::outputCGI (parameters, genomeLengths, finalResults, fileName);

  //report output as matrix
  if (parameters.matrixOutput)
    cgi::outputPhylip (parameters, genomeLengths, finalResults, fileName);
}
