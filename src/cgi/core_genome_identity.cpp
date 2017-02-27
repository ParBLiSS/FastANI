/**
 * @file    core_genome_identity.cpp
 * @ingroup src
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>

//Own includes
#include "map_parameters.hpp"
#include "base_types.hpp"
#include "parseCmdArgs.hpp"
#include "winSketch.hpp"
#include "computeMap.hpp"
#include "computeCoreIdentity.hpp" 

//External includes
#include "argvparser.hpp"

int main(int argc, char** argv)
{
  CommandLineProcessing::ArgvParser cmd;
  using namespace std::placeholders;  // for _1, _2, _3...

  //Setup command line options
  skch::initCmdParser(cmd);

  //Parse command line arguements   
  skch::Parameters parameters;        //sketching and mapping parameters

  skch::parseandSave(argc, argv, cmd, parameters);   

  //Redirect mapping output to null fs, using file name for CGI output
  std::string fileName = parameters.outFileName;

#ifdef DEBUG
  parameters.outFileName = parameters.outFileName + ".map";
#else
  parameters.outFileName = "/dev/null";
#endif

  auto t0 = skch::Time::now();

  //Build the sketch for reference
  skch::Sketch referSketch(parameters);

  std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;
  std::cerr << "INFO, skch::main, Time spent sketching the reference : " << timeRefSketch.count() << " sec" << std::endl;

  //Map the sequences in query file
  t0 = skch::Time::now();

  skch::MappingResultsVector_t mapResults;
  uint64_t totalQueryFragments;

  auto fn = std::bind(skch::Map::insertL2ResultsToVec, std::ref(mapResults), _1);
  skch::Map mapper = skch::Map(parameters, referSketch, totalQueryFragments, fn);

  std::chrono::duration<double> timeMapQuery = skch::Time::now() - t0;
  std::cerr << "INFO, skch::main, Time spent mapping the query : " << timeMapQuery.count() << " sec" << std::endl;

  std::cerr << "INFO, skch::main, mapping results saved in : " << parameters.outFileName << std::endl;

  t0 = skch::Time::now();

  cgi::computeCGI(parameters, mapResults, referSketch, totalQueryFragments, fileName);

  std::chrono::duration<double> timeCGI = skch::Time::now() - t0;
  std::cerr << "INFO, skch::main, Time spent post mapping : " << timeCGI.count() << " sec" << std::endl;
}
