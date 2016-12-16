/**
 * @file    map_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_CONFIG_HPP 
#define SKETCH_CONFIG_HPP

#include <vector>

//Switch to enable timing of L1 and L2 stages for each read
//Timings are reported in a file
#define ENABLE_TIME_PROFILE_L1_L2 0

namespace skch
{
  /**
   * @brief   configuration parameters for building sketch
   *          expected to be initialized using command line arguments
   */
  struct Parameters
  {
    int kmerSize;                                     //kmer size for sketching
    int windowSize;                                   //window size used for sketching 
    int minReadLength;                                //minimum read length which code maps
    int alphabetSize;                                 //alphabet size
    uint64_t referenceSize;                           //Approximate reference size
    float percentageIdentity;                         //user defined threshold for good similarity
    double p_value;                                   //user defined threshold for p value
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string outFileName;                          //output file name
    bool reportAll;                                   //Report all alignments if this is true
    //const static int L2slideJump = 5;                 //parameter to jump read during L2 evaluation, 
                                                      //  1 is the most accurate

  };
}

#endif
