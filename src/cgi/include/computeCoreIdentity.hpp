/**
 * @file    computeCoreIdentity.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef CGI_IDENTITY_HPP 
#define CGI_IDENTITY_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>

//Own includes
#include "base_types.hpp"
#include "cgid_types.hpp"

//External includes
#include "prettyprint.hpp"

namespace cgi
{
  /**
   * @brief                       Use reference sketch's sequence to file (genome) mapping 
   *                              and revise reference ids to genome id
   * @param[in/out] shortResults
   */
  void reviseRefIdToGenomeId(std::vector<MappingResult_CGI> &shortResults, skch::Sketch &refSketch)
  {
    for(auto &e : shortResults)
    {
      auto referenceSequenceId = e.genomeId;
      auto upperRangeIter = std::upper_bound(refSketch.sequencesByFileInfo.begin(), 
          refSketch.sequencesByFileInfo.end(),
          referenceSequenceId);

      auto genomeId = std::distance(refSketch.sequencesByFileInfo.begin(), upperRangeIter);
      e.genomeId = genomeId;
    }
  }

  /**
   * @brief                 compute and report AAI/ANI 
   * @param[in] parameters  algorithm parameters
   * @param[in] results     mapping results
   * @param[in] refSketch   reference sketch
   * @param[in] fileName    file name where results will be reported
   */
  void computeCGI(skch::Parameters &parameters,
      skch::MappingResultsVector_t &results,
      skch::Sketch &refSketch,
      std::string &fileName)
  {

    //Vector to save relevant fields from mapping results
    std::vector<MappingResult_CGI> shortResults;

    shortResults.reserve(results.size());

    ///Parse map results and save fields which we need
    // reference id (R), query id (Q), estimated identity (I)
    for(auto &e : results)
    {
      shortResults.push_back(MappingResult_CGI{
          e.refSeqId,
          e.querySeqId,
          e.nucIdentity
          });
    }

    //Sort the vector shortResults
    std::sort(shortResults.begin(), shortResults.end());

    /*
     * NOTE: We assume single file contains the sequences for single genome
     * We revise reference sequence id to genome (or file) id
     */
    reviseRefIdToGenomeId(shortResults, refSketch);

    //We need best identity match for each genome, query pair
    std::vector<MappingResult_CGI> singleQueryRefResults;

    //Code below fetches best identity match for each genome, query pair
    for(auto &e : shortResults)
    {
      if(singleQueryRefResults.empty())
        singleQueryRefResults.push_back(e);

      else if ( !(
            e.genomeId == singleQueryRefResults.back().genomeId && 
            e.querySeqId == singleQueryRefResults.back().querySeqId))
        singleQueryRefResults.push_back(e);

      else
        singleQueryRefResults.back().nucIdentity = e.nucIdentity;
    }

    //Final output vector of ANI/AAI computation
    std::vector<cgi::CGI_Results> CGI_ResultsVector;


    //Do average for ANI/AAI computation 
    for(auto it = singleQueryRefResults.begin(); it != singleQueryRefResults.end();)
    {
      skch::seqno_t currentRefId = it->genomeId;

      //Bucket by genome id
      auto rangeEndIter = std::find_if(it, singleQueryRefResults.end(), [&](const MappingResult_CGI& e) 
          { 
            return e.genomeId != currentRefId; 
          } );

      float sumIdentity = 0.0;

      for(auto it2 = it; it2 != rangeEndIter; it2++)
      {
        sumIdentity += it2->nucIdentity;
      }

      //Save the result 
      CGI_Results currentResult;
      currentResult.genomeId = currentRefId;
      currentResult.countSeq = std::distance(it, rangeEndIter);
      currentResult.identity = sumIdentity/currentResult.countSeq;

      CGI_ResultsVector.push_back(currentResult);

      //Advance the iterator it
      it = rangeEndIter;
    }

    //Sort in decreasing order of matches
    std::sort(CGI_ResultsVector.rbegin(), CGI_ResultsVector.rend());

    std::ofstream outstrm(fileName);

    //Report results
    for(auto &e : CGI_ResultsVector)
    {
      outstrm << parameters.refSequences[e.genomeId]
        << " " << e.identity 
        << " " << e.countSeq
        << "\n";
    }

    std::ofstream outstrm2(fileName + ".map.best");

    //Report all mappings that contribute to core-genome identity estimate
    for(auto &e : singleQueryRefResults)
    {
      outstrm2 << parameters.refSequences[e.genomeId]
        << " " << e.querySeqId 
        << " " << e.nucIdentity
        << "\n";
    }
  }

}

#endif
