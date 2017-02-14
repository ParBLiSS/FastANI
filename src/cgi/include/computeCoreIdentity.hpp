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
      auto referenceSequenceId = e.refSequenceId;
      auto upperRangeIter = std::upper_bound(refSketch.sequencesByFileInfo.begin(), 
          refSketch.sequencesByFileInfo.end(),
          referenceSequenceId);

      e.genomeId = std::distance(refSketch.sequencesByFileInfo.begin(), upperRangeIter);
    }
  }

  /**
   * @brief       compute mean of the nucleotide identites while discarding outliers
   * @details     Discard values if they are outside [1st quartile - 1.5*IQR, 3rd 
   *              quartile + 1.5*IQR] range
   * @param[in]   mappings_2way_beg
   * @param[in]   mappings_2way_end
   * @param[out]  sumIdentity
   * @param[out]  fragment_count
   */
  template <typename Iter>
    void tuckeyCorrection(Iter &mappings_2way_beg, Iter &mappings_2way_end, 
        float &sumIdentity, int &fragment_count)
    {
      std::sort(mappings_2way_beg, mappings_2way_end, cmp_identity);

      auto totalElementCount = std::distance(mappings_2way_beg, mappings_2way_end);
      auto Q1 = (mappings_2way_beg + totalElementCount / 4)->nucIdentity;
      auto Q3 = (mappings_2way_beg + totalElementCount * 3/4)->nucIdentity;
      auto IPR = Q3 - Q1; 

      auto allowedRange = std::make_pair(Q1 - 1.5 * IPR, Q3 + 1.5 * IPR); 

      fragment_count= 0;
      sumIdentity = 0;

      for(auto it = mappings_2way_beg; it != mappings_2way_end; it++)
      {
        if(it->nucIdentity >= allowedRange.first && it->nucIdentity <= allowedRange.second)
        {
          sumIdentity += it->nucIdentity;
          fragment_count += 1;
        }
        else
        {
          //INVALIDATE 
          it->nucIdentity = 0;
        }
      }
    }

  /**
   * @brief                       Compute N50 statistics and total bp length of all reference genomes
   * @param[out]  refLenStats
   * @param[in]   refSketch
   */
  template <typename VEC>
    void computeRefLenStatistics(VEC &refLenStats, 
        skch::Sketch &refSketch) 
    {
      skch::seqno_t seqBeginOffset = 0;
      for(auto &e : refSketch.sequencesByFileInfo)
      {
        skch::seqno_t seqEndOffset = e;

        //Compute length statistics of sequences saved in metadata [seqBeginOffset, seqEndOffset)
        {
          uint64_t length, N50;

          std::vector<uint64_t> lens;
          for(skch::seqno_t i = seqBeginOffset; i < seqEndOffset; i++)
          {
            lens.push_back(refSketch.metadata[i].len);
          }

          //Sum of length of all the contigs in this genome
          length = std::accumulate(lens.begin(), lens.end(), 0, plus<uint64_t>());

          //Sort length values in descending order
          std::sort(lens.rbegin(), lens.rend());

          uint64_t partialSum = 0;
          for(auto &f : lens)
          {
            partialSum += f;

            if(partialSum >= length/2)
            {
              N50 = f;
              break;
            }
          }

          refLenStats.emplace_back(N50, length);
        }//Next reference genome

        seqBeginOffset = e;
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
      shortResults.emplace_back(MappingResult_CGI{
          e.refSeqId,
          0,                  //this value will be revised to genome id
          e.querySeqId,
          e.refStartPos/(parameters.minReadLength),
          e.nucIdentity
          });
    }

    /*
     * NOTE: We assume single file contains the sequences for single genome
     * We revise reference sequence id to genome (or file) id
     */
    reviseRefIdToGenomeId(shortResults, refSketch);

    //Compute assembly statistics of all the reference genomes
    std::vector<std::pair< uint64_t, uint64_t >> refLenStats;
    computeRefLenStatistics(refLenStats, refSketch); 

    //We need best reciprocal identity match for each genome, query pair
    std::vector<MappingResult_CGI> mappings_1way;
    std::vector<MappingResult_CGI> mappings_2way;

    ///1. Code below fetches best identity match for each genome, query pair
    //For each query sequence, best match in the reference is preserved
    {
      //Sort the vector shortResults
      std::sort(shortResults.begin(), shortResults.end(), cmp_query_bucket);

      for(auto &e : shortResults)
      {
        if(mappings_1way.empty())
          mappings_1way.push_back(e);

        else if ( !(
              e.genomeId == mappings_1way.back().genomeId && 
              e.querySeqId == mappings_1way.back().querySeqId))
        {
          mappings_1way.emplace_back(e);
        }
        else
        {
          mappings_1way.back() = e;
        }
      }
    }

    ///2. Now, we compute 2-way ANI
    //For each mapped region, and within a reference bin bucket, single best query mapping is preserved
    {
      std::sort(mappings_1way.begin(), mappings_1way.end(), cmp_refbin_bucket);

      for(auto &e : mappings_1way)
      {
        if(mappings_2way.empty())
          mappings_2way.push_back(e);

        else if ( !(
              e.refSequenceId == mappings_2way.back().refSequenceId && 
              e.mapRefPosBin == mappings_2way.back().mapRefPosBin))
        {
          mappings_2way.emplace_back(e);
        }
        else
        {
          mappings_2way.back() = e;
        }
      }
    }

    //Final output vector of ANI/AAI computation
    std::vector<cgi::CGI_Results> CGI_ResultsVector;


    //Do average for ANI/AAI computation 
    //mappings_2way should be sorted by genomeId 

    for(auto it = mappings_2way.begin(); it != mappings_2way.end();)
    {
      skch::seqno_t currentGenomeId = it->genomeId;

      //Bucket by genome id
      auto rangeEndIter = std::find_if(it, mappings_2way.end(), [&](const MappingResult_CGI& e) 
          { 
            return e.genomeId != currentGenomeId; 
          } );

      float sumIdentity;
      int countFrag;

      tuckeyCorrection(it, rangeEndIter, sumIdentity, countFrag);

      //Save the result 
      CGI_Results currentResult;
      currentResult.genomeId = currentGenomeId;
      currentResult.countSeq = countFrag;
      currentResult.identity = sumIdentity/countFrag;

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
        << " " << refLenStats[e.genomeId].first
        << " " << refLenStats[e.genomeId].second
        << "\n";
    }

#ifdef DEBUG
    std::ofstream outstrm2(fileName + ".map.best");

    //Report all mappings that contribute to core-genome identity estimate
    for(auto &e : mappings_2way)
    {
      if(e.nucIdentity != 0) 
        outstrm2 << parameters.refSequences[e.genomeId]
          << " " << e.querySeqId 
          << " " << e.mapRefPosBin
          << " " << e.nucIdentity
          << "\n";
    }
#endif

  }

}

#endif
