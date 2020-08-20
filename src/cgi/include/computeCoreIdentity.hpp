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
#include <omp.h>
#include <zlib.h>  

//Own includes
#include "map/include/base_types.hpp"
#include "cgi/include/cgid_types.hpp"

//External includes
#include "common/kseq.h"
#include "common/prettyprint.hpp"

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
   * @brief                       compute genome lengths in reference and query genome set
   * @param[out] genomeLengths
   */
  void computeGenomeLengths(skch::Parameters &parameters, std::unordered_map <std::string, uint64_t> &genomeLengths) 
  { 
    for(auto &e : parameters.querySequences)
    {
      //Open the file using kseq
      gzFile fp = gzopen(e.c_str(), "r");
      kseq_t *seq = kseq_init(fp);
      int l; uint64_t genomeLen = 0;

      while ((l = kseq_read(seq)) >= 0) {
        if (l >= parameters.minReadLength) {
          uint64_t _l_ = (((uint64_t)strlen(seq->seq.s)) / parameters.minReadLength) * parameters.minReadLength;
          genomeLen = genomeLen + _l_;
        }
      }

      genomeLengths[e] = genomeLen;

      kseq_destroy(seq);  
      gzclose(fp); //close the file handler 
    }

    for(auto &e : parameters.refSequences)
    {
      if( genomeLengths.find(e) == genomeLengths.end() )
      {
        //Open the file using kseq
        gzFile fp = gzopen(e.c_str(), "r");
        kseq_t *seq = kseq_init(fp);
        int l; uint64_t genomeLen = 0;

      while ((l = kseq_read(seq)) >= 0) {
        if (l >= parameters.minReadLength) {
          uint64_t _l_ = (((uint64_t)strlen(seq->seq.s)) / parameters.minReadLength) * parameters.minReadLength;
          genomeLen = genomeLen + _l_;
        }
      }

        genomeLengths[e] = genomeLen;

        kseq_destroy(seq);  
        gzclose(fp); //close the file handler 
      }
    }
  }

  /**
   * @brief                             output blast tabular mappings for visualization 
   * @param[in]   parameters            algorithm parameters
   * @param[in]   results               bidirectional mappings
   * @param[in]   mapper                mapper object used for mapping
   * @param[in]   refSketch             reference sketch
   * @param[in]   queryFileNo           query genome is parameters.querySequences[queryFileNo]
   * @param[in]   fileName              file name where results will be reported
   */
  void outputVisualizationFile(skch::Parameters &parameters,
      std::vector<MappingResult_CGI> &mappings_2way,
      skch::Map &mapper,
      skch::Sketch &refSketch,
      uint64_t queryFileNo,
      std::string &fileName)
  {
    std::ofstream outstrm(fileName + ".visual", std::ios::app);

    //Shift offsets for converting from local (to contig) to global (to genome)
    std::vector <skch::offset_t> queryOffsetAdder (mapper.metadata.size());
    std::vector <skch::offset_t> refOffsetAdder (refSketch.metadata.size());

    for(int i = 0; i < mapper.metadata.size(); i++)
    {
      if(i == 0)
        queryOffsetAdder[i] = 0;
      else
        queryOffsetAdder[i] = queryOffsetAdder[i-1] + mapper.metadata[i-1].len;
    }

    for(int i = 0; i < refSketch.metadata.size(); i++)
    {
      if(i == 0)
        refOffsetAdder[i] = 0;
      else
        refOffsetAdder[i] = refOffsetAdder[i-1] + refSketch.metadata[i-1].len;
    }

    //Report all mappings that contribute to core-genome identity estimate
    //Format the output to blast tabular way (outfmt 6)
    for(auto &e : mappings_2way)
    {
      outstrm << parameters.querySequences[queryFileNo]
        << "\t" << parameters.refSequences[e.genomeId]
        << "\t" << e.nucIdentity
        << "\t" << "NA"
        << "\t" << "NA"
        << "\t" << "NA"
        << "\t" << e.queryStartPos                                    + queryOffsetAdder[e.querySeqId] 
        << "\t" << e.queryStartPos + parameters.minReadLength - 1     + queryOffsetAdder[e.querySeqId]
        << "\t" << e.refStartPos                                      + refOffsetAdder[e.refSequenceId]
        << "\t" << e.refStartPos + parameters.minReadLength - 1       + refOffsetAdder[e.refSequenceId]
        << "\t" << "NA"
        << "\t" << "NA"
        << "\n";
    }
  }

  /**
   * @brief                             compute and report AAI/ANI 
   * @param[in]   parameters            algorithm parameters
   * @param[in]   results               mapping results
   * @param[in]   mapper                mapper object used for mapping
   * @param[in]   refSketch             reference sketch
   * @param[in]   totalQueryFragments   count of total sequence fragments in query genome
   * @param[in]   queryFileNo           query genome is parameters.querySequences[queryFileNo]
   * @param[in]   fileName              file name where results will be reported
   * @param[out]  CGI_ResultsVector     FastANI results
   */
  void computeCGI(skch::Parameters &parameters,
      skch::MappingResultsVector_t &results,
      skch::Map &mapper,
      skch::Sketch &refSketch,
      uint64_t totalQueryFragments,
      uint64_t queryFileNo,
      std::string &fileName,
      std::vector<cgi::CGI_Results> &CGI_ResultsVector
      )
  {
    //Vector to save relevant fields from mapping results
    std::vector<MappingResult_CGI> shortResults;

    shortResults.reserve(results.size());

    // Note to self: For debugging any issue, it is often useful to print
    // shortResults, mappings_1way and mappings_2way vectors

    ///Parse map results and save fields which we need
    // reference id (R), query id (Q), estimated identity (I)
    for(auto &e : results)
    {
      shortResults.emplace_back(MappingResult_CGI{
          e.refSeqId,
          0,                  //this value will be revised to genome id
          e.querySeqId,
          e.refStartPos,
          e.queryStartPos,
          e.refStartPos/(parameters.minReadLength - 20),
          e.nucIdentity
          });
    }

    /*
     * NOTE: We assume single file contains the sequences for single genome
     * We revise reference sequence id to genome (or file) id
     */
    reviseRefIdToGenomeId(shortResults, refSketch);


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

    {
      if(parameters.visualize)
      {
        outputVisualizationFile(parameters, mappings_2way, mapper, refSketch, queryFileNo, fileName);
      }
    }

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

      float sumIdentity = 0.0;

      for(auto it2 = it; it2 != rangeEndIter; it2++)
      {
        sumIdentity += it2->nucIdentity;
      }

      //Save the result 
      CGI_Results currentResult;

      currentResult.qryGenomeId = queryFileNo;
      currentResult.refGenomeId = currentGenomeId;
      currentResult.countSeq = std::distance(it, rangeEndIter);
      currentResult.totalQueryFragments = totalQueryFragments;
      currentResult.identity = sumIdentity/currentResult.countSeq;

      CGI_ResultsVector.push_back(currentResult);

      //Advance the iterator it
      it = rangeEndIter;
    }
  }

  /**
   * @brief                             output FastANI results to file
   * @param[in]   parameters            algorithm parameters
   * @param[in]   genomeLengths
   * @param[in]   CGI_ResultsVector     results
   * @param[in]   fileName              file name where results will be reported
   */
  void outputCGI(skch::Parameters &parameters,
      std::unordered_map <std::string, uint64_t> &genomeLengths,
      std::vector<cgi::CGI_Results> &CGI_ResultsVector,
      std::string &fileName)
  {
    //sort result by identity
    std::sort(CGI_ResultsVector.rbegin(), CGI_ResultsVector.rend());

    std::ofstream outstrm(fileName);

    //Report results
    for(auto &e : CGI_ResultsVector)
    {
      std::string qryGenome = parameters.querySequences[e.qryGenomeId];
      std::string refGenome = parameters.refSequences[e.refGenomeId];

      assert(genomeLengths.find(qryGenome) != genomeLengths.end());
      assert(genomeLengths.find(refGenome) != genomeLengths.end());

      uint64_t queryGenomeLength = genomeLengths[qryGenome];
      uint64_t refGenomeLength = genomeLengths[refGenome]; 
      uint64_t minGenomeLength = std::min(queryGenomeLength, refGenomeLength);
      uint64_t sharedLength = e.countSeq * parameters.minReadLength;

      //Checking if shared genome is above a certain fraction of genome length
      if(sharedLength >= minGenomeLength * parameters.minFraction)
      {
        outstrm << qryGenome
          << "\t" << refGenome
          << "\t" << e.identity 
          << "\t" << e.countSeq
          << "\t" << e.totalQueryFragments
          << "\n";
      }
    }

    outstrm.close();
  }

  /**
   * @brief                             output FastANI results as lower triangular matrix
   * @param[in]   parameters            algorithm parameters
   * @param[in]   genomeLengths
   * @param[in]   CGI_ResultsVector     results
   * @param[in]   fileName              file name where results will be reported
   */
  void outputPhylip(skch::Parameters &parameters,
      std::unordered_map <std::string, uint64_t> &genomeLengths,
      std::vector<cgi::CGI_Results> &CGI_ResultsVector,
      std::string &fileName)
  {
    std::unordered_map <std::string, int> genome2Int;    // name of genome -> integer
    std::unordered_map <int, std::string> genome2Int_rev; // integer -> name of genome

    //Assign unique index to the set of query and reference genomes
    for(auto &e : parameters.querySequences)
    {
      auto id = genome2Int.size();
      if( genome2Int.find(e) == genome2Int.end() )
      {
        genome2Int [e] = id;
        genome2Int_rev [id] = e;
      }
    }

    for(auto &e : parameters.refSequences)
    {
      auto id = genome2Int.size();
      if( genome2Int.find(e) == genome2Int.end() )
      {
        genome2Int [e] = id;
        genome2Int_rev [id] = e;
      }
    }

    int totalGenomes = genome2Int.size();

    //create a square 2-d matrix
    std::vector< std::vector<float> > fastANI_matrix (totalGenomes,  std::vector<float> (totalGenomes, 0.0));

    //transform FastANI results into 3-tuples
    for(auto &e : CGI_ResultsVector)
    {
      std::string qryGenome = parameters.querySequences[e.qryGenomeId];
      std::string refGenome = parameters.refSequences[e.refGenomeId];

      assert(genomeLengths.find(qryGenome) != genomeLengths.end());
      assert(genomeLengths.find(refGenome) != genomeLengths.end());

      uint64_t queryGenomeLength = genomeLengths[qryGenome];
      uint64_t refGenomeLength = genomeLengths[refGenome]; 
      uint64_t minGenomeLength = std::min(queryGenomeLength, refGenomeLength);
      uint64_t sharedLength = e.countSeq * parameters.minReadLength;

      //Checking if shared genome is above a certain fraction of genome length
      if(sharedLength >= minGenomeLength * parameters.minFraction)
      {
        int qGenome = genome2Int [ qryGenome ];
        int rGenome = genome2Int [ refGenome ];

        if (qGenome != rGenome)   //ignore if both genomes are same
        {
          if (qGenome > rGenome)
          {
            if (fastANI_matrix[qGenome][rGenome] > 0)
              fastANI_matrix[qGenome][rGenome] = (fastANI_matrix[qGenome][rGenome] + e.identity)/2;
            else
              fastANI_matrix[qGenome][rGenome] = e.identity;
          }
          else
          {
            if (fastANI_matrix[rGenome][qGenome] > 0)
              fastANI_matrix[rGenome][qGenome] = (fastANI_matrix[rGenome][qGenome] + e.identity)/2;
            else
              fastANI_matrix[rGenome][qGenome] = e.identity;
          }
        }
      }
    }

    std::ofstream outstrm(fileName + ".matrix");

    outstrm << totalGenomes << "\n";

    //Report matrix
    for (int i = 0; i < totalGenomes; i++)
    {
      //output genome name
      outstrm << genome2Int_rev[i];

      for (int j = 0; j < i; j++)
      {
        //output ani values
        //average if computed twice
        std::string val = fastANI_matrix[i][j] > 0.0 ? std::to_string (fastANI_matrix[i][j]) : "NA";
        outstrm << "\t" << val; 
      }
      outstrm << "\n";
    }

    outstrm.close();
  }

  /**
   * @brief                         generate multiple parameter objects from one
   * @details                       purpose it to divide the list of reference genomes
   *                                into as many buckets as there are threads 
   * @param[in]   parameters
   * @param[out]  parameters_split
   */
  void splitReferenceGenomes(skch::Parameters &parameters,
      std::vector <skch::Parameters> &parameters_split)
  {
    for (int i = 0; i < parameters.threads; i++)
    {
      parameters_split[i] = parameters;

      //update the reference genomes list
      parameters_split[i].refSequences.clear();

      //assign ref. genome to threads in round-robin fashion
      for (int j = 0; j < parameters.refSequences.size(); j++)
      {
        if (j % parameters.threads == i)
          parameters_split[i].refSequences.push_back (parameters.refSequences[j]);
      }
    }
  }

  /**
   * @brief                             update thread local reference genome ids to global ids
   * @param[in/out] CGI_ResultsVector
   */
  void correctRefGenomeIds (std::vector<cgi::CGI_Results> &CGI_ResultsVector)
  {
    int tid = omp_get_thread_num();
    int thread_count = omp_get_num_threads(); 
    
    for (auto &e : CGI_ResultsVector)
      e.refGenomeId = e.refGenomeId * thread_count +  tid;
  }
}

#endif
