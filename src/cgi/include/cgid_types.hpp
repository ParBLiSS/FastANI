/**
 * @file    cgid_types.hpp
 * @brief   Critical type defintions for ANI algorithm
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef CGID_TYPES_HPP 
#define CGID_TYPES_HPP

#include <tuple>
#include "map/include/base_types.hpp"

namespace cgi
{
  //Saves mapping results from Mashmap
  struct MappingResult_CGI
  {
    skch::seqno_t refSequenceId;        //id of the reference contig
    skch::seqno_t genomeId;             //id of a genome
    skch::seqno_t querySeqId;           //name of query sequence
    skch::offset_t refStartPos;         //start position of the mapping on reference
    skch::offset_t queryStartPos;       //start position of the query for this mapping
    skch::offset_t mapRefPosBin;        //bin of mapped region on the reference sequence
    float nucIdentity;                  //calculated identity
  };

  /**
   * @brief     functor for comparing cgi mapping results by nucleotide identity, 
   *            while bucketing with genome id and query sequence id
   */
  struct compareMappingResult_withQuerySeqBucket
  {
    bool operator() (const MappingResult_CGI &x, const MappingResult_CGI &y)
    {
      return std::tie(x.genomeId, x.querySeqId, x.nucIdentity, x.refSequenceId, x.refStartPos) 
        < std::tie(y.genomeId, y.querySeqId, y.nucIdentity, y.refSequenceId, y.refStartPos);
      //Added ref. id and pos also to make sort output deterministic [issue #57]
    }
  } cmp_query_bucket;

  /**
   * @brief     functor for comparing cgi mapping results by nucleotide identity, 
   *            while bucketing with contig id and reference mapped region bin
   */
  struct compareMappingResult_withRefBinBucket
  {
    bool operator() (const MappingResult_CGI &x, const MappingResult_CGI &y)
    {
      //Note that when bucketing based on mapping position, reference sequence id should be used
      return std::tie(x.refSequenceId, x.mapRefPosBin, x.nucIdentity) 
        < std::tie(y.refSequenceId, y.mapRefPosBin, y.nucIdentity);
    }
  } cmp_refbin_bucket;

  /**
   * @brief     functor for comparing cgi mapping results by nucleotide identity, 
   */
  struct compareMappingResult_withIdentity
  {
    bool operator() (const MappingResult_CGI &x, const MappingResult_CGI &y)
    {
      //Note that when bucketing based on mapping position, reference sequence id should be used
      return x.nucIdentity < y.nucIdentity;
    }
  } cmp_identity;

  //Final format to save CGI results
  struct CGI_Results
  {
    skch::seqno_t refGenomeId;
    skch::seqno_t qryGenomeId;
    skch::seqno_t countSeq;
    skch::seqno_t totalQueryFragments;
    float identity;

    bool operator <(const CGI_Results& x) const {
      return std::tie(x.qryGenomeId, identity) 
        < std::tie(qryGenomeId, x.identity);
    }
  };
}

#endif
