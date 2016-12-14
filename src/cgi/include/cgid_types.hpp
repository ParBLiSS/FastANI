/**
 * @file    cgid_types.hpp
 * @brief   Critical type defintions for ANI algorithm
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef CGID_TYPES_HPP 
#define CGID_TYPES_HPP

#include <tuple>
#include "base_types.hpp"

namespace cgi
{
  //Saves mapping results from Mashmap
  struct MappingResult_CGI
  {
    skch::seqno_t genomeId;           //internal file id of the reference genome
    skch::seqno_t querySeqId;         //name of query sequence
    skch::offset_t mapRefPosBin;       //bin of mapped region on the reference sequence
    float nucIdentity;                //calculated identity
  };

  /**
   * @brief     functor for comparing cgi mapping results by nucleotide identity, 
   *            while bucketing with genome id and query sequence id
   */
  struct compareMappingResult_withQuerySeqBucket
  {
    bool operator() (const MappingResult_CGI &x, const MappingResult_CGI &y)
    {
      return std::tie(x.genomeId, x.querySeqId, x.nucIdentity) 
        < std::tie(y.genomeId, y.querySeqId, y.nucIdentity);
    }
  } cmp_query_bucket;

  /**
   * @brief     functor for comparing cgi mapping results by nucleotide identity, 
   *            while bucketing with genome id and reference mapped region bin
   */
  struct compareMappingResult_withRefBinBucket
  {
    bool operator() (const MappingResult_CGI &x, const MappingResult_CGI &y)
    {
      return std::tie(x.genomeId, x.mapRefPosBin, x.nucIdentity) 
        < std::tie(y.genomeId, y.mapRefPosBin, y.nucIdentity);
    }
  } cmp_refbin_bucket;


  //Final format to save CGI results
  struct CGI_Results
  {
    skch::seqno_t genomeId;
    skch::seqno_t countSeq;
    float identity;

    //Default comparison is by identity
    bool operator <(const CGI_Results& x) const {
      return identity < x.identity;
    }
  };

}

#endif
