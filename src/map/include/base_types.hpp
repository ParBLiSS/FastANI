/**
 * @file    base_types.hpp
 * @brief   Critical type defintions for mapping algorithm
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef BASE_TYPES_MAP_HPP 
#define BASE_TYPES_MAP_HPP

#include <tuple>

namespace skch
{
  typedef uint32_t hash_t;    //hash type
  typedef int offset_t;       //position within sequence
  typedef int seqno_t;        //sequence counter in file
  typedef int32_t strand_t;   //sequence strand 

  //C++ timer
  typedef std::chrono::high_resolution_clock Time;

  //Information about each minimizer
  struct MinimizerInfo
  {
    hash_t hash;                              //hash value
    seqno_t seqId;                            //sequence or contig id
    offset_t wpos;                            //First (left-most) window position when the minimizer is saved
    strand_t strand;                          //strand information

    //Lexographical less than comparison
    bool operator <(const MinimizerInfo& x) {
      return std::tie(hash, seqId, wpos, strand) 
        < std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    //Lexographical equality comparison
    bool operator ==(const MinimizerInfo& x) {
      return std::tie(hash, seqId, wpos, strand) 
        == std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    bool operator !=(const MinimizerInfo& x) {
      return std::tie(hash, seqId, wpos, strand) 
        != std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    static bool equalityByHash(const MinimizerInfo& x, const MinimizerInfo& y) {
      return x.hash == y.hash;
    }

    static bool lessByHash(const MinimizerInfo& x, const MinimizerInfo& y) {
      return x.hash < y.hash;
    }

  };

  //Type for map value type used for
  //L1 stage lookup index
  struct MinimizerMetaData
  {
    seqno_t seqId;          //sequence or contig id
    offset_t wpos;          //window position (left-most window)
    strand_t strand;        //strand information

    bool operator <(const MinimizerMetaData& x) const {
      return std::tie(seqId, wpos, strand) 
        < std::tie(x.seqId, x.wpos, x.strand);
    }
  };

  typedef hash_t MinimizerMapKeyType;
  typedef std::vector<MinimizerMetaData> MinimizerMapValueType;

  //Metadata recording for contigs in the reference DB
  struct ContigInfo
  {
    std::string name;       //Name of the sequence
    offset_t len;           //Length of the sequence
  };

  //Label tags for strand information
  enum strnd : strand_t
  {
    FWD = 1,  
    REV = -1
  };  

  //Information about query sequence during L1/L2 mapping
  template <typename KSEQ, typename MinimizerVec>
    struct QueryMetaData
    {
      KSEQ seq;                           //query sequence object pointer (kseq library) 
      seqno_t seqCounter;                 //query sequence counter
      offset_t len;                       //length of this query sequence
      int sketchSize;                     //sketch size
      MinimizerVec minimizerTableQuery;   //Vector of minimizers in the query 
    };

  //Final mapping result
  struct MappingResult
  {
    offset_t queryLen;                //length of the query sequence
    offset_t refStartPos;             //start position of the mapping on reference
    offset_t refEndPos;               //end pos
    seqno_t refSeqId;                 //internal sequence id of the reference contig
    seqno_t querySeqId;               //internal sequence id of the query sequence
    float nucIdentity;                //calculated identity
    float nucIdentityUpperBound;      //upper bound on identity (90% C.I.)
    int sketchSize;                   //sketch size
    int conservedSketches;            //count of conserved sketches
    strand_t strand;                  //strand
    float mappedRegionComplexity;     //estimated entropy in the mapped region on reference
    std::string queryName;            //name of query sequence
  };

  typedef std::vector<MappingResult> MappingResultsVector_t;
}

#endif
