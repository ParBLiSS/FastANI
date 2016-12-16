/**
 * @file    computeMap.hpp
 * @brief   implments the sequence mapping logic
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_MAP_HPP 
#define SKETCH_MAP_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <zlib.h>  

//Own includes
#include "base_types.hpp"
#include "map_parameters.hpp"
#include "commonFunc.hpp"
#include "winSketch.hpp"
#include "map_stats.hpp"
#include "slidingMap.hpp"

//External includes

namespace skch
{
  /**
   * @class     skch::Map
   * @brief     L1 and L2 mapping stages
   */
  class Map
  {
    public:

      //Type for Stage L1's predicted candidate location
      struct L1_candidateLocus_t 
      {
        seqno_t seqId;                    //sequence id where read is mapped

        /* read could be mapped with its begin location
         * from [rangeStartPos, rangeEndPos]
         */
        offset_t rangeStartPos;
        offset_t rangeEndPos;  
      };

      //Type for Stage L2's predicted mapping coordinate within each L1 candidate
      struct L2_mapLocus_t 
      {
        seqno_t seqId;                    //sequence id where read is mapped
        offset_t optimalStartPos;         //optimal start mapping position 
        int sharedSketchSize;             //count of shared sketch elements
        int uniqueRefHashes;              //count of unique reference hashes in the mapping region
        strand_t strand;                  //mapping strand
      };

    private:

      //algorithm parameters
      const skch::Parameters &param;

      //reference sketch
      const skch::Sketch &refSketch;

      //Container type for saving read sketches during L1 and L2 both
      typedef Sketch::MI_Type MinVec_Type;

      typedef Sketch::MIIter_t MIIter_t;

      //Custom function for post processing the results, by default does nothing 
      typedef std::function< void(const MappingResult&) > PostProcessResultsFn_t;
      PostProcessResultsFn_t processMappingResults;

    public:

      /**
       * @brief                 constructor
       * @param[in] p           algorithm parameters
       * @param[in] refSketch   reference sketch
       * @param[in] f           optional user defined custom function to post process the reported mapping results
       */
      Map(const skch::Parameters &p, const skch::Sketch &refsketch,
          PostProcessResultsFn_t f = nullptr) :
        param(p),
        refSketch(refsketch),
        processMappingResults(f)
    {
      this->mapQuery();
    }

    private:

      /**
       * @brief   parse over sequences in query file and map each on the reference
       */
      void mapQuery()
      {
        //Count of reads mapped by us
        //Some reads are dropped because of short length
        seqno_t totalReadsPickedForMapping = 0;
        seqno_t totalReadsMapped = 0;
        seqno_t seqCounter = 0;

        std::ofstream outstrm(param.outFileName);

        for(const auto &fileName : param.querySequences)
        {
          //Open the file using kseq
          FILE *file = fopen(fileName.c_str(), "r");
          gzFile fp = gzdopen(fileno(file), "r");
          kseq_t *seq = kseq_init(fp);

#ifdef DEBUG
          std::cout << "INFO, skch::Map::mapQuery, mapping reads in " << fileName << std::endl;
#endif

          //size of sequence
          offset_t len;

          while ((len = kseq_read(seq)) >= 0) 
          {
            //Is the read too short?
            if(len < param.windowSize || len < param.kmerSize || len < param.minReadLength)
            {
              seqCounter++;

#ifdef DEBUG
              std::cout << "WARNING, skch::Map::mapQuery, read is not long enough for mapping" << std::endl;
#endif

              continue;
            }
            else 
            {
              QueryMetaData <decltype(seq), MinVec_Type> Q;

              Q.seq = seq;
              Q.seqCounter = seqCounter;
              Q.len = len; 

              //Map this sequence
              bool mappingReported = mapSingleQuerySeq(Q, outstrm);
              if(mappingReported)
                totalReadsMapped++;

              seqCounter++;
              totalReadsPickedForMapping++;
            }
          }

          //Close the input file
          kseq_destroy(seq);  
          gzclose(fp);  
        }

        std::cout << "INFO, skch::Map::mapQuery, [count of mapped reads, reads qualified for mapping, total input reads] = [" << totalReadsMapped << ", " << totalReadsPickedForMapping << ", " << seqCounter << "]" << std::endl;

      }

      /**
       * @brief               map the parsed query sequence (L1 and L2 mapping)
       * @param[in] Q         metadata about query sequence
       * @param[in] outstrm   outstream stream where mappings will be reported
       */
      template<typename Q_Info>
        inline bool mapSingleQuerySeq(Q_Info &Q, std::ofstream &outstrm)
        {
#if ENABLE_TIME_PROFILE_L1_L2
          auto t0 = skch::Time::now();
#endif

#if ENABLE_TIME_PROFILE_L1_L2
          auto t1 = skch::Time::now();
#endif
          //L1 Mapping
          std::vector<L1_candidateLocus_t> l1Mappings; 
          doL1Mapping(Q, l1Mappings);

#if ENABLE_TIME_PROFILE_L1_L2
          std::chrono::duration<double> timeSpentL1 = skch::Time::now() - t1;
          t1 = skch::Time::now();
#endif

          //L2 Mapping
          MappingResultsVector_t l2Mappings;
          bool mappingReported = doL2Mapping(Q, l1Mappings, l2Mappings);

          //Write mapping results to file
          reportL2Mappings(l2Mappings, outstrm);


#if ENABLE_TIME_PROFILE_L1_L2
          {
            std::chrono::duration<double> timeSpentL2 = skch::Time::now() - t1;
            std::chrono::duration<double> timeSpentMappingRead = skch::Time::now() - t0;
            int countL1Candidates = l1Mappings.size();

            std::cerr << Q.seq->name.s << " " << Q.len
              << " " << countL1Candidates 
              << " " << timeSpentL1.count() 
              << " " << timeSpentL2.count()
              << " " << timeSpentMappingRead.count()
              << "\n";
          }
#endif

          return mappingReported;
        }

      /**
       * @brief       Find candidate regions for a read using level 1 (seed-hits) mapping
       * @details     The count of hits that should occur within a region on the reference is 
       *              determined by the threshold similarity
       *              The resulting start and end target offsets on reference is (are) an 
       *              overestimate of the mapped region. Computing better bounds is left for
       *              the following L2 stage.
       * @param[in]   Q                         query sequence details 
       * @param[out]  l1Mappings                all the read mapping locations
       */
      template <typename Q_Info, typename Vec>
        void doL1Mapping(Q_Info &Q, Vec &l1Mappings)
        {
          //Vector of positions of all the hits 
          std::vector<MinimizerMetaData> seedHitsL1;

          ///1. Compute the minimizers

          CommonFunc::addMinimizers(Q.minimizerTableQuery, Q.seq, param.kmerSize, param.windowSize, param.alphabetSize);

#ifdef DEBUG
          std::cout << "INFO, skch::Map:doL1Mapping, read id " << Q.seqCounter << ", minimizer count = " << Q.minimizerTableQuery.size() << "\n";
#endif

          ///2. Find the hits in the reference, pick 's' unique minimizers as seeds, 

          std::sort(Q.minimizerTableQuery.begin(), Q.minimizerTableQuery.end(), MinimizerInfo::lessByHash);

          //note : unique preserves the original relative order of elements 
          auto uniqEndIter = std::unique(Q.minimizerTableQuery.begin(), Q.minimizerTableQuery.end(), MinimizerInfo::equalityByHash);

          //This is the sketch size for estimating jaccard
          Q.sketchSize = std::distance(Q.minimizerTableQuery.begin(), uniqEndIter);

          int totalMinimizersPicked = 0;

          for(auto it = Q.minimizerTableQuery.begin(); it != uniqEndIter; it++)
          {
            //Check if hash value exists in the reference lookup index
            auto seedFind = refSketch.minimizerPosLookupIndex.find(it->hash);

            if(seedFind != refSketch.minimizerPosLookupIndex.end())
            {
              auto hitPositionList = seedFind->second;

              //Save the positions (Ignore high frequency hits)
              if(hitPositionList.size() < refSketch.getFreqThreshold())
              {
                seedHitsL1.insert(seedHitsL1.end(), hitPositionList.begin(), hitPositionList.end());
              }

            }
          }

          int minimumHits = Stat::estimateMinimumHitsRelaxed(Q.sketchSize, param.kmerSize, param.percentageIdentity);

          this->computeL1CandidateRegions(Q, seedHitsL1, minimumHits, l1Mappings);

#ifdef DEBUG
          std::cout << "INFO, skch::Map:doL1Mapping, read id " << Q.seqCounter << ", Count of L1 hits in the reference = " << seedHitsL1.size() << ", minimum hits required for a candidate = " << minimumHits << ", Count of L1 candidate regions = " << l1Mappings.size() << "\n";

          for(auto &e : l1Mappings)
            std::cout << "INFO, skch::Map:doL1Mapping, read id " << Q.seqCounter << ", L1 candidate : [" << this->refSketch.metadata[std::get<0>(e)].name << " " << this->refSketch.metadata[std::get<0>(e)].len << " " << std::get<1>(e) << " " << std::get<2>(e) << "]\n";
#endif

        }

      /**
       * @brief                     Helper function to doL1Mapping()
       * @param[in]   Q             query
       * @param[in]   seedHitsL1    minimizer hits in the reference
       * @param[in]   minimumHits   estimated minimum hits required for significant match
       * @param[out]  l1Mappings    all the read mapping locations
       */
      template <typename Q_Info, typename Vec1, typename Vec2>
        void computeL1CandidateRegions(Q_Info &Q, Vec1 &seedHitsL1, int minimumHits, Vec2 &l1Mappings)
        {
          if(minimumHits < 1)
            minimumHits = 1;

          //Sort all the hit positions
          std::sort(seedHitsL1.begin(), seedHitsL1.end());

          for(auto it = seedHitsL1.begin(); it != seedHitsL1.end(); it++)
          {
            if(std::distance(it, seedHitsL1.end()) >= minimumHits)
            {
              auto it2 = it + minimumHits -1;
              //[it .. it2] are 'minimumHits' consecutive hits 

              //Check if consecutive hits are close enough
              //NOTE: hits may span more than a read length for a valid match, as we keep window positions 
              //      for each minimizer
              if(it2->seqId == it->seqId && it2->wpos - it->wpos < Q.len + param.windowSize)
              {
                //Save <1st pos --- 2nd pos>
                L1_candidateLocus_t candidate{it->seqId, 
                    std::max(0, it2->wpos - Q.len + 1), it->wpos};

                //Check if this candidate overlaps with last inserted one
                auto lst = l1Mappings.end(); lst--;

                //match seq_no and see if this candidate begins before last element ends
                if( l1Mappings.size() > 0 
                    && candidate.seqId == lst->seqId 
                    && lst->rangeEndPos >= candidate.rangeStartPos)
                {
                  //Push the end pos of last candidate locus further out
                  lst->rangeEndPos = std::max(candidate.rangeEndPos, lst->rangeEndPos);
                }
                else
                  l1Mappings.push_back(candidate);
              }
            }
          }
        }

      /**
       * @brief                                 Revise L1 candidate regions to more precise locations
       * @param[in]   Q                         query sequence information
       * @param[in]   l1Mappings                candidate regions for query sequence found at L1
       * @param[out]  l2Mappings                Mapping results in the L2 stage
       * @return      T/F                       true if atleast 1 mapping region is proposed
       */
      template <typename Q_Info, typename VecIn, typename VecOut>
        bool doL2Mapping(Q_Info &Q, VecIn &l1Mappings, VecOut &l2Mappings)
        {
          bool mappingReported = false;

          ///2. Walk the read over the candidate regions and compute the jaccard similarity with minimum s sketches
          for(auto &candidateLocus: l1Mappings)
          {
            L2_mapLocus_t l2;
            computeL2MappedRegions(Q, candidateLocus, l2);

            //Compute mash distance using calculated jaccard
            float mash_dist = Stat::j2md(1.0 * l2.sharedSketchSize/Q.sketchSize, param.kmerSize);

            //Compute lower bound to mash distance within 90% confidence interval
            float mash_dist_lower_bound = Stat::md_lower_bound(mash_dist, Q.sketchSize, param.kmerSize, 0.9);

            float nucIdentity = 100 * (1 - mash_dist);
            float nucIdentityUpperBound = 100 * (1 - mash_dist_lower_bound);

            //Compute reference region complexity
            float actualDensity = l2.uniqueRefHashes * 1.0 / Q.len;
            float expectedDensity = 2.0 / param.windowSize;
            float referenceDNAComplexity = actualDensity/expectedDensity;

            //Report the alignment
            if(nucIdentityUpperBound >= param.percentageIdentity && referenceDNAComplexity >= 0.75)
            {
              MappingResult res;

              //Save the output
              {
                res.queryLen = Q.len;
                res.refStartPos = l2.optimalStartPos ;
                res.refEndPos = l2.optimalStartPos + Q.len - 1;
                res.refSeqId = l2.seqId;
                res.querySeqId = Q.seqCounter;
                res.nucIdentity = nucIdentity;
                res.nucIdentityUpperBound = nucIdentityUpperBound;
                res.sketchSize = Q.sketchSize;
                res.conservedSketches = l2.sharedSketchSize;
                res.strand = l2.strand;
                res.mappedRegionComplexity = referenceDNAComplexity;
                res.queryName = Q.seq->name.s; 

                l2Mappings.push_back(res);
              }

              mappingReported = true;
            }
          }

          return mappingReported;
        }

      /**
       * @brief                                 Find optimal mapping within an L1 candidate
       * @param[in]   Q                         query sequence information
       * @param[in]   candidateLocus            L1 candidate location
       * @param[out]  l2_out                    L2 mapping inside L1 candidate 
       */
      template <typename Q_Info>
        void computeL2MappedRegions(Q_Info &Q, 
            L1_candidateLocus_t &candidateLocus, 
            L2_mapLocus_t &l2_out)
        {
          /// 1. Create an array of minimizers with wpos [START, END), continuous in the position space
          offset_t START = candidateLocus.rangeStartPos; 
          offset_t END = candidateLocus.rangeEndPos + Q.len;

          // This is the vector of minimizers we need for computing Jaccard estimates
          MinVec_Type allMinimizersInRange;

          {
            //Look up within the reference minimizerIndex
            MIIter_t minimizerIndexRangeStart = this->refSketch.searchIndex(candidateLocus.seqId, candidateLocus.rangeStartPos);

            for(auto it = minimizerIndexRangeStart; it != this->refSketch.getMinimizerIndexEnd() && it->seqId == candidateLocus.seqId && it->wpos < END; it++)
            {
              if(allMinimizersInRange.size() > 0)
              {
                auto last_min = allMinimizersInRange.back();

                //last_min minimizer stays minimum for [it->wpos - last_min.wpos - 1] additional positions
                allMinimizersInRange.insert(allMinimizersInRange.end(), it->wpos - last_min.wpos - 1, last_min);
              }

              //Now insert the new minimizer
              allMinimizersInRange.emplace_back(*it);
            }
          }

          ///2. Walk the read over windows in 'allMinimizersInRange'
          {
            l2_out.sharedSketchSize = 0;
            l2_out.seqId = candidateLocus.seqId;

            /**
             * The count of winnowing windows in a sequence of a length L is
             *    L - windowSize - kmerSize
             */
            auto countMinimizerWindows = Q.len - (param.windowSize-1) - (param.kmerSize-1);

            //Define map such that it contains only the query minimizers
            SlideMapper<Q_Info> slidemap(Q);


            for(int i = 0; i + countMinimizerWindows < allMinimizersInRange.size(); i++)
            {
              //Consider minimizers in 'allMinimizersInRange' in the range [i, j)
              int j = i + countMinimizerWindows;

              //First super-window or we see new first hash value or new last hash value, we should compute then
              if(i == 0 || allMinimizersInRange[i] != allMinimizersInRange[i-1] || allMinimizersInRange[j-1]  != allMinimizersInRange[j-2])
              {

                //Push reference minimizers into map
                if(i == 0) 
                {
                  slidemap.insert_ref(allMinimizersInRange.begin() + i, allMinimizersInRange.begin() + j);
                }
                else 
                {
                  //Minimizer to delete?
                  if(allMinimizersInRange[i] != allMinimizersInRange[i-1])
                    slidemap.delete_ref( allMinimizersInRange[i-1] );

                  //New minimizer to insert?
                  if(allMinimizersInRange[j-1]  != allMinimizersInRange[j-2])
                    slidemap.insert_ref( allMinimizersInRange[j-1] );
                }

                int currentSharedMinimizers, strandVotes, uniqueRefHashes;

                //Compute the count of shared sketch elements as well as the strand
                slidemap.computeSharedMinimizers(currentSharedMinimizers, strandVotes, uniqueRefHashes);

                //Is this sliding window the best we have so far?
                if(currentSharedMinimizers > l2_out.sharedSketchSize)
                {
                  l2_out.sharedSketchSize = currentSharedMinimizers;
                  l2_out.strand = strandVotes > 0 ? strnd::FWD : strnd::REV;
                  l2_out.uniqueRefHashes = uniqueRefHashes;
                  l2_out.optimalStartPos = allMinimizersInRange[i].wpos;
                }

              }//End of if condition
            }//End of loop for sliding the read
          }//End of the phase 2 (computing optimal position)
        }

      /**
       * @brief                     Report the final L2 mappings to output stream
       * @param[in]   l2Mappings    mapping results
       * @param[in]   outstrm       file output stream object
       */
      void reportL2Mappings(MappingResultsVector_t &l2Mappings, 
          std::ofstream &outstrm)
      {
        float bestNucIdentity = 0;

        //Scan through the mappings to check best identity mapping
        for(auto &e : l2Mappings)
        {
          if(e.nucIdentity > bestNucIdentity)
            bestNucIdentity = e.nucIdentity;
        }

        //Print the results
        for(auto &e : l2Mappings)
        {
          //Report top 1% mappings (unless reportAll flag is true, in which case we report all)
          if(param.reportAll == true || e.nucIdentity >= bestNucIdentity - 1.0)
          {
            outstrm << e.queryName 
              << " " << e.queryLen 
              << " " << "0"
              << " " << e.queryLen - 1 
              << " " << (e.strand == strnd::FWD ? "+" : "-") 
              << " " << this->refSketch.metadata[e.refSeqId].name
              << " " << this->refSketch.metadata[e.refSeqId].len
              << " " << e.refStartPos 
              << " " << e.refEndPos
              << " " << e.nucIdentity;

            //Print some additional statistics
            outstrm << " " << e.conservedSketches 
              << " " << e.sketchSize 
              << " " << e.nucIdentityUpperBound
              << " " << e.mappedRegionComplexity;

            outstrm << "\n";

            //User defined processing of the results
            if(processMappingResults != nullptr)
              processMappingResults(e);
          }
        }
      }

    public:

      /**
       * @brief     An optional utility function to save the 
       *            reported results by the L2 stage into a vector
       */
      static void insertL2ResultsToVec(MappingResultsVector_t &v, const MappingResult &reportedL2Result)
      {
        v.push_back(reportedL2Result);
      }
  };

}

#endif
