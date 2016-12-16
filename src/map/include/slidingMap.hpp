/**
 * @file    slidingMap.hpp
 * @brief   implements ordered map to compute Jaccard
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SLIDING_MAP_HPP 
#define SLIDING_MAP_HPP

#include <vector>
#include <algorithm>
#include <map>

//Own includes
#include "base_types.hpp"

//External includes

namespace skch
{
  /**
   * @class     skch::SlideMapper
   * @brief     L1 and L2 mapping stages
   */
  template <typename Q_Info>
    class SlideMapper
    {

      private:

        //Metadata for the minimizers saved in sliding ordered map during L2 stage
        struct slidingMapContainerValueType
        {
          offset_t wposQ;                   //wpos and strand of minimizers in the query
          strand_t strandQ;
          offset_t wposR;                   //wpos and strand of minimizers in the reference
          strand_t strandR;
        };

        //Container type for saving read sketches during L1 and L2 both
        typedef Sketch::MI_Type MinVec_Type;

        typedef Sketch::MIIter_t MIIter_t;

        //reference to query's metadata
        const Q_Info &Q;

        //Define a Not available position marker
        static const offset_t NAPos = std::numeric_limits<offset_t>::max();

      public:

        //Ordered map to save unique sketch elements, and associated value as 
        //a pair of its occurrence in the query and the reference
        std::map< hash_t, slidingMapContainerValueType > slidingWindowMinhashes;

        //Delete default constructor
        SlideMapper() = delete;

        /**
         * @brief                 constructor
         * @param[in]   Q         query meta data
         */
        SlideMapper(Q_Info &Q_) :
          Q(Q_)
        {
          this->init();
        }

      private:

        /**
         * @brief       Fills map with minimum 's' minimizers in the query
         */
        inline void init()
        {
          //Range of sketch in query
          //Assuming unique query minimizers were placed at the start during L1 mapping
          auto uniqEndIter = std::next(Q.minimizerTableQuery.begin(), Q.sketchSize);

          for(auto it = Q.minimizerTableQuery.begin(); it != uniqEndIter; it++)
          {
            this->slidingWindowMinhashes.emplace_hint(slidingWindowMinhashes.end(), it->hash, slidingMapContainerValueType{it->wpos, it->strand, NAPos, 0});
          }
        }

      public:

        /**
         * @brief               insert a minimizer from the reference sequence into the map
         * @param[in]   m       reference minimizer to insert
         */
        inline void insert_ref(const MinimizerInfo &m)
        {
          hash_t hashVal = m.hash;

          //if hash doesn't exist in the map, add to it
          if(slidingWindowMinhashes.find(hashVal) == slidingWindowMinhashes.end())
            slidingWindowMinhashes[hashVal] = slidingMapContainerValueType{this->NAPos, 0, m.wpos, m.strand};   //add the hash to window
          else
          {
            //if hash already exists in the map, just revise it
            //Note that, if hash exists in the map from reference as well, 
            //           we still need to revise the wposR to keep the right-most entry
            //           When we delete a duplicate entry from the map, we should not remove this one.
            slidingWindowMinhashes[hashVal].wposR = m.wpos;
            slidingWindowMinhashes[hashVal].strandR = m.strand;
          }
        }

        /**
         * @brief               insert a range of minimizers from the reference sequence into the map
         * @param[in]   begin   begin iterator
         * @param[in]   end     end iterator
         */
        inline void insert_ref(MIIter_t begin, MIIter_t end)
        {
          for(auto it = begin; it != end; it++)
            this->insert_ref(*it);
        }

        /**
         * @brief               delete a minimizer from the reference sequence from the map
         * @param[in]   m       reference minimizer to remove
         */
        inline void delete_ref(const MinimizerInfo &m)
        {
          hash_t hashVal = m.hash;
          
          assert(this->slidingWindowMinhashes.find(hashVal) != this->slidingWindowMinhashes.end());

          //This hashVal may exist with different wpos from reference, do nothing in that case
          
          if(this->slidingWindowMinhashes[hashVal].wposR == m.wpos)
          {
            if(this->slidingWindowMinhashes[hashVal].wposQ == NAPos)
              this->slidingWindowMinhashes.erase(hashVal);              //Remove the entry from the map
            else
              this->slidingWindowMinhashes[hashVal].wposR = NAPos;      //Just mark the reference hash absent
          }
        }

        /**
         * @brief                                     compute shared sketch elements between 
         *                                            query and the reference window
         * @param[out]    currentSharedMinimizers     #shared minimizers among the smallest s  
         * @param[out]    strandVotes                 #consensus strand vote among the shared minimizers
         * @param[out]    uniqueRefHashes             #unique minimizers from reference in the complete slidingMap 
         */
        inline void computeSharedMinimizers(int &currentSharedMinimizers, int &strandVotes, int &uniqueRefHashes)
        {
          int uniqueHashes = 0;
          currentSharedMinimizers = strandVotes = uniqueRefHashes = 0;

          //Iterate over map
          for (auto it = this->slidingWindowMinhashes.cbegin(); it != this->slidingWindowMinhashes.cend(); ++it)
          {
            uniqueHashes++;

            //Fetch the value in map
            auto m = it->second;

            //Check if minimizer is shared (among s unique sketches)
            if(uniqueHashes <= Q.sketchSize && m.wposQ != this->NAPos && m.wposR != this->NAPos)
            {
              currentSharedMinimizers += 1;
              strandVotes += m.strandQ * m.strandR; //Assuming FWD=1, BWD=-1
            }

            //Check if minimizer occurs comes from the reference
            if(m.wposR != this->NAPos)
              uniqueRefHashes++;
          }
        }

    };
}

#endif
