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

        //Ordered map to save unique sketch elements, and associated value as 
        //a pair of its occurrence in the query and the reference
        typedef std::map< hash_t, slidingMapContainerValueType > MapType;
        MapType slidingWindowMinhashes;

        //Iterator pointing to the smallest 's'th element in the map
        typename MapType::iterator pivot;

        //Label status while inserting reference minimizer
        enum IN : int
        {
          //reference minimizer is inserted into map as a new map entry
          UNIQ = 1,

          //reference minimizer is coupled with a query minimizer, 
          //previously, there was no ref. minimizer at this entry
          CPLD = 2, 

          //ref. minimizer just revises the hash position of already 
          //existing reference minimizer
          REV = 3
        };  

        //Label status while deleting reference minimizer
        enum OUT : int
        {
          //entry in the map is deleted
          DEL = 1,

          //just the reference minimizer is updated to null
          UPD = 2,

          //Nothing changed in the map
          NOOP = 3
        };

      public:

        //Count of shared sketch elements between query and the reference
        //Updated after insert or delete operation on map 
        int sharedSketchElements;


        //Delete default constructor
        SlideMapper() = delete;

        /**
         * @brief                 constructor
         * @param[in]   Q         query meta data
         */
        SlideMapper(Q_Info &Q_) :
          Q(Q_),
          pivot(this->slidingWindowMinhashes.end()),
          sharedSketchElements(0)
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

          //Insert query sketch elements to map
          for(auto it = Q.minimizerTableQuery.begin(); it != uniqEndIter; it++)
          {
            this->slidingWindowMinhashes.emplace_hint(slidingWindowMinhashes.end(), it->hash, slidingMapContainerValueType{it->wpos, it->strand, NAPos, 0});
          }

          //Point pivot to last element in the map
          this->pivot = std::next(this->slidingWindowMinhashes.begin(), Q.sketchSize - 1);

          //Current count of shared sketch elements is zero
          this->sharedSketchElements = 0;
        }

      public:

        /**
         * @brief               insert a minimizer from the reference sequence into the map
         * @param[in]   m       reference minimizer to insert
         */
        inline void insert_ref(MIIter_t m)
        {
          hash_t hashVal = m->hash;
          int status;

          //if hash doesn't exist in the map, add to it
          if(slidingWindowMinhashes.find(hashVal) == slidingWindowMinhashes.end())
          {
            slidingWindowMinhashes[hashVal] = slidingMapContainerValueType{this->NAPos, 0, m->wpos, m->strand};   //add the hash to window
            status = IN::UNIQ;
          }
          else
          {
            status = (slidingWindowMinhashes[hashVal].wposR == NAPos) ? IN::CPLD 
              : IN::REV;

            //if hash already exists in the map, just revise it
            slidingWindowMinhashes[hashVal].wposR = m->wpos;
            slidingWindowMinhashes[hashVal].strandR = m->strand;
          }

          updateCountersAfterInsert(status, m);

          assert(this->sharedSketchElements >= 0);
          assert(this->sharedSketchElements <= Q.sketchSize);
        }

        /**
         * @brief               delete a minimizer from the reference sequence from the map
         * @param[in]   m       reference minimizer to remove
         */
        inline void delete_ref(MIIter_t m)
        {
          hash_t hashVal = m->hash;
          int status;
          bool pivotDeleteCase = false;
          
          assert(this->slidingWindowMinhashes.find(hashVal) != this->slidingWindowMinhashes.end());

          //This hashVal may exist with different wpos from 
          //reference, do nothing in that case
          
          if(this->slidingWindowMinhashes[hashVal].wposR == m->wpos)
          {
            if(this->slidingWindowMinhashes[hashVal].wposQ == NAPos)
            {
              //Handle pivot deletion as a separate case
              if(this->slidingWindowMinhashes.find(hashVal) == pivot)
              {
                pivot++;

                if( (this->pivot->second).wposQ != NAPos && (this->pivot->second).wposR != NAPos)
                  this->sharedSketchElements += 1;

                pivotDeleteCase = true;
              }

              this->slidingWindowMinhashes.erase(hashVal);              //Remove the entry from the map
              status = OUT::DEL; 
            }
            else
            {
              this->slidingWindowMinhashes[hashVal].wposR = NAPos;      //Just mark the reference hash absent
              status = OUT::UPD; 
            }
          }
          else
          {
            status = OUT::NOOP;
          }

          if(!pivotDeleteCase) updateCountersAfterDelete(status, m);

          assert(this->sharedSketchElements >= 0);
          assert(this->sharedSketchElements <= Q.sketchSize);
        }

        /**
         * @brief               insert a range of minimizers from the reference sequence into the map
         * @param[in]   begin   begin iterator
         * @param[in]   end     end iterator
         */
        inline void insert_ref(MIIter_t begin, MIIter_t end)
        {
          for(auto it = begin; it != end; it++)
            this->insert_ref(it);
        }

        /**
         * @brief       compute strand consensus and unique reference hashes
         * @param[out]  strandVotes
         * @param[out]  uniqueRefHashes
         */
        inline void computeStatistics(int &strandVotes, int &uniqueRefHashes)
        {
          int uniqueHashes = 0;
          strandVotes = uniqueRefHashes = 0;

          //Iterate over map
          for (auto it = this->slidingWindowMinhashes.cbegin(); it != this->slidingWindowMinhashes.cend(); ++it)
          {
            uniqueHashes++;

            //Fetch the value in map
            auto m = it->second;

            if(uniqueHashes <= Q.sketchSize && m.wposQ != this->NAPos && m.wposR != this->NAPos)
            {
              strandVotes += m.strandQ * m.strandR; //Assuming FWD=1, BWD=-1
            }

            //Check if minimizer occurs comes from the reference
            if(m.wposR != this->NAPos)
              uniqueRefHashes++;
          }
        }

      private:

        /**
         * @brief             logic to update internal counters after insert to map
         * @param[in] status  insert status
         * @param[in] m       reference minimizer that was inserted
         */
        inline void updateCountersAfterInsert(int status, MIIter_t m)
        {
          //Revise internal counters
          if(m->hash <= this->pivot->first)
          {
            if(status == IN::CPLD)
            {
              //Increase count of shared sketch elements by 1
              this->sharedSketchElements += 1;
            }
            else if(status == IN::UNIQ)
            {
              //Pivot needs to be decremented
              if( (this->pivot->second).wposQ != NAPos && (this->pivot->second).wposR != NAPos)
                this->sharedSketchElements -= 1;

              std::advance(this->pivot, -1);
            }
            else if(status == IN::REV)
            {
              //Do nothing
            }
          }
        }

        /**
         * @brief             logic to update internal counters after delete to map
         * @param[in] status  insert status
         * @param[in] m       reference minimizer that was inserted
         */
        inline void updateCountersAfterDelete(int status, MIIter_t m)
        {
          //Revise internal counters
          if(m->hash <= this->pivot->first)
          {
            if(status == OUT::UPD)
            {
              //Decrease count of shared sketch elements by 1
              this->sharedSketchElements -= 1;
            }
            else if(status == OUT::DEL)
            {
              //Pivot needs to be advanced
              std::advance(this->pivot, 1);

              if( (this->pivot->second).wposQ != NAPos && (this->pivot->second).wposR != NAPos)
                this->sharedSketchElements += 1;
            }
            else if(status == OUT::NOOP)
            {
              //Do nothing
            }
          }
        }

    };
}

#endif
