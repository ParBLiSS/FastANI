/**
 * @file    MIIteratorL2.hpp
 * @brief   implements iterator over minimizer index to process L2 mapping stage
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef INDEX_ITERATOR_L2_HPP 
#define INDEX_ITERATOR_L2_HPP

#include <vector>
#include <algorithm>
#include <map>

//Own includes
#include "base_types.hpp"

//External includes

namespace skch
{
  /**
   * @class     skch::MIIteratorL2
   * @brief     L1 and L2 mapping stages
   */
  class MIIteratorL2
  {
    private:
      //Container type for saving read sketches during L1 and L2 both
      typedef Sketch::MI_Type MinVec_Type;

      typedef Sketch::MIIter_t MIIter_t;

      //Current begin position of a super-window
      //each next operation advances this value
      offset_t sw_pos;

      //Count of minimizer windows in a 'super-window'
      offset_t countMinimizerWindows;

    public:

      //Delete default constructor
      MIIteratorL2() = delete;

      //Current super-window range is bounded by [beg, end) on minimizer index
      MIIter_t sw_beg;
      MIIter_t sw_end;

      /**
       * @brief     constructor for MIIteratorL2 class
       * @details   sets the position and iterator values for the 
       *            first super-window
       */
      MIIteratorL2( MIIter_t firstSuperWindowRangeStart, 
          MIIter_t firstSuperWindowRangeEnd,
          offset_t countMinimizerWindows_ ) 
        :
        sw_beg (firstSuperWindowRangeStart),
        sw_end (firstSuperWindowRangeEnd),
        countMinimizerWindows (countMinimizerWindows_)
      {
        //Search for the end iterator of the first super-window over index
        this->sw_pos = this->sw_beg->wpos;
      }

      /**
       * @brief     advances the current super-window by a minimum shift
       *            possible on the minimizer index
       * @details   begin position of super-window is advanced by atleast 1 position
       *            Either or both 'beg' and 'end' iterators get right shifted
       * @NOTE      Calling function is expected to respect the bound of the reference index,
       *            this function will mis-behave if reached outside the index range
       */
      inline void next()
      {
        offset_t beginPos = this->sw_pos;
        offset_t lastPos = this->sw_pos + this->countMinimizerWindows - 1;

        // Always, range [beg, end) represents the minimizers in the 
        // current super-window
        assert( (this->sw_beg+1)->wpos - beginPos > 0 );
        assert( (this->sw_end  )->wpos - lastPos > 0 );

        offset_t advanceBy = std::min( (this->sw_beg+1)->wpos - beginPos, (this->sw_end)->wpos - lastPos);

        //Advance current super-window
        this->sw_pos += advanceBy;

        //Advance 'beg' and 'end' iterators

        if(advanceBy == (this->sw_beg + 1)->wpos - beginPos)
          this->sw_beg++;

        if(advanceBy == (this->sw_end)->wpos - lastPos)
          this->sw_end++;
      }
  };
}

#endif
