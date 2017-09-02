/**
 * @file    map_stats.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef MAP_STATS_HPP 
#define MAP_STATS_HPP

#include <vector>
#include <algorithm>
#include <deque>
#include <cmath>

#ifdef USE_BOOST
    #include <boost/math/distributions/binomial.hpp>
    using namespace::boost::math;
#else
    #include <gsl/gsl_cdf.h>
#endif

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"

//External includes
#include "common/murmur3.h"
#include "common/kseq.h"
#include "common/prettyprint.hpp"

namespace skch
{
  /**
   * @namespace skch::Stat
   * @brief     Implements utility functions that involve statistical computation
   */
  namespace Stat
  {
    /**
     * @brief         jaccard estimate to mash distance
     * @param[in] j   jaccard estimate
     * @param[in] k   kmer size 
     * @return        mash distance [0.0 - 1.0]
     */
    inline float j2md(float j, int k)
    {
      if(j == 0)
        return 1.0; //jaccard estimate 0 -> 1.0 mash distance

      if(j == 1)
        return 0.0; //jaccard estimate 1 -> 0.0 mash distance

      float mash_dist = (-1.0 / k) * log(2.0 * j/(1+j) );
      return mash_dist;
    }

    /**
     * @brief         mash distance to jaccard estimate
     * @param[in] d   mash distance [0.0 - 1.0]
     * @param[in] k   kmer size 
     * @return        jaccard estimate 
     */
    inline float md2j(float d, int k)
    {
      float jaccard = 1.0 / (2.0 * exp( k*d ) - 1.0);
      return jaccard;
    }

    /**
     * @brief               Given a distance d, compute the lower bound on d within required confidence interval 
     * @details             If a given match has distance d in the L2 stage, we compare its lower distance bound 
     *                      against the assumed cutoff to decide its significance. This makes the mapping algorithm
     *                      more sensitive to true alignments
     * @param[in]   d       calculated mash distance
     * @param[in]   s       sketch size
     * @param[in]   k       kmer size
     * @param[in]   ci      confidence interval [0-1], example 0.9 implies 90% confidence interval
     * @return              computed lower bound on d within 'ci' confidence interval
     */
    inline float md_lower_bound(float d, int s, int k, float ci)
    {
      //One side interval probability
      float q2 = (1.0 - ci)/2;

      //Computing count of sketches using confidence interval
#ifdef USE_BOOST
      
      //Inverse binomial 
      int x = quantile(complement(binomial(s, md2j(d,k)), q2));

#else   
      //GSL 
      int x = std::max( int(ceil(s * md2j(d,k))), 1 );    //Begin search from jaccard * s
      while(x <= s)
      {
        //probability of having x or more shared sketches
        double cdf_complement = gsl_cdf_binomial_Q(x-1, md2j(d,k), s);

        if (cdf_complement < q2)
        {
          x--;  //Last guess was right
          break;
        }

        x++;
      }
#endif

      float jaccard = float(x) / s;
      float low_d = j2md(jaccard, k);
      return low_d; 
    }

    /**
     * @brief                 Estimate minimum number of shared sketches to achieve the desired identity
     * @param[in] s           sketch size
     * @param[in] k           kmer size
     * @param[in] identity    percentage identity [0-100]
     * @return                minimum count of hits
     */
    inline int estimateMinimumHits(int s, int k, float perc_identity)
    {
      //Compute the estimate
      float mash_dist = 1.0 - perc_identity/100.0;
      float jaccard = md2j(mash_dist, k);

      //function to convert jaccard to min hits
      //Atleast these many minimizers should match for achieving the required jaccard identity
      int minimumSharedMinimizers = ceil (1.0 * s * jaccard); 

      return minimumSharedMinimizers;
    }

    /**
     * @brief                 Estimate minimum number of shared sketches 
     *                        s.t. upper bound identity is >= desired identity
     *                        Upper bound is computed using the 90% confidence interval
     * @param[in] s           sketch size
     * @param[in] k           kmer size
     * @param[in] identity    percentage identity [0-100]
     * @return                count of min. shared minimizers
     */
    inline int estimateMinimumHitsRelaxed(int s, int k, float perc_identity)
    {
      // The desired value has be between [0, min  s.t. identity >= perc_identity]
      auto searchRange = std::pair<int, int>( estimateMinimumHits(s, k, perc_identity) , 0);

      int minimumSharedMinimizers_relaxed = searchRange.first;

      for(int i = searchRange.first ; i >= searchRange.second; i--)
      {
        float jaccard = 1.0 * i/s;
        float d = j2md(jaccard, k);

        float d_lower = md_lower_bound(d, s, k, 0.9);

        //Upper bound identity
        float id_upper = 100.0 * (1.0 - d_lower);

        //Check if it satisfies the criteria
        if(id_upper >= perc_identity)
          minimumSharedMinimizers_relaxed = i;
        else
          break;    //Stop the search
      }

      return minimumSharedMinimizers_relaxed;
    }

    /**
     * @brief                     calculate p-value for a given alignment identity, sketch size..
     * @param[in] s               sketch size
     * @param[in] k               kmer size
     * @param[in] alphabetSize    alphabet size
     * @param[in] identity        mapping identity cut-off
     * @param[in] lengthQuery     query length
     * @param[in] lengthReference reference length
     * @return                    p-value
     */
    inline double estimate_pvalue (int s, int k, int alphabetSize, 
        float identity, 
        int lengthQuery, uint64_t lengthReference)
    {
      //total space size of k-mers
      double kmerSpace = pow(alphabetSize, k);

      //probability of a kmer match by random in |query| sized sequence 
      double pX, pY; 
      pX = pY = 1. / (1. + kmerSpace / lengthQuery);

      //Jaccard similarity of two random given sequences
      double r = pX * pY / (pX + pY - pX * pY);

      int x = estimateMinimumHitsRelaxed(s, k, identity);

      //P (x or more minimizers match)
      double cdf_complement;
      if(x == 0)
      {
        cdf_complement = 1.0;
      }
      else
      {
#ifdef USE_BOOST
      cdf_complement = cdf(complement(binomial(s, r), x-1));
#else
      cdf_complement =  gsl_cdf_binomial_Q(x-1, r, s);
#endif
      }

      double pVal = lengthReference * cdf_complement;

      return pVal;
    }

    /**
     * @brief                     calculate minimum window size for sketching that satisfies
     *                            the given p-value threshold
     * @param[in] pValue_cutoff   cut off p-value threshold
     * @param[in] k               kmer size
     * @param[in] alphabetSize    alphabet size
     * @param[in] identity        mapping identity cut-off
     * @param[in] lengthQuery     query length
     * @param[in] lengthReference reference length
     * @return                    optimal window size for sketching
     */
    inline int recommendedWindowSize(double pValue_cutoff,
        int k, int alphabetSize,
        float identity,
        int lengthQuery, uint64_t lengthReference)
    {
      //Push all the sketch values that we should try out in a vector
      //{1, 2, 5, 10, 20, 30...}
      std::vector<int> potentialSketchValues{1,2,5};
      for(int i = 10; i < lengthQuery; i+= 10)
        potentialSketchValues.push_back(i);

      int optimalSketchSize;

      for(auto &e : potentialSketchValues)
      {
        //Compute pvalue
        double pVal = estimate_pvalue(e, k, alphabetSize, identity, lengthQuery, lengthReference);

        //Check if pvalue is <= cutoff
        if(pVal <= pValue_cutoff)
        {
          optimalSketchSize = e;
          break;
        }
      }
      
      int w =  2.0 * lengthQuery/optimalSketchSize;

      // 1 <= w <= lengthQuery
      return std::min( std::max(w,1), lengthQuery);
    }
  }
}

#endif
