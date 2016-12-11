/**
 * @file    commonFunc.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMMON_FUNC_HPP 
#define COMMON_FUNC_HPP

#include <vector>
#include <algorithm>
#include <deque>
#include <cmath>
#include <fstream>

//Own includes
#include "map_parameters.hpp"

//External includes
#include "murmur3.h"
#include "kseq.h"
#include "prettyprint.hpp"

namespace skch
{
  /**
   * @namespace skch::CommonFunc
   * @brief     Implements frequently used common functions
   */
  namespace CommonFunc
  {
    //seed for murmerhash
    const int seed = 42;

    /**
     * @brief   reverse complement of kmer (borrowed from mash)
     */
    inline void reverseComplement(const char * src, char * dest, int length) 
    {
      for ( int i = 0; i < length; i++ )
      {    
        char base = src[i];

        switch ( base )
        {    
          case 'A': base = 'T'; break;
          case 'C': base = 'G'; break;
          case 'G': base = 'C'; break;
          case 'T': base = 'A'; break;
          default: break;
        }    

        dest[length - i - 1] = base;
      }    
    }

    template <typename KSEQ>
      inline void makeUpperCase(KSEQ *seq)
      {
        for ( int i = 0; i < seq->seq.l; i++ )
        {
          if (seq->seq.s[i] > 96 && seq->seq.s[i] < 123)
          {
            seq->seq.s[i] -= 32;
          }
        }
      }

    /**
     * @brief   hashing kmer string (borrowed from mash)
     */
    inline hash_t getHash(const char * seq, int length)
    {
      char data[16];
      MurmurHash3_x64_128(seq, length, seed, data);

      hash_t hash;

      hash = *((hash_t *)data);

      return hash;
    }

    /**
     * @brief       compute winnowed minimizers from a given sequence and add to the index
     * @param[out]  minimizerIndex  minimizer table storing minimizers and their position as we compute them
     * @param[in]   seq             kseq fasta/q parser
     * @param[in]   kmerSize
     * @param[in]   windowSize
     * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
     */
    template <typename T, typename KSEQ>
      inline void addMinimizers(std::vector<T> &minimizerIndex, KSEQ *seq, int kmerSize, 
          int windowSize,
          int alphabetSize,
          seqno_t seqCounter)
      {
        /**
         * Double-ended queue (saves minimum at front end)
         * Saves pair of the minimizer and the position of hashed kmer in the sequence
         * Position of kmer is required to discard kmers that fall out of current window
         */
        std::deque< std::pair<MinimizerInfo, offset_t> > Q;

        makeUpperCase(seq);

        //length of the sequencd
        offset_t len = seq->seq.l;

        //Compute reverse complement of seq
        char *seqRev = new char[len];

        if(alphabetSize == 4) //not protein
          CommonFunc::reverseComplement(seq->seq.s, seqRev, len);

        for(offset_t i = 0; i < len - kmerSize + 1; i++)
        {
          //The serial number of current sliding window
          //First valid window appears when i = windowSize - 1
          offset_t currentWindowId = i - windowSize + 1;

          //Hash kmers
          hash_t hashFwd = CommonFunc::getHash(seq->seq.s + i, kmerSize); 
          hash_t hashBwd;

          if(alphabetSize == 4)
            hashBwd = CommonFunc::getHash(seqRev + len - i - kmerSize, kmerSize);
          else  //proteins
            hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

          //Consider non-symmetric kmers only
          if(hashBwd != hashFwd)
          {
            //Take minimum value of kmer and its reverse complement
            hash_t currentKmer = std::min(hashFwd, hashBwd);

            //Check the strand of this minimizer hash value
            auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

            //If front minimum is not in the current window, remove it
            while(!Q.empty() && Q.front().second <=  i - windowSize)
              Q.pop_front();

            //Hashes less than equal to currentKmer are not required
            //Remove them from Q (back)
            while(!Q.empty() && Q.back().first.hash >= currentKmer) 
              Q.pop_back();

            //Push currentKmer and position to back of the queue
            //0 indicates the dummy window # (will be updated later)
            Q.push_back( std::make_pair(
                  MinimizerInfo{currentKmer, seqCounter, 0, currentStrand},
                  i)); 

            //Select the minimizer from Q and put into index
            if(currentWindowId >= 0)
            {
              //We save the minimizer if we are seeing it for first time
              if(minimizerIndex.empty() || minimizerIndex.back() != Q.front().first)
              {
                //Update the window position in this minimizer
                //This step also ensures we don't re-insert the same minimizer again
                Q.front().first.wpos = currentWindowId;     
                minimizerIndex.push_back(Q.front().first);
              }
            }
          }
        }

#ifdef DEBUG
        std::cout << "INFO, skch::CommonFunc::addMinimizers, inserted following minimizers: ";
        for(auto &e: minimizerIndex) std::cout << "\n" << e;
        std::cout << std::endl;
#endif

        delete [] seqRev;
      }

    /**
     * @brief       overloaded function for case where seq. counter does not matter
     */
    template <typename T, typename KSEQ>
      inline void addMinimizers(std::vector<T> &minimizerIndex, KSEQ *seq, int kmerSize,
          int windowSize, 
          int alphabetSize)
      {
        addMinimizers(minimizerIndex, seq, kmerSize, windowSize, alphabetSize, 0);
      }

   /**
     * @brief           Functor for comparing tuples by single index layer
     * @tparam layer    Tuple's index which is used for comparison
     * @tparam op       comparator, default as std::less
     */
    template <size_t layer, template<typename> class op = std::less>
      struct TpleComp
      {
        //Compare two tuples using their values
        template<typename T>
          bool operator() (T const &t1, T const &t2)
          {
            return op<typename std::tuple_element<layer, T>::type>() (std::get<layer>(t1), std::get<layer>(t2));
          }
      };

    /**
     * @brief                   computes the total size of reference in bytes
     * @param[in] refSequences  vector of reference files
     * @return                  total size
     */
    inline uint64_t getReferenceSize(const std::vector<std::string> &refSequences)
    {
      uint64_t count = 0;

      for(auto &f : refSequences)
      {
        //Open the file as binary, and set the position to end
        std::ifstream in(f, std::ifstream::ate | std::ifstream::binary);

        //the position of the current character
        count += (uint64_t)(in.tellg());
      }

      return count;
    }
  }
}

#endif
