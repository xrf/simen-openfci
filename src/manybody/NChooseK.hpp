#ifndef _NCHOOSEK_HPP_
#define _NCHOOSEK_HPP_

//
// Copyright (c) 2008 Simen Kvaal
//
// This file is part of OpenFCI.
//
// OpenFCI is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// OpenFCI is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with OpenFCI. If not, see <http://www.gnu.org/licenses/>.
//

/**
 * \file NChooseK.hpp
 * \author Simen Kvaal
 * \date 6-13-08
 *
 * \brief NChooseK is a class for computing subsets of K elements out of N.
 * 
 * We simply interpret the subsets as bit patterns: Out of N bits, K can be lit.
 * NChooseK<T> does this using integers of type T as bit patterns.
 *
 * The very similar class NChooseKBitset<N> uses std::bitset<N> instead. This should
 * really be preferred, since in general it will be just as efficient but much more general.
 *
 */


#include <vector>
#include <cassert>
#include <climits>

using namespace std;

/// \brief Class that helps with computing subsets of K elements
/// from a collection of N elements.
template<class T>
class NChooseK
{
  protected:

    /// \brief Recursive helper function for choose.
    ///
    /// In the code, notice from the commented out code,
    /// that it is easy to introduce a constraint on there
    /// choices made. This can be useful ...
    void generate(int N, int K, std::vector<T>& words, size_t& c, size_t& calls, T prefix = T(0));


  public:
    int N; ///< Number of elements in set.
    int K; ///< Number of elements to choose.
    std::vector<T> subset; ///< Holds the subset after choose() is called.

    /// \brief Default constructor
    NChooseK()
    {
      N = 0;
      K = 0;
      subset.resize(0);
    }

    /// \brief Find all subsets
    void choose()
    {

      // Make sure that typename T holds N bits.
      assert(sizeof(T) * CHAR_BIT >= (size_t)N);

      subset.resize(0); // Empty contents, if any
      subset.resize(N); // Allocate memory; tentative size...

      size_t c = 0; // Counts the number of states during recursion
      size_t calls = 0; // Counts the number of calls during recursion

      //
      // Call the recursive generation function.
      //
      // Result:
      //   subset = vector of created words
      //   c = number of words created
      //   calls = calls to recursive funcion. not used ...
      //
      generate(N, K, subset, c, calls);

      //
      // Erase surplus elements.
      //
      subset.resize(c);
    }
    


};

template<class T>
void NChooseK<T>::generate(int N, int K, vector<T>& words, size_t& c, size_t& calls, T prefix)
{

  ++calls;
  T word = 0;

  if (N == K)
  {

    // N == K, only one possible word 111....1, N=K bits set.

    // Add word/state if it fulfills constraint.
    word = prefix;
    word |= T(1) << N;
    --word;
    //if (constraint(word)) 
    {
       if (c == words.size())
	 words.resize(2*words.size()); 
       words[c] = word;
       c++;
    }
  }
  else
    if (K == 0)
    {
      // N == K, only one possible word 000....0, N bits set.
 
      // Add word/state if it fulfills constraint.
      word = prefix;
      //if (constraint(word)) 
      {
 	if (c == words.size())
 	  words.resize(2*words.size()); 
 	words[c] = word;
 	c++;
      }
    }
    else
      if (N > 0)
      {
	// add 0 to prefix, and process N-1 levels, still K particles.
	generate(N-1,K,words,c,calls,prefix);
	// add 1 to prefix, and process N-1 levels, K-1 particles.
	generate(N-1,K-1,words,c,calls,prefix + ( T(1) << (N-1U) )  ) ;
      }

}


#endif
