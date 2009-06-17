#ifndef _NCHOOSEK_BITSET_HPP_
#define _NCHOOSEK_BITSET_HPP_
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
 * \file NChooseKBitset.hpp
 * \author Simen Kvaal
 * \date 6-13-08
 *
 * \brief NChooseKBitset is a class for computing subsets of K
 * elements out of N, using a bitset to represent a subset. A set bit
 * represents the presence of an element in the set. The class
 * interface is somewhat sloppy, and *should* be tightened.
 *
 * The class allows for a constraint in the choice of subsets
 * (bitsets) actually stored. The feature must be turned on with
 * setUseConstraint(true). The default constraint is trivial --
 * everything is stored. To use a differenc constraint, derive a
 * subclass and re-implement constraint().
 *
 * If this is not flexible enough, use the template version instead,
 * where the template argument is a class that has the constraint
 * function.
 *
 */


#include <vector>
#include <cassert>
#include "bitsethack.hpp"

using namespace std;

/// \brief Class that helps with computing subsets of K elements
/// from a collection of N elements.
template<int N0>
class NChooseKBitset
{

  protected:
    bool use_constraint;      ///< Indicate whether we wish to use a constraint or not.

    /// \brief Recursive helper function for choose.
    ///
    /// In the code, notice from the commented out code,
    /// that it is easy to introduce a constraint on the
    /// choices made. This can be useful ...
    ///
    void generate(int N, int K, std::vector<std::bitset<N0> >& words, size_t& c, size_t& calls, std::bitset<N0> prefix = 0);


  public:
    int N; ///< Number of elements in set.
    int K; ///< Number of elements to choose.
    std::vector<bitset<N0> > subset; ///< Holds the subset after choose() is called.

    /// \brief Default constructor
    NChooseKBitset()
    {
      N = 0;
      K = 0;
      subset.resize(0);
      use_constraint = false;
    }

    /// \brief Turn on/off usage of constraint function in choose() and generate().
    /// \param x On or off?
    void setUseConstraint(bool x)
    {
      use_constraint = x;
    }

    /// \brief Find all subsets
    virtual void choose()
    {

      subset.resize(0); // Empty contents, if any

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
    
    /// \brief A constraint function for selecting interesting subset of subset generated. Called by choose().
    ///  
    /// Re-implement in derived class in order to choose a different interesting subset. This default
    /// implementation is trivial -- returns true. 
    ///
    /// \param x Candidate to check against constraint
    virtual bool constraint(const bitset<N0>& x)
    {
      return true;
    }



};

template<int N0>
void NChooseKBitset<N0>::generate(int N, int K, std::vector<std::bitset<N0> >& words, size_t& c, size_t& calls, bitset<N0> prefix)
{

  ++calls;
  std::bitset<N0> word = 0;

  if (N == K)
  {

    // N == K, only one possible word 111....1, N=K bits set.

    // Add word/state if it fulfills constraint.
    word = prefix;
    for (size_t k = 0; k<(size_t)N; ++k)
      word.set(k, 1);
    //words.push_back(word);
    //c++;
    if ( (!use_constraint) || (use_constraint && constraint(word)) )
    {
      words.push_back(word);
      c++;
    }
  }
  else
    if (K == 0)
    {
      // N == K, only one possible word 000....0, N bits set.
 
      // Add word/state if it fulfills constraint.
      word = prefix;
      //words.push_back(prefix);
      //c++;
      if ( (!use_constraint) || (use_constraint && constraint(word)) )
      {
	words.push_back(word);
  	c++;
      }
    }
    else
      if (N > 0)
      {
	// add 0 to prefix, and process N-1 levels, still K particles.
	generate(N-1,K,words,c,calls,prefix);
	// add 1 to prefix, and process N-1 levels, K-1 particles.
	prefix.set(N-1, 1);
	generate(N-1,K-1,words,c,calls,prefix) ;
      }

}

/// \brief Template version of a NChooseK-class with constraint.
///
/// The template argument must implement the constraint function.
///
template<class Caller, int N0>
class NChooseKBitsetExternal : public NChooseKBitset<N0>
{
  private:
    Caller& c_obj;   ///< Reference to object that has the actual constraint function.

  public:
    /// \brief Default constructor
    /// \param c   Reference to object of type Caller that has the constraint() function.
    NChooseKBitsetExternal(Caller& c) : NChooseKBitset<N0>(), c_obj(c)  {  }

    /// \brief Constraint function
    /// \param x Candidate to check agains constraint.
    virtual bool constraint(const bitset<N0>& x)
    {
      return c_obj.constraint(x);
    }

};

#endif
