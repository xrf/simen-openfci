#ifndef _BITSET_SLATER_HPP_
#define _BITSET_SLATER_HPP_

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

#include "bitsethack.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>


namespace manybody
{


/**
 *
 * \file BitsetSlater.hpp
 * \author Simen Kvaal
 * \date 9-9-08 (prev. 8-19-08)
 *
 * \brief Definition of class template BitsetSlater<int N>; a Slater determinant
 * based on the std::bitset<N> class. 
 *
 */

  /// \brief Class definition of BitsetSlater.
  ///
  /// This class encapsulates the occupation number formalism version
  /// of a Slater determinant.  It basically \em is a std::bitset<N>
  /// instance with additional functionality.
  ///
  /// Let A be the number of particles in a system, i.e., in a Slater
  /// determinant.  In the occupation number formalism, instead of
  /// specifying A distinct orbitals (integer indices), one specifies
  /// single integer (huge in general!) with A bits set, leaving the
  /// others blank. So, if 3 particles resides in orbitals 4, 5, and
  /// 10, the occupation number formalism specification of the slater
  /// determinant is phi = |...00010000110000>.
  ///
  /// Notice, that while the oroginal Slater determinant specification
  /// was not unique (it had 3! = 6 different representations with
  /// different signs), the occupation number formalism specifies a
  /// \em unique integer with 3 bits set. Conversely, all integers
  /// with A bits set is a legal A-body Slater determinant.
  ///
  /// When we speak of 'bits' in the documentation here, we speak of
  /// orbitals ('bit numbers') and particles in orbitals ('set bits').
  ///
  /// <b>Implementation note:</b> We chose to have the std::bitset<N>
  /// instance as a private member for clarity, but it would be better
  /// to derive it. Most importantly, since std::bitset<N> does not
  /// have an operator<() defined, and this is needed in \em some
  /// applications, but not essential, we had to hack the original
  /// header file.  Moreover, the bitset_base class defined in
  /// <bitset> is private to bitset, so we cannot access it
  /// directly... In my opininon, the STL should define operator<()
  /// etc for bitsets.
  ///
  template<int N>
  class BitsetSlater 
  {
    private:
      std::bitset<N> bits;   ///< Sole private member.

    public:
      /// \brief Set a bit (place particle)
      /// \param alpha    Position of bit (orbital number)
      void set(size_t alpha)
      {
	bits[alpha] = 1;
      }
      /// \brief Clear a bit (remove particle)
      /// \param alpha    Position of bit (orbital number)
      void unset(size_t alpha)
      {
	bits[alpha] = 0;
      }

      /// \brief Read a bit (check for particle).
      /// \param alpha    Position of bit (orbital number)
      int get(size_t alpha) const
      {
	return bits[alpha];
      }

      /// \brief Compute the memory in bytes used by the BitsetSlater object.
      /// \return Number of bytes
      size_t mem() const
      {
	return sizeof(std::bitset<N>);
      }

      /// \brief Count the number of bits set, i.e., the number of particles.
      /// \return Number of bits set (particles)
      size_t count() const
      {
	return bits.count();
      }

      /// \brief Count bits set (particles) in the range [a, b) (i.e., exclusive b).
      ///
      /// \param a     First endpoint
      /// \param b     Last endpoint
      /// \return      Number of bits (particles)
      size_t count_between(size_t a, size_t b)
      {
	size_t sum = 0;
	for (size_t k = a; k < b; ++k)
	  if (bits[k]) 
	    sum++;
	return sum;
      }

      /// \brief Count the number of even-numbered and odd-numbered bits set.
      ///
      /// This function is typically used whenever even spins are 'up'
      /// and odd are 'down'.
      ///
      /// \param even Output: number of even-numbered bits set
      /// (starting with 0) 
      /// \param odd Output: number of odd-numbered
      /// bits set (starting with 1)
      void count_evenodd(int& even, int& odd) const 
      {
	even = 0;
	odd = 0;
	for (size_t k = 0; k<N; ++k)
	  if (k%2)
	    odd += bits[k];
	  else
	    even += bits[k];
	
      }

      /// \brief Analyze the difference between two Slater dets, i.e.,
      /// how many and which occupied orbitals are different in the two.
      ///
      /// \param phi BitsetSlater to compare with
      /// \param alpha output; indices of bits in *this different from phi
      /// \param beta output; indices of bits in phi different from *this
      /// \param p output; bit numbers of the bits in alpha, i.e.,
      /// number of transpositions needed to bring creation operator
      /// to front.
      /// \param q output; bit numbers of the bits in beta
      /// \param n input; maximum number of differing bits to detect
      /// \return 0, 1, 2, 3, ... n, (or n+1 for >n; then alpha, beta, p and q are not defined) bits different.
      int analyze_difference(const BitsetSlater& phi, 
			       size_t alpha[], size_t beta[], 
			       size_t p[], size_t q[], int n)
      {

	std::bitset<N> x = bits ^ phi.bits; // differing bits.
	int ndiff = (x.count() >> 1);
	
	// if more than n bits differ, return n+1.
	if (ndiff > n)
	  return n+1;

	// if no bits differ, just exit.
	if (ndiff == 0)
	  return 0;

	// this is for work.
	std::bitset<N> the_bit;

	// if one bit differ, find it (them)
	for (int num = 0; num < ndiff; ++num)
	{
	  the_bit = bits & x; // has exactly 1, 2 or 3 bit(s) lit.
	  alpha[num] = the_bit._Find_first(); // find the first bit's index
	  the_bit.set(); // now, set all bits
	  the_bit >>= N - alpha[num]; // shift left, so that bits 0...(alpha - 1) are 1's, the rest zeros
	  the_bit &= bits;  // & with original pattern, and ...
	  p[num] = the_bit.count(); // ... count to get position of the differing bit!
	  
	  // rinse and repeat ... (for phi)
	  the_bit = phi.bits & x; 
	  beta[num] = the_bit._Find_first();
	  the_bit.set();
	  the_bit >>= N - beta[num];
	  the_bit &= phi.bits;
	  q[num] = the_bit.count();

	  // Prepare for next bit pair
	  x.set(alpha[num], 0);
	  x.set(beta[num], 0);
	}

	return ndiff;

      }



      /// \brief "Create a particle" at orbital alpha. Returns sign change of operation.
      ///
      /// \param alpha   Position of particle to create
      /// \return -1 if odd permutation is reqiured, 1 if even permutation, 0 if particle already exists.
      int create(size_t alpha)
      {
	// Retun 0 if alredy set.
	if (bits[alpha])
	  return 0;

	// Count number of bits set before alpha.
	size_t num_before = count_between(0, alpha);
	// Compute sign.
	int sign = 1;
	if (num_before % 2 == 1)
	  sign = -1;

	// Finally, set bit.
	bits[alpha] = 1;

	return sign;
      }


      /// \brief "Destroy a particle" at orbital alpha. Returns sign change of operation.
      ///
      /// \param alpha   Position of particle to destroy
      /// \return -1 if odd permutation is reqiured, 1 if even permutation, 0 if particle didn't already exist.
      int annihilate(size_t alpha)
      {
	// Retun 0 if not set.
	if (!bits[alpha])
	  return 0;

	// Count number of bits set before alpha.
	size_t num_before = count_between(0, alpha);
	// Compute sign.
	int sign = 1;
	if (num_before % 2 == 1)
	  sign = -1;

	// Finally, clear bit.
	bits[alpha] = 0;

	return sign;
      }



      /// \brief Compute binary number representation of Slater determinant
      ///
      /// \param K    Number of bits to use; defaults to all bits.
      /// \return     String sequence of 0's and 1's
      std::string to_binary(int K = N) const
      {
      
	return bits.to_string().substr(N-K,N);
      }

      /// \brief Compute string representation of Slater determinant;
      /// each occupied orbital becomes an integer with a space after.
      /// \return String, example: 01001 --> "0 3 "
      std::string to_string() const
      {
	std::stringstream ss;
	for (size_t k = 0; k<N; ++k)
	  if (bits[k])
	    ss << k << " ";

	return ss.str();
    
      }

      /// \brief Compute a vector containing the positions of all the set bits, i.e.,
      /// occupied orbitals. 
      ///
      /// \param ones Reference to std::vector<size_t> where the
      /// result is stored. Resized automatically.
      void find_ones(std::vector<size_t>& ones)
      {
	ones.resize(count());
	size_t pos = 0;
	for (size_t k = 0; k<N; ++k)
	  if (bits[k])
	    ones[pos++] = k;

      }

      /// \brief Write the slater determinant to cerr. 
      void print() const
      {
	std::cerr << to_string() << std::endl;
      }

      /// \brief Comparison operator
      bool operator==(const BitsetSlater<N>& c2) const
      {
	return bits == c2.bits;
      }


#ifdef BITSET_HACK
      /// \brief Comparison operator. Requires the hacked version of std::bitset.
      bool
      operator<(const BitsetSlater<N>& x) const
      { return bits < x.bits; }

      /// \brief Compute the last set bit, i.e., the particle in the
      /// orbital with the highest alpha, i.e., the MSB of the bitset.
      size_t last_index() const
      {
	for (int k = N-1; k>=0; --k)
	  if (bits[k])
	    return k;

	return 0;
      }

#endif

      /// \brief Bitwise XOR of bit patterns
      void bitwise_xor(const BitsetSlater<N>& c)
      {
	bits ^= c.bits;
      }
      /// \brief Bitwise AND of bit patterns
      void bitwise_and(const BitsetSlater<N>& c)
      {
	bits &= c.bits;
      }

      /// \brief Copy from a std::bitset<N>r.
      BitsetSlater<N>& operator=(const std::bitset<N>& x)
      {
	bits = x;
	return *this;
      }
    
  };


  /// \brief Bitwise XOR of two BitsetSlater instances
  /// \param c1       First argument
  /// \param c2       Second argument
  /// \param result   Result
  template<int N>
  inline void bitwise_xor(const BitsetSlater<N>& c1, const BitsetSlater<N>& c2, BitsetSlater<N>& result)
  {

    result = c1;
    result.bitwise_xor(c2);

  }

  /// \brief Bitwise AND of two BitsetSlater instances
  /// \param c1       First argument
  /// \param c2       Second argument
  /// \param result   Result
  template<int N>
  inline void bitwise_and(const BitsetSlater<N>& c1, const BitsetSlater<N>& c2, BitsetSlater<N>& result)
  {

    result = c1;
    result.bitwise_and(c2);

  }

} // end of namespace manybody


#endif
