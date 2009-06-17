#ifndef _CSF_HPP_
#define _CSF_HPP_

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

#include <vector>
#include <algorithm>
#include <map>

#include "Slater.hpp"
#include "linalg.hpp"


namespace manybody
{



/**
 * \file CSF.hpp
 * \author Simen Kvaal
 * \date 9-15-2008, update 10-7-2008
 *
 * \brief Declaration of various classes and functions for compuation
 * of configurational state functions (CSFs).
 *
 * Configurational State Functions (CSFs) constitute a common basis
 * for Sz and S^2 spin operators.  They are constructed by taking a
 * proper linear combination of Slater determinants whose orbitals are
 * a single-particle basis, diagonal in Sz. We assume that two
 * consecutive orbitals in a Slater determinants (i.e., consecutive
 * bits) is the same \em spatial orbital with spin up in the even bit,
 * and down in the odd bit. (First bit is the LSB, number zero, which
 * is "spin up".)
 *
 * The funcitonality of the class CsfMachine is basically to take a
 * given basis of Slater determinants and compute the proper linear
 * combinations that produce the S^2-basis functions with a given
 * eigenvalue. More on this in the class documentation.
 *
 * For simple applications, the global function computeCSFBlocks() can
 * be used, and the use of CsfMachine can be skipped altogether. For
 * applications where spin is undefined, consider
 * computeTrivialCSFBlocks().
 *
 */


  /// \brief Compute Clebsh Gordan coefficients for n spins up and m
  /// spins down.
  ///
  /// The result is the columns of CG, and s; a vector of the
  /// corresponding total spin of each column in CG. Thus, CG(:,j) is
  /// the vector of coefficients for an S^2-eigenfunction of eigenvalue
  /// s(j)*(s(j)+1).
  ///
  /// \param basis Upon exit: a vector of unsigned ints (bit patterns)
  /// holding the spin configurations used as basis. In total
  /// binom(n+m,n) bit patterns.
  /// \param CG Upon exit: Matrix of Clebsch-Gordan
  /// coefficients. Unitary matrix. Each column is the coefficients
  /// needed to construct a S^2-eigenfunction.
  /// \param s Upon exit: Vector of S^2-eigenvalues of the coumns of CG.
  /// \param n Input: Number of spins up (ones in the bit patterns in
  /// basis)
  /// \param m Input: Number of spins down (zeroes in the bit patterns
  /// in basis)
  void computeCGCoeffs(dense_matrix& CG, dense_vector& s, std::vector<unsigned int>& basis, int n, int m);



  /// \brief This structure is used internally in CsfMachine. 
  ///
  /// equivalence_class represents all the Slater determinants in the
  /// basis that has the same singly and doubly occupied levels.  When
  /// this struct is used, the basis is already ordered according to
  /// equivalence class, so that the only information needed is the end
  /// points of the class and the spin configurations.
  typedef struct
  {
      size_t begin;                ///< First SD
      size_t end;                  ///< Last SD + 1
      short int up;                ///< Number of spins up
      short int down;              ///< Number of spins down
      short int n_doubly;          ///< Number of doubly occupied levels.
      std::vector<unsigned int> spins;  ///< Spin configuration of relevant particles.
  } equivalence_class;


  /// \brief This struct holds a set of CSFs in the sense of: Coefficients 
  /// of linear combinations of SDs that give an eigenfunction of S^2. (The 
  /// eigenvalue is not specified in this structure.)
  typedef struct
  {
      size_t sd_begin;    ///< First SD
      size_t sd_end;      ///< Last SD + 1
      size_t n_csf;       ///< Number of CSFs in the block
      dense_matrix coeffs; ///< CG coeffs for states; (sd_end-sd_begin) x n_csf.
  } csf_block;

 
  /// \brief Compute trivial CSFs from a vector of SDs. 
  ///
  /// Stand-alone function for treating spinless problems, where every 
  /// SD constitute a S^2-eigenfunction. Thus, no need for computing CG coefficients
  /// and so forth. The result of running this function is vector of csf_blocks
  /// of same length as the basis. Each block has exactly one CSF, mathematically identical
  /// to a SD in the basis. So, <tt>sd_basis[k]</tt> is mathematically identical to <tt>CSF[k]</tt>.
  /// 
  /// \param sd_basis    Input: Vector of Slater determinants.
  /// \param CSF         Output: Vector of CSF blocks, each block has 1 CSF.
  inline void computeTrivialCSFBlocks(const std::vector<Slater>& sd_basis, std::vector<csf_block>& CSF)
  {

    CSF.resize(0);
    csf_block temp;
    dense_matrix just1(1,1);
    just1(1,1) = 1;
    for (size_t k = 0; k<sd_basis.size(); ++k)
    {
      temp.sd_begin = k;
      temp.sd_end = k+1;
      temp.n_csf = 1;
      temp.coeffs = just1;
      CSF.push_back(temp);
    }

  }


 
  /// \brief Short-cut function that takes a Slater determinant basis (a std::vector<Slater>)
  /// And produces \em all CSF blocks of a given S value.
  ///
  /// While CsfMachine only handles a single Sz value, this function allows several Sz values,
  /// by picking out the individual constant Sz subspaces and building CSF bases for each.
  ///
  /// \param sd_basis    Basis of Slater determinants.
  /// \param csf_basis   CSF block vector to store result.
  /// \param S           The value for the total spin.
  void computeCSFBlocks(std::vector<Slater>& sd_basis, std::vector<csf_block>& csf_basis, int S);



  /// \brief This class computes CSF basis functions.
  ///
  /// Let a <tt>sd_basis</tt> be a std::vector<Slater>, where all elements have
  /// the same Sz eigenvalue (number of odd bits set - number of even bits set)/2. Let a desired
  /// S^2-eigenvalue s/2 be given. Then, this class produces a std::vector<csf_block>
  /// containing all the common eigenvalues of Sz and S^2.
  ///
  /// Example usage:
  /// \code
  /// // Example goes here.
  /// \endcode
  ///
  class CsfMachine
  {
    private:
      int sz;                                    ///< The value for 2*Sz for the basis.
      int nup;                                   ///< Number of up spins
      int ndown;                                 ///< Number of down spins
      std::vector<Slater> sd_basis;                   ///< Slater determinant basis.
      std::vector<int> class_number;                  ///< Class number for each SD.
      std::vector<equivalence_class> eq_classes;      ///< Holds the equivalence classes
      std::vector<dense_matrix> cg_coeffs;            ///< Holds the needed Clebsh Gordan coefficients.
      std::vector<std::vector<unsigned int> > cg_basis;    ///< Holds the spin basis used for computing CG coeffs.
      std::vector<int> cg_n;                          ///< not really very useful, but holds the number of spins for each cg block.
      std::vector<dense_vector> s_values;             ///< Holds the S^2 eigenvalues for each CG coeff matrix.

    protected:

      /// \brief Compute the singly occupied levels of phi.
      ///
      /// A level is two consecutive bits, one even bit and one odd, i.e., a pair of orbitals.
      /// This function checks all levels (pairs of orbitals) and decides which are singly occupied.
      /// It returns a new slater with each bit set indicating a singly occupied orbital.
      ///
      /// Example: if <tt>phi = |00110110></tt>, then the two first levels (rightmost two pairs of bits) 
      /// are singly occupied, so that <tt>singly = |0011></tt>. (Rightmost bit is first bit.)
      ///
      /// \param phi     Slater determinant to work on
      /// \param singly  "Condensed" Slater determinants (no spin) of singly occupied levels. Really only used as a bit pattern.
      void findSinglyOccupied(const Slater& phi, Slater& singly)
      {
	// Get the largest index, and round down to nearest even number.
	orbital_t alphamax = ((phi.last_index()>>1)<<1) + 2;
	// Loop though all indices, and check for singly occupied positions.
	for (orbital_t alpha = 0; alpha<alphamax; alpha+=2)
	{
	  bool x1 = phi.get(alpha) ? true : false;
	  bool x2 = phi.get(alpha+1) ? true : false;
	  if ( (x1 && (!x2)) || (x2 && (!x1)) )
	  {
	    singly.set(alpha >> 1);
	  }
	}

      }

      /// \brief Compute the doublyy occupied levels of phi.
      ///
      /// A level is two consecutive bits, one even bit and one odd, i.e., a pair of orbitals.
      /// This function checks all levels (pairs of orbitals) and decides which are doubly occupied.
      /// It returns a new slater with each bit set indicating a doubly occupied orbital.
      ///
      /// Example: if <tt>phi = |00110110></tt>, then the third level (rightmost two pairs of bits) 
      /// is doubly occupied, so that <tt>doubly = |0100></tt>. (Rightmost bit is first bit.)
      ///
      /// \param phi     Slater determinant to work on
      /// \param doubly  "Condensed" Slater determinants (no spin) of doubly occupied levels. Really only used as a bit pattern.
      void findDoublyOccupied(const Slater& phi, Slater& doubly)
      {
	// Get the largest index, and round down to nearest even number.
	orbital_t alphamax = ((phi.last_index()>>1)<<1) + 2;
	// Loop though all indices, and check for doubly occupied positions.
	for (orbital_t alpha = 0; alpha<alphamax; alpha+=2)
	{
	  bool x1 = phi.get(alpha) ? true : false;
	  bool x2 = phi.get(alpha+1) ? true : false;
	  if ( x1 && x2 )
	  {
	    doubly.set(alpha >> 1);
	  }
	}

      }

      /// \brief Compute the Sz, <tt>nup</tt> and <tt>ndown</tt> values for <tt>sd_basis[0]</tt>,
      /// and also check that all Slater determinants have these values.
      void computeSz();

    public:

      /// \brief Assign a SD basis set to work on.
      ///
      /// It is assumed, that it is a complete basis for
      /// a given value for Sz, so that all S^2 eigenfunctions can be created.
      ///
      /// \param basis    a reference to a std::vector<Slater>. This is copied.
      void setSDBasis(std::vector<Slater>& basis) { 
	sd_basis = basis; 
	computeSz();
      }

      /// \brief Get the SD basis. 
      ///
      /// This is useful, since in general <tt>sd_basis</tt> is modified by permutation.
      ///
      /// \return Reference to <tt>sd_basis</tt>.
      std::vector<Slater>& getSDBasis() { return sd_basis; }

      /// \brief Compute eqivalence classes for all SDs, and rearrange
      /// basis according to equivalence class.
      ///
      /// The equivalence class is defined as all the Slater determinants with with the same set
      /// of singly and doubly occupied levels. This is useful since S^2 becomes block diagonal
      /// (acts only within each equivalence class) after rearrangement. See also getSDBasis().
      void rearrangeBasis();

      /// \brief Compute blocks of CSFs, assuming that the rearrangeBasis() is already
      /// called.
      ///
      /// The function loops through all equivalence classes, exctract configurations of 
      /// spins in the singly occupied levels, and takes proper linear combinations
      /// to produce S^2 eigenfunctions. (S^2 does not couple the doubly occupied levels.)
      ///
      /// The block structure is stored in block_structure, such that
      /// <tt>[ block_structure[i], block_structure[i+1] )</tt> is the range
      /// of SDs belonging to the block with class number
      /// <tt>eq_class[k].first</tt>.
      ///
      /// <tt>cg_coeffs</tt>, which is a std::vector<dense_matrix> ends up with
      /// holding the Clebsh-Gordan coefficients within each equivalence class. (The coefficients
      /// of the S^2 eigenfunctions), and <tt>s_values</tt> is a dense_vector of
      /// the corresponding S^2 eigenvalues.
      void findBlocks();

      /// \brief Print some debugging information to cerr.
      void print();

      /// \brief Return the number of equivalence class blocks.
      int numBlocks() { return eq_classes.size();  }

      /// \brief Get the sub-CSF block with total spin s, in the given eqivalence
      /// class. Note; it may be empty! Callers responsibility to check.
      ///
      /// The function getAllBlocks() is probably more useful in general.
      ///
      /// \param CSF       Output: CG coefficients of the sub-block
      /// \param S         an integer being 2*s, where s is the desired eigenvalue of S^2.
      /// \param eq_index  Which equivalence class to scan
      void getBlock(csf_block& CSF, int S, size_t eq_index);
    
      /// \brief Create a std::vector<csf_block> containing all the blocks
      /// belonging to a given s value. No empty blocks are added.
      /// (In effect extracts the S^2-basis with eigenvalue s.)
      void getAllBlocks(std::vector<csf_block>& CSF, int S);

      /// Similar to getAllBlocks(), but creates only trivial blocks from scratch,
      /// with no assumptions on sd_basis. It simply creates one block per SD,
      /// so that each SD \em is a CSF. Useful for e.g. spinless applications. See also the global
      /// function computeTrivialCSFBlocks(), which is perhaps more appropriate to use.
      void getAllBlocksTrivial(std::vector<csf_block>& CSF);
    
  };



} // end of namespace manybody


#endif // _CSF_HPP_
