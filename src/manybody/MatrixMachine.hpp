#ifndef _MATRIX_MACHINE_HPP_
#define _MATRIX_MACHINE_HPP_

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


#include "CSF.hpp"
#include "sparse.hpp"
#include "tools.hpp"
#include "linalg.hpp"

/** \file   MatrixMachine.hpp
 *  \author Simen Kvaal
 *  \date   9-12-08, update 10-7-2008
 *
 * \brief   Implementation of class template MatrixMachine.
 *
 * The MatrixMachine class is a templated class that handles the
 * building of matrices in of operators in a given basis of CSFs.
 *
 * At the moment it handles both 1, 2 and 3-body matrices.
 *
 * The template argument is a class that must implement the
 * computation of individual matrix elements, i.e., <a|X|b>, 
 * <a1 a2|Y|b1 b2>, and <a1 a2 a3|Z|b1 b2 b3>. At the time of writing,
 * these are \em not assumed to be \em anti-symmetrized. The function in
 * question must be compatible with the following prototype: 
 *
 *
 * \code
 * double element(size_t* a, size_t* b); 
 * \endcode
 *
 * The pointers assigned are C-style arrays with number of elements equal to rank. For example, the
 * function should compute <a[0] a[1] a[2] | Z | b[0] b[1] b[2]> in the 3-body case.
 *
 * For 2-body matrices (and 1-body matrices as a trivial case), the condition
 * <a1 a2|Y|b1 b2> = <a2 a1|Y|b2 b1> is common. This results in simpler computations,
 * and can be indicated with setSymmetryAssumption(bool y). If y = true, then the matrix
 * elements are assumed to be symmetric in this way.
 *
 * For 1-body matrices, the symmetry condition is vacuous, and for 3-body matriced the \em non-symmetric
 * case is \em not implemented, so setSymmetryAssumption() will have no effect.
 *
 *
 * Typical use:
 * 
 * \code
 * MatrixMachine m<elementHandler>; // elementHandler.element(size_t*,size_t*) returns a double.
 * m.setRank(2);
 * m.setSymmetryAssumption( false );
 * m.setSDBasis(sd_basis);
 * m.setCSFBlocks(csf_blocks);
 * m.buildMatrix(V); // V is a simple_sparse::SparseMatrixCrs<double> object
 * \endcode
 *
 */

namespace manybody {


  /// \brief The MatrixMachine class definition.
  template<class Caller>
  class MatrixMachine
  {
  
    private:
      /// Pointer to the Slater determinant basis
      std::vector<Slater>* sd_basis;
      /// Pointer to the CSF-block vector
      std::vector<csf_block>* CSF;
      /// The object that has the matrix element computation functions, assigned in the constructor.
      Caller& c_obj;
      /// The type of matrix: one-body, two-body, or three-body.
      int rank; 
      /// Number of particles. Assigned during buildMatrix(), where it is detected from sd_basis.
      int particles;
      /// If true, then progress messages during building will be output.
      bool prog_messages;
      /// True if symmetry of two-particle operators is assumed; gives shorter calculations in that case.
      /// (For one-body, the assumtion is vacuous, and for three-body only the symmetric case is implemented.)
      bool symmetry_assumption;
    public:
      /// \brief Default constructor: Erases pointers, etc.
      /// 
      /// \param c_obj1   Reference to object of type Caller that holds 
      /// the matrix element computation function.
      MatrixMachine(Caller& c_obj1) : c_obj(c_obj1)
      {
	symmetry_assumption = true;
	prog_messages = false;
	sd_basis = 0;
	CSF = 0;
	rank = 2;
      }

      /// \brief Assign symmetry assumption.
      ///
      /// \param s    true if matrix elements are assumed to be symmetric; 
      /// see main documentation of MatrixMachine.hpp.
      void setSymmetryAssumption(bool s) { symmetry_assumption = s; }

      /// \brief Assign rank.
      ///
      /// \param r    rank of matrix: 1, 2 or 3 are legal values
      void setRank(int r) { 
	rank = r; 
	assert((r>=1) && (r<=3)); 
      }

      /// \brief Assign a SD basis.
      ///
      /// \param basis    Reference to std::vector<Slater> instance. A pointer is extracted.
      void setSDBasis(std::vector<Slater>& basis)
      {
	sd_basis = &basis;
      }

      /// \brief Assign a vector of CSF blocks.
      ///
      /// \param CSFVec   Reference to std::vector<csf_block> instance. A pointer is extracted.
      void setCSFBlocks(std::vector<csf_block>& CSFVec)
      {
	CSF = &CSFVec;
      }

      /// \brief Build a matrix. 
      ///
      /// \param A            Reference to matrix where we store the result
      /// \param is_diagonal  Set to true if we know that the matrix is diagonal. Saves a lot of time. 
      void buildMatrix(simple_sparse::SparseMatrixCrs<double>& A, bool is_diagonal = false);

      /// \brief This function computes a single element of the matrix 
      /// A computed in buildMatrix(), rank = 1 case.
      ///
      /// \param sd_row Index of row Slater determinant.
      /// \param sd_col Index of column Slater determinant.
      /// \param diff   Output: Reference to integer which holds the number of orbitals different for this matrix elements.
      /// \return       Matrix element <sd_row|A|sd_col>
      double singleElement1(size_t sd_row, size_t sd_col, int& diff);
      /// \brief Two-body matrix element, see singleElement1().
      double singleElement2(size_t sd_row, size_t sd_col, int& diff);
      /// \brief Two-body matrix element in the <em>asymmetric case</em>, see singleElement1().
      double singleElementAsymmetric2(size_t sd_row, size_t sd_col, int& diff);
      /// \brief Three-body matrix element, see singleElement1().
      double singleElement3(size_t sd_row, size_t sd_col, int& diff);



      /// \brief Given a dense matrix \f$A\f$, perform
      /// \f$U_1AU_2^\dag\f$. Not interesting for user...
      ///
      ///
      /// Used in buildMatrix() to compute matrix
      /// elements in CSF basis. Typically, U1 and U2 will be
      /// Clebsh-Gordan coefficient matrices, see class CsfMachine.
      ///
      /// (The user will never use this ...)
      ///
      /// \param A Inner matrix; block of matrix in SD basis
      /// \param U1 CG coeffs
      /// \param U2 CG coeffs.
      void blockTransform(const dense_matrix& A, const dense_matrix& U1, const dense_matrix& U2, dense_matrix& A2)
      {
	const int n1 = U1.numCols();
	const int n2 = U2.numCols();
	const int m1 = U1.numRows();
	const int m2 = U2.numRows();
      
      
	for (int j = 1; j<= m1; ++j)
	  for (int k = 1; k <= m2; ++k)
	  {
	    A2(j,k) = 0.0;
	  
	    for (int p = 1; p <= n1; ++p)
	      for (int q = 1; q <= n2; ++q)
	      {
		A2(j,k) += U1(j,p)*A(p,q)*U2(k,q);
	      }
	  
	  }
      
      }
    
    
    
      /// \brief Turn on / off progress messages.
      ///
      /// \param yesno    Messages or not, Sir/Madam?
      void messages(bool yesno) { prog_messages = yesno; }

  };



  template<class Caller>
  inline
  double MatrixMachine<Caller>::singleElement1(size_t sd_row, size_t sd_col, int& diff)
  {
  
    // Fetch Slater determinants.

    Slater& phi_row = (*sd_basis)[sd_row];
    Slater& phi_col = (*sd_basis)[sd_col];

    // alpha and beta are used for pinpointing different orbitals,
    // p and q for their positions.
    orbital_t alpha[2], beta[2], p[2], q[2];
    diff = phi_row.analyze_difference(phi_col, alpha, beta, p, q, 1);

    // bail out if too many occupied orbitals are different
    if (diff > 1)
      return 0.0;

    // The number of orbitals.
    orbital_t n_orbitals = std::max(phi_col.last_index(), phi_row.last_index()) + 1;


    // elm holds the result.
    double elm = 0.0;

    switch (diff)
    {
      case 0: // diagonal matrix element
	for (alpha[0] = 0; alpha[0]<(orbital_t)n_orbitals; ++alpha[0])
	  if (phi_col.get(alpha[0]))
	    elm += c_obj.element(alpha,alpha);
      
	return elm;
      
	break;
      case 1: // singly excited matrix element
      
	elm += m1pow<int>(p[0]+q[0]) * c_obj.element(alpha, beta);
      
	return elm;
	break;
      default:
	return 0.0;
	break;
    }
  

  }




  template<class Caller>
  inline
  double MatrixMachine<Caller>::singleElement2(size_t sd_row, size_t sd_col, int& diff)
  {
  
    // Fetch Slater determinants.

    Slater& phi_row = (*sd_basis)[sd_row];
    Slater& phi_col = (*sd_basis)[sd_col];

    // alpha and beta are used for pinpointing different orbitals,
    // p and q for their positions.
    orbital_t alpha[2], beta[2], p[2], q[2];
    diff = phi_row.analyze_difference(phi_col, alpha, beta, p, q, 2);

    // bail out if too many occupied orbitals are different
    if (diff > 2)
      return 0.0;

    // The number of orbitals.
    orbital_t n_orbitals = std::max(phi_col.last_index(), phi_row.last_index()) + 1;

    // Used for holding common orbitals
    Slater z;

    // elm holds the result.
    double elm = 0.0;
  
    // Various helper things.  
    size_t temp;
    int sign;

    switch (diff)
    {
      case 0:

	for (alpha[0] = 0; alpha[0]<(orbital_t)n_orbitals; ++alpha[0])
	{
	  if (phi_col.get(alpha[0]))
	  {
	    for (alpha[1] = alpha[0]+1; alpha[1]<(orbital_t)n_orbitals; ++alpha[1])
	    {
	      if (phi_col.get(alpha[1]))
	      {
		beta[0] = alpha[1]; beta[1] = alpha[0];
		elm += c_obj.element(alpha, alpha) - c_obj.element(alpha, beta);
	      }
	    }
	  }
	}
	  
	return elm;
	break; // irrelevant...
      case 1:

	  
	// Create a Slater determinant with only the common bits
	// of c and c2 set.
	z = phi_row;
	z.unset(alpha[0]);

	for (int betaz=0; betaz<(int)n_orbitals; ++betaz)
	{
	  //if (z.bit(beta)) // (common) bit at position beta !
	  if (z.get(betaz)) // (common) bit at position beta !
	  {
	    //alpha[0] = b;
	    alpha[1] = betaz;
	    //beta[0] = a;
	    beta[1] = betaz;
	    int sign = p[0] + q[0];
	    elm += m1pow<int>(sign)*c_obj.element(alpha,beta);
	    size_t temp = beta[0];
	    beta[0] = beta[1]; beta[1] = temp;
	    elm -= m1pow<int>(sign)*c_obj.element(alpha,beta);
	    beta[0] = beta[1]; beta[1] = temp;

	  }
	}


	return elm;


	break;
      case 2:

	sign = m1pow<int>(p[0]+p[1]+q[0]+q[1]);
  
	elm = c_obj.element(alpha, beta);
	temp = beta[0];
	beta[0] = beta[1]; beta[1] = temp;
	elm -= c_obj.element(alpha, beta);
	elm *= sign;
	//alpha[0] = alpha[1]; alpha[1] = temp;
	  
	return elm;

	break;
      default:
	return 0.0;
	break;
    } 

  }




  template<class Caller>
  inline
  double MatrixMachine<Caller>::singleElementAsymmetric2(size_t sd_row, size_t sd_col, int& diff)
  {
  
    // Fetch Slater determinants.

    Slater& phi_row = (*sd_basis)[sd_row];
    Slater& phi_col = (*sd_basis)[sd_col];

    // alpha and beta are used for pinpointing different orbitals,
    // p and q for their positions.
    orbital_t alpha[2], beta[2], p[2], q[2];
    diff = phi_row.analyze_difference(phi_col, alpha, beta, p, q, 2);

    // bail out if too many occupied orbitals are different
    if (diff > 2)
      return 0.0;

    // The number of orbitals.
    orbital_t n_orbitals = std::max(phi_col.last_index(), phi_row.last_index()) + 1;

    // Used for holding common orbitals
    Slater z;

    // elm holds the result.
    double elm = 0.0;
  
    // Various helper things.  
    //size_t temp;
    int sign;

    // flipped indices.
    orbital_t alphaflip[2], betaflip[2];
    //alphaflip[0] = alpha[1]; alphaflip[1] = alpha[0];
    //betaflip[0] = beta[1]; betaflip[1] = beta[0];

    switch (diff)
    {
      case 0:

	for (alpha[0] = 0; alpha[0]<(orbital_t)n_orbitals; ++alpha[0])
	{
	  if (phi_col.get(alpha[0]))
	  {
	    for (alpha[1] = alpha[0]+1; alpha[1]<(orbital_t)n_orbitals; ++alpha[1])
	    {
	      if (phi_col.get(alpha[1]))
	      {
		alphaflip[0] = alpha[1]; alphaflip[1] = alpha[0];
	      
		elm += 0.5*c_obj.element(alpha, alpha) - 0.5*c_obj.element(alpha, alphaflip);
		elm += -0.5*c_obj.element(alphaflip, alpha) + 0.5*c_obj.element(alphaflip, alphaflip);
	      }
	    }
	  }
	}
	  
	return elm;
	break; // irrelevant...
      case 1:

	  
	// Create a Slater determinant with only the common bits
	// of c and c2 set.
	z = phi_row;
	z.unset(alpha[0]);

	for (int betaz=0; betaz<(int)n_orbitals; ++betaz)
	{
	  //if (z.bit(beta)) // (common) bit at position beta !
	  if (z.get(betaz)) // (common) bit at position beta !
	  {
	    //alpha[0] = b;
	    alpha[1] = betaz;
	    //beta[0] = a;
	    beta[1] = betaz;
	    int sign = p[0] + q[0];


	    alphaflip[0] = alpha[1]; alphaflip[1] = alpha[0];
	    betaflip[0] = beta[1]; betaflip[1] = beta[0];

	    elm += m1pow<int>(sign)*0.5*c_obj.element(alpha,beta);

	    elm -= m1pow<int>(sign)*0.5*c_obj.element(alpha,betaflip);

	    elm -= m1pow<int>(sign)*0.5*c_obj.element(alphaflip,beta);

	    elm += m1pow<int>(sign)*0.5*c_obj.element(alphaflip,betaflip);


	  }
	}


	return elm;


	break;
      case 2:

	sign = m1pow<int>(p[0]+p[1]+q[0]+q[1]);

	alphaflip[0] = alpha[1]; alphaflip[1] = alpha[0];
	betaflip[0] = beta[1]; betaflip[1] = beta[0];
  
	elm = 0.5*c_obj.element(alpha, beta);
	elm -= 0.5*c_obj.element(alpha, betaflip);
	elm -= 0.5*c_obj.element(alphaflip, beta);
	elm += 0.5*c_obj.element(alphaflip, betaflip);

	elm *= sign;
	  
	return elm;

	break;
      default:
	return 0.0;
	break;
    } 

  }






  template<class Caller>
  inline
  double MatrixMachine<Caller>::singleElement3(size_t sd_row, size_t sd_col, int& diff)
  {
  
    // Fetch Slater determinants.

    Slater& phi_row = (*sd_basis)[sd_row];
    Slater& phi_col = (*sd_basis)[sd_col];

    // alpha and beta are used for pinpointing different orbitals,
    // p and q for their positions.
    orbital_t alpha[3], beta[3], p[3], q[3];
    diff = phi_row.analyze_difference(phi_col, alpha, beta, p, q, 3);

    // bail out if too many occupied orbitals are different
    if (diff > 3)
      return 0.0;

    // The number of orbitals.
    orbital_t n_orbitals = std::max(phi_col.last_index(), phi_row.last_index()) + 1;

    // Used for holding common orbitals
    Slater z;

    // elm holds the result.
    double elm = 0.0;
  

    size_t betaperm[3];

    switch (diff)
    {
      case 0:

	// the matrix element is a sum over all triples of
	// occupied orbitals.
	for (alpha[2] = 2; alpha[2] < n_orbitals; ++alpha[2])
	  if (phi_col.get(alpha[2]))
	    for (alpha[1] = 1; alpha[1] < alpha[2]; ++alpha[1])
	      if (phi_col.get(alpha[1]))
		for (alpha[0] = 0; alpha[0] < alpha[1]; ++alpha[0])
		  if (phi_col.get(alpha[0]))
		  {
		    beta[0] = alpha[0];
		    beta[1] = alpha[1];
		    beta[2] = alpha[2];
		    elm += 6*c_obj.element(alpha, beta);

		    beta[1] = alpha[0];
		    beta[2] = alpha[1];
		    beta[0] = alpha[2];
		    elm += 6*c_obj.element(alpha, beta);

		    beta[2] = alpha[0];
		    beta[0] = alpha[1];
		    beta[1] = alpha[2];
		    elm += 6*c_obj.element(alpha, beta);

		    beta[1] = alpha[0];
		    beta[0] = alpha[1];
		    beta[2] = alpha[2];
		    elm -= 6*c_obj.element(alpha, beta);

		    beta[0] = alpha[0];
		    beta[2] = alpha[1];
		    beta[1] = alpha[2];
		    elm -= 6*c_obj.element(alpha, beta);

		    beta[2] = alpha[0];
		    beta[1] = alpha[1];
		    beta[0] = alpha[2];
		    elm -= 6*c_obj.element(alpha, beta);

		  }

	  
	return elm;
	break; 

      case 1:

	// compute the common orbitals in phi_row and phi_col

	z = phi_row;
	z.unset(alpha[0]);


	for (alpha[2] = 1; alpha[2] < n_orbitals; ++alpha[2])
	  if (z.get(alpha[2]))
	    for (alpha[1] = 0; alpha[1] < alpha[2]; ++alpha[1])
	      if (z.get(alpha[1]))
	      {

		betaperm[0] = beta[0];
		betaperm[1] = alpha[1];
		betaperm[2] = alpha[2];
		elm += 6*c_obj.element(alpha, betaperm);

		betaperm[1] = beta[0];
		betaperm[2] = alpha[1];
		betaperm[0] = alpha[2];
		elm += 6*c_obj.element(alpha, betaperm);

		betaperm[2] = beta[0];
		betaperm[0] = alpha[1];
		betaperm[1] = alpha[2];
		elm += 6*c_obj.element(alpha, betaperm);

		betaperm[1] = beta[0];
		betaperm[0] = alpha[1];
		betaperm[2] = alpha[2];
		elm -= 6*c_obj.element(alpha, betaperm);

		betaperm[0] = beta[0];
		betaperm[2] = alpha[1];
		betaperm[1] = alpha[2];
		elm -= 6*c_obj.element(alpha, betaperm);

		betaperm[2] = beta[0];
		betaperm[1] = alpha[1];
		betaperm[0] = alpha[2];
		elm -= 6*c_obj.element(alpha, betaperm);


	      }

	elm *= m1pow(p[0]+q[0]);

	  
	return elm;
	break;

      case 2:

	// compute the common orbitals in phi_row and phi_col

	z = phi_row;
	z.unset(alpha[0]);
	z.unset(alpha[1]);

	for (alpha[2] = 0; alpha[2] < n_orbitals; ++alpha[2])
	  if (z.get(alpha[2]))
	  {
	    betaperm[0] = beta[0];
	    betaperm[1] = beta[1];
	    betaperm[2] = alpha[2];
	    elm += 6*c_obj.element(alpha, betaperm);

	    betaperm[1] = beta[0];
	    betaperm[2] = beta[1];
	    betaperm[0] = alpha[2];
	    elm += 6*c_obj.element(alpha, betaperm);

	    betaperm[2] = beta[0];
	    betaperm[0] = beta[1];
	    betaperm[1] = alpha[2];
	    elm += 6*c_obj.element(alpha, betaperm);

	    betaperm[1] = beta[0];
	    betaperm[0] = beta[1];
	    betaperm[2] = alpha[2];
	    elm -= 6*c_obj.element(alpha, betaperm);

	    betaperm[0] = beta[0];
	    betaperm[2] = beta[1];
	    betaperm[1] = alpha[2];
	    elm -= 6*c_obj.element(alpha, betaperm);

	    betaperm[2] = beta[0];
	    betaperm[1] = beta[1];
	    betaperm[0] = alpha[2];
	    elm -= 6*c_obj.element(alpha, betaperm);


	  }
      

	elm *= m1pow(p[0]+p[1]+q[0]+q[1]); // permutation sign

	return elm;
	break;

      case 3:

	// exactly 3 orbitals are different in phi_row and phi_col.
	// alpha/p belongs to phi_row, beta/q to phi_col.

	// there are a total of 36 contributions to the matrix element, some
	// are equal by symmetry assumption on tree-body operator.

	size_t beta0[3];
	beta0[0] = beta[0];
	beta0[1] = beta[1];
	beta0[2] = beta[2];

	elm += 6*c_obj.element(alpha, beta);
      
	beta[1] = beta0[0];
	beta[0] = beta0[1];
	beta[2] = beta0[2];

	elm -= 6*c_obj.element(alpha, beta);

	beta[0] = beta0[0];
	beta[2] = beta0[1];
	beta[1] = beta0[2];

	elm -= 6*c_obj.element(alpha, beta);

	beta[2] = beta0[0];
	beta[1] = beta0[1];
	beta[0] = beta0[2];

	elm -= 6*c_obj.element(alpha, beta);

	beta[1] = beta0[0];
	beta[2] = beta0[1];
	beta[0] = beta0[2];

	elm += 6*c_obj.element(alpha, beta);

	beta[2] = beta0[0];
	beta[0] = beta0[1];
	beta[1] = beta0[2];

	elm += 6*c_obj.element(alpha, beta);

	elm *= m1pow(p[0]+q[0]+p[1]+q[1]+p[2]+q[2]);

	return elm;
	break;

      default:

	return 0.0;
	break;
    } 

  }






  template<class Caller>
  inline void MatrixMachine<Caller>::buildMatrix(simple_sparse::SparseMatrixCrs<double>& A, bool is_diagonal)
  {

    using namespace std;

    // block indices for rows/columns.
    size_t b1, b2;
    // number of blocks
    size_t nblocks = CSF->size();

    // row and column indices for total matrix.
    size_t row = 0, col = 0;
    // boundary of the block to be processed.
    size_t row0 = 0, col0 = 0;


    // compute number of particles in basis
    particles = (*sd_basis)[0].count();

    // compute matrix dimension
    size_t mat_dim = 0;
    for (size_t j = 0; j<CSF->size(); ++j)
      mat_dim += (*CSF)[j].n_csf;
    cerr << "Matrix will be " << mat_dim << " dimensional." << endl;
    cerr << "There are " << nblocks << " CSF blocks." << endl;

//   simple_sparse::sparse_matrix<size_t, double> Atemp;
//   Atemp.setNrows(mat_dim);
//   Atemp.setNcols(mat_dim);
    A.setNrows(mat_dim);
    A.setNcols(mat_dim);
    A.row_ptr.resize(mat_dim+1);
    fill(A.row_ptr.begin(), A.row_ptr.end(), 0.0);
    A.data.resize(0);
    A.col.resize(0);

    // Temporary space for a few rows.
    size_t temp_rows = 0;
    // compute maximum number of rows in a block
    for (b1 = 0; b1 < nblocks; ++b1)
      temp_rows = std::max(temp_rows, (*CSF)[b1].n_csf); 
    dense_matrix_int col_temp(temp_rows, mat_dim);
    std::vector<size_t> nelms_temp(temp_rows); // number of elements found in a row
    dense_matrix data_temp(temp_rows, mat_dim);
    size_t nelms_tot = 0;

    dense_matrix Ablock;

    // this vector holds the indices of the nonzero blocks for a given b1.
    std::vector<size_t> nonzero_blocks;
    // see of a sparsity pattern is assigned. If not, initialize nonzero_blocks
    // with a full block structure.
    nonzero_blocks.resize(nblocks);
    for (size_t t = 0; t<nblocks; ++t)
      nonzero_blocks[t] = t;

    // loop over row blocks.
    for (b1 = 0; b1 < nblocks; ++b1)
    {
      csf_block& block1 = (*CSF)[b1];


      if (block1.n_csf == 0)
	continue;

      if (prog_messages)
	cerr << move_to(0) << "( Row block #" << b1 << " of " << nblocks << ", nnz so far: " << A.data.size() << " ) ";

      // initialize temporary variables for matrix.
      fill(nelms_temp.begin(), nelms_temp.end(), 0);


      // start a new row, set block boundary = 0.
      col0 = 0;

      // loop over col blocks.
      for (b2 = 0; b2 < nblocks; ++b2)
      {
	bool process_block = true;
	if (is_diagonal && (b1 != b2))
	  process_block = false;

	csf_block& block2 = (*CSF)[b2];
	if (process_block)
	{

       
	  if (block2.n_csf == 0)
	    continue;

	  // block1 and block2 are now csf_blocks for rows/columns, resp.
	  // compute corresponding matrix block in Slater determinant basis

	  size_t sd1, sd2;
	  size_t dim1 = block1.sd_end - block1.sd_begin;
	  size_t dim2 = block2.sd_end - block2.sd_begin;
	  Ablock.resize(dim1, dim2);
	  int block_nnz = 0; // number of matrix elements not identical to zero for this block.
	  for (sd1 = 0; sd1<dim1; ++sd1)
	    for (sd2 = 0; sd2<dim2; ++sd2)
	    {
	      int diff = 0;

	      switch (rank)
	      {
		case 1:
		  Ablock(sd1+1, sd2+1) = 
		    singleElement1(sd1+block1.sd_begin, sd2+block2.sd_begin, diff);
		  break;
		case 2:
		  if (symmetry_assumption)
		    Ablock(sd1+1, sd2+1) = 
		      singleElement2(sd1+block1.sd_begin, sd2+block2.sd_begin, diff);
		  else
		    Ablock(sd1+1, sd2+1) = 
		      singleElementAsymmetric2(sd1+block1.sd_begin, sd2+block2.sd_begin, diff);
		  break;
		case 3:
		  Ablock(sd1+1, sd2+1) = 
		    singleElement3(sd1+block1.sd_begin, sd2+block2.sd_begin, diff);
		  break;
		default:
		  Ablock(sd1+1, sd2+1) = 0;
		  break;
	      }

	      if (diff <= rank)
		block_nnz++;
	    }

	  if (block_nnz) // add block to matrix if some nonzeroes were found.
	  {

	    dense_matrix Ablock3(block1.n_csf, block2.n_csf);
	    blockTransform(Ablock, block1.coeffs, block2.coeffs, Ablock3);
	  
	    // update the total matrix A.
	    //cerr << Ablock3 << endl;
	    for (row = 0; row < block1.n_csf; ++row)
	    {
	      // add new row to matrix.
	      for (col = 0; col < block2.n_csf; ++col)
	      {
		if (Ablock3(row+1,col+1)  != 0.0)
		{
		  col_temp(row+1, nelms_temp[row]+1) = col + col0;
		  data_temp(row+1, nelms_temp[row]+1) = Ablock3(row+1,col+1);
		  nelms_temp[row]++;
		  nelms_tot++;
		  //Atemp[row+row0][col+col0] = Ablock3(row+1,col+1);
		}
	      }
	    }
	  
	  } // if (block_nnz)

	} // if (process_block) ...      
	// increment block boundary in total matrix
	col0 += block2.n_csf;

      } // b2 loop

      // add temporary matrix data to full matrix A.
      //A.data.resize(nelms_tot);
      //A.col.resize(nelms_tot);
      for (row = 0; row < block1.n_csf; ++row)
      {
	A.row_ptr[row0+row+1] = A.row_ptr[row0+row] + nelms_temp[row];
	//size_t start = A.row_ptr[row0+row];
	for (size_t k = 1; k<=nelms_temp[row]; ++k)
	{
	  A.data.push_back(data_temp(row+1,k));
	  A.col.push_back(col_temp(row+1,k));
// 	A.data[start+k-1] = data_temp(row+1, k);
// 	A.col[start+k-1] = col_temp(row+1,k);
	}
      }


      // increment block boundary in total matrix
      row0 += block1.n_csf;

    } // b1 loop


    //A.convertFrom(Atemp);

    if (prog_messages)
      cerr << endl;
    cerr << "Density of matrix : " << A.data.size() / sqr((double)A.getNrows()) << endl;

  }

} // namespace manybody


#endif // _MATRIX_MACHINE_HPP_
