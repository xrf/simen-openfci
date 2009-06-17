#ifndef _SIMPLE_LINALG_HPP_
#define _SIMPLE_LINALG_HPP_

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


#include "simple_matrix.hpp"
#include "simple_vector.hpp"
#include "lpp/lapack.hh"

namespace simple_dense {

  /// \brief Compute eigenvalues and/or eigenvectors of a matrix A, known to be symmetric.
  /// It uses only upper triangular part for computations.
  ///
  /// \param A      Matrix to diagonalize; output eigenvectors if compute_evec = true
  /// \param lambda Output vector of eigenvalues
  /// \param compute_evec  True if eigenvectors should be computed. They are stored in A if true.
  /// \return 0 if success, -i of this wrapper has a bug, +i if algorithm failed to converge.
  inline
  long int symeig(simple_matrix<double>& A, simple_vector<double>& lambda, bool compute_evec)
  {

    std::string jobz = compute_evec ? "v" : "n";
    long int n = A.numRows();
    long int lda = A.getLd();
    long int info;

    assert(n == A.numCols());

    simple_matrix<double> temp(n,n);
    // Backup A if !compute_evec, since syev destroys it.
    if (!compute_evec)
    {
      temp = A;
    }


    lpp::syev(jobz.c_str(), "u", &n, &(A(1,1)), 
    &lda, &(lambda(1)), &info);

    // restore A of !compute_evec.
    if (!compute_evec)
    {
      A = temp;
    }

    return info;
      

  }




  /// \brief Solve a system of tridiagonal linear equations Ax = b.
  ///
  /// Calls the LAPACK function lpp:gtsv().
  ///
  /// \param A Coefficient matrix. Unaffected.
  /// \param b In: Right-hand-side. Out: Solution.
  /// \return result of gtsv. == 0 if ok.
  inline
  long int solveTridiagonal(const simple_matrix<double>& A, simple_vector<double>& b)
  {

    long int n = A.numRows();
    long int nrhs = 1;
    long int ldb = b.length();
    long int info;



    //simple_vector<double> dl(n-1), d(n), du(n-1); // diagonals
    simple_vector<double> dl(n), d(n), du(n); // diagonals


    for (index_t j = 1; j <= n; ++j)
      d(j) = A(j,j);

    for (index_t j = 1; j < n; ++j)
    {
      dl(j) = A(j+1, j);
      du(j) = A(j, j+1);
    }

    lpp::gtsv(&n, &nrhs, &dl(1), &d(1), &du(1), &b(1), &ldb, &info);

    return info;

  }




  /// \brief Dump a dense_matrix to a stream. The output is a Matlab script.
  ///
  /// \param A        matrix to dump
  /// \param mstream  matlab script stream
  /// \param name     name of matrix
  /// \param debug    true if some fprintf statements are to be included in the script.
  template<class T>
  inline
  void dense_matrix_dump(simple_matrix<T>& A, std::ostream& mstream, std::string name, const bool debug = true)
  {
    index_t n = A.numRows();
    index_t m = A.numCols();

    //mstream.precision(16);
    //mstream << scientific;

    if (debug)
      mstream << "fprintf('Loading dense matrix " << name << " ...\\n');" << std::endl;

    mstream << name << " = zeros(" << n << ", " << m << ");" << std::endl;
    for (index_t col = 1; col <= m; ++col)
      for (index_t row = 1; row <= n; ++row)
	mstream << name << "(" << row << ", " << col << ") = " << A(row, col) << ";" << std::endl;
  
  }

  /// \brief Dump a dense_vector to a stream. The output is a Matlab script.
  ///
  /// \param A        vector to dump
  /// \param mstream  matlab script stream
  /// \param name     name of matrix
  /// \param debug    true if some fprintf statements are to be included in the script.
  template<class T>
  inline
  void dense_vector_dump(simple_vector<T>& A, std::ostream& mstream, std::string name, const bool debug = true)
  {
    int n = A.length();

    //mstream.precision(16);
    //mstream << scientific;
    if (debug)
      mstream << "fprintf('Loading dense vector " << name << " ...\\n');" << std::endl;

    mstream << name << " = zeros(" << n << ", " << 1 << ");" << std::endl;
    for (index_t row = 1; row <= n; ++row)
      mstream << name << "(" << row << ") = " << A(row) << ";" << std::endl;
  
  }


  /// \brief Analog to MATLAB's linspace function.
  ///
  /// Computes a vector v of length n, linearly interpolating between a and b.
  /// Thus, v(1) == a, v(n) == b, and v(j+1)-v(j) is constant.
  ///
  /// \param v      reference to destination vector
  /// \param a      becomes v(1)
  /// \param b      becomes v(n)
  /// \param n      number of points
  template<class T>
  inline
  void linspace(simple_vector<T>& v, T a, T b, index_t n)  
  {
    assert(n>0);
    
    v.resize(n);
    
    if (n==1)
    {
      v(1) = a;
      return;
    }
    
    for (index_t k = 1; k<=n; ++k)
    {
      v(k) = a + (b-a)*(T)(k-1)/(n-1);
    }
    
  }

  /// \brief Compute the n by n identity matrix
  ///
  /// \param I     Reference to destination matrix
  /// \param n     dimension of matrix
  template<class T>
  inline
  void identity_matrix(simple_matrix<T>& I, index_t n)  
  {
    I.resize(n,n);
    I = T(0);
    for (index_t k=1; k<=n; ++k)
      I(k,k) = T(1);
  }

  /// \brief Add matrices: C = a*A + b*B.
  /// \param A First matrix
  /// \param a First coeff
  /// \param B Second matrix
  /// \param b Second coeff
  /// \param C Destination
  template<class T>
  inline void addMatrices(simple_matrix<T>& C, const simple_matrix<T>& A, const T& a,
			  const simple_matrix<T>& B, const T& b)
  {

   
    C = A;
    C *= a;
    for (index_t col = 1; col <= A.numCols(); ++col)
      for (index_t row = 1; row <= A.numRows(); ++row)
	C(row,col) += b*B(row, col);

  }

  /// \brief Multiply matrices: C = A*B.
  /// \param A First matrix
  /// \param B Second matrix
  /// \param C Destination
  /// \param transB  True if B is to be transposed.
  template<class T>
  inline void multiplyMatrices(simple_matrix<T>& C, const simple_matrix<T>& A,
			       const simple_matrix<T>& B, bool transB = false)
  {

    if (!transB)
    {
      assert(A.numCols() == B.numRows());

      C.resize(A.numRows(), B.numCols());
      C = 0;

      for (index_t col = 1; col <= C.numCols(); ++col)
	for (index_t row = 1; row <= C.numRows(); ++row)
	  for (index_t k = 1; k <= A.numCols(); ++k)
	    C(row,col) += A(row,k)*B(k,col);
    }
    else
    {
      assert(A.numCols() == B.numCols());

      C.resize(A.numRows(), B.numRows());
      C = 0;

      for (index_t col = 1; col <= C.numCols(); ++col)
	for (index_t row = 1; row <= C.numRows(); ++row)
	  for (index_t k = 1; k <= A.numCols(); ++k)
	    C(row,col) += A(row,k)*B(col,k);
    }

  }


} // namespace simple_dense



#endif // _SIMPLE_LINALG_HPP_
