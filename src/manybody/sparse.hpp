#ifndef _SIMPLE_SPARSE_HPP_
#define _SIMPLE_SPARSE_HPP_

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

#include <map>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>


/**
 * \file sparse.hpp
 * \author Simen Kvaal
 * \date Early 2008.
 *
 * \brief This is a simple implementation of two
 * very simple sparse matrix class templates.
 * There is little bounds-checking and so on.
 *
 * SparseMatrix<I, T>
 *    - I should be an integer type; typically size_t.
 *    - T is the atomic numerical type, typically double.
 *
 * This class is derived from std::map<I, std::map<I, T> >.
 * In other words: The sparsity is implemented through use
 * of associative arrays. This is actually not so silly: It
 * is very easy to implement matrix manipulation algorithms
 * in this picture. Moreover, it is very simple to convert
 * this matrix to a more efficient sparse matrix, such as
 * SparseMatrixCrs<T> also implemented here.
 *
 * SparseMatrixCrs<T>
 *    - T should be the atomic numerical unit; typically
 *      double.
 *
 * This class is a very simple implementation of compressed
 * row storage sparse matrices.
 *
 *
 * Typically, one would use SparseMatrix<I, T> to build matrices
 * elementwise, compute matrix-matrix products and so on,
 * and one would then convert to SparseMatrixCrs<T> before
 * performing a lot of matrix-vector products, which are
 * a factor 10 faster with SparseMatrixCrs<T>; actually
 * in my tests equal in performance to Matlab sparse matrices.
 *
 */



namespace simple_sparse
{

  ///
  /// \brief A simple sparse vector class.
  ///
  template<class I, class T>
  class SparseVector : public std::map<I, T>
  {
    public:
      /// compute the length of the vector, i.e., the largest key + 1.
      I length() const { 
	typename SparseVector<I, T>::const_iterator it = this->end();
	--it;
	return it->first + 1; 
      }
  };



  
  /// \brief Define a very simple sparse matrix class,
  /// as a "sparse vector of sparse vectors"
  template<class I, class T>
  class SparseMatrix : public std::map<I, SparseVector<I, T> >
  { 
    protected:
      /// Number of rows
      I nrows;
      /// Number of columns
      I ncols;
      /// Set to true if nrows and ncols are properly set.
      bool structure_computed;
    public:
      
      /// Default constructor
      SparseMatrix<I, T>() : std::map<I, SparseVector<I, T> >(), 
			      nrows(0), ncols(0), structure_computed(false)
			      { }
      
      /// Compute number of rows and cols.
      void computeStructure() 
      {
	// rows
	typename SparseMatrix<I, T>::const_iterator it = this->end();
	--it;
	nrows = it->first + 1;
	
	// columns
	ncols = 0;
	for (it = this->begin(); it!=this->end(); ++it)
	  if ( ncols < it->second.length() )
	    ncols = it->second.length();
	
	
	structure_computed = true;
	
      }
      /// Set number of rows.
      void setNrows(I n) { 
	nrows = n;
	if ((nrows>0) && (ncols>0))
	  structure_computed = true;
      }
      /// Set number of columns.
      void setNcols(I n) { 
	ncols = n;
	if ((nrows>0) && (ncols>0))
	  structure_computed = true;
      }
      
      /// Return number of rows
      I getNrows() const { 
	assert(structure_computed);
	//if (!structure_computed)
	//  computeStructure();
	
	return nrows;
      }
      /// Return number of columns.
      I getNcols() const {
	assert(structure_computed);
	//if (!structure_computed)
	//  computeStructure();
	
	return ncols;
      }
      
      /// Compute number of nonzero entries.
      I getNnz() {
	I count = 0;
	
	typename SparseMatrix<I, T>::const_iterator it;
	for (it=this->begin(); it!=this->end(); ++it)
	  count += it->second.size();
	
	return count;
	
      }

      /// Create matrix from a simple_sparse::SparseMatrixCrs<T>
      /// matrix.
      template<class Matrix>
      void convertFrom(Matrix& src)
      {
	
	// set the number of columns/rows.
	setNrows( src.getNrows() );
	setNcols( src.getNcols() );
	// erase everything stored, if anything.
	this->clear();
	
	I row = 0, n = src.getNrows();
	
	// traverse data and copy.
	for (row=0; row < n; ++row)
	  for (I k=src.row_ptr[row]; k<src.row_ptr[row+1]; ++k)
	    (*this)[row][src.col[k]] = src.data[k];

	
	

      }

      /// Computes the transpose of
      /// the current matrix and stores
      /// in dst.
      void saveTranspose(SparseMatrix<I, T>& dst)
      {
	typename SparseMatrix<I, T>::iterator it;

	// for each row...
	for (it = this->begin(); it != this->end(); ++it)
	{
	  typename SparseMatrix<I, T>::value_type::second_type::iterator it2;

	  // for each element in the row...
	  for (it2 = it->second.begin(); it2 != it->second.end(); ++it2)
	    dst[it2->first][it->first] = it2->second;
	}

	dst.setNrows(getNcols());
	dst.setNcols(getNrows());

      }
      

  };
  


  ///
  /// \brief A compressed row storage sparse matrix
  ///
  template<class T>
  class SparseMatrixCrs
  {
    protected:
      /// Number of rows
      std::size_t nrows;
      /// Number of columns
      std::size_t ncols;
    public:
      /// The value type. Typically double.
      typedef T value_type;
      /// Holds pointers (into col and data) 
      /// to the first element in each row.
      std::vector<std::size_t> row_ptr; 
      /// Column number array
      std::vector<std::size_t> col;
      /// Element array
      std::vector<T> data;


      /// Default constructor.
      SparseMatrixCrs<T>() : nrows(0), ncols(0) { }
      
      /// Set number of rows.
      void setNrows(std::size_t n) { nrows = n; }
      /// Set number of columns.
      void setNcols(std::size_t n) { ncols = n; }
      /// Get number of rows.
      std::size_t getNrows() { return nrows; }
      /// Get number of columns.
      std::size_t getNcols() { return ncols; }

      
      /// Build the sparse matrix from a
      /// SparseMatrix object.
      template<class Matrix>
      void convertFrom(Matrix& src)
      {
	
	// Compute the numer of nonzeroes.
	std::size_t nnz = src.getNnz();
	
	// Resize data/index arrays accordingly.
	row_ptr.resize( src.getNrows() + 1 );
	col.resize(nnz);
	data.resize(nnz);
	// Update nrows and ncols.
	nrows = src.getNrows();
	ncols = src.getNcols();

	row_ptr[0] = 0; // first element. first row.

	// Loop through all rows of src.
	for (std::size_t row = 0; row<src.getNrows(); ++row)
	{
	  
	  // number of elements in row k is row_ptr[k+1]-row_ptr[k].
	  row_ptr[row+1] = row_ptr[row] + src[row].size();
	  
	  // fill the data/column arrays with info from
	  // current row.
	  std::size_t element_index = row_ptr[row];
	  typename Matrix::value_type::second_type::iterator it;
	  for (it = src[row].begin(); it != src[row].end(); ++it)
	  {
	    col[element_index] = it->first;
	    data[element_index] = it->second;
	    element_index++;
	  }
	}
	
      }

      /// Matrix-vector product, using C-style arrays.
      /// \param u C-array operand
      /// \param v C-array result
      void matrix_vector_product(T* u, T* v)
      {
	
	//std::cerr << "optimized ... " << std::endl;
	std::cerr << "*";
	
	// make sure that the vector is ok.
	//assert(M.getNcols() == u.size());
	
	typename std::size_t row = 0, n = this->getNrows();
	
	// resize target vector
	//v.resize(n);
	
	// compute product.
	for (row=0; row < n; ++row)
	{
	  v[row] = 0;
	  for (std::size_t k=this->row_ptr[row]; k<this->row_ptr[row+1]; ++k)
	  {
	    v[row] += u[this->col[k]]*this->data[k];
	  }
	}
    
      }
   
      /// Used memory; approximately, in kb.
      size_t mem()
      {
	return ((row_ptr.size()+col.size())>>10)*sizeof(size_t) + (data.size()>>10)*sizeof(T);
      }

      /// Read a specific matrix element. Slow, but functional.
      double element(size_t row, size_t the_col)
      {
	// Loop through all matrix elements k in the row and seek a match for col[k].
	for (size_t k = row_ptr[row]; k < row_ptr[row+1]; ++k)
	  if (col[k] == the_col)
	    return data[k];
	      
	return 0.0;

      }
   
      
  };
  




  ///
  /// Compute matrix-matrix product: M3 = M1*M2.
  ///
  template<class Matrix>
  void matrix_matrix_product(Matrix& M1, Matrix& M2, Matrix& M3 )
  {

    //if (!left_transpose)
    {

      // bail out if matrices not compatible for product
      assert(M1.getNcols() == M2.getNrows());
      
      // delete all data in M3.
      M3.clear();
      M3.setNrows(M1.getNrows());
      M3.setNcols(M2.getNcols());
      
      
      typename Matrix::const_iterator it1;
      typename Matrix::value_type::second_type::const_iterator it2;
      
      // loop through the rows of M1.
      // it1->first = row number
      // it1->second = row vector
      for (it1=M1.begin(); it1 != M1.end(); ++it1)
      {
	// loop through the elements of the row vector
	// it2->first = column number
	// it2->second = M1[row][col]
	for (it2=it1->second.begin(); it2 != it1->second.end(); ++it2)
	{
	  
	  // select the row of M2 corresponding to the current column in M1.
	  typename Matrix::value_type::second_type& row(M2[it2->first]);
	  
	  // if M2 has an element matching that of the current element in M1,
	  // then update M3.
	  typename Matrix::value_type::second_type::const_iterator it3;
	  for (it3 = row.begin(); it3 != row.end(); ++it3)
	  {
	    //if (it3->first == it2->first)
	    M3[it1->first][it3->first] += it2->second * it3->second;
	  }
	  
	}
	
      }
    } 

  }


  ///
  /// compute v = M*u, M is a SparseMatrixCrs,
  ///
  /// WARNING: no checking is done on the vectors.
  /// => only needed operation on vectors is [].
  template<class T, class Vector> 
  void matrix_vector_product(SparseMatrixCrs<T>& M, Vector& u, Vector& v)
  {
    
    //std::cerr << "optimized ... " << std::endl;

    // make sure that the vector is ok.
    //assert(M.getNcols() == u.size());
    
    typename std::size_t row = 0, n = M.getNrows();
    
    // resize target vector
    //v.resize(n);
    
    // compute product.
    for (row=0; row < n; ++row)
    {
      v[row] = T(0.0);
      for (std::size_t k=M.row_ptr[row]; k<M.row_ptr[row+1]; ++k)
      {
	v[row] += u[M.col[k]]*M.data[k];
      }
    }
    
  }


  ///
  /// compute v = M*u, M is a SparseMatrixCrs,
  ///
  /// WARNING: no checking is done on the vectors.
  /// => only needed operation on vectors is [].
  ///
  /// v is NOT initialized!
  template<class T, class Vector> 
  void matrix_vector_product_noinit(SparseMatrixCrs<T>& M, Vector& u, Vector& v)
  {
    
    
    typename std::size_t row = 0, n = M.getNrows();
    
    // compute product.
    for (row=0; row < n; ++row)
    {
      //v[row] = T(0.0);
      for (std::size_t k=M.row_ptr[row]; k<M.row_ptr[row+1]; ++k)
      {
	v[row] += u[M.col[k]]*M.data[k];
      }
    }
    
  }

  
  ///
  /// Compute v = M*u
  /// The vector must support [], resize(), size().
  /// The matrix is supposed to be SparseMatrix<I, T>
  ///
  template<class Matrix, class Vector> 
  void matrix_vector_product(Matrix& M, Vector& u, Vector& v)
  {

    std::cerr << "non-optimized ... " << std::endl;

    assert(M.getNcols() == u.size());
    
    v.resize(M.getNrows());
    
    typename Matrix::iterator row_it;
    
    for (row_it = M.begin(); row_it != M.end(); ++row_it)
    {
      v[row_it->first] = 0;
      
      typename Matrix::value_type::second_type::iterator col_it;
      
      for (col_it = row_it->second.begin(); col_it != row_it->second.end(); ++col_it)
	v[row_it->first] += col_it->second*u[col_it->first];
    }
    
  } 
  


  /// Dump *dense* matrix to Matlab-compliant dense matrix
  template<class Matrix>
  void matlab_dump2(Matrix& M, std::ostream& o, const char *name)
  {

    // statement that initializes the sparse matrix.
    o << name << " = zeros(" << M.size1() << ", " 
      << M.size2() << ");" << std::endl;

    o.precision(16);
    o << std::scientific;

    // loop through all the elements row by row and write them

    for (size_t row = 0; row < M.size1(); ++row)
      for (size_t col = 0; col < M.size2(); ++col)
	o << name << "(" << (row+1) << "," << (col+1) << ") = "
//	  << 1 << ";" << std::endl;
//	  << real(M(row,col)) << "+" << imag(M(row,col)) << "j;" << std::endl;
	  << M(row,col) << ";" << std::endl;
    
    
    

  }
  

  /// Dump matrix to Matlab-compliant sparse matrix
  /// Template specialized to SparseMatrixCrs<T>.
  template<class T>
  void matlab_dump(SparseMatrixCrs<T>& M, std::ostream& o, const char *name)
  {
    // statement that initializes the sparse matrix.
    o << name << " = sparse(" << M.getNrows() << ", " 
      << M.getNcols() << ");" << std::endl;

    o.precision(16);
    o << std::scientific;

    // loop through all the elements row by row and write them

    for (size_t row = 0; row < M.getNrows(); ++row)
    {
      for (size_t k=M.row_ptr[row]; k<M.row_ptr[row+1]; ++k)
      {
	o << name << "(" << (row+1) << "," << (M.col[k]+1) << ") = "
//	  << real(M.data[k]) << "+" << imag(M.data[k]) << "j;" << std::endl;
	  << M.data[k] << ";" << std::endl;
      }
    }
    

  }


  /// Dump matrix to Matlab-compliant sparse matrix
  /// Template specialized to SparseMatrixCrs<T>.
  template<class I,class T>
  void matlab_dump(SparseMatrix<I,T>& M, std::ostream& o, const char *name)
  {
    // statement that initializes the sparse matrix.
    o << name << " = sparse(" << M.getNrows() << ", " 
      << M.getNcols() << ");" << std::endl;

    o.precision(16);
    o << std::scientific;

    // loop through all the elements row by row and write them

    typename SparseMatrix<I,T>::const_iterator It;

    for (It = M.begin(); It != M.end(); ++It)
    {
      const SparseVector<I,T>& row = It->second;

      typename SparseVector<I,T>::const_iterator It2;

      for (It2 = row.begin(); It2 != row.end(); ++It2)
      {
	o << name << "(" << (It->first + 1) << "," << (It2->first+1) << ") = "
//	  << real(M.data[k]) << "+" << imag(M.data[k]) << "j;" << std::endl;
	  << It2->second << ";" << std::endl;
      }
    }
    

  }




  /// Save a vector to a binary file.
  /// The template parameter Vector must support size() and [].
  /// It must also define value_type as the type of ... the values!
  /// The function writes as many values to the file as possible.
  /// \sa load_vector
  template<class Vector>
  long save_vector(const Vector& v, std::string filename)
  {
    std::ofstream f(filename.c_str());
    long written = 0;
    
    if (!f.is_open())
    {
      return 0; // error! could not open file.
    }

    for (size_t k = 0; k<v.size(); ++k)
    {

      long bytes = sizeof(Vector::value_type);
      char *data = (char*)&v[k];

      f.write(data,bytes);

      if (f.bad())
      {
	return written; // error! could not write.
      }
      ++written;

    }

    
    f.close();


    return written; // success.

  }

  /// Load a vector from a binary file.
  /// The template parameter Vector must support resize() and [].
  /// It must also define value_type as the type of ... the values!
  /// The function reads as many values from the file as possible
  /// and stores them sequentially in the vector.
  /// \sa save_vector
  template<class Vector>
  long load_vector(Vector& v, std::string filename)
  {
    long read = 0;
    std::ifstream f(filename.c_str(), std::ios::binary|std::ios::ate);

    if (f.is_open())
    {
      long size = f.tellg();
      size_t elems = size / sizeof(Vector::value_type);

      v.resize(elems);
      
      f.seekg(0, std::ios::beg);

      for (size_t k=0; k<elems; ++k)
      {
	char *data = (char*)&v[k];
	size_t bytes = sizeof(Vector::value_type);
	f.read (data, bytes);
	++read;
      }

      f.close();
    
    }

    return read;

  }


  /// Save a sparse matrix to a text file.
  /// First line:
  ///    <Nrows> <Ncols>
  /// Each afterwards line contains: 
  ///    <row> <col> <element>
  /// \sa load_SparseMatrix
  /// Note: This function should be reprogrammed. It is very crude,
  /// and does not heed proper template parameter handling.
  template<class Matrix>
  int save_SparseMatrix(Matrix& A, const char* fname)
  {

    std::ofstream f(fname);
    if (!f.is_open())
    {
      return 0; // error! could not open file.
    }


    f << A.getNrows() << " " << A.getNcols() << std::endl;

    typename Matrix::const_iterator row_it;
    


    // Loop through all rows.
    for (row_it = A.begin(); row_it != A.end(); ++row_it)
    {
      // row_it->first = row number.
      // row_it->second = row data; a std::map<I, T>.
      
      typename Matrix::value_type::second_type::const_iterator col_it;
      
      // Loop through all elements in row.
      for (col_it = row_it->second.begin(); 
	   col_it != row_it->second.end(); ++col_it)
	f << row_it->first << " " << col_it->first << " " 
	  << col_it->second << std::endl;

    }
    
    // close file.
    f.close();


    return 1; // success
  }



  /// Load a sparse matrix from a text file.
  /// First line:
  ///    <Nrows> <Ncols>
  /// Each afterwards line contains: 
  ///    <row> <col> <element>
  /// \sa save_SparseMatrix
  /// Note: This function should be reprogrammed. It is very crude,
  /// and does not heed proper template parameter handling.
  template<class Matrix>
  int load_SparseMatrix(Matrix& A, const char* fname)
  {

    //std::cerr << "Reading SparseMatrix data from '" << fname << "'." << std::endl;

    std::ifstream f(fname);
    if (!f.is_open())
    {
      //std::cerr << "Some error occured." << std::endl;
      return 0; // error! could not open file.
    }

    size_t rows, cols;
    f >> rows >> cols;
    A.clear();
    A.setNcols(cols);
    A.setNrows(rows);


    //Matrix::key_type row;
    size_t row, col;
    //Matrix::value_type::second_type::key_type col;
    double value;
    //Matrix::value_type::second_type::value_type value;
    
    size_t old_row = 0;
    while (!f.eof())
    {
      
      old_row = row;
      f >> row >> col >> value;
      A[row][col] = value;
//       if ((row > 0) && (old_row < row))
// 	cerr << "Reading row " << row << endl;
		     
    }


    return 1; // success

    

  }

  /// Read a dense matrix from a text file.
  /// First line:
  ///    <Nrows> <Ncols>
  /// Each afterwards line contains: 
  ///    <row> <col> <element>
  /// \sa save_SparseMatrix load_SparseMatrix
  /// Note: This function should be reprogrammed. It is very crude,
  /// and does not heed proper template parameter handling.
  template<class Matrix>
  int load_dense_matrix(Matrix& A, const char* fname)
  {

    //std::cerr << "Reading SparseMatrix data from '" << fname << "'." << std::endl;

    std::ifstream f(fname);
    if (!f.is_open())
    {
      //std::cerr << "Some error occured." << std::endl;
      return 0; // error! could not open file.
    }

    size_t rows, cols;
    f >> rows >> cols;
    A.resize(rows, cols);
    A.clear();
    
    


    //Matrix::key_type row;
    size_t row, col;
    //Matrix::value_type::second_type::key_type col;
    double value;
    //Matrix::value_type::second_type::value_type value;
    
    size_t old_row = 0;
    while (!f.eof())
    {
      
      old_row = row;
      f >> row >> col >> value;
      A(row,col) = value;
//       if ((row > 0) && (old_row < row))
// 	cerr << "Reading row " << row << endl;
		     
    }


    return 1; // success

    

  }



  /// Convert a SparseMatrix<size_t, double> to a dense_matrix.
  /// A very simple function.
  inline
  int convert_to_dense_matrix(SparseMatrix<size_t, double>& A, std::vector<std::vector<double> >& B)
  {
    
    B.resize(A.getNrows());

    for (size_t row = 0; row < B.size(); ++row)
    {
      B[row].resize(A.getNcols());
      for (size_t col = 0; col < B[0].size(); ++col)
	B[row][col] = A[row][col];

    }

    return 0;
  }




  /// C = A*B,
  /// where A and B are symmetric matrices.
  /// This is then very efficient.
  template<class I, class T>
  void SparseMatrix_mult_symm(SparseMatrix<I,T>& A, 
			       SparseMatrix<I,T>& B,
			       SparseMatrix<I,T>& C)


  {
    // make sure all matrices are of equal dimension.
    I n = A.getNrows();
    assert(A.getNcols() == B.getNrows());
    assert(A.getNrows() == A.getNcols());
    assert(B.getNrows() == B.getNcols());

    C.clear();

    I row, col;

    for (row=0; row<n; ++row)
    {
      const SparseVector<I,T>& u = A[row];
      for (col=0; col<n; ++col)
      {
	SparseVector<I,T>& v = B[col];

	T sum(0);
	typename SparseVector<I,T>::const_iterator it;
	for (it = u.begin(); it != u.end(); ++it)
	  if (v.count(it->first))
	    sum += v[it->first]*it->second;

	if (sum != 0)
	  C[row][col] = sum;
	
      }
      
    }
    
  }
  
  
}



#endif // _SPARSE_HPP_
