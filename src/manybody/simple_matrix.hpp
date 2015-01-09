#ifndef _SIMPLE_MATRIX_HPP_
#define _SIMPLE_MATRIX_HPP_

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
#include <cassert>
#include <ostream>

#include "simple_vector.hpp"

namespace simple_dense {


  /// \brief Class definition of simple_matrix<T>: a simple implementation
  /// of column major dense matrices. Supports arbitrary leading dimension size.
  template<class T>
  class simple_matrix
  {
    private:
      std::vector<T> data_;   ///< Data.
      index_t nrows_;             ///< Number of rows
      index_t ncols_;             ///< Number of columns
      index_t ld_;                ///< Leading dimension size

      range active_rows_;         ///< Range of active rows  
      range active_cols_;         ///< Range of active cols
      bool use_block_;         ///< True if a block has been selected.


      typedef typename std::vector<T>::const_iterator const_data_iterator;    
      typedef typename std::vector<T>::iterator data_iterator;    
      

    public:
      /// Default constructor. Simple initialization.
      simple_matrix() : data_(0), nrows_(0), ncols_(0), ld_(0),
			active_rows_(), active_cols_(), use_block_(false) { }

      /// Constructor that allocates memory for matrix.
      /// \param nrows   Number of rows in matrix
      /// \param ncols   Number of columns in matrix
      /// \param ld Leading dimension size. Default is equivalent to
      /// using ld = nrows. Must be >= nrows.
      simple_matrix(index_t nrows, index_t ncols, index_t ld = 0)
	: active_rows_(), active_cols_(), use_block_(false)
      {
	this->resize(nrows, ncols, ld);
      }

      /// Element access.
      /// \param row  The row
      /// \param col  The column
      T& operator()(index_t row, index_t col)
      {
	return *getElmIterator(row,col);
      }

      /// Element access, const version
      /// \param row  The row
      /// \param col  The column
      T operator()(index_t row, index_t col) const
      {
	return *getElmIterator(row,col);
      }

      /// \brief Get pointer (iterator) to element, const version
      /// \param row row number
      /// \param col column number
      const_data_iterator getElmIterator(index_t row, index_t col) const
      {
	assert( (row>=1) && (row <= nrows_) );
	assert( (col>=1) && (col <= ncols_) );
	const index_t index = ld_*(col-1) + (row-1);
	assert( static_cast<typename std::vector<T>::size_type>(index) < data_.size() );
	return data_.begin() + index;
	
      }

      /// \brief Get pointer (iterator) to element.
      /// \param row row number
      /// \param col column number
      data_iterator getElmIterator(index_t row, index_t col)
      {
	assert( (row>=1) && (row <= nrows_) );
	assert( (col>=1) && (col <= ncols_) );
	const index_t index = ld_*(col-1) + (row-1);
	assert( static_cast<typename std::vector<T>::size_type>(index) < data_.size() );
	return data_.begin() + index;
	
      }


      /// Get number of rows
      index_t numRows() const { return nrows_; }

      /// Get number of columns
      index_t numCols() const { return ncols_; }

      /// Get leading dimension size
      index_t getLd() const { return ld_; }

      /// Resize matrix
      /// \param nrows   Number of rows in matrix
      /// \param ncols   Number of columns in matrix
      /// \param ld      Leading dimension. Default is = nrows. Must be >= nrows.
      void resize(index_t nrows, index_t ncols, index_t ld = 0)
      {
	//assert(nrows > 0);
	//assert(ncols > 0);

	nrows_ = nrows;
	ncols_ = ncols;

	if (ld < nrows)
	  ld_ = nrows_;
	else
	  ld_ = ld;

	if ( data_.size() != static_cast<typename std::vector<T>::size_type>(ld_ * ncols_) )
	  data_.resize(ld_ * ncols_);
      }


      /// \brief Fill matrix with a value.
      /// \param value   Value to fill with
      void fill(const T value)
      {
	std::fill(data_.begin(), data_.end(), value);
      }

      /// \brief Assignment operator that fills with value
      /// \param value   Value to fill with
      simple_matrix<T>& operator=(const T value)
      {
	this->fill(value);
	return *this;
      }

      /// \brief Multiply by constant
      /// \note Block version not implemented yet.
      /// \param factor factor to multiply with
      simple_matrix<T>& operator*=(const T factor)
      {
	for (index_t col = 1; col <= ncols_; ++col)
	{
	  data_iterator first = getElmIterator(1, col);
	  for (index_t k = 0; k < nrows_; ++k)
	    *(first + k) *= factor;
	} 
	return *this;

      }
      



      // ******* Functions for obtaining blocks *********

      /// \brief Obtain a sub-block by setting
      /// active ranges
      simple_matrix<T>& operator()(range rows, range cols)
      {

	// Expand ranges if all-ranges.
	rows.expand_all(1, nrows_);
	cols.expand_all(1, ncols_);
	
	// Check for legality of ranges.
	assert(rows.first_ >= 0);
	assert(cols.first_ >= 0);
	assert(rows.last_ <= nrows_);
	assert(cols.last_ <= ncols_);

	// Set active ranges.
	active_rows_ = rows;
	active_cols_ = cols;

	// Indicate that block has ben selected.
	use_block_ = true;

	// Exit.
	return *this;

      }


      /// \brief Copy constructor
      simple_matrix<T>& operator=(const simple_matrix<T>& src)
      {

	this->nrows_ = src.nrows_;
	this->ncols_ = src.ncols_;
	this->data_ = src.data_;
	this->ld_ = src.ld_;
	this->use_block_ = src.use_block_;
	this->active_cols_ = src.active_cols_;
	this->active_rows_ = src.active_rows_;
	return *this;

      }

      /// \brief Copy from simple_matrix<T>
      /// \param src  Matrix to copy
      simple_matrix<T>& operator=(simple_matrix<T>& src)
      {

	// Copy entire matrix. Note that data_ is changed.
	if ((!this->use_block_) && (!src.use_block_))
	{
	  this->nrows_ = src.nrows_;
	  this->ncols_ = src.ncols_;
	  this->data_ = src.data_;
	  this->ld_ = src.ld_;
	}
	// Copy a block. Note that data_ is changed.
	else if ((!this->use_block_) && (src.use_block_))
	{
	  
	  // Resize matrix to accomodate data.
	  this->resize(src.active_rows_.length(), src.active_cols_.length());

	  // Fill columns with data.
	  for (index_t k = 1; k <= ncols_; ++k)
	  {
	    index_t src_col = src.active_cols_.first_ + k - 1;
	    this->fillColumn(k, src.getElmIterator(src.active_rows_.first_, src_col));
	  }

	  // Erase active rows/cols in src.
	  src.use_block_ = false;

	}
	// Fill a block with entire matrix src. Note that data_ is *not* changed
	else if ((this->use_block_) && (!src.use_block_))
	{
	  
	  
	  // Make sure sizes are compatible.
	  assert(active_cols_.length() == src.numCols());
	  assert(active_rows_.length() == src.numRows());


	  // Copy data.
	  for (index_t k = active_cols_.first_; k <= active_cols_.last_; ++k)
	  {
	    index_t src_col = 1 + k - active_cols_.first_;
	    index_t n = active_rows_.last_ - active_rows_.first_ + 1;
	    copy2(getElmIterator(active_rows_.first_, k), 
	    src.getElmIterator(1, src_col), n);
	  }


	  // Reset block.
	  this->use_block_ = false;

	}
	// Fill a block with block from src.
	else
	{

	  // Make sure sizes are compatible.
	  assert(active_cols_.length() == src.active_cols_.length());
	  assert(active_rows_.length() == src.active_rows_.length());


	  // Copy data.
	  for (index_t k = active_cols_.first_; k <= active_cols_.last_; ++k)
	  {
	    index_t src_col = src.active_cols_.first_ + k - active_cols_.first_;
	    index_t n = active_rows_.last_ - active_rows_.first_ + 1;
	    copy2(getElmIterator(active_rows_.first_, k), 
	    src.getElmIterator(src.active_rows_.first_, src_col), n);
	  }

	  // Reset blocks.
	  this->use_block_ = false;
	  src.use_block_ = false;

	}
	
	
	return *this;

	

      }


      /// \brief Fill column with data
      /// \brief k  Column number
      /// \brief src   Data source
      void fillColumn(index_t k,  const_data_iterator src)
      {
	fillColumnSection(k, src, nrows_);
      }

      /// \brief Fill column section with data
      /// \brief k  Column number
      /// \brief src   Data source
      void fillColumnSection(index_t k,  const_data_iterator src, index_t n)
      {
	//T* dest = &(*this)(k, 1);
	data_iterator dest = getElmIterator(1,k);
	std::copy(src, src + n, dest);
      }
      
      /// \brief Fill section of data data
      /// \brief k  Column number
      /// \brief src   Data source
      void fillDataSection(data_iterator dest,  const_data_iterator src, index_t n)
      {
	std::copy(src, src + n, dest);
      }

      

  };



  

// ***************** GLOBAL FUNCTIONS *********************

  /// Ostream support for simple_matrix<T>
  template<class T>
  std::ostream& operator << (std::ostream& os, const simple_matrix<T>& A)
  {
    os << "size = [" << A.numRows() << ", " << A.numCols() << "], leading = " << A.getLd() << std::endl;
    index_t j, k;
    for (j = 1; j <= A.numRows(); ++j)
    {
      for (k = 1; k <= A.numCols(); ++k)
	os << A(j,k) << "  ";
      os << std::endl;
    }
    return os;
  }



} // namespace simple_dense

#endif // _SIMPLE_DENSE_HPP_
