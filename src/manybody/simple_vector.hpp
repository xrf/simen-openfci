#ifndef _SIMPLE_VECTOR_HPP_
#define _SIMPLE_VECTOR_HPP_

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

#include <ostream>
#include <cassert>
#include <vector>

namespace simple_dense {



  /// \brief Alternative copy function
  template<class I, class CI>
  void copy2(I dest, CI src, size_t n)
  {
    std::copy(src, src + n, dest);
  }

  /// \brief Type to use for indices.
  typedef long int index_t ;

  /// \brief This struct defines a range of indices. 
  ///
  /// Use the constructor
  /// to construct a range, e.g., range(5,6). Use range() to construct
  /// the default range, which is everything.
  struct range
  {
    public:
      index_t first_;   ///< First index of range
      index_t last_;    ///< Last index of range
      bool all_;    ///< True if all is to be seleced, regardless of first_ and last_ settings
      /// Constructur, constructs a range.
      range(index_t first, index_t last)
      {
	first_ = first;
	last_ = last;
	all_ = false;

      }
      /// \brief Default constructor. Constructs default range, which is all.
      range() : first_(0), last_(0), all_(true) { }

      /// \brief Expand an all-range to the limits given
      void expand_all(index_t a, index_t b)
      {
	if (all_)
	{
	  first_ = a; 
	  last_ = b;
	}
      }

      /// \brief Return length of interval
      index_t length() const { return last_ - first_ + 1; }


  };


  /// \brief Class definition of simple_vector<T>: a simple implementation
  /// of a (column) vector.
  template<class T>
  class simple_vector
  {
    private:
      std::vector<T> data_;      ///< Data vector.
      index_t length_;           ///< Length of vector.

      typedef typename std::vector<T>::const_iterator const_data_iterator;    
      typedef typename std::vector<T>::iterator data_iterator;    

    public:
      /// \brief Default constructor
      simple_vector() : length_(0), data_(0) { }

      /// \brief Constructur that allocates memory (if OWNS = true)
      /// \brief length     Length of vector
      simple_vector(index_t length) : data_(0), length_(0)
      {
	this->resize(length);
      }

      /// \brief Resize vector.
      /// \param length   new length of vector
      void resize(index_t length)
      {
	//assert(length > 0);

	if (length != length_)
	{
	  length_ = length;
	  data_.resize(length_);
	}

      }

      /// \brief Get vector length.
      index_t length() const
      {
	return length_;
      }

      /// \brief Get pointer to element (const).
      /// \param j  Index of element
      const_data_iterator getElmIterator(index_t j) const
      {
	assert((j >= 1) && (j <= length_));
	return data_.begin() + j - 1;
      }

      /// \brief Get pointer to element.
      /// \param j  Index of element
      data_iterator getElmIterator(index_t j)
      {
	assert((j >= 1) && (j <= length_));
	return data_.begin() + j - 1;
      }

      /// \brief Get an element reference.
      /// \param j  Index of element
      T& operator()(index_t j)
      {
	return *getElmIterator(j);
      }

      /// \brief Get an element (const version).
      /// \param j  Index of element.
      T operator()(index_t j) const
      {
	return *getElmIterator(j);
      }

      /// \brief Fill with a value
      /// \param value  Value to fill with
      void fill(const T value)
      {
	std::fill(data_.begin(), data_.end(), value);
      }

      /// \brief Assigment operator that fills with value
      simple_vector& operator=(const T value)
      {
	this->fill(value);
	return(*this);
      }


      /// \brief Scale with scalar
      /// \param factor  Scaling factor
      simple_vector& operator*=(const T factor)
      {
	for (index_t j = 1; j <= length_; ++j)
	  data_[j] *= factor;

	return(*this);
      }


      /// \brief Copy from simple_vector<T>.
      /// \param sec   Vector to copy.
      simple_vector<T>& operator=(const simple_vector<T>& src)
      {
	this->resize(src.length());
	copy2(data_.begin(), src.data_.begin(), length_);
      }

      
  };


// ******************** GLOBAL FUNCTIONS *************************

  /// \brief Ostream support for simple_vector<T, OWNS>
  template<class T>
  std::ostream& operator<<(std::ostream& os, const simple_vector<T>& A)
  {
    os << "length = [" << A.length() << "]" << std::endl;

    index_t j;
    for (j = 1; j <= A.length(); ++j)
      os << A(j) << "  ";
    os << std::endl;
    return os;

    
  }



} // namespace simple_linalg





#endif // _SIMPLE_VECTOR_HPP_
