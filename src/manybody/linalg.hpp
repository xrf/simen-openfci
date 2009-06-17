#ifndef _LINALG_HPP_
#define _LINALG_HPP_

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
 *
 * \file linalg.hpp
 * \author Simen Kvaal
 * \date 10-11-2008
 *
 * \brief Typedefs and inclusions for dense linear algebra.
 *
 **/

#include "simple_vector.hpp"
#include "simple_matrix.hpp"
#include "simple_linalg.hpp"
#include <lpp/lapack.hh>

/// \brief Dense matrix
typedef simple_dense::simple_matrix<double> dense_matrix;
/// \brief Dense vector
typedef simple_dense::simple_vector<double> dense_vector;
/// \brief Dense matrix of ints
typedef simple_dense::simple_matrix<int> dense_matrix_int;



#endif // _LINALG_HPP_
