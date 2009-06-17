#ifndef _SLATER_HPP_
#define _SLATER_HPP_

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

#include "BitsetSlater.hpp"

/**
 * \file Slater.hpp
 * \author Simen Kvaal
 * \date 9-16-08
 *
 * \brief Some typedefs and definitions.
 *
 **/

namespace manybody {

#ifndef SLATER_BITS
/// The number of bits to use for a Slater determinant.
#define SLATER_BITS 512
#endif

/// \brief Define orbital_t as the indices used by the Slater class.
/// (Originally I played with different types of Slater classes using different
/// index types.)
typedef size_t orbital_t;
/// \brief We define a Slater thus:
typedef BitsetSlater<SLATER_BITS> Slater;

}

#endif
