#ifndef _QDOT_MATRIX_ELEMENTS_
#define _QDOT_MATRIX_ELEMENTS_

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


#include "../manybody/Slater.hpp"

/**
 * \file QdotMatrixElements.hpp
 * \author Simen Kvaal
 * \date 9-16-08 (Original date: 6-14-08)
 *
 * \brief Class definition of MatrixElements: matrix elements for 2D quantum dots. 
 *
 * This class is the "top level" class in the matrix element generator class for
 * the quantum dot problem. While QdotInteraction produces the spin-free interaction and
 * handles some intricate indexing things. QdotInteraction also builds effective interactions,
 * if that is desired. At the bottom is RadialPotential, which is attached to QdotInteraction; see
 * QdotInteraction for details.
 *
 *
 */

namespace quantumdot {


/// \brief Class declaration of QdotMatrixElements: Computes single matrix elements for MatrixMachine.
///
/// This class has (a pointer to) an instance of QdotInteraction,
/// in addition to a single particle space (fetched from a QdotHilbertSpace).
/// It implements the functions needed for a manybody::MatrixMachine to compute
/// the interaction matrix, as well as the HO matrix.
///
/// This class is the "top level" class in the matrix element generator class for
/// the quantum dot problem. While QdotInteraction produces the spin-free interaction and
/// handles some intricate indexing things. QdotInteraction also builds effective interactions,
/// if that is desired. At the bottom is RadialPotential, which is attached to QdotInteraction; see
/// QdotInteraction for details.
///
class QdotMatrixElements
{
  private:
    QdotInteraction* interaction;          ///< Pointer to a QdotInteraction instance
    std::vector<quantum_numbers> sp_basis; ///< Vector of SP states.
    int which_matrix;                      ///< Indicates if HO matrix (0) or interaction (1) is to be computed.

  public:
    /// \brief Default constructor. 
    QdotMatrixElements()
    {
      interaction = (QdotInteraction*)0;
    }

    /// \brief Assign single particle basis.
    /// \param sp Single particle basis instance
    void setSPBasis(std::vector<quantum_numbers>& sp)
    {
      sp_basis = sp;
    }

    /// \brief Assign interaction. Address is taken and stored in a pointer.
    /// Object is not copied, because it can be a lot of information.
    /// \param i Reference to interaction
    void setInteraction(QdotInteraction& i)
    {
      interaction = &i;
    }
    
    /// \brief Set which matrix the "outside world" (i.e., a MatrixMachine instance) tries
    /// to calculate.
    ///
    /// Legal values are:
    ///  - 0 == HO matrix (which is diagonal)
    ///  - 1 == interaction matrix.
    /// \param which Which matrix, 0 or 1.
    void setWhichMatrix(int which)
    {
      which_matrix = which;
      assert((which == 0) || (which == 1));
    }

    /// \brief Function that computes the matrix element.
    /// 
    /// This is called by manybody::MatrixMachine to compute a single matrix element
    /// See the documentation of MatrixMachine.hpp for details.
    double element(const size_t* alpha, const size_t* beta)
    {
      switch (which_matrix)
      {
	case 0: // HO matix element
	  if (alpha[0] == beta[0])
	    return (double)( sp_basis[alpha[0]].N + 1 );
	  else
	    return 0;

	case 1: // interaction matrix element
	  if ((sp_basis[alpha[0]].sigma == sp_basis[beta[0]].sigma) // conservation of spin
	      && (sp_basis[alpha[1]].sigma == sp_basis[beta[1]].sigma))	// conservation of spin
 	    return interaction->singleElement(sp_basis[alpha[0]].N, sp_basis[alpha[0]].m,
 					      sp_basis[alpha[1]].N, sp_basis[alpha[1]].m,
 					      sp_basis[beta[0]].N, sp_basis[beta[0]].m,  
 					      sp_basis[beta[1]].N, sp_basis[beta[1]].m);
	  else
	    return 0.0;


	default:
	  return 0.0;
      }

    }


};



} // namespace quantumdot

#endif // _MATRIX_ELEMENTS_
