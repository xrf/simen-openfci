#ifndef _RADIAL_POTENTIAL_HPP_
#define _RADIAL_POTENTIAL_HPP_

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


#include <cstdlib>
#include <vector>
#include "linalg.hpp"

namespace quantumdot
{

  /// \brief This class encapsulates a quite general interaction
  /// potential in d=2 dimensions.
  ///
  /// The potential has the form:
  /// \f[   U(r) = \lambda r^{-\alpha} p(r) e^{-\beta r^2}   \f]
  /// The potential can be integrated exactly 
  /// using generalized half-range Hermite quadrature. One may choose to use 
  /// standard Hermite quadrature as well, but this requires alpha = 1  AND even p
  /// in order to be sensible.
  ///
  /// Use setIntegrationMethod() to set the integration method:
  /// STANDARD_HERMITE_QUAD or GENERALIZED_HERMITE_QUAD The integrator
  /// will automatically detect if alpha = 1 is true and p is even, in
  /// addition to whether the basis type is compatible. (If one uses
  /// GENERALIZED_HERMITE_BASIS, only GENERALIZED_HERMITE_QUAD may be
  /// used.)
  ///
  ///
  class RadialPotential
  {
    public:
      /// Enum type for selection of quadrature rules
      enum quad_enum { STANDARD_HERMITE_QUAD, GENERALIZED_HERMITE_QUAD };
      /// Enum type for basis type
      enum basis_enum { GENERALIZED_LAGUERRE_BASIS, GENERALIZED_HERMITE_BASIS };
     

    private:
      double lambda;         ///< Potential strength
      double alpha;          ///< Parameter
      double beta;           ///< Parameter
      std::vector<double> p; ///< Parameter (polynomial, p[0] = highest order coeff, p[n+1] = constant coeff.
      bool p_iseven;         ///< Parameter, true if p is even polynomial. (Each coeff is treated as r^2j-coeff instead of r^j-coeff)
      quad_enum quad_type;   ///< Which integration technique to use. STANDARD_HERMITE or GENERALIZED_HERMITE.
      basis_enum basis_type; ///< Which basis polynomials to use. GENERALIZED_HERMITE_BASIS or GENERALIZED_LAGUERRE_BASIS.

    public:
      /// \brief Default constructor. Constructs standard Coulomb potential
      RadialPotential()
      {
	// Set default potential parameters: Coulomb potential, lambda = 1.
	alpha = 1;
	beta = 0;
	p.resize(1);
	p_iseven = true;
	p[0] = 1.0;
	lambda = 1;
	quad_type = STANDARD_HERMITE_QUAD;
	basis_type = GENERALIZED_LAGUERRE_BASIS;
      }

      /// \brief Set integration type. Possible are STANDARD_HERMITE_QUAD and GENERALIZED_HERMITE_QUAD.
      void setIntegrationMethod(quad_enum x) { quad_type = x; }

      /// \brief Set basis type. Possible are GENERALIZED_LAGUERRE_BASIS and GENERALIZED_HERMITE_BASIS.
      void setBasisType(basis_enum x) { basis_type = x; }

      /// Set whether p is even or odd.
      void setPIsEven(bool x) { p_iseven = x; }

      /// Set p.
      void setP(std::vector<double> poly) 
      { 
	p = poly; 
      }

      /// Set alpha
      void setAlpha(double x) { alpha = x; }
      /// Set beta
      void setBeta(double x) { beta = x; }
      /// Set lambda
      void setLambda(double x) { lambda = x; }


      /// \brief Compute a matrix.
      /// \param m   Angular momentum quantum number
      /// \param C   Reference to destination matrix
      /// \param n_max  Compute elements for 0 <= n, n' <= n_max.
      void computeMatrix(dense_matrix& C, int n_max, int m = 0);

      /// \brief Compute effective interaction matrix.
      ///
      /// This is always with respect to Laguerre polynomial basis.
      /// The matrix is built in Generalized Hermite basis first
      /// in order to extract the exact effective interaction.
      /// The *result* is in Laguerre basis.
      /// (Integration method and basis type are not affected after function is finished.)
      ///
      /// \param Ceff    Result matrix.
      /// \param n_max   Ceff will be (n_max+1) times (n_max + 1)
      /// \param m       Angular momentum, default is 0.
      void computeEffectiveMatrix(dense_matrix& Ceff, int n_max, int m = 0);

      /// \brief Print brief info to stream.
      std::ostream& operator<<(std::ostream& os)
      {
	using namespace std;
	os << "Potential:" << endl;
	os << "  r^(-"<<alpha<<")*p(r)*exp(-r*"<<beta<<")"<<endl;
	return os;
      }


  };



  /// \brief Compute the harmonic oscillator matrix at
  /// given angular momentum in generalized half-range Hermite
  /// basis. Used for computing effective interaction with RadialPotential::computeEffectiveMatrix().
  ///
  /// \param H0   Destination matrix
  /// \param N    Dimension of H0 will be N+1.
  /// \param m    Angular momentum.
  void computeHOMatrixGenHer(dense_matrix& H0, int N, int m);

  /// Compute transition matrix from generalized Half-range Hermite basis
  /// to Laguerre basis (with angular momentum). Used in RadialPotential::computeEffectiveMatrix().
  /// \param U   Destination matrix
  /// \param N_lag  Use N_lag + 1 Laguerre polynomials
  /// \param N_her  Use H_her + 1 Hermite polynomials
  /// \param m  Angular momentum.
  void computeLagGenHerTrafo(dense_matrix& U, int N_lag, int N_her, int m);

  


}





#endif // _RADIAL_POTENTIAL_HPP_
