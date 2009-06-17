#ifndef _GAUSS_TOOLS_HPP_
#define _GAUSS_TOOLS_HPP_

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

#include "linalg.hpp"

/**
 * \file gauss_tools.hpp
 * \author Simen Kvaal
 * \date 6-13-08, update 10-7-2008
 *
 * \brief Tools for Gaussian quadrature and orthogonal polynomials.
 *
 */

namespace gauss {


  /// \brief Compute \em normalized Hermite polynomials:
  /// P_n(x) = (1/sqrt[2^n*n!*sqrt(pi)])*H_n(x),
  ///
  /// Nice LaTex formula:
  /// \f[ P_n(x) = [2^n n! \sqrt{\pi}]^{-1/2} H_n(x) \f]
  ///
  /// H_n(x) are the regular Hermite polynnomials; 
  /// they are orhogonal wth weight w(x) = exp(-x^2).
  ///
  /// Recursion formula: P_0(x) = pi^(-1/4)
  /// P_1(x) = sqrt(2/sqrt(pi))*x
  /// P_(n+1)(x) = sqrt(2/(n+1))*x*P_n(x) - sqrt(n/(n+1))*P_(n-1)(x).
  ///
  /// Nice LaTeX formulas:
  /// \f[ P_0(x) = \pi^{-1/4} \f]
  /// \f[ P_1(x) = \sqrt{2} x P_0 \f]
  /// \f[ P_{n+1}(x) = \sqrt{\frac{2}{n+1}} x P_n(x) - \sqrt{\frac{n}{n+1}} P_{n-1}(x) \f]
  ///
  /// N = order of highest polynomial.
  /// x = points at which we evaluate P_n(x), 0<=n<=N.
  /// Result: P(i,n+1) = P_n(x_i); a matrix.
  ///
  void computeHermitePolys(dense_vector& x, dense_matrix& P, int N);


  /// \brief Compute generalized Laguerre polynomials:
  /// P^a_n(x).
  ///
  /// Recursion formula: P_0(x) = 1
  /// P_1(x) = 1+alpha-x
  /// P_(n+1)(x) = (1/(n+1)*[(2n+1+alpha-x)P_n(x) - (n+alpha)P_(n-1)(x)].
  ///
  /// LaTeX formulas:
  /// \f[ P_0(x) = 1 \f]
  /// \f[ P_1(x) = 1 + \alpha - x \f]
  /// \f[ (n+1)P_{n+1}(x) = (2n + 1 + \alpha - x)P_n(x) - (n+\alpha)P_{n-1}(x) \f]
  ///
  /// P is resized automatically. x is not altered.
  ///
  /// \param N order of highest polynomial.
  /// \param x points at which we evaluate P_n(x), 0<=n<=N.
  /// \param P Output: P(i,n+1) = P_n(x_i); a matrix.
  void computeLaguerrePolys(dense_vector& x, dense_matrix& P, int N, double alpha);

  /// \brief Compute \em normalized generalized Laguerre polynomials.
  ///
  /// The normalization constant is sqrt(Gamma(n+1)/Gamma(n+alpha+1)),
  /// which can be used to rederive the reecursion formula. This is
  /// then implemented here.
  /// 
  /// See also computeLaguerrePolys() for other details.
  ///
  /// P is resized automatically. x is not altered.
  ///
  /// \param N order of highest polynomial.
  /// \param x points at which we evaluate P_n(x), 0<=n<=N.
  /// \param P Output: P(i,n+1) = P_n(x_i); a matrix.
  void computeNormalizedLaguerrePolys(dense_vector& x, dense_matrix& P, int N, double alpha);


  /// \brief Compute Gauss-Hermite quadrature rules using Golub-Welsch algorithm.
  ///
  /// x becomes the abscissa, w the weights, n is the number of points.
  /// The vectors are resized automatically.
  ///
  /// \param x Output: abscissa
  /// \param w Output: weights
  /// \param n Number of quadrature points to generate.
  void computeGaussHermite(dense_vector& x, dense_vector& w, long int n);

  /// \brief Compute Gauss-Laguerre quadrature rules using Golub-Welsch algorithm.
  ///
  /// NOTE: This could be extended to generalized Laguerre polynomials ...
  ///
  /// x becomes the abscissa, w the weights, n is the number of points.
  /// The vectors are resized automatically.
  ///
  /// \param x Output: abscissa
  /// \param w Output: weights
  /// \param n Number of quadrature points to generate.
  void computeGaussLaguerre(dense_vector& x, dense_vector& w, long int n);


  /// \brief Compute generalized galf-range Gauss-Hermite quadrature rules.
  //
  /// x becomes the abscissa, w the weights, n is the number of
  /// points.  The vectors are resized automatically.
  ///
  /// Adapted from matlab script written earlier.  Note: The matlab
  /// script generates N+1 points; we generate N points.
  ///
  /// References: 
  /// * Ball, JS: "Half-range generalized Hermite polynomials and the
  ///   related Gaussian quadratures", SIAM J. Numer. Anal. <b>40</b> (2003), 
  ///   pp. 2311-2317
  /// * Golub, GH and Welsch, JH: "Calculation of Gauss Quadrature Rules",
  ///   Math. Comp. <b>23</b> (1969), 221 
  ///
  ///
  /// The weight function for these polynomials are
  ///     w(x) = x^gam * exp(-x^2), 0<x<+infinity,
  /// and are useful for radial Schroedinger equations
  /// using harmonic oscillator basis.
  //
  /// \param x Output: abscissa
  /// \param w Output: weights
  /// \param N Number of quadrature points to generate.
  /// \param gam Parameter
  void computeGenHalfGaussHermite(dense_vector& x, dense_vector& w, double gam, long int N);

  /// \brief Compute recursion coefficients of the monic, generalized
  /// half-range Hermite polynomials.
  ///
  /// They obey the recursion formula:
  ///
  /// p[n+1] = (x - alpha(n+1))*p[n] - beta(n+1)*p[n-1], n >= 0.
  ///
  /// LaTeX:
  /// \f[  p_{n+1}(x) = (x-\alpha_{n+1})p_n(x) - \beta_{n+1} p_{n-1}(x), \quad n \geq 0    \f]
  ///
  /// N+1 coefficients are stored in alpha and beta. The vectors are
  /// resized by the function call.
  ///
  /// References: 
  /// * Ball, JS: "Half-range generalized Hermite polynomials and the
  ///   related Gaussian quadratures", SIAM J. Numer. Anal. *40* (2003), 
  ///   pp. 2311-2317
  /// * Golub, GH and Welsch, JH: "Calculation of Gauss Quadrature Rules",
  ///   Math. Comp. *23* (1969), 221
  ///
  ///
  /// The weight function for these polynomials are
  ///     w(x) = x^gam * exp(-x^2), 0<x<+infinity,
  /// and are useful for radial Schroedinger equations
  /// of perturbed HOs.
  ///
  /// \param alpha Reference to vector which will hold half of the coefficients
  /// \param beta Reference to vector which will hold half of the coefficients  
  /// \param gam Parameter
  /// \param N Number of polynomials to create.
  void computeGenHalfHermiteCoeffs(dense_vector& alpha, dense_vector& beta, double gam, int N);

  /// \brief Compute coefficients of recursion of \em normalized
  /// generalized half-range Hermite polynomials.
  /// 
  /// p[n+1] = (a(n+1)*x - b(n+1))*p[n] - c(n+1)*p[n-1], n >= 0.
  ///
  /// Returns 3 * (N+1) coefficients.
  ///
  /// \param a Reference to vector which will hold coefficients
  /// \param b Reference to vector which will hold coefficients
  /// \param c Reference to vector which will hold coefficients
  /// \param gam Parameter
  /// \param N Number of poly-coeffs to generate
  void computeGenHalfHermiteCoeffsNormalized(dense_vector& a, dense_vector& b, dense_vector& c, double gam, int N);

  /// \brief Compute normalized generalized half-range Hermite
  /// polynomials of order 0..N evaluated at x.
  ///
  /// If diff == true, then the derivatives are also computed, and
  /// stored in P, after the polynomials themselves.
  ///
  /// P(:,j+1) is the polynomial number j, P(:,(N+1) + j+1) its derivative.
  ///
  /// \param x Points at which to evaluate
  /// \param P Matrix which will hold the polynomials.
  /// \param gam Parameter.
  /// \param N Generate N+1 polynomials
  /// \param diff Generate derivatives or not?
  void computeNormalizedGenHalfHermitePolys(dense_vector& x, dense_matrix& P, double gam, int N, bool diff);


} // namespace gauss


#endif // _GAUSS_TOOLS_HPP_
