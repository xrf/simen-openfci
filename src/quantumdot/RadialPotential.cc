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

#include "RadialPotential.hpp"
#include "gauss_tools.hpp"
#include "tools.hpp"

using namespace gauss;

namespace quantumdot 
{


  void RadialPotential::computeMatrix(dense_matrix& C, int n_max, int m)
  {
    
    int j, k;

    // Compute order of Gauss quadrature
    int deg_p = p.size() - 1;
    int gauss_order;
    if (basis_type == GENERALIZED_LAGUERRE_BASIS)
      gauss_order = 2*n_max + abs(m) + deg_p/2 + 1;
    else // GENERALIZED_HERMITE_BASIS
    {
      // Only GENERALIZED_HERMITE_QUAD makes sense.
      assert(quad_type == GENERALIZED_HERMITE_QUAD);
      gauss_order = n_max + deg_p/2 + 2;
    }
    // Compute weights and abscissa.
    dense_vector y, w;
    if (quad_type == GENERALIZED_HERMITE_QUAD)
      computeGenHalfGaussHermite(y, w, 1.0-alpha, gauss_order);
    else // quad_type == STANDARD_HERMITE_QUAD.
    {
      // this quadrature is only legal if alpha = 1.0 and if p is even!
      assert(alpha == 1.0);
      assert(p_iseven);
      computeGaussHermite(y, w, gauss_order);

      // scale weights to accomodate for double range.
      for (k = 1; k <= gauss_order; ++k)
      {
	w(k) *= 0.5;
      }

      // use reflection symmetry to eliminate about half of the points.
      //  ... not impemented ...

    }


    // Compute r and r^2 as function of y.
    double a = sqrt(2*beta + 1);
    dense_vector r(gauss_order), r2(gauss_order);
    for (k = 1; k <= gauss_order; ++k)
    {
      r(k) = y(k)/a;
      r2(k) = r(k)*r(k);
    }

    // Evaluate basis functions at r.
    dense_matrix Poly(gauss_order, n_max + 1);
    if (basis_type == GENERALIZED_LAGUERRE_BASIS)
    {
      computeNormalizedLaguerrePolys(r2, Poly, n_max, abs(m));
      Poly *= sqrt(2); // normalization ...
      // Multiply with r^(2|m|).
      for (j = 1; j <= gauss_order; ++j)
	for (k = 1; k <= n_max + 1; ++k)
	  Poly(j,k) *= m1pow(k-1)*pow(r(j), abs(m));
    }
    else
    {

      computeNormalizedGenHalfHermitePolys(r, Poly, 1.0,  n_max, false);

    }

    // Evaluate polynomial at sqrt(2)*r or 2*r^2.
    dense_vector p_eval(gauss_order);
    for (k = 1; k <=gauss_order; ++k)
    {
      p_eval(k) = 0;
      double arg = 1.0;
      for (j = p.size() - 1; j >= 0; --j)
      {
	p_eval(k) += arg*p[j];
	arg *= p_iseven ? 2*r2(k) : sqrt(2)*r(k);
      }

    }
    
    // Compute matrix elements.
    C.resize(n_max + 1, n_max + 1);
    double prefactor = lambda * pow(0.5, alpha*.5) * pow(a, alpha - 2);

    for (j = 1; j <= n_max + 1; ++j)
    {
      for (k = j; k <= n_max + 1; ++k)
      {
	double sum = 0;
	
	int ell;
	for (ell = 1; ell <= gauss_order; ++ell)
	{
	  sum +=  Poly(ell, j)*Poly(ell, k)*p_eval(ell)*w(ell);
	}
	C(j,k) = prefactor*sum;
	C(k,j) = C(j,k);

      }
    }


  }




  void RadialPotential::computeEffectiveMatrix(dense_matrix& Ueff, int n_max, int m)
  {

    // Dimension of "full" matrix, whose eigenpairs will be exact.
    // N_extra and N_min are determined from empirical observations.
    // m = 0 is worst case (large N_min) and N_extra is choosen valid for lambda <= 20, at least.
    const int N_extra = 20;
    const int N_min = 20;
    long int dim = std::max(2*n_max + abs(m) + 1, N_min) + N_extra;

    // Dimension of effective matrix Ueff.
    long int dim2 = n_max + 1;


    // Build exact problem
    dense_matrix T, V, H;
    computeHOMatrixGenHer(T, dim-1, abs(m));

    // Save quad/basis settings.
    basis_enum basis_type_temp = basis_type;
    quad_enum quad_type_temp = quad_type;

    // Integrate potential matrix.
    this->setBasisType(GENERALIZED_HERMITE_BASIS);
    this->setIntegrationMethod(GENERALIZED_HERMITE_QUAD);
    this->computeMatrix(V, dim-1);

    // Restore settings.
    this->setBasisType(basis_type_temp);
    this->setIntegrationMethod(quad_type_temp);

    // Compute total Hamiltonian
    //H = T + V;
    simple_dense::addMatrices(H, T, 1.0, V, 1.0);


    // Diagonalize.
    // The eigenvalues are automatically sorted, which we want.
    dense_matrix U = H;
    dense_vector E(dim);
    long int info;
    lpp::syev("v","u", &dim, &U(1,1), &dim, &E(1), &info);


    // Compute transformation coefficients from g.h. Hermite basis
    // to Laguerre basis.
    dense_matrix ToL;
    //computeLagGenHerTrafo(ToL, n_max, 2*n_max + abs(m), abs(m)); 
    computeLagGenHerTrafo(ToL, n_max, dim - 1, abs(m)); 


    // Compute projection matrix P.
    dense_matrix P(dim,dim);
  
    for (int j=1; j<=dim; ++j)
      for (int k=1; k<=dim; ++k)
	for (int n = 0; n<=n_max; ++n)
	  if ((j-1 <= 2*n + abs(m)) && (k-1 <= 2*n + abs(m)))
	    P(j,k) += ToL(n+1,j)*ToL(n+1,k);
    
    
    
    // Compute matrix of desired eigenvectors in U.
    dense_matrix PUP;
    //PUP = P*U;
    simple_dense::multiplyMatrices(PUP, P, U);
    dense_matrix temp;
    temp = PUP(simple_dense::range(1, dim), simple_dense::range(1,n_max+1));
    PUP = temp;
    

    // Compute economy size SVD of PUP.
    // Result: PUP = X*diag(sigma)*Yt.
    dense_matrix X(dim,dim2), Yt(dim2, dim2);
    dense_vector sigma(dim2);
    lpp::gesvd("s", "s", &dim, &dim2, &PUP(1,1), &dim, &sigma(1), 
	       &X(1,1), &dim, &Yt(1,1), &dim2, &info);
    
    // Compute new eigenvector set.
    dense_matrix Q;
    //Q = X*Yt;
    simple_dense::multiplyMatrices(Q, X, Yt);
    temp = Q(simple_dense::range(1, dim), simple_dense::range(1, dim2));
    simple_dense::multiplyMatrices(Q, ToL, temp);


    // Compute effective Hamiltonian
    dense_matrix Heff(dim2,dim2);
    dense_matrix Eeff(dim2,dim2);
    Eeff = 0;
    for (int j=1; j<=dim2; ++j)
    {
      Eeff(j,j) = E(j);
    }
    simple_dense::multiplyMatrices(Ueff, Q, Eeff);
    simple_dense::multiplyMatrices(Heff, Ueff, Q, true);
    Ueff = Heff;

    for (int j = 1; j<=dim2; ++j)
      Ueff(j,j) -= 2*(j-1) + abs(m) + 1;
    
    
  }




  void computeHOMatrixGenHer(dense_matrix& H0, int N, int m)
  {
    // New version of the code uses gam = 1, instead of
    // the "natural" gam = 2|m|+1.
    // It requires more basis functions (scales linearly with |m|)
    // but is much more stable.

    int M = abs(m);
    int gauss_points = N + M + 1;

    dense_vector r, w, r2, w2;
  
    computeGenHalfGaussHermite(r, w, 1.0, gauss_points);
    computeGenHalfGaussHermite(r2, w2, 0.0, gauss_points);
    dense_matrix P, P2;
    computeNormalizedGenHalfHermitePolys(r, P, 1.0, N, true);
    computeNormalizedGenHalfHermitePolys(r2, P2, 1.0, N, false);
    int K = N + 1;

    H0.resize(N+1,N+1);
    int n;
    for (n=0; n<=N; ++n)
    {
      int n2;
      for (n2 = 0; n2<=N; ++n2)
      {
	double temp = 0;
	for (int k=1; k<=gauss_points; ++k)
	{
	  temp +=  ( 0.5*P(k,n+1+K)*P(k,n2+1+K) + P(k,n+1)*P(k,n2+1)*sqr(r(k)) - 
		     0.5*(P(k,n+1)*P(k,n2+1+K) + P(k,n2+1)*P(k,n+1+K))*r(k) ) * w(k);//        P(k,n+1) * P(k,n2+1) * w(k);
	  temp += 0.5*M*M * P2(k,n+1)*P2(k,n2+1)*w2(k)/r2(k);
	}
	H0(n+1,n2+1) = temp ;
      }
    }



  }
  

  void computeLagGenHerTrafo(dense_matrix& U, int N_lag, int N_her, int m)
  {
    int n, j, k;
    int M = abs(m);


    // Compute Gauss-g.h. Hermite quadrature rules needed
    // for inner products
    int gauss_points = std::max(2*N_lag + M, N_her) + 1;
    dense_vector r, w, r2;
    computeGenHalfGaussHermite(r, w, 1.0, gauss_points);

 
    // Compute polynomials.
    dense_matrix P, L;
    computeNormalizedGenHalfHermitePolys(r, P, 1.0, N_her, false);

    r2.resize(gauss_points);
    for (j = 1; j<= gauss_points; ++j)
      r2(j) = sqr(r(j));

    computeNormalizedLaguerrePolys(r2, L, N_lag, M);


    // Compute transformation matrix
    U.resize(N_lag+1,N_her+1);
    U = 0;
    for (n = 0; n<= N_lag; ++n)
      for (k = 0; k<= N_her; ++k)
	for (j=1; j<= gauss_points; ++j)
	  U(n+1,k+1) += m1pow(n)*sqrt(2)*pow(r(j),M)*L(j, n+1)*P(j, k+1)*w(j);


    // Finished!

  }


}
