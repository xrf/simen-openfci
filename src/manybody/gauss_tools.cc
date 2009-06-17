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

#include "gauss_tools.hpp"
#include "tools.hpp"

namespace gauss {


  void computeHermitePolys(dense_vector& x, dense_matrix& P, int N)
  {
    // allocate matrix
    P.resize(x.length(), N+1);

    int nx = x.length();

    assert(N>=0);

    int n,j;

    // Generate first two polynomials.
    for (j=1; j<=nx; ++j)
    {
      P(j,1) = pow(M_PI,-0.25);
      if (N>0)
	P(j,2) = sqrt(2)*pow(M_PI,-0.25)*x(j);

    }

    // Generate the rest by recursion
    for (n = 1; n<N; ++n)
    {
      for (j=1; j<=nx; ++j)
	P(j,n+2) = sqrt(2.0/(n+1.0))*x(j)*P(j,n+1) - sqrt(n/(n+1.0))*P(j,n);
    }
  

  }



  void computeLaguerrePolys(dense_vector& x, dense_matrix& P, int N, double alpha)
  {
    // allocate matrix
    P.resize(x.length(), N+1);

    int nx = x.length();

    assert(N>=0);

    int n,j;

    // Generate first two polynomials.
    for (j=1; j<=nx; ++j)
    {
      P(j,1) = 1.0;
      if (N>0)
	P(j,2) = -x(j) + alpha + 1.0;

    }

    // Generate the rest by recursion
    for (n = 1; n<N; ++n)
    {
      for (j=1; j<=nx; ++j)
	P(j,n+2) = 1.0/(n+1)*((2*n+1+alpha-x(j))*P(j,n+1) - (n+alpha)*P(j,n));
    }
  

  }

  void computeNormalizedLaguerrePolys(dense_vector& x, dense_matrix& P, int N, double alpha)
  {
    // allocate matrix
    P.resize(x.length(), N+1);

    int nx = x.length();

    assert(N>=0);

    int n,j;

    // Generate first two polynomials.
    for (j=1; j<=nx; ++j)
    {
      //P(j,1) = 1.0;
      P(j,1) = 1.0/exp(0.5*gamma(alpha+1));
      if (N>0)
	//P(j,2) = -x(j) + alpha + 1.0;
	P(j,2) = (-x(j) + alpha + 1.0)/exp(0.5*gamma(alpha+2));

    }

    // Generate the rest by recursion
    for (n = 1; n<N; ++n)
    {
    
      double c1 = sqrt((n+alpha+1)*(n+1));
      double c2 = sqrt(n*(n+alpha));
      for (j=1; j<=nx; ++j)
      {
	P(j,n+2) = (2*n+1+alpha-x(j))*P(j,n+1)/c1 - (c2/c1)*P(j,n);
      }


    }
  

  }


  void computeGaussHermite(dense_vector& x, dense_vector& w, long int n)
  {

  
    x.resize(n);
    w.resize(n);

    //dense_vector d(n-1);
    dense_matrix J(n,n);

    int j, k;

    for (j=1; j<=n; ++j)
      for (k=1; k<=n; ++k)
	J(j,k) = 0;

    for (j = 1; j<=(n-1); ++j)
    {
      //d[j] = sqrt(.5*(j+1));
      J(j+1,j) = sqrt(.5*(j));
      J(j,j+1) = sqrt(.5*(j));
    }

    // diagonalize it.
    long int info = 0;
    lpp::syev("v","u", &n, &J(1,1), &n, &x(1), &info);

    for (j=1; j<=n; ++j)
    {
      w(j) = sqrt(M_PI)*J(1,j)*J(1,j);
    }


  }

  void computeGaussLaguerre(dense_vector& x, dense_vector& w, long int n)
  {

  
    x.resize(n);
    w.resize(n);

    dense_matrix J(n,n);

    int j, k;

    for (j=1; j<=n; ++j)
      for (k=1; k<=n; ++k)
	J(j,k) = 0;

    for (j = 1; j<=(n-1); ++j)
    {
      J(j,j) = 2*j-1;
      J(j+1,j) = -j;
      J(j,j+1) = -j;
    }
    J(n,n) = 2*n-1;

    // diagonalize it.
    long int info = 0;
    lpp::syev("v","u", &n, &J(1,1), &n, &x(1), &info);

    for (j=1; j<=n; ++j)
    {
      w(j) = J(1,j)*J(1,j);
    }


  }



  void computeGenHalfGaussHermite(dense_vector& x, dense_vector& w, double gam, long int N)
  {

    // Decrement N so N points instead of N+1 points are generated.
    // NB: This is different from the matlab implementation.

    N--;

    dense_vector alpha, beta;


    computeGenHalfHermiteCoeffs(alpha, beta, gam, std::max(1L,N));

    // compute the Golub-Welsch matrix J.

    dense_matrix J(N+1,N+1);
    J = 0;
    for (int n = 0; n<=N; ++n)
    {
      J(n+1,n+1) = alpha(n+1);
      if (n>0)
	J(n+1,n) = sqrt(beta(n+1));
    
      if (n<N)
	J(n+1,n+2) = sqrt(beta(n+2));
    
    }
  
    // compute abscissa and weights.

    long int info = 0;
    x.resize(N+1);
    w.resize(N+1);
    long int N1 = N+1;
    lpp::syev("v","u", &N1, &J(1,1), &N1, &x(1), &info);

    for (int j=1; j<=N1; ++j)
    {
      w(j) = J(1,j)*J(1,j) * exp(gamma((gam+1)/2))/2;
    }



  }


  void computeGenHalfHermiteCoeffs(dense_vector& alpha, dense_vector& beta, double gam, int N)
  {

    // Initially, we do min(N,8) evaluations by recursion.
    int Nrec = std::min( (long int)N, (long int)8 ); 

    dense_vector g(N+1);
    g = 0;

    // Compute first Nrec coefficients with direct recursion

    double alpha0 = exp(gamma(gam/2+1))/exp(gamma((gam+1)/2));

    g(1) =  -alpha0*alpha0 + ((1+gam)/2 - (1+gam/2)/6);

    double g0 = -gam/12;

    int n;

    for (n = 1; n<=Nrec-1; ++n)
    {
      double gp, gc;
      double Y = 2.0*n + gam;
      if (n == 1)
      {
	gp = g0;
	gc = g(n);
      }
      else
      {
	gp = g(n-1);
	gc = g(n);
      }
      double p1 = (sqr(Y/6 - gc) - sqr(gam/4));
      double p2 = (Y/12 + gc);
      double p3 = (Y-1)/3 -(gc + gp);
      //cerr << sqr(p1/p2)/p3 << endl;
      double gn = (Y+1)/3 - sqr(p1/p2)/p3 - gc;
      g(n+1) = gn;
      //cerr << "g(" << n+1 << ") == " << gn << endl;

    } // end of direct recursion

    // Compute the remaining N - Nrec using asyptotic formula

    double C0 = (double)1/36 - sqr(gam)/8;
    double C1 = (double)23/432 - (double)11/48*sqr(gam) + (double)3/32*pow(gam,4);
    double C2 = (double)1189/2592 - (double)409/192*sqr(gam) 
      + (double)75/64*pow(gam,4) + (double)9/64*pow(gam,6);
    double C3 = (double)196057/20736 - (double)153559/3456*sqr(gam) 
      + (double)7111/256*pow(gam,4) + (double)639/128*pow(gam,6) + (double)135/512*pow(gam,8);


    for (n = Nrec+1; n<=N+1; ++n)
    {
      double Y = 2*n + gam;
      g(n) = C0/Y + C1/pow(Y,3) + C2/pow(Y,5) + C3/pow(Y,7);
    }

    // Finally, use Newton iteration to make the coefficients converge properly.

    int iter;
    int maxiter = 2; // After 2 iterations, it *is* fully converged, in all my experiments!!
    dense_vector F(N);
    dense_matrix J(N,N);
    for (iter = 1; iter <= maxiter; ++iter)
    {
      F = 0; // clear the function
      J = 0; // clear the Jacobian

      // Fill function/Jacobian
      for (int j = 1; j <= N; ++j)
      {
	double Y = 2*j + gam;
	double gn = g(j+1);
	double gc = g(j);
	double gp;
	if (j==1)
	  gp = g0;
	else
	  gp = g(j-1);
        
	double p1 = ((Y+1)/3 - gn - gc);
	double p2 = ((Y-1)/3 - gc - gp);
	double p3 = (Y/12 + gc);
	double p4 = (sqr(Y/6 - gc) - sqr(gam)/16);
	F(j) = p1 * p2 * sqr(p3) - sqr(p4);
        
	J(j,j) = -1*p2*sqr(p3) -p1*sqr(p3) + 2*p3*p1*p2 + 4*p4*(Y/6 - gc);
	if (j>1)
	  J(j,j-1) =-p1*sqr(p3);
      
	if (j<N)
	  J(j,j+1) = -p2*sqr(p3);
      
        
     
      } // end of fill func/Jacobian

      // Okay, now we do the actual Newton step.

      simple_dense::solveTridiagonal(J, F);

      // (note that g(N+1) is held fixed.)
      for (int j=1; j<=N; ++j)
	g(j) = g(j) - F(j);
    
    } // end of Newton iterations.
    
    // now we have the coefficients g(n), n = 0..N.
    // compute the alpha and beta coeffs.
    
    alpha.resize(N+1);
    beta.resize(N+1); 
    alpha = 0;
    beta = 0;
    alpha(1) = alpha0;
    beta(1) = 0;
  
    for (n=1; n<=N; ++n)
    {
      beta(n+1) = g(n) + (n+gam/2)/6;
      alpha(n+1) = sqrt((2*n+gam+1)/3 - g(n+1) - g(n));
    }
  


  }



  void computeGenHalfHermiteCoeffsNormalized(dense_vector& a, dense_vector& b, dense_vector& c, double gam, int N)
  {

    // Note: one more coefficient is needed for the computation
    // of normalized coefficients.
    dense_vector alpha(N+2), beta(N+2);
    computeGenHalfHermiteCoeffs(alpha, beta, gam, N+1);
  
    a.resize(N+1);
    b.resize(N+1);
    c.resize(N+1);
    for (int n=1; n<=N+1; ++n)
    {
      a(n) = 1/sqrt(beta(n+1));
      b(n) = a(n)*alpha(n);
      c(n) = sqrt(beta(n)/beta(n+1));
    }

  }



  void computeNormalizedGenHalfHermitePolys(dense_vector& x, dense_matrix& P, double gam, int N, bool diff)
  {
    // allocate matrix
    if (diff)
      P.resize(x.length(), 2*(N+1));
    else
      P.resize(x.length(), (N+1));

    int K = N+1; // index for first derivative, i.e. P(:,K+n) == dP(:,n)

    // zero out everything.

    P = 0;

    int nx = x.length();

    assert(N>=0);

    // compute recursion coefficients.
    dense_vector a, b, c;
    computeGenHalfHermiteCoeffsNormalized(a, b, c, gam, N);

    int n,j;

    // Generate first polynomial -- a constant.
    for (j=1; j<=nx; ++j)
      P(j,1) = pow(0.5*exp(gamma((gam+1)/2)), -0.5);

    // Generate the rest by recursion
    for (n = 0; n<N; ++n)
    {
      if (n>0)
      {
	for (j=1;j<=nx; ++j)
	{
	  P(j,n+2) = (a(n+1)*x(j) - b(n+1))*P(j,n+1) - c(n+1)*P(j,n);
	  if (diff)
	    P(j,K + n+2) = (a(n+1)*x(j) - b(n+1))*P(j,K + n+1) - c(n+1)*P(j,K+n) + a(n+1)*P(j,n+1);
	}
      }
      else
      {
	for (j=1;j<=nx; ++j)
	{
	  P(j,n+2) = (a(n+1)*x(j) - b(n+1))*P(j,n+1);
	  if (diff)
	    P(j,K + n+2) = (a(n+1)*x(j) - b(n+1))*P(j,K + n+1) + a(n+1)*P(j,n+1);

	}
      }

    }
  

  }


} // namespace gauss
