#ifndef _ARPACK_HPP
#define _ARPACK_HPP


//
// Copyright (c) 2009 Simen Kvaal
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
 * \file arpack.hpp
 * \date 6-17-09
 * \author Simen Kvaal

 * \brief Defines a class symmetric_eigensolver; a wrapper around the 
 * Fortran77 ARPACK library.
 *
 * It relies on (slightly modified) files culled from the arpack++ library
 * (being public domain), being basically prototypes for function calls to
 * central ARPACK subroutines.
 *
 */

#include "saupp.h"
#include "seupp.h"
#include <cassert>
#include <iostream>

namespace arpack
{
  
  /// \brief Simple wrapper for symmetric, standard eigenvalue problem.
  ///
  /// Currently only supports finding a few of the lowest eigenvalues;
  /// other parts of the spectrum will be added later.
  ///
  template<class matrix_t>
  class symmetric_eigensolver
  {
    private:
      /// \brief Define type for operation function.
      typedef void (matrix_t::* operation_func_t)(double*, double*);

      bool debug;                 ///< Turn on/off messages.

      matrix_t *mat;              ///< Pointer to class instance that holds the matrix-vector operation function.
      operation_func_t matrix_op; ///< Function that evaluates matrix-vector operation required for ARPACK.
      
      ARint n;                    ///< Dimension of problem
      ARint nev;                  ///< Number of eigenvalues to search for
      ARint ncv;                  ///< Number of Lanczos vectors generated at each iteration
      double tol;                 ///< Tolerance. Set to 0 to use machine precision

      ARint ido;                  ///< Reverse communication job indicator
      char bmat;                  ///< 'I' = standard eigenvalue problem, 'G' = generalized
      char which[3];              ///< Part of spectrum. Should be set to 'SM' in this simple implementation

      std::vector<double> resid;       ///< Residual vector; for reverse communication
      std::vector<double> V;           ///< Matrix (as a vector); for reverse communication
      ARint ldv;                  ///< Leading dimension of V
      std::vector<ARint> iparam;       ///< For reverse communication
      std::vector<ARint> ipntr;        ///< For reverse communication
      std::vector<double> workd;       ///< For reverse communication
      std::vector<double> workl;       ///< For reverse communication
      ARint lworkl;               ///< Length of workl
      ARint info;                 ///< Status flags from ARPACK calls

      ARint maxit;                ///< Max number of iterations

      ARint nconv;                ///< Number of eigenvalues that converged
      std::vector<double> d;           ///< Eigenvalues
      std::vector<double> Z;           ///< Eigenvectors

    public:
      /// \brief Attach an object instance and matrix-vector operation function.
      void set_matrix_op(matrix_t *the_mat, void (matrix_t::* the_matrix_op)(double*, double*))
      {
	mat = the_mat;
	matrix_op = the_matrix_op;
      }

      /// \brief Set parameters for problem.
      void set_parameters(ARint the_n, ARint the_nev,  double the_tol, ARint the_ncv = 0)
      {
	assert(the_n > 0);
	assert(the_nev > 0);
	assert(the_nev < the_n);
	
	n = the_n;
	nev = the_nev;
	tol = the_tol;
	if (the_ncv == 0)
	{
	  ncv = min(2*nev+1,n);
	}
	else
	  ncv = the_ncv;
	
      }
      

      /// \brief Default constructor.
      symmetric_eigensolver() 
      { 
	maxit = 1000;
	debug = true;
      }

      /// \brief Constructor that assigns object ant matrix-vector operation
      symmetric_eigensolver(matrix_t *the_mat, void (matrix_t::* the_matrix_op)(double*, double*))  
      {	
	debug = true;
	maxit = 1000;
	set_matrix_op(the_mat, the_matrix_op);
      }

      /// \brief Turn on/off debugging
      void set_debug(bool d) { debug = d; }



      /// \brief Solve eigenvalue problem.
      ARint solve()
      {
	

	bool finished;
	
	finished = false;

	ido = 0;
	bmat = 'I';
	strcpy(which,"SM");
	tol = 1e-8;
	resid.resize(n); // initial vector if info != 0.
	info = 0; // use random initial vector.
	V.resize(ncv*n+1);
	ldv = n;
	iparam.resize(12); // set them!
	iparam[1] = 1;
	iparam[3] = maxit;
	iparam[4] = 1;
	iparam[7] = 1;
	ipntr.resize(12);
	
	workd.resize(3*n + 1);
	lworkl = ncv*(ncv+8);
	workl.resize(lworkl + 1);
	
	
	int it = 1;

	while (finished == false)
	{
	  it++; // increase iteration counter
	  if (debug)
	  {
	    std::cerr << "Iteration number " << it << " commencing." << endl;
	  }
	  saupp(ido, bmat, n, which, nev, tol, &resid[0],
	  ncv, &V[0], ldv, &iparam[0], &ipntr[0], &workd[0], &workl[0], lworkl, info);
	  
	  if (debug)
	    std::cerr << "(ido, info) == " << ido << ", " << info << std::endl;
	  
	  if ((ido == -1) || (ido == 1))
	  {
	    (mat->*matrix_op)(&workd[ipntr[1]], &workd[ipntr[2]]);
	  }
	  else
	    finished = true;
	  
	  
	}
	
	if (info != 0)
	{
	  std::cerr << "Warning! Exit status of ARPACK was info == " << info << "." << std::endl;
	}
	

	nconv = iparam[5];
	bool rvec = true;
	double sigma = 0;

	if (debug)
	{
	  std::cerr << "Found nconv == " << nconv << " out of nev == " << nev << " eigenvalues." << std::endl;
	}
	
	d.resize(n);
	Z.resize(nconv*n);
	
	seupp(rvec, 'A', &d[0], &Z[0],
  	  n, sigma, bmat, n,
	  which, nev, tol, &resid[0],
	  ncv, &V[0], ldv, &iparam[0],
	  &ipntr[0], &workd[0], &workl[0],
	  lworkl, info);
	
	
	return nconv;
	
      }
      

      /// \brief Store eigenvalues
      void eigenvalues(std::vector<double>& evals)
      {
	evals.resize(nconv);
	for (size_t k = 0; k<(size_t)nconv; ++k)
	  evals[k] = d[k];
      
      }

      /// \brief Get a single eigenvalue
      double  get_eigenvalue(size_t k)
      {
	return d[k];
      }

      /// \brief Get single component k of eigenvector j.
      double get_eigenvector(size_t j, size_t k)
      {
	return Z[n*j + k];
      }

  };
  

} // namespace arpack


#endif
