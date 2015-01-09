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

#include "QdotFci.hpp"

//#include "arpack++/arssym.h"
#include "../arpack/arpack.hpp"

#include "MatrixMachine.hpp"
#include "timer.hpp"


using namespace std;

namespace quantumdot {

  void QdotFci::buildHamiltonian()
  {

    // Do some initial matlab output.
    if (matlab.is_open())
    {
      matlab << "%" << endl << "% Output from FCI_CODE:" << endl << "%" << endl;
      matlab << "A = " << A << ";" << endl;
      matlab << "R = " << R << ";" << endl;
      matlab << "S = " << S << ";" << endl;
      matlab << "M = " << M << ";" << endl;
      matlab << "Sz = " << Sz << ";" << endl;
      matlab << "use_energy_cut = " << use_energy_cut << ";" << endl;
      matlab << "lambda = " << format(lambda) << ";" << endl;
    }

    // Initialize the Hilbert space machine
    hilbert.setParticles(A);
    hilbert.setSpinValues(Sz,S);
    hilbert.setHilberSpaceParameters(R, use_energy_cut); // This also generates SP basis.
    hilbert.setAngmom(M);

    // Generate Slater determinant basis and check that its functions
    // obey constant M and Sz, has A particles.
    hilbert.generateSDBasis();
    hilbert.printBasis(cerr);
    //hilbert.generateSDBasisBruteForce();

    if (hilbert.getSDBasis().size() == 0)
    {
      cerr << "No states found! Exiting." << endl;
      exit(0);
    }

    hilbert.checkBasis();


    // Build CSFs.
    states.setSDBasis( hilbert.getSDBasis() );

    cerr << "Building CSF basis ..." << endl;

    states.rearrangeBasis();

    // Re-attach the Hilbert space, for completeness
    hilbert.setSDBasis( states.getSDBasis() );
    //hilbert.printBasis(cerr);

    states.findBlocks();
  
    // set lambda in inter.
    double alpha = pow(1 + 0.25*B*B, 0.25);
    double lambda_star = lambda / alpha;
    if (B != 0.0)
    {
      cerr << "Correction for B-field gives: lambda_star = " << format(lambda_star) << endl;
      cerr << "Energy unit is now: " << pow(1 + 0.25*B*B, 0.5) << endl;
    }
    inter.setLambda(lambda_star);
  
    inter.setRadialPotential(potential);

    inter.setR( R );
    if (!use_veff)
      inter.buildInteractionComBlocks(); // just compute Coulomb ...
    else
    {
      if (use_energy_cut)
	inter.buildEffectiveInteractionComBlocks(1);
      else
	inter.buildEffectiveInteractionComBlocks(2);
    }

    if (precompute_lab_frame)
      inter.precomputeLabFrameMatrix();

    // Attach the single particle basis and the interaction
    mat.setSPBasis( hilbert.getSPBasis() );
    mat.setInteraction( inter );

    // Create a matrix-builder
    manybody::MatrixMachine<QdotMatrixElements> m_builder(mat);
    m_builder.setSDBasis( hilbert.getSDBasis() );

    // Obtain all CSF blocks for the currens value for S.
    vector<manybody::csf_block> csf;
    //states.getAllBlocks(csf, (double)0.5*S );
    states.getAllBlocks(csf, S );
    m_builder.setCSFBlocks( csf );

    timer t;


    // Build the matrices.
    cerr << "Building H0 ..." << endl;
    m_builder.messages(false);
    mat.setWhichMatrix(0);
    m_builder.setRank(1);
    t.start();
    m_builder.buildMatrix(H0, true);
    t.stop();
    cerr << "That took " << t() << " seconds." << endl << endl;
    double t1 = t();

    cerr << "Building interaction H1 ..." << endl;
    if (matrix_messages)
      m_builder.messages(true);
    mat.setWhichMatrix(1);
    m_builder.setRank(2);
    t.restart();
    m_builder.buildMatrix(H1);
    t.stop();
    cerr << "That took " << t() << " seconds." << endl << endl;

    if (matlab.is_open())
    {
      matlab << "T_build = [" << format(t1) << ", " << format(t()) << "];" << endl;
      matlab << "nonzeroes = " << H1.data.size() << ";" << endl;
      matlab << "dim = " << (H1.row_ptr.size()-1) << ";" << endl;
    }

  }



  void QdotFci::diagonalize(int nev)
  {
    size_t dim = H1.getNrows();
  
    // we cannot compute more eigenvalues than dim-1 using ARPACK.
    int nev2 = min((int)dim-1, nev);

    if (dim > 1)
    {
    
//      // Create a ARPACK++ object instance.
//      ARSymStdEig<double, QdotFci> 
//	eigensolver(dim, nev2, this, 
//		    &QdotFci::applyHamiltonian, "SM");

      // Create an interface to ARPACK
      arpack::symmetric_eigensolver<QdotFci> 
	eigensolver(this, &QdotFci::applyHamiltonian);

      eigensolver.set_parameters(dim, nev2, 0.0);
    
//      eigensolver.Trace(5, 0, 1, 1, 1, 0, 0, 0, 0);
      cerr << "Trying to find " << nev2 << " (matrix dim = " << dim << ") eigenpairs ..." << endl;
    
      timer t;
      t.start();
//      size_t nconv = eigensolver.FindEigenvectors();
      size_t nconv = eigensolver.solve();
      t.stop();
    
      cerr << "ARPACK found " << nconv << " eigenpairs in " << t() << " s. " << endl;
      if (matlab.is_open())
      {
	matlab << "T_diag = " << format(t()) << ";" << endl;
      }
    
      // Bail out if we didn't find everything.
      assert(nconv == (size_t)nev2);
    
      // Allocate space for result.
      eigenvalues.resize(nconv);
      eigenvectors.resize(dim, nconv);
    
      //
      // Save eigenvalues
      //
      for (size_t k = 0; k<nconv; ++k)
      {
	eigenvalues(k+1) = eigensolver.get_eigenvalue(k);
	for (size_t j = 0; j<dim; ++j)
	  eigenvectors(j+1, k+1) = eigensolver.get_eigenvector(k,j);
      }
    
//       for (size_t k = 0; k<nconv; ++k)
//       {
// 	eigenvalues(k+1) = eigensolver.Eigenvalue(k);
// 	for (size_t j = 0; j<dim; ++j)
// 	  eigenvectors(j+1, k+1) = eigensolver.Eigenvector(k,j);
//       }

    }
    else
    {
      // dim == 1.

      eigenvalues.resize(1);
      eigenvectors.resize(1,1);
      eigenvectors(1,1) = 1;
      eigenvalues(1) = H0.data[0] + H1.data[0];

    }

  }


  void QdotFci::postProcess()
  {

    double omega = pow(1 + 0.25*B*B, 0.5);

    // Multiply eigenvalues with omega and add angular momentum
    // dependent shift.
    for (size_t k = 1; k<=static_cast<size_t>(eigenvalues.length()); ++k)
      eigenvalues(k) = eigenvalues(k)*omega - 0.5*B*M;


  }


} // namespace quantumdot
