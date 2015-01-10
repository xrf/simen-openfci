#ifndef _QDOT_FCI_HPP_
#define _QDOT_FCI_HPP_

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

#include <fstream>
#include <cassert>
#include <iostream>

#include "QdotHilbertSpace.hpp"
#include "QdotInteraction.hpp"
#include "QdotMatrixElements.hpp"
#include "RadialPotential.hpp"

#include "../manybody/CSF.hpp"
#include "../manybody/tools.hpp"


/**
 * \file QdotFci.hpp
 * \date 9-17-08
 * \author Simen Kvaal

 * \brief Defines a class QdotFci that sews together all concepts
 * to create a FCI solver for the quantum dot system in two dimensions.
 *
 */

namespace quantumdot {


/// \brief FCI solver for quantum dot system
///
/// QdotFci sews together the different classes in the namespace
/// quantumdot to create a complete solver.  The class collaboration
/// diagram shows how the different classes work together.
///
/// The class is used in the program "qdot".
///
/// Here is a brief overview of the roles each class play:
/// * QdotHilbertSpace: Encapsulates Hilbert space of Slater determinants and single-particle space,
///   i.e., orbitals and their structure.
/// * CsfMachine: Computes basis of CSFs.
/// * MatrixMachine: Takes QdotMatrixElements and a CSF basis and computes matrices
/// * QdotMatrixElements: Delivers matrix elements to MatrixMachine; has a QdotInteraction
/// * QdotInteraction: Delivers lab frame matrix elements to QdotMatrixElements; has a RadialPotential
/// * RadialPotential: Delivers COM frame matrix elements to QdotInteraction.
///
/// That's about it ...
///
/// For diagonalization, QdotFci uses ARPACK++. It supplies a function applyHamiltonian() to
/// compute phi' = H*phi, which is needed for the Arnoldi algorithm (Laczos, actually, since H is symmetric.)
///
/// After diagonalization, eigenvalues and eigenvectors may be saved to a Matlab script
/// for further processing.
///
/// This class also supports a magnetic field, see setMagneticField(). This is equivalent to
/// chaning the energy scale and the interaction strength lambda, in addition to a splitting
/// in the eigenvalues, proportional to the spin projection Sz and the magnetic field.
/// This is also output when saveEigenvalues() is called.
///
  class QdotFci
  {
    private:
      int R; ///< Largest shell index
      int A; ///< Number of particles
      int M; ///< Total angular momentum
      int Sz;///< Total spin projection (*2)
      int S; ///< Total spin (*2)

      double lambda; ///< Interaction strength.

      double B;      ///< Magnetic field
      double g_factor; ///< Dimless effective g-factor. g_factor := (g^*)*(m^*)/m_e.

      QdotHilbertSpace hilbert;        ///< Hilbert space / basis of Slater determinants.
      QdotInteraction inter;           ///< Interactions are ocmputed with this
      manybody::CsfMachine states;     ///< Creates CSFs
      QdotMatrixElements mat;          ///< Matrix element computer
      RadialPotential potential;       ///< Potential function interaction.

      bool use_energy_cut;             ///< True if energy cut space is desired over direct product space.
      bool use_veff;                   ///< True if we want effective interaction instead of bare interaction
      bool matrix_messages;            ///< True if messages during matrix build is desired.
      bool precompute_lab_frame;       ///< True if we want to precompute lab frame interaction

      dense_vector eigenvalues;        ///< Holds the eigenvalues computed in diagonalize()
      dense_matrix eigenvectors;       ///< Holds the eigenvectors computed in diagonalize()
    
      simple_sparse::SparseMatrixCrs<double> H0; ///< Holds the one-body Hamiltonian
      simple_sparse::SparseMatrixCrs<double> H1; ///< Holds the two-body Hamiltonian

      std::string filename;            ///< Filename for matlab output.
      std::ofstream matlab;            ///< Stream for matlab output.

    public:
      /// \brief Default constructor.
      QdotFci()
      {
	filename = "";
      }
      /// \brief Open MATLAB output file.
      /// \param filename  File name of MATLAB script.
      void setMatlabOutput(const std::string filename)
      {
	matlab.open(filename.c_str());
      }
      /// \brief Set number of particles
      /// \param the_A   The number of particles, must be greater than 0.
      void setParticles(int the_A) { A = the_A; assert(A > 0); }
      /// \brief Set spin projection values Sz and S.
      /// \param the_Sz   2 x the value for Sz; must be even (odd) if A is even (odd)
      /// \param the_S    2 x the value for S; must be even (odd) if A is even (odd)
      void setSpinValues(int the_Sz, int the_S)
      {
	S = the_S; 
	assert(S <= A);
	assert(S >= 0);
	Sz = the_Sz;
	assert((abs(Sz) % 2) == (S % 2));
	assert(abs(Sz) <= S);
      }
      /// \brief Set total angular momentum
      /// \param the_M The total angular momentum, an integer
      void setAngularMomentum(int the_M)
      {
	M = the_M;
      }
      /// \brief Set max shell index R. Number of shells = R+1.
      /// \param the_R  The thell index.
      void setMaxShell(int the_R)
      {
	R = the_R;
	assert(R>=0);
      }
      /// \brief Set interaction strength lambda
      /// \param the_lambda The interaction srength
      void setLambda(double the_lambda)
      {
	lambda = the_lambda;
      }
      /// \brief Set magnetic field B and effective g-factor.
      /// \param the_B     The magnetic field
      /// \param the_g_factor   The g-factor
      void setMagneticField(double the_B, double the_g_factor = 0.0)
      {
	B = the_B;
	g_factor = the_g_factor;
      }

      /// \brief Set whether we want energy cut space or direct product space
      void setUseEnergyCut(bool use)  {  use_energy_cut = use;  }

      /// \brief Set whether we want effective interaction or bare interaction
      void setUseVeff(bool use)  {  use_veff = use;  }

      /// \brief Turn on/off messages for matrix progress
      void setMatrixMessages(bool onoff) { matrix_messages = onoff; }

      /// \brief Turn on/off precompute lab frame
      void setPrecomputeLabFrame(bool onoff) { precompute_lab_frame = onoff; }

      /// \brief Attach a radial potential object.
      void setRadialPotential(const RadialPotential& u) { potential = u; }
      /// \brief Build complete Hamiltonian matrix.
      void buildHamiltonian();

      /// \brief Save Hamiltonian to a MATLAB script, to std::cerr otherwise
      void saveHamiltonian()
      {
	if (matlab.is_open())
	{
	  simple_sparse::matlab_dump(H0, matlab, "H0");
	  simple_sparse::matlab_dump(H1, matlab, "H1");
	}
	else
	{
	  simple_sparse::matlab_dump(H0, std::cerr, "H0");
	  simple_sparse::matlab_dump(H1, std::cerr, "H1");
	}
      }
      /// Save eigenvalues to matlab script
      void saveEigenvalues()
      {
	using namespace std;

	vector<double> shifts;
	for (int Sz2 = -S; Sz2 <= S; Sz2 += 2)
	  shifts.push_back(Sz2*g_factor*B);

	if (matlab.is_open())
	  matlab << "E_spin = [ ";
	else
	  cerr << "E_spin = [ ";
	for (size_t k = 0; k<shifts.size(); ++k)
	{
	  if (matlab.is_open())
	    matlab << format(shifts[k]) << " ";
	  else
	    cerr << format(shifts[k]) << " ";
	}
	if (matlab.is_open())
	  matlab << "];" << endl;
	else
	  cerr << "];" << endl;

	if (matlab.is_open())
	  matlab << "E = zeros(" << eigenvalues.length() << ", 1);" << endl;
	for (int k = 1; k<=eigenvalues.length(); ++k)
	{
	  if (matlab.is_open())
	    matlab << "E("<<k<<") = " << format(eigenvalues(k)) << ";" << endl;
	  else
	    cerr << "E("<<k<<") = " << format(eigenvalues(k)) << ";" << endl;
	
	}
      }
      /// \brief Save eigenvectors to matlab script, to std::cerr otherwise
      void saveEigenvectors(int number)
      {
	using namespace std;

	int j, k;

	int N = min(number, (int)eigenvectors.numCols());

	if (matlab.is_open())
	  matlab << "U = zeros(" << eigenvectors.numRows() << ", " << N << ");" << endl;
	else
	  cerr << "U = zeros(" << eigenvectors.numRows() << ", " << N << ");" << endl;

	for (j = 1; j<=eigenvectors.numRows(); ++j)
	{
	  if (matlab.is_open())
	    matlab << "U("<<j<<",:) = [";
	  else
	    cerr << "U("<<j<<",:) = [";

	  for (k = 1; k<=N; ++k)
	  {
	    if (matlab.is_open())
	    {
	      matlab << format(eigenvectors(j,k));
	      if (k<N) 
		matlab << ", ";
	      else
		matlab << "];" << endl;
	    }
	    else
	    {
	      cerr << format(eigenvectors(j,k));
	      if (k<N) 
		cerr << ", ";
	      else
		cerr << "];" << endl;
	    }
	  }
	}
      }

      /// \brief Find a few eigenvalues and eigenvectors of the Hamiltonian.
      /// \param nev Number of eigenvalues to compute.
      void diagonalize(int nev);

      /// \brief Post process eigenvalues to accomodate B-field.
      void postProcess();


    protected:
      /// \brief Apply Hamiltonian H = H0 + H1 to vector u, and store in v.
      ///
      /// This function is called by ARPACK.
      ///
      /// \param u Pointer to first element of u.
      /// \param v Pointer to first element of v.
      void applyHamiltonian(double *u, double *v)
      {
	size_t k;
	const size_t dim = H1.getNrows();

	// Initialize result vector v.
	for (k=0; k<dim; ++k)
	  v[k] = 0;

	// Result: v = H1_opt * u
	simple_sparse::matrix_vector_product_noinit(H1, u, v);

	// Result: v = (H0_opt + lambda * H1_opt)*u
	simple_sparse::matrix_vector_product_noinit(H0, u, v);


      }


  };

} // namespace quantumdot


#endif // _QDOT_FCI_HPP_
