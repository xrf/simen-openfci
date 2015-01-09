#include <iostream>
#include <vector>

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

#include "NChooseKBitset.hpp"
#include "ConfigFile.h"
#include "CSF.hpp"
#include "MatrixMachine.hpp"


/**
 *
 * \file pairing.cc
 * \author Simen Kvaal
 * \date 10-7-2008
 *
 * \brief Program demonstratning basic functionality of the core functions in namespace manybody by solving the simple pair model.
 *
 * This small program builds the Hamiltonian of a simple pairing model. The intention
 * is to demonstrate the features used in the quantum dot application "qdot" but without
 * the need for the rather complex framework specific to that model.
 *
 * Consider a single-particle space consisting of L levels where spin-1/2 fermions "live".
 * Each level is equally spaced from each other -- in fact it can be viewed as a 1D harmonic
 * oscillator truncated at L shells. Anyway, the Hamiltonian has the form
 * 
 * \f[ H = \sum_{j,\sigma} j a^\dag_{(j,\sigma)} a_{j,\sigma} + g
 *     \sum_{j,k} a^\dag_{j,+} a^\dag_{k,-} a_{k,-} a_{j,+} + f
 *     \sum_{j,k,\ell} a^\dag_{j,+} a^\dag_{k,-} a_{\ell,-} a_{\ell,+}
 *     + \mathrm{h.c}, \f]
 *
 * and it is easily seen that the interaction conserves S^2 and
 * S_z. The interaction is split into two parts: the part proportional
 * to g is called the pair interaction, and that only \em pairs of
 * particles <em>at the same level</em> interact. The part proportional to f is called
 * a particle-hole interaction, and it breaks pairs.
 *
 *
 * The program constructs Hilbert space blocks of constant Sz and/or constant S^2.
 * Different settings can be played with. The program does no diagonalization, but writes
 * the matrices to a Matlab script.
 *
 * The program reads a configuration file (see ConfigFile) like this one:
 *
 * \code
 * g = 0.5        # pair interaction strength
 * f = 0.001      # particle-hole interaction strength
 * levels = 8     # number of levels in model
 * P_levels = 3   # number of levels in model space
 * particles = 4  # number of particles
 * Sz = 0         # spin projection * 2
 * S = 0          # total spin * 2
 *
 * Sz_projection = yes   # project onto Sz subspace?
 * S2_projection = yes   # project onto S^2 subspace?
 * \endcode
 **/



using namespace std;
using namespace manybody;
using namespace simple_sparse;

/// \brief Class for a simple pairing model.
class Pairing
{
  private:
    int particles;            ///< Number of particles.
    int Sz;                   ///< Spin projection.
    int S;                    ///< Total spin
    bool Sz_projection;       ///< Project onto Sz subspace?
    bool S2_projection;       ///< Project onto S^2 subspace?
    int levels;               ///< Number of levels in model
    int P_levels;         ///< Number of levels in P-space.
    vector<Slater> sd_basis;  ///< Basis of Slater determinants.
    vector<csf_block> csf_basis; ///< CSF functions based on Slater determinants to be used as basis.

    
    double g;                 ///< Coupling constant og pair interaction U
    SparseMatrixCrs<double> Jz; ///< Spin projection operator
    SparseMatrixCrs<double> H0; ///< One-body Hamiltonian
    SparseMatrixCrs<double> U;  ///< Pair interaction matrix
    SparseMatrixCrs<double> V;  ///< Particle-hole interaction matrix
    
  public:
    /// \brief Build space of Slater determinants.
    void buildHilbertSpace();

    /// \brief Find the P-space and write it to cout.
    void findPSpace();

    /// \brief Build Hamiltonian matrices H0 and H1.
    void buildMatrices();

    /// \brief Set number of levels in model.
    /// \param ell The number of levels
    void setLevels(int ell)
    {
      levels = ell;
      assert(2*(size_t)levels <= SLATER_BITS);
    }
    /// \brief Set number of levels in P-space.
    /// \param ell The number of levels
    void setPLevels(int ell)
    {
      P_levels = ell;
    }

    /// \brief Set spin parameters.
    /// \param sz  The spin projection value
    /// \param s   The eigenvalue for S^2
    void setSpinParameters(int sz, int s)
    {
      Sz = sz;
      S = s;
    }
    /// \brief Select subspace to project onto: Sz, S^2 or both or none 
    /// (i.e., use all of Hilbert space)!
    /// \param Sz  true if Sz-subspace is to be chosen
    /// \param S2  true if S^2-subspace is to be chose
    void setProjections(bool Sz, bool S2)
    {
      Sz_projection = Sz;
      S2_projection = S2;
    }

    /// \brief Set number of particles
    /// \param n Number of particles
    void setParticles(int n)
    {
      particles = n;
    }

    /// \brief Return true if the argument is in Sz-subspace.
    bool constraint(const bitset<SLATER_BITS>& x)
    {
      // Count number of even bits (up spins) and odd bits (down spins).
      // Spin projection = up - down.

      int up, down;
      Slater phi;
      phi = x;
      phi.count_evenodd(up, down);

      return up - down == Sz;
    }


};

/// \brief Matrix elements of pairing interaction
class PairModelMatrixElements
{
  public:
    /// To select which matrix to compute elements for.
    enum Which { ONE_BODY, PAIR_INTERACTION, PARTICLE_HOLE_INTERACTION, JZ };

  private:
    Which which;    ///< Which matrix shall we compute?

  public:

    /// Select which matrix to compute.
    /// \param w See PairModelMatrixElements::Which
    void setWhich(Which w) { which = w; }

    double element(size_t* a, size_t* b)
    {
      // Compute quantum numbers: pi, qi = level, si, ti = spin; 0 = up, 1 = down
      const int p0 = a[0]/2; const int p1 = a[1]/2;
      const int s0 = a[0]%2; const int s1 = a[1]%2;

      const int q0 = b[0]/2; const int q1 = b[1]/2;
      const int t0 = b[0]%2; const int t1 = b[1]%2;

      double elm = 0.0; // for temp computations

      switch (which) // which (switch) ? :-)
      {
	case JZ: // Jz matrix elements
	  return ((p0 == q0) && (s0 == t0)) ? 1.0 - 2.0*s0 : 0.0;
	  break;

	case ONE_BODY: // One-body Hamiltonian
	  return ((p0 == q0) && (s0 == t0)) ? (double)p0 : 0.0;
	  break;

	case PAIR_INTERACTION: // pair interaction
	  if (  (p0 == p1) && (q0 == q1) 
		&& (s0 == 0) && (s1 == 1) 
		&& (t0 == 0) && (t1 == 1)  )
	    return 1.0;
	  else
	    return 0.0;
	  break;

	case PARTICLE_HOLE_INTERACTION: // particle-hole interaction
	  elm = 0.0;
	  if (  (p0 == p1)  
		&& (s0 == 0) && (s1 == 1) 
		&& (t0 == 0) && (t1 == 1)  )
	    elm += 1.0;
	  if (  (q0 == q1)  
		&& (s0 == 0) && (s1 == 1) 
		&& (t0 == 0) && (t1 == 1)  )
	    elm += 1.0;

	  return elm;
	  break;
	default:
	  return 0.0;

      } // switch (which)
    }

};

void Pairing::buildMatrices()
{
  // Create a matrix element object and a matrix builder
  PairModelMatrixElements mat_elems;
  MatrixMachine<PairModelMatrixElements> mat_machine(mat_elems);


  // Insert basis & CSF blocks into matrix machine object
  mat_machine.setSDBasis(sd_basis);
  mat_machine.setCSFBlocks(csf_basis);
  mat_machine.setSymmetryAssumption(false);

  // Build Jz matrix.
  cerr << "Building Jz ..." << endl;
  mat_machine.setRank(1);
  mat_elems.setWhich(PairModelMatrixElements::JZ);
  mat_machine.buildMatrix(Jz, true);

  // Build H0 matrix
  cerr << "Building H0 ..." << endl;
  mat_elems.setWhich(PairModelMatrixElements::ONE_BODY);
  mat_machine.buildMatrix(H0, true);

  // Build pair interaction
  cerr << "Building U ..." << endl;
  mat_machine.setRank(2);
  mat_elems.setWhich(PairModelMatrixElements::PAIR_INTERACTION);
  mat_machine.buildMatrix(U);

  // Build particle-hole interaction
  cerr << "Building V ..." << endl;
  mat_elems.setWhich(PairModelMatrixElements::PARTICLE_HOLE_INTERACTION);
  mat_machine.buildMatrix(V);

  matlab_dump(Jz, cout, "Jz");
  matlab_dump(H0, cout, "H0");
  matlab_dump(U, cout, "U");
  matlab_dump(V, cout, "V");

}

void Pairing::buildHilbertSpace()
{
  // Erase current basis.
  sd_basis.clear();

  cerr << "Building Hilbert space... ";
  if (Sz_projection)
    cerr << "(for a single value of Sz, i.e., Sz = " << Sz << ")" << endl;
  else
    cerr << "(all Slater determinants)" << endl;

  // Get a basis-building object.
  NChooseKBitsetExternal<Pairing,SLATER_BITS> generator(*this);
  generator.setUseConstraint(Sz_projection);
  generator.N = 2*levels;
  generator.K = particles;
  generator.choose();

  // Copy the basis.
  size_t k;
  cerr << "Here is the Slater determinant basis: " << endl;
  for (k=0; k<generator.subset.size(); ++k)
  {
    Slater phi;
    phi = generator.subset[k];    sd_basis.push_back(phi);
  }


  // Create CSF basis
  if (!S2_projection)
    computeTrivialCSFBlocks(sd_basis, csf_basis);
  else
  {
    computeCSFBlocks(sd_basis, csf_basis, S);

    for (size_t j = 0; j < csf_basis.size(); ++j)
      cerr << "CSF block " << j << ": " << endl << csf_basis[j].coeffs << endl;
  }

  // Copy the basis, also print it.
  cerr << "Here is the Slater determinant basis (after possible rearrangement by CsfMachine): " << endl;
  for (k=0; k<generator.subset.size(); ++k)
  {
    cerr << "phi["<<k<<"]  =  | " << sd_basis[k].to_binary(2*levels) << " >" << endl;
  }


}


void Pairing::findPSpace()
{
  // Go through the CSF blocks and identify which are in P-space.
  size_t n_csf = csf_basis.size();
  vector<bool> in_P_space(n_csf);
  size_t dim = 0;

  for (size_t k = 0; k < n_csf; ++k)
  {
    // Increment dimension.
    dim += csf_basis[k].n_csf;
    // Take a Slater determinant in the block.
    Slater phi = sd_basis[ csf_basis[k].sd_begin ];
    // Check if a particle is above P-space level.
    int highest_particle = phi.last_index() / 2;
    in_P_space[k] = (highest_particle < P_levels);
    //cerr << "Is CSF block " << k << " in P-space ? " << yesno(in_P_space[k]) << endl;
  }

  cout << "P = zeros(" << dim << ", 1);" << endl;
  dim = 0;
  for (size_t k = 0; k < n_csf; ++k)
  {
    for (size_t j = 0; j < csf_basis[k].n_csf; ++j)
      cout << "P(" << (dim + j+1) << ") = " << in_P_space[k] << "; ";
    dim += csf_basis[k].n_csf;
  }
  cout << endl;

}


int main(int argc, char **argv)
{

  print_copyright(cerr);

  string configname = "pairing.ini";
  if (argc > 1)
    configname = argv[1];


  try
  {
    ConfigFile config(configname);
    config.read<double>("g", 0.25);
    config.read<double>("f", 0.25);
    int levels = config.read<int>("levels", 10);
    int P_levels = config.read<int>("P_levels", 3);
    int particles = config.read<int>("particles", 4);
    int Sz = config.read<int>("Sz", 0);
    int S = config.read<int>("S", 0);
    bool Sz_projection = config.read<bool>("Sz_projection", true);
    bool S2_projection = config.read<bool>("S2_projection", true);

    Pairing solver;
    solver.setLevels(levels);
    solver.setPLevels(P_levels);
    solver.setParticles(particles);
    solver.setSpinParameters(Sz,S);
    solver.setProjections(Sz_projection, S2_projection);
    solver.buildHilbertSpace();
    solver.findPSpace();
    solver.buildMatrices();


  }
  catch( ConfigFile::file_not_found& e ) {
    cerr << "The configuration '" << e.filename << "' was not found. Aborting.";
    cerr << endl << endl;
    return 1;
  }




}
