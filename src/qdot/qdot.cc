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

#include <iostream>

#include "ConfigFile.h"

#include "NChooseK.hpp"
#include "CSF.hpp"

#include "QdotFci.hpp"
#include "RadialPotential.hpp"

/**
 * \file qdot.cc
 * \date 10-7-2008
 * \author Simen Kvaal
 *
 * \brief Program to find eigenvalues and eigenvectors of 2D parabolic
 * quantum dots.
 *
 *
 * This program finds eigenpairs of the parabolic quantum dot
 * Hamiltonian.
 *
 * It reads a configuration file given as argument ('default.ini' is
 * default), and solves the corresponding problem. None of the
 * settings in the configuration file are mandatory -- sensible
 * defaults are always given.
 *
 *
 * The ConfigFile library by R. Wagner is used; see ConfigFile.hpp for
 * details.
 *
 * The configuration file is on the following form:
 *
 * \code
 * A = 2          # number of particles
 * S = 0          # total spin eigenvalue * 2
 * M = 0          # total angular momentum quantum number
 * R = 20         # model space size parameter
 *
 * lambda = 1.0   # interaction strength parameter
 *
 * # Use the below commands to set up a more general potential than the
 * # default Coulomb potential: U(r) = r^(-alpha)*p(r)*exp(-beta*r^2)
 *
 * # pot_p = 1.0 2.0 -0.2 (...)     # Polynomial; last coeff is constant coeff.
 * # pot_p_iseven = yes             # Poly = 1.0*r^4 + 2.0*r^2 - 0.2.
 * # pot_alpha = -1.0               
 * # pot_beta = 0.23
 *
 * # Additional support for magnetic field. Default is B = 0.
 * # Adding a B-field is treated by re-scaling the interaction strength
 * # lambda into lambda^*, and the harmonic oscillator energy scale.
 * # The eigenvalues are split according to spin projection Sz proportional
 * # to g_factor.
 * #
 * # One may safely ignore the B-field here, and use the default value, which
 * # doesn't change a thing.
 *
 * # B = 0                  # field strength.
 * # g_factor = 2.0*0.0067  # effective g-factor.
 *
 *
 *
 * nev = 20       # number of eigenpairs to find (lowest eigenvalues are searched for)
 *
 * use_energy_cut = yes    # energy cut model space or direct product model space
 * use_veff = no           # effective interaction? (only sensible if use_energy_cut = yes)
 * matrix_messages = yes   # gives some indication of matrix building progress
 * precompute_lab_frame = no   # yes: precompute interaction in lab frame instead of
 *                             # computing COM transformation "on the fly".
 *                             # smart if A>3, not smart if A<=3.
 *
 * matlab_output = script.m    # filename for output in Matlab script format
 *
 * # what to save in matlab script:
 *
 * save_eigenvalues = yes
 * save_eigenvectors = yes
 * save_eigenvectors_n = 2    # number of eigenvectors to save; they are often big!
 * save_matrices = no         # saves Hamiltonian if yes.
 *
 *
 * \endcode
 *
 **/



using namespace std;
using namespace quantumdot;

/// \brief Input stream operator overloading for std::vector<T>
template<class T>
istream& operator>>(istream& is, std::vector<T>& x)
{
  T val;

  x.resize(0);

  while (!is.eof())
  {
    is >> val;
    x.push_back(val);
  }
  return is;
  
}

template<class T>
ostream& operator<<(ostream& os, std::vector<T>& x)
{
  for (size_t k = 0; k < x.size(); ++k)
    os << x[k] << " ";
  return os;
}

/// \brief Set up radial potential from options in config.
/// \param potential Destination potential
/// \param config ConfigFile instant with configuration
void initPotential(RadialPotential& potential, const ConfigFile& config)
{

  double alpha = 0;
  std::vector<double> p, p_default;
  bool p_iseven;

  p_default.push_back(1.0);


  // Read parameters from config.
  potential.setAlpha( config.read<double>("pot_alpha", 1.0) );
  potential.setBeta( config.read<double>("pot_beta", 0.0) );
  potential.setP( config.read<vector<double> >("pot_p", p_default) );
  p_iseven = config.read<bool>("pot_p_iseven", true);
  potential.setPIsEven(p_iseven);

  // Set the basis type to Laguerre (Fock-Darwin)
  potential.setBasisType(RadialPotential::GENERALIZED_LAGUERRE_BASIS);
 
  // Set integration method to standard Hermite if applicable, otherwise
  // to generalized half-range Hermite quad.
  if ((alpha == 1.0) && (p_iseven))
    potential.setIntegrationMethod(RadialPotential::STANDARD_HERMITE_QUAD);
  else
    potential.setIntegrationMethod(RadialPotential::GENERALIZED_HERMITE_QUAD);

		     
}


int main(int argc, char **argv)
{

  // Get name of configuration file from command line.
  // If not given, 'default.ini' will be used.
  string configname = "default.ini";
  if (argc > 1)
  {
    configname = argv[1];
  }

  try
  {
    ConfigFile config(configname);
    int A = config.read<int>("A", 2);
    int R = config.read<int>("R", 20);
    int S = config.read<int>("S", A);
    int Sz = config.read<int>("Sz", S);
    int M = config.read<int>("M", 0);
    double B = config.read<double>("B", 0.0);
    double g_factor = config.read<double>("g_factor", 2.0*0.067);
    double lambda = config.read<double>("lambda", 1.0);
    bool use_energy_cut = config.read<bool>("use_energy_cut", true);
    bool matrix_messages = config.read<bool>("matrix_messages", true);
    
    bool precompute_lab_frame = config.read<bool>("precompute_lab_frame", false);
    int nev = config.read<int>("nev", 20);
    
    bool save_matrices = config.read<bool>("save_matrices", false);
    bool save_eigenvalues = config.read<bool>("save_eigenvalues", true);
    bool save_eigenvectors = config.read<bool>("save_eigenvectors", false);
    int save_eigenvectors_n = config.read<int>("save_eigenvectors_n", nev);
    //bool compute_block_sparsity = config.read<bool>("compute_block_sparsity", false);
    bool use_veff = config.read<bool>("use_veff", false);
    
    string filename = config.read<string>("matlab_output", "");
    
    RadialPotential potential;
    initPotential(potential, config);

    QdotFci solver;
  
    if (filename != "")
      solver.setMatlabOutput(filename);
  
    solver.setParticles(A);
    solver.setSpinValues(Sz,S);
    solver.setAngularMomentum(M);
    solver.setMaxShell(R);
    solver.setLambda(lambda);
    solver.setMagneticField(B, g_factor);
    solver.setUseEnergyCut(use_energy_cut);
    solver.setUseVeff(use_veff);
    solver.setMatrixMessages(matrix_messages);
    solver.setPrecomputeLabFrame(precompute_lab_frame);
    solver.setRadialPotential(potential);
  
  
    solver.buildHamiltonian();   
  
    if (save_matrices)
      solver.saveHamiltonian();
  
    solver.diagonalize(nev);
    solver.postProcess();

    cerr << "Saving results ..." << endl;
  
    if (save_eigenvalues)
      solver.saveEigenvalues();
  
    if (save_eigenvectors)
      solver.saveEigenvectors(save_eigenvectors_n);
  
    cerr << "Finished." << endl;

  }
  catch( ConfigFile::file_not_found& e ) {
    cerr << "The configuration '" << e.filename << "' was not found. Aborting.";
    cerr << endl << endl;
    return 1;
  }
  
  return 0;
}


