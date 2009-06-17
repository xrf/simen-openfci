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

#include "QdotInteraction.hpp"
#include "RadialPotential.hpp"
#include "tools.hpp"

#include <iostream>

using namespace std;
using namespace quantumdot;

/**
 * \file radial.cc
 * \date 10-7-2008
 * \author Simen Kvaal
 *
 * \brief Simple program to solve the radial equation for the two-body problem.
 * 
 * Play around with it, test different potentials and stuff (see RadialPotential).
 *
 *
 **/

int main()
{
  
  print_copyright(cerr);

  // Get some data from user
  
  int N = 10;
  int m = 0;
  cerr << "Enter N (matrix dim): ";
  cin >> N;
  cerr << "Enter m (angular momentum): ";
  cin >> m;

//   cerr << "Enter omega: ";
//   double omega;
//   cin >> omega;

  // Declare matrices.
  dense_matrix T; // HO matrix
  dense_matrix V, Veff; // Potential and effective potential
  
  // Compute HO matrix.
  T.resize(N+1,N+1);
  T = 0;
  for (int k = 1; k <=N+1; ++k)
    T(k,k) = 2*(k-1) + abs(m) + 1;

  // Declare a potential; it is a pure Coulomb potential by default.
  RadialPotential coulomb;
  coulomb.setLambda(1.0);

  // Set integration method.
  coulomb.setIntegrationMethod(RadialPotential::GENERALIZED_HERMITE_QUAD);

  // Compute potential matrix and effective potential matrix.
  coulomb.computeMatrix(V, N, m);
  coulomb.computeEffectiveMatrix(Veff, N, m);


  // Write to Matlab script on standard out.

  cout.precision(16);
  cout << scientific;
  dense_matrix_dump(T, cout, "H0");
  dense_matrix_dump(V, cout, "V");
  dense_matrix_dump(Veff, cout, "Veff");

 
  
}
