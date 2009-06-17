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

#include "tools.hpp"
#include "Slater.hpp"
#include <iostream>
#include <cassert>
#include <vector>

using namespace manybody;
using namespace std;


/**
 *
 * \file slater_demo.cc
 * \author Simen Kvaal
 * \date 10-7-2008
 *
 * \brief Program demonstratning basic functionality of the Slater class.
 *
 * Interactively, two Slater determinants are built by using Slater::create() and Slater::annihilate().
 * Then, the program analyzes the differing occupied orbitals. This procedure is at the heart of
 * MatrixMachine, where matrix elements of arbitrary operators are implemented.
 *
 * Play around with it...
 *
 * 
 **/




int main()
{


  int N = 20;
  assert(N <= SLATER_BITS);
  vector<Slater> psi(2);

  print_copyright(cerr);

  cerr << "Interactive demo of slater determinants." << endl;
  cerr << endl;
  cerr << "I will now ask you to enter pairs on the form 'c <alpha>' or 'd <alpha>'." << endl;
  cerr << "The former instructs me to create a particle in orbital alpha, the latter to destroy." << endl;
  cerr << "For example, if you enter 'c 3' (without the quotes), I will create a particle in alpha=3." << endl;
  cerr << "Enter 'x' to finish." << endl;
  cerr << endl;


  for (size_t k = 0; k < 2; ++k)
  {
    Slater& phi = psi[k];
    
    cerr << "Let's build Slater determinant number " << k+1 << " !" << endl;

    cerr << "Initial vacuum: 1 * | " <<  phi.to_binary(N) << " > " << endl;
    
    int sign = 1;
    
    while (true)
    {
      size_t alpha;
      string s;
      cerr << "Enter a pair: ";
      cin >> s;
      if (s == "x")
	break;
      cin >> alpha;
      
      if (s == "c")
      {
	cerr << "Operating with a+_" << alpha << ": " << endl;
	sign *= phi.create(alpha);
      }
      else if (s == "d")
      {
	cerr << "Operating with a_" << alpha << ": " << endl;
	sign *= phi.annihilate(alpha);
      }
      else
      {
	cerr << "I didn't understand that, sorry. :-(" << endl;
	break;
      }
      if (sign == 0)
      {
	cerr << "Zero vector!" << endl;
	break;
      }
      cerr << "Current Slater determinant: " << sign << " * | " << phi.to_binary(N) << " > "<< endl;
    }
   

 
  }

  cerr << endl << endl;
  cerr << "Let's analyze the difference between the Slater determinants." << endl;
  cerr << "Slater number 1 has "<< psi[0].count() << " particles. " << endl;
  cerr << "Slater number 2 has "<< psi[1].count() << " particles. " << endl;

  int max;
  cerr << "Enter max number of different occupied orbitals to search for: ";
  cin >> max;

  size_t alpha[max], beta[max], p[max], q[max];
  
  int d = psi[0].analyze_difference(psi[1], alpha, beta, p, q, max);

  if (d > max)
    cerr << "I found MORE than " << max << " different bits." << endl;
  else
  {
    cerr << "I found " << d << " different occupied orbitals." << endl;
    for (size_t k = 0; k < d; ++k)
    {
      cerr << "alpha[" << k << "] = " << alpha[k] << ", p[" << k << "] = " << p[k] << endl;
      cerr << "beta[" << k << "]  = " << beta[k]  << ", q[" << k << "] = " << q[k] << endl;
    }

  }
}
