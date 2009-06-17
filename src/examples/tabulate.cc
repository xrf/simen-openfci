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
#include <iomanip>

#include "ConfigFile.h"

#include "QdotInteraction.hpp"

#include "tools.hpp"

using namespace std;
using namespace quantumdot;


/**
 * \file tabulate.cc
 * \date 10-7-2008
 * \author Simen Kvaal
 *
 * \brief Program for tabulation of interaction matrix elements.
 *
 * This small program is provided mostly as an example of how one may use
 * the OpenFCI library. The matrix elements tabulated are the <a b | U | c d>
 * matrix elements, and these could be input for other methods, e.g., coupled cluster
 * methods. Indeed, such a method could be implemented using parts of the OpenFCI library.
 * 
 * The program reads a configuration file (see ConfigFile) on the following form:
 * \code
 * # Tabulate Coulomb matrix elements for SP energies 2*n + |m| <= R.
 * R = 20
 * lambda = 1.00
 * use_veff = true
 * \endcode
 *
 *
 **/




/// \brief Tabulate Coulomb matrix elements for R+1 shells.
/// 
/// Output is written to std::cout, one line at a time,
/// each line has one matrix element. The line has the format:
/// \verbatim n1 m1 n2 m2 n3 m3 n4 m4 element  \endverbatim
/// All values 2nj+|mj| <= R are processed, but only nonzeroes are written.
///
/// \param lambda  The interaction strength
/// \param R   The value for R
/// \param use_veff   Set to true if effective interaction is to be tabulated instead of bare interaction.
void tabulateCoulomb(double lambda, int R, bool use_veff);


void tabulateCoulomb(double lambda, int R, bool use_veff)
{
  QdotInteraction inter;
  cerr << "Building COM to LAB transformation ..." << endl;
  inter.setR(R);
  inter.setLambda(lambda);
  cerr << "Computing COM blocks of interaction ..." << endl;
  if (use_veff)
    inter.buildEffectiveInteractionComBlocks(1);
  else
    inter.buildInteractionComBlocks();

  // Tabulate values!

  int N1, N2, N3, N4;
  int n1, n2, n3, n4;
  int m1, m2, m3, m4;
  double value;

  const int w = 3;

  for (N1 = 0; N1 <=R; ++N1)
    for (m1 = -N1; m1 <= N1; m1 += 2)
      for (N2 = 0; N2 <=R-N1; ++N2)
	for (m2 = -N2; m2 <= N2; m2 += 2)

	  for (N3 = N1; N3 <=R; ++N3)
	    for (m3 = (N3==N1 ? m1 : -N3); m3 <= N3; m3 += 2)
	      for (N4 = N2; N4 <=R-N3; ++N4)
		for (m4 = (N4==N2 ? m2 : -N4); m4 <= N4; m4 += 2)
		{
		  n1 = (N1 - abs(m1)) / 2;
		  n2 = (N2 - abs(m2)) / 2;
		  n3 = (N3 - abs(m3)) / 2;
		  n4 = (N4 - abs(m4)) / 2;
		  value = inter.singleElement(N1,m1, N2,m2, N3,m3, N4,m4);
		  if (value != 0.0)
		  {
		    cout << setw(w) << n1 << " " << setw(w) << m1 << " " << setw(w) << n2 << " " << setw(w) << m2 << " " 
			 << setw(w) << n3 << " " << setw(w) << m3 << " " << setw(w) << n4 << " " << setw(w) << m4 << " " 
			 << "    " << format(value) << endl;
		  }
		}

}





int main(int argc, char **argv)
{

  print_copyright(cerr);

  string configname = "tabulate.ini";
  if (argc > 1)
    configname = argv[1];


  try
  {
    ConfigFile config(configname);
    double lambda = config.read<double>("lambda", 1.0);
    int R = config.read<int>("R", 10);
    bool use_veff = config.read<bool>("use_veff", false);

    tabulateCoulomb(lambda, R, use_veff);


  }
  catch( ConfigFile::file_not_found& e ) {
    cerr << "The configuration '" << e.filename << "' was not found. Aborting.";
    cerr << endl << endl;
    return 1;
  }


}

