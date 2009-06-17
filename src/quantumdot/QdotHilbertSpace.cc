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

#include "QdotHilbertSpace.hpp"

#include "timer.hpp"
#include "NChooseK.hpp"
#include "NChooseKBitset.hpp"
#include "tools.hpp"

namespace quantumdot {

  using namespace manybody;

  void QdotHilbertSpace::buildSPBasis()
  {
    qnumbers.resize(0);

    // Loop through all possible quantum numbers and
    // Add to vector.

    // Note that alpha *even* is a spin *up*.

    int N, m, sigma, alpha = 0, n;
    for (N = 0; N<=R; ++N)
      for (m = -N; m<=N; m+=2)
	for (sigma = 0; sigma<=1; ++sigma)
	{
	  alpha = (N*(N+2)+m) + sigma;
	  n = (N - abs(m)) >> 1;

	  quantum_numbers x;
	  x.N = N;
	  x.n = n;
	  x.m = m;
	  x.sigma = 1 - sigma;
	  x.alpha = alpha;

	  qnumbers.push_back(x);

	  // Make sure that alpha has the correct value.
	  assert(alpha == (int)qnumbers.size() - 1);

	}


    n_orbitals = alpha + 1;

    cerr << "Number of orbitals in single particle basis = " << n_orbitals << endl;

    assert(n_orbitals <= SLATER_BITS);
  }



  void QdotHilbertSpace::generateSDBasisBruteForce()
  {

    // first, generate all bit patterns.

    cerr << "Generating binom("<<n_orbitals<<","<<A<<") = " << round(exp(lnbinom(n_orbitals,A))) << " slater determinants." << endl;

    NChooseKBitsetExternal<QdotHilbertSpace, SLATER_BITS> patterns(*this);
    patterns.setUseConstraint(true);
    patterns.N = n_orbitals;
    patterns.K = A;
    patterns.choose();

    sd_basis.resize(0);
    for (size_t k = 0; k < patterns.subset.size(); ++k)
    {
      Slater phi;
      phi = patterns.subset[k];
      sd_basis.push_back(phi);
    }

    cerr << "I found " <<sd_basis.size() << " states." << endl;

  }



  void QdotHilbertSpace::generateSDBasis()
  {

    cerr << "Generating Slater determinants ..." << endl;

    timer t, t2;

    // Count discarded states.
    int discarded = 0;
    int j, k; // counters

    vector<int> m(A); // angmoms for particles.
    int mmax = R; // maximum value for angmom. 
    int msum; // always the sum of the A-1 first angoms.
    bool angmom_finished = false; // indicates if all angular moms have been stepped through
    //in most cases, it actually is less.

    int dim = 0; // counts dimension of created space

    // reset space.
    sd_basis.resize(0);

    NChooseK<unsigned int> spin_chooser;
    spin_chooser.N = A; // # spins in total
    spin_chooser.K = (A + Sz)/2; // # spins up
    spin_chooser.choose(); // generate spin configurations.


    // initialize angmoms
    msum = 0;
    for (j=0; j<A-1; ++j)
    {
      m[j] = -mmax;
      msum += m[j];
    }
    m[A-1] = M - msum;

    t2.start();

    // Do a loop that covers all possible
    // distributions of the m[j].
    while (!angmom_finished)
    {

      // increment angmoms until we have a legal value
      bool non_incr = false;
      bool bad_combo = false;
      while ( (abs(m[A-1]) > mmax) || bad_combo)
      {
	vector<int> temp;
	temp = m;
	temp.resize(A-1);
	angmom_finished = general_inc(temp, -mmax, mmax+1);
	m = temp;
	m.resize(A);
	msum = 0;
	for (j=0; j<A-1; ++j)
	  msum += m[j];
	m[A-1] = M - msum;

	if (angmom_finished)
	  break;

	// see if m is increasing. must be increasing
	// so we do not double-count a configuration.
	bad_combo = false;
	for (j = 0; j<A-1; ++j)
	  if (m[j]>m[j+1])
	  {
	    bad_combo = true;
	    break;
	  }

	if (use_energy_cut)
	{
	  int sum = 0;
	  for (j = 0; j<A-1; ++j)
	    sum += abs(m[j]);

	  // see if the minimum energy taken by this angmom config
	  // is too large for the energy cut subspace. this is 
	  // very likely, in general.
	  if (sum > R)
	    bad_combo = true;
		  

	}
     
      }
    
      // If we have tried all possible combinations, break off.
      if (angmom_finished)
	break;

      // create a vector containing partitioning
      // of spin configuration into m-sections.
      vector<int> part;
      part.push_back(0);
      part.push_back(1);
      for (j=1; j<A; ++j)
	if (m[j] == m[j-1])
	{
	  part[part.size()-1]++;
	}
	else
	{
	  part.push_back(part[part.size()-1] + 1);
	}
    
      // now, for example part = {1, 1, 3, 5}, 

      // Create all QUNIQUE spin configurations, by
      // pruning away all that produce duplicated, i.e.,
      // "permutations within one m-value"
      spin_chooser.choose();
    
      for (j = 0; j<(int)part.size()-1; ++j)
      {
	vector<unsigned int>& subset = spin_chooser.subset;
	// remove all spin configs that differ from 
	// unique by a permutation within the j'th m-value
	// section.
	unsigned int mask0 = (1U << A)-1;
	unsigned int mask1 = (1U << part[j+1]) - 1;
	unsigned int mask2 = (1U << part[j]) - 1;
	unsigned int mask = (mask1 ^ mask2) ^ mask0;
	size_t p = 1, q = 0;
	while (q < subset.size())
	{
	  p = q+1;
	  unsigned int unique = subset[q];
	  while (p < subset.size())
	  {
	    if ((subset[p] & mask) == (unique & mask))
	    {
	      subset.erase(subset.begin() + p);
	    }
	    else
	      p++;
	  }
	  q++;
	}
      }

      // we now have an ok set of angmoms.
      // process the configuration.

      // Loop though all possible *spin (projection)* configurations (computed above)
      // and create all possible slater determinants for the
      // current m-config and the spin config.
      for (size_t s_index = 0; s_index < spin_chooser.subset.size(); ++s_index)
      {
     
	// spin configuration:
	unsigned int spin_conf = spin_chooser.subset[s_index];
	
     

	// We now define a "column" as an (m, sigma) pair; sigma = 0,1. Every
	// particle resides in a column. Each column is independent
	// of the other, so the total set of particle configurations
	// becomes the *direct product* of all the possible sub-configurations
	// defined within the columns.
	// The number of *orbitals* in one column is floor((R-|m|)/2) + 1.
	// The *number of columns* is 2*(2*R+1).

	// We now count the particles in each column.
	// index for column (m,sigma) is 2*(R+m) + sigma. (No absolute
	// value on m there.) We also compute the minimum shell energy Nmin
	// of the total configuration.
	const int col_dim = 2*(2*R + 1);
	vector<int> column_count(col_dim);
	fill(column_count.begin(), column_count.end(), 0);
      
	int Nmin = 0;
	bool fit = true;
	for (j=0; j<A; ++j)
	{
	  int sigma = ((spin_conf >> j) & 1);
	  int i = 2*(R+m[j]) + sigma; // col index
	  column_count[i]++; // add particle to column
	  Nmin += abs(m[j]); // update smallest energy for a slater determinant.
	  // check if the state "fits"; that there are not more
	  // particles in a column than it can hold.
	  if (column_count[i] > ((int)R-abs(m[j]))/2 + 1)
	  {
	    fit = false;
	    break;
	  }
	}
      
      
	// Everything fits (at least one particle config in each column).
	// Generate direct products of column-configs.
	if (fit)
	{
	
	  typedef unsigned long long int the_int; // representing a column config.
	  assert(sizeof(the_int)*CHAR_BIT >= (size_t)R);

	  NChooseK<the_int> helper; // builds column states.
	  vector<vector<the_int> > states(col_dim); // a vector of column state sets.
	  //vector<int> levels;
	
	  // Generate column configurations: result is states[j] for each j.
	  for (j = 0; j<col_dim; ++j)
	  {
	    int sigma = j%2;
	    int mj = (j - sigma)/2 - R;
	    helper.N = (R - abs(mj))/2 + 1; // number of levels
	    helper.K = column_count[j];
	    if (helper.K)
	    {
	      helper.choose();
	      states[j] = helper.subset;
	    }
	  }

	  // Initialize multi-index thing:
	  vector<int> imin(col_dim); // minimum value for index
	  vector<int> imax(col_dim); // maximum value for index
	  vector<int> ivec(col_dim); // index
	  for (k = 0; k<(int)col_dim; ++k)
	  {
	    imin[k] = 0;
	    imax[k] = max(states[k].size(), (size_t)1);
	    ivec[k] = 0;
	  }
	
	  // Loop through all multi-indices and build DP states.
	  bool finished = false;
	  while (!finished)
	  {

	    // build complete state phi.
	    Slater phi; // the Slater determinant
	    int Nsum = 0; // total HO energy of state
	    // loop over each column...
	    for (k=0; k<(int)col_dim; ++k)
	    {

	      // if there aren't any configs in the column, skip
	      if (states[k].size() > 0)
	      {
		// Compute sigma and m for the current column...
		int sigma = k % 2;
		int mj = (k - sigma)/2 - R;
	      
		// fetch bit pattern for sub-state
		the_int substate = states[k][ivec[k]];
	      
	      
		// map bits into a A-body state
		int levels = (R - abs(mj))/2 + 1;
		for (j=0; j<levels; ++j)
		{
		  int N2 = 2*j + abs(mj); // the energy of the site
		  int alpha = (N2*(N2+2) + mj) + (1-sigma); // the orbital index
		  //if (bit(substate,j))
		  if ((substate >> j) & 1) // if bit is set, particle is there.
		  {
		    Nsum += N2;
		    if ((Nsum > R) && (use_energy_cut))
		    {
		      discarded++;
		      break;
		    }
		    phi.set(alpha);
		  }
		} // end of level-in-column loop
	      
	      } // if (states[k].size() > 0)
	    } // end of column loop
	  
	    // The state is finished. Add state to basis list.
	    if (((Nsum <= R) && (use_energy_cut)) || (!use_energy_cut))
	    {
	      sd_basis.push_back(phi);
	      dim++;

	    }
	  
	    finished = more_general_inc(ivec,imin,imax);
	  } // while (!finished), i.e., product state loop
	
	
	
	} // if fit.
      
      
      
      } // spin-state loop


      // make sure at least one step forward is taken
      m[A-1] = mmax + 1;

    } // angmom loop

    t.stop();
    t2.stop();

    // sort the basis.
    //std::sort(basis.begin(), basis.end());

    cerr << "I found " << dim << " states." << endl;
    cerr << "I discarded " << discarded << " states." << endl;
    //cerr << "Original method found " << Hilb.getBasis().size() << " states." << endl;
  
  
  
  }



  void QdotHilbertSpace::checkBasis()
  {

    cerr << "Chekcing basis functions against constraints ..." << endl;

    size_t dim = sd_basis.size();
    size_t k;
    bool ok = true;

    for (k=0; k<dim; ++k)
    {
      int nup, ndown;
      sd_basis[k].count_evenodd(nup, ndown);

      orbital_t last;
      last = sd_basis[k].last_index();

      if (nup-ndown != Sz)
      {
	cerr << "Basis function #" << k << " violates Sz = " << Sz << ", since " 
	     << "up-down == " << (nup-ndown) << endl;
	ok = false;
      }
      if (nup+ndown != A)
      {
	cerr << "Basis function #" << k << " has " << (nup+ndown) << " particles!" << endl;
	ok = false;
      }
    
      if ((int)last >= n_orbitals)
      {
	cerr << "Basis function #" << k << " has an alpha outside model space! " << endl;
	ok = false;
      }

      int m = 0;

      for (orbital_t alpha = 0; alpha<(size_t)n_orbitals; ++alpha)
      {
	if (sd_basis[k].get(alpha))
	{
	  m += qnumbers[alpha].m;
	}
      }
      if (m != M)
      {
	cerr << "Basis function #" << k << " has angmom == " << m << " != " << M << endl;
	ok = false;
      }


      assert(ok);

    }

  


  }


  void QdotHilbertSpace::printBasis(std::ostream& os)
  {
    size_t k;
    size_t dim = sd_basis.size();

    for (k = 0; k<dim; ++k)
    {
      os << "basis("<<k<<"+1) = struct('bin','" << sd_basis[k].to_binary(n_orbitals) << "','occ',[" << sd_basis[k].to_string() <<"]);" << endl;
    }

  }


} // namespace quantumdot
