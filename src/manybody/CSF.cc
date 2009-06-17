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

#include <sstream>
#include "CSF.hpp"
#include "tools.hpp"
#include "lpp/lapack.hh"
#include "NChooseK.hpp"


namespace manybody {


  void computeCGCoeffs(dense_matrix& CG, dense_vector& s, vector<unsigned int>& basis, int n, int m)
  {
    // Create spin states.
    NChooseK<unsigned int> chooser;
    chooser.N = n + m;
    chooser.K = n;
    chooser.choose();
    basis = chooser.subset; // store basis.
    long int dim = basis.size();

    // build S^2 matrix
    CG.resize(dim, dim);
    s.resize(dim);
  
    double sz = 0.5*(n-m);
    double diag_elem = n + sz*sz - sz;
    long int j, k;
    for (j = 1; j<=dim; ++j)
    {
      CG(j,j) = diag_elem;
      for (k = j+1; k<=dim; ++k)
      {
	unsigned int x = basis[j-1] ^ basis[k-1];
	if (count_bits(x) == 2)
	{
	  CG(j,k) = 1;
	  CG(k,j) = 1;
	}
      }
    }

    // diagonalize
    long int info = 0;
    s.resize(dim);
    lpp::syev("v", "u", &dim, &CG(1,1), &dim, &s(1), &info);

    // convert s*(s+1) to s.
    for (j = 1; j<=dim; ++j)
      s(j) = sqrt(0.25 + s(j)) - 0.5;

  

  }




  void CsfMachine::computeSz()
  {
    size_t k;
    size_t dim = sd_basis.size();
    sd_basis[0].count_evenodd(nup, ndown);
    sz = nup - ndown;

    for (k = 0; k < dim; ++k)
    {
      int nup2, ndown2;
      sd_basis[k].count_evenodd(nup2, ndown2);
      if ((nup != nup2) || (ndown != ndown2))
      {
	sd_basis[0].print();
	sd_basis[k].print();
      }
      assert(nup == nup2);
      assert(ndown == ndown2);
    }

  }


  void CsfMachine::rearrangeBasis()
  {
    size_t dim = sd_basis.size();
    size_t k, k2;

    // Compute all the equivalence classes.

    vector<pair<Slater, Slater> > eq_class(dim);
    class_number.resize(dim);

    cerr << "Computing singly and doubly occupied orbitals..." << endl;
    for (k = 0; k<dim; ++k)
    {
      findSinglyOccupied(sd_basis[k], eq_class[k].first);
      findDoublyOccupied(sd_basis[k], eq_class[k].second);
      class_number[k] = 0;
    }

    // Group equivalence classes.

    cerr << "Grouping equivalence classes (rearranges Slater dets) ..." << endl;

    int current_class = 0;
    for (k = 0; k<dim; ++k)
    {
      // Create a new class if no class is assigned.
      if (class_number[k] == 0)
      {
	current_class++;
	class_number[k] = current_class;
    
	// Loop though the rest and assign class number if class signature
	// is equal.
	for (k2 = k+1; k2<dim; ++k2)
	{
	  if ((eq_class[k].first == eq_class[k2].first)
	      && (eq_class[k].second == eq_class[k2].second))
	  {
	    class_number[k2] = current_class;
	  }
	
	}
      }
    }

    // Rearrange basis based on equivalence class.

    //cerr << "Rearranging basis ..." << endl;

    // Copy class indices and corresponding SD index.
    vector<pair<int, size_t> > class_sorter(dim);
    for (k=0; k<dim; ++k)
      class_sorter[k] = make_pair(class_number[k], k);
  
    // Sort the vector (based on class index)
    sort(class_sorter.begin(), class_sorter.end());

    // Rearrange sd_basis and class_number.
    vector<Slater> temp2(dim);
    vector<int> temp(dim);
    temp2 = sd_basis;
    temp = class_number;

    for (size_t k=0; k<dim; ++k)
    {
      class_number[k] = temp[class_sorter[k].second]; 
      sd_basis[k] = temp2[class_sorter[k].second];
    }  


  }



  void CsfMachine::findBlocks()
  {
    size_t dim = sd_basis.size();
    size_t row;
  
    // Upon entry:
    // sd_basis is sorted accoding to class number,
    // which is defined by class_number[k] for sd_basis[k].

    // Resize the eq_classes vector. The number of elements is
    // the number of equivalence classes, which is the last element
    // of class_number.
    eq_classes.resize(class_number[class_number.size() - 1]);
    equivalence_class zero_class;
    zero_class.begin = 0;
    zero_class.end = 0;
    zero_class.up = 0;
    zero_class.down = 0;
    zero_class.n_doubly = 0;
    fill(eq_classes.begin(), eq_classes.end(), zero_class);
  

    for (row = 0; row<dim; ++row)
    {
      // Get the index.
      int c = class_number[row] - 1;
      assert(c < (int)eq_classes.size());

      // Update the corresponding eq_class element.
      if (eq_classes[c].end == eq_classes[c].begin)
      {
	eq_classes[c].begin = row;
	eq_classes[c].end = row+1;
      }
      else
      {
	eq_classes[c].end = row + 1;
      }

    }


    // Compute up and down spins.
    int n_doubly_max = 0;
    for (row = 0; row < eq_classes.size(); ++row)
    {

      Slater doubly, singly;
      findDoublyOccupied(sd_basis[eq_classes[row].begin], doubly);
      findSinglyOccupied(sd_basis[eq_classes[row].begin], singly);
      int n_doubly = doubly.count();
      eq_classes[row].up = nup - n_doubly;
      eq_classes[row].down = ndown - n_doubly;
      eq_classes[row].n_doubly = n_doubly;
      n_doubly_max = max(n_doubly, n_doubly_max);

      eq_classes[row].spins.resize(eq_classes[row].end - eq_classes[row].begin);
      fill(eq_classes[row].spins.begin(), eq_classes[row].spins.end(), 0u);
      vector<size_t> positions;
      singly.find_ones(positions); // where are the singly occupied orbitals?
      for (size_t k = eq_classes[row].begin; k<eq_classes[row].end; ++k)
      {
	for (size_t ell = 0; ell < positions.size(); ++ell)
	{
	  if (sd_basis[k].get(2*positions[ell]))  // a spin up here?
	    eq_classes[row].spins[k - eq_classes[row].begin] |= (1u << ell);
	}
      }

      assert(eq_classes[row].up >= 0);
      assert(eq_classes[row].down >= 0);
      assert(eq_classes[row].end > eq_classes[row].begin);
    
    
    }


    // compute clebsh gordan coefficients.
    cg_coeffs.resize(0);
    s_values.resize(0);
    cg_basis.resize(0);
    cg_n.resize(0);
    for (int n_doubly = 0; n_doubly<=n_doubly_max; ++n_doubly)
    {
      dense_matrix CG;
      dense_vector s;
      vector<unsigned int> basis(0);
      // Special treatment for zero singly occupied
      if ((nup - n_doubly == 0) && (ndown - n_doubly == 0))
      {
	CG.resize(1,1);
	CG(1,1) = 1.0;
	s.resize(1);
	s(1) = 0.0; // Always S = 0 then!
	basis.push_back(0u);
      }
      else
	computeCGCoeffs(CG, s, basis, nup - n_doubly, ndown - n_doubly) ;

      cg_coeffs.push_back(CG);
      s_values.push_back(s);
      cg_basis.push_back(basis);
      cg_n.push_back(nup + ndown - 2*n_doubly);
    }			    


  }



  void CsfMachine::print()
  {

    size_t k;

    for (k = 0; k<eq_classes.size(); ++k)
    {
      cerr << "Eq. class #" << k+1 << ": ";
      cerr << "  [" << eq_classes[k].begin << ", " << eq_classes[k].end << "), up/down = " << eq_classes[k].up << "/" << eq_classes[k].down  << endl;;
      int relevant_spins_down =  eq_classes[k].down;
      int relevant_spins_up = eq_classes[k].up;
      NChooseK<unsigned int> blabla;
      blabla.N = relevant_spins_up + relevant_spins_down;
      blabla.K = relevant_spins_up;
      if (blabla.N > 0)
	blabla.choose();
      int correct_block_size;
      if (blabla.N > 0)
	correct_block_size = blabla.subset.size();
      else
	correct_block_size = 1;
      assert(correct_block_size == (int)eq_classes[k].end - (int)eq_classes[k].begin);
    }

    for (k=0; k<cg_coeffs.size(); ++k)
    {
      cerr << cg_coeffs[k] << endl << s_values[k] << endl;
      for (size_t ell = 0; ell<cg_basis[k].size(); ++ell)
	cerr << to_binary(cg_basis[k][ell], cg_n[k]) << " ";
      cerr << endl;
    }


  }


  void CsfMachine::getBlock(csf_block& CSF, int S, size_t eq_index)
  {
    // Assign eng points of block of Slater determinants
    CSF.sd_begin = eq_classes[eq_index].begin;
    CSF.sd_end = eq_classes[eq_index].end;
    size_t dim = CSF.sd_end - CSF.sd_begin;
    CSF.n_csf = 0;

    // Which CG matrix do we need?
    int cg_index = eq_classes[eq_index].n_doubly;
    assert(cg_coeffs[cg_index].numRows() == (int)dim);

    // Deduce which CSFs to include
    vector<size_t> which(0);
    for (size_t k = 1; k<=dim; ++k)
    {
      if (S == (int)(2*s_values[cg_index](k) + 0.5 ))
	which.push_back(k);
    }

    CSF.n_csf = which.size();

    // Allocate matrix elements.
    if (CSF.n_csf > 0)
    {
      CSF.coeffs.resize(CSF.n_csf, dim);
    }

    // Save matrix elements. 
    for (size_t k = 1; k<=CSF.n_csf; ++k)
    {
      size_t k2 = which[k-1];
      for (size_t j = 1; j<=dim; ++j)
      {
	for (size_t ell = 0; ell<cg_basis[cg_index].size(); ++ell)
	  if (cg_basis[cg_index][ell] == eq_classes[eq_index].spins[j-1])
	  {
	    CSF.coeffs(k, j) = cg_coeffs[cg_index](ell+1, k2);
	    break;
	  }
      }
    }


  }

  void CsfMachine::getAllBlocks(vector<csf_block>& CSF, int S)
  {
    CSF.resize(0);
    csf_block temp;
    for (size_t k = 0; k<(size_t)numBlocks(); ++k)
    {
      getBlock(temp, S, k);
      if (temp.n_csf > 0)
	CSF.push_back(temp);
    }

  }

  void CsfMachine::getAllBlocksTrivial(vector<csf_block>& CSF)
  {
    // Create a vector with trivial CSFs, i.e., one CSF per Slater determinant,
    // with 1x1 coefficient matrix equal to [1].

    CSF.resize(0);
    csf_block temp;
    for (size_t k = 0; k<sd_basis.size(); ++k)
    {
      temp.sd_begin = k;
      temp.sd_end = k+1;
      temp.n_csf = 1;
      temp.coeffs.resize(1,1);
      temp.coeffs(1,1) = 1.0;
      CSF.push_back(temp);
    }

  }


  void computeCSFBlocks(std::vector<Slater>& sd_basis, std::vector<csf_block>& csf_basis, int S)
  {
    // Detect constant Sz subspaces.

    // int --> Sz, vector<Slater> --> Slater dets with spin projection Sz.
    map<int, vector<Slater> > sz_spaces;

    // Build the individual Sz bases.
    for (size_t k = 0; k < sd_basis.size(); ++k)
    {
      int up, down;
      sd_basis[k].count_evenodd(up, down);
      int Sz = up - down;
      sz_spaces[Sz].push_back(sd_basis[k]);
    }


    // Now: For each value of Sz encountered, sz_spaces[Sz] is the corresponding Slaters.

    sd_basis.clear(); // erase all vectors; we're going to insert them again in different order.
    vector<csf_block> one_csf_basis; // will hold the csf_basis for a single value of Sz.

    // loop through all Sz values ...
    map<int, vector<Slater> >::iterator it;
    for (it = sz_spaces.begin(); it != sz_spaces.end(); ++it)
    {
      cerr << "Current Sz = " << it->first << " with " << it->second.size() << " Slater determinants." << endl;
      // Find the CSF basis for Sz = it->first.
      CsfMachine csf_machine;
      csf_machine.setSDBasis(it->second);
      csf_machine.rearrangeBasis();
      csf_machine.findBlocks();
      one_csf_basis.clear();
      csf_machine.getAllBlocks(one_csf_basis, S);
      it->second = csf_machine.getSDBasis(); // the basis is modified!
      
      // when adding the CSF blocks to the global csf_basis,
      // the Slater determinant indices are wrong, since they do not
      // know of the Slater determinants of the *previous* Sz spaces inserted
      // into sd_basis.
      size_t basis_offset = sd_basis.size();

      // insert all basis functions for the Sz value into sd_basis
      for (size_t k = 0; k < it->second.size(); ++k)
	sd_basis.push_back(it->second[k]);

      // Modify the CSF blocks by offseting the begining and endpoint 
      // of the sd-block used.
      // Also add one_csf_basis's elements to the global csf_basis vector.
      for (size_t k = 0; k < one_csf_basis.size(); ++k)
      {
	one_csf_basis[k].sd_begin += basis_offset;
	one_csf_basis[k].sd_end += basis_offset;
	csf_basis.push_back(one_csf_basis[k]);
      }
      

    }

  }



}
