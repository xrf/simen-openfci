#ifndef _QDOT_INTERACTION_HPP_
#define _QDOT_INTERACTION_HPP_

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


#include <cassert>

#include "linalg.hpp"
#include "sparse.hpp"
#include "RadialPotential.hpp"

/**
 * \file QdotInteraction.hpp
 * 
 * \brief Class definition of quantum dot interaction matrix elements, plus various helper
 * functions related to this.
 *
 * QdotInteraction takes care of building lab frame matrix elements <a b | U | c d> 
 * given a radial potential U(r_12), encapsulated in RadialPotential.
 *
 *
 *
 */


namespace quantumdot {


/// \brief Internal struct used in QdotInteraction.
typedef struct
{
    int N1, m1, N2, m2;
} orb_pair;

/// \brief Computes lab frame matrix elements <a b|V|c d> in HO basis of
/// the Coulomb interaction (or another radial potential).
///
/// By default, the interaction matrix is given in COM coordinates,
/// then transformed the lab frame on the fly using the Talmi matrices. This
/// is the fastest option for 2 or 3-body problems.
///
/// The user \em may also precompute lab frame matrix elements. This is especially
/// handy if many-body calculations (as opposed to 2 or 3-body) are performed.
/// Then, precomputeLabFrameMatrix() should be called before use. This
/// is really not useful for 2-body problems, since the number of matrix elements created
/// will be much higher than the resulting Hamiltonian matrix in general! For 3-body problems
/// this is not so \em asymptotically, but in \em practice it is.
///
/// Most likely, an instance of QdotInteraction will be fed to a QdotMatrixElements instance.
/// The latter computes the matrix elements <(a sigma) (b tau) | V | (c sigma') (d tau')>, i.e.,
/// it adds spin degrees of freedom.
///
/// QdotInteraction has an instance of RadialPotential which defines the actual potential.
/// By default, this is a Coulomb potential.
///
/// QdotInteraction also can compute the effective interaction in the lab frame.
///
///
/// Typical use:
/// \code
/// QdotInteraction m;
/// m.setR(5);
/// m.setRadialPotential( u ); // assign a RadialPotential instance, if not Coulomb is wanted.
/// m.buildInteractionBlocks(); // builds Coulomb interaction.
/// if (i_want_to_precompute_in_lab_frame)
///   m.precoputeLabFrameMatrix();
///
/// // Ready to use! Call singleElement(N1,m1,N2,m2,N3,m3,N4,m4) to calculate!
/// \endcode
class QdotInteraction
{

  private:
    /// Maximum shell index. (Number of shells = R + 1.)
    int R;

    /// Holds matrix elements of interaction in relative coordinates. C[|m|] belongs to angular momentum m, with the bare interaction.
    std::vector<dense_matrix> C;

    /// For the effective interaction we need *more* blocks, for varying COM energies.
    /// This vector helps with this: C[C_index[delta_R] + abs(m)] is the coulomb block
    /// for COM energy delta_R and relative coordinate angular momentum m. For the standard Coulomb, C_index[i] = 0 for all,
    /// so that we regain C[abs(m)].
    std::vector<int> C_index;

    /// Holds matrix elements of LAB to COM coordinates, or rather
    /// building blocks of this. (Also called Talmi matrices.)
    std::vector<dense_matrix> T;

    /// Strength of interaction.
    double lambda;

    /// Vector of sparse matrices holding a precomputed LAB-frame interaction, 
    /// if desired. U[M] is the matrix elements belonging to M = m1 + m2 = m3 + m4.
    std::vector<simple_sparse::SparseMatrix<size_t, double> > U;

    /// Helps fetching row and column numbers for U[M].
    std::vector<std::map<int, size_t> > state_map;

    /// True if we have precomputed the lab frame interaction, false otherwise.
    bool precomputed;

    /// The number of orbitals in the single-particle basis, needed for the precomputed LAB frame interaction.
    int n_orbitals;

    /// Interaction potential.
    RadialPotential potential;

  protected:
    /// \brief Build the Talmi matrices up to 2*R shells. 
    /// (We need 2*R blocks if we use
    /// direct product Hilbert space; only R otherwise. But this is insignificant in the greater scheme of things.)
    void buildTalmiBlocks();


    /// \brief Map quantum number pair to orbital number.
    /// Used for precomputing lab frame interaction. Harmonic oscillator energy
    /// is given by E = N + 1, and the usual radial quantum number n is defined through N = 2*n + abs(m).
    /// \param N    Shell number quantum number
    /// \param m    Angular momentum quantum number.
    int orbitalMap2(int N, int m)
    {
      return (N*(N+2)+m)/2;
    }
    
  public:
    /// \brief Default constructor.
    QdotInteraction()
    {
      lambda = 1.0;
      R = -1;
      precomputed = false;
    }

    /// \brief Assign a RadialPotential object. (The default object is a Coulomb potential.)
    void setRadialPotential(const RadialPotential& u) { potential = u; }

    /// \brief Build Coulomb matrix blocks in COM frame. This can easily be modified to
    /// produce a different interaction.
    void buildInteractionComBlocks();

    /// \brief Build *effective* Coulomb matrix blocks.
    /// \param g Use R*g shells for effective Coulomb blocks.
    void buildEffectiveInteractionComBlocks(int g);

    /// \brief Set the interaction strength.
    /// \param ell The value for lambda.
    void setLambda(double ell) { lambda = ell; }

    /// \brief Set the maximum shell number R.
    /// \param the_R The value for R.
    void setR(int the_R)
    {
      R = the_R;
      assert(R >= 0);
      buildTalmiBlocks();
    }
    
    /// \brief Compute the matrix element <a1 a2|U|b1 b2> of
    /// the interaction. 
    ///Here, a1 = (N1, m1), a2 = (N2, m2), b1 = (N1pr, m2pr) and
    /// b2 = (N2pr, m2pr). 
    double singleElement(int N1, int m1, int N2, int m2, int N1pr, int m1pr, int N2pr, int m2pr);

    /// \brief Analytic computation of LAB frame elements. Note: The order of arguments is not wrong!
    /// This is included only for completeness, and it is really slow and inaccurate for many shells.
    /// For many-particle problems (>3) it can be useful since we are then restricted to only a 
    /// few shells.
    ///
    /// \b NOTE: In the code in QdotInteraction.cc, it would be better to take the log of f0, f1 etc and use
    /// exp(f0+....) instead of f0*f1...
    /// Should implement this; if I haven't done it when you read this, do so
    /// yourself. :-)
    ///
    double singleElementAnalytic(int N1, int m1, int N2, int m2, int N4, int m4, int N3, int m3);

    /// \brief Precopute lab frame version of the interaction.
    void precomputeLabFrameMatrix();


};


/// \brief Compute the "Talmi matrix"; transforms lab -> com for two particles in 1D HO.
///
/// \param TN   Reference to destination matrix
/// \param N    Index of matrix
void computeTalmiMatrix(dense_matrix& TN, int N);



} // namespace quantumdot


#endif 
