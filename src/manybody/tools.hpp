#ifndef _TOOLS_HPP_
#define _TOOLS_HPP_

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

/**
 *
 * \file tools.hpp
 * \author Simen Kvaal
 * \date 6-13-08
 * 
 * \brief Some mathematical and other tools.
 *
 */

#include <string>
#include <sstream>
#include <cmath>
#include <vector>

/// \brief Simple function to format a double
/// with a given precision. Useful if we don't want to mess
/// with, e.g., cout's settings and produce predictable double-output.
///
/// \param x    number to write
/// \param prec precision
inline
std::string format(double x, int prec = 16)
{
  std::stringstream ss;
  ss.precision(prec);
  ss << std::right;
  ss.setf(std::ios::fixed,std::ios::floatfield);
  ss << x;
  return ss.str();
}

/// \brief Print copyright notice.
inline void print_copyright(std::ostream& os)
{
  using namespace std;

  os << "Copyright (c) 2008 Simen Kvaal" << endl
     << endl
//     << "This executable is part of OpenFCI." << endl
//     << endl
//      << "This is free software: you can redistribute it and/or modify" << endl
//      << "it under the terms of the GNU General Public License as published by" << endl
//      << "the Free Software Foundation, either version 3 of the License, or" << endl
//      << "(at your option) any later version." << endl
//      << endl
     << "This program is distributed in the hope that it will be useful," << endl
     << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl
     << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl
     << "GNU General Public License for more details." << endl
//     << endl
//      << "You should have received a copy of the GNU General Public License" << endl
//      << "along with OpenFCI. If not, see <http://www.gnu.org/licenses/>." << endl
      << endl;

}


/// \brief Move to position (using ANSI escape seq.)
/// \param N   position on line, should be >=0.
inline
std::string move_to(int N)
{
  std::stringstream s;
  s << "\033[" << N << "G";
  return s.str();
}



/// \brief Compute ln(n!)
/// \param n number to take factorial of
inline double lnfact(int n)
{
  double result = 0;
  for (int j=2; j<=n; ++j)
    result += log(j);
  return result;
}


/// \brief Compute n!. 
/// \param n number to take factorial of
inline double fact(int n)
{
  return exp(lnfact(n));
}
 
/// \brief Compute log of factorial, ln(n!/m!/(n-m)!)
/// \param n a number
/// \param m a number
inline double lnbinom(int n,int m)
{
  return lnfact(n) - lnfact(m) - lnfact(n-m);
}

/// \brief Compute "n ober m" =  n!/(m!(n-m)!).  
/// \param n a number
/// \param m a number
inline double binom(int n,int m)
{
  return exp(lnbinom(n,m));
}

/// \brief Return string "yes" if true, "no" if false
template<typename T>
std::string yesno(const T& x)
{
  return x ? "yes" : "no";
}



/// \brief Compute pow(-1,i). Used a lot!
template<typename T>
inline int m1pow(const T& i)
{
  if (i % T(2) == T(1))
    return -1;
  else
    return 1;
}



/// \brief Compute square of number; used a lot!
template<typename T>
inline T sqr(const T& x)
{
  return x*x;
}

// /// Round to nearest integer.
// inline
// long int round(double a) {
//   return (long int)(a + 0.5);
// }

/// \brief Count bits in a word. See http://graphics.stanford.edu/~seander/bithacks.html
template<class T>
inline size_t count_bits(const T& x)
{
  T v(x);
  
  v = v - ((v >> 1U) & (T)~(T)0U/3U);                           // temp
  v = (v & (T)~(T)0U/15U*3U) + ((v >> 2U) & (T)~(T)0U/15U*3U);  // temp
  v = (v + (v >> 4)) & (T)~(T)0U/255U*15U;                      // temp
  
  return (T)(v * ((T)~(T)0U/255U)) >> (unsigned int)(sizeof(v) - 1U) * CHAR_BIT; // count
  
}

/// \brief Compute a binary string version of number
/// \param a    Number to convert
/// \param n    Number of bits to use; n = 0 (default) means all bits.
template<typename T>
std::string to_binary(const T& a, const size_t& n = 0)
{
  size_t bits = n;
  if (bits == 0)
    bits = CHAR_BIT * sizeof(T);

  std::stringstream ss;

  for (int k = bits-1; k>=0; --k)
  {
    if ((a >> k) & T(1))
      ss << "1";
    else 
      ss << "0";
  }
  
  return ss.str();

}


/// \brief Simple recursive function to step through multi-indices.
///
/// The template parameter Vector must implement Vector::size() and Vector::operator[]().
///
/// 
/// \param alpha     Multi-index
/// \param min       Minimum value of each index
/// \param max       Maximum value + 1 of each index
/// \param digit     Which digit to increase. User should \em not use this; used for recursion.
/// \return true if alpha = (max-1,max-1,....,max-1) upon entry, such that it wraps around to (min,min,...min); false otherwise.
template<class Vector>
inline bool general_inc(Vector& alpha, int min, int max, size_t digit = 0)
{
  alpha[digit]++;
  if (alpha[digit] == max)
  {
    alpha[digit] = min;
    if (digit+1 < alpha.size())
      return general_inc(alpha, min, max, digit+1);
    else
    {
      return true;
    }
  }
  return false;
}

/// \brief Simple recursive function to step through multi-indices. More general version of
/// general_inc(). Here, the max and min values are different for each index in the multi-index.
///
/// The template parameter Vector must implement Vector::size() and Vector::operator[]().
///
/// 
/// \param alpha     Multi-index
/// \param min       Minimum value of each index
/// \param max       Maximum value + 1 of each index
/// \param digit     Which digit to increase. User should \em not use this; used for recursion.
/// \return true if alpha = (max[0]-1,max[1]-1,....,max[n]-1) upon entry, such that it wraps around to (min[0],min[1]...,min[n]); false otherwise.
inline bool more_general_inc(std::vector<int>& alpha, 
                             const std::vector<int>& min, 
                             const std::vector<int>& max, int digit = 0)
{
  alpha[digit]++;
  if (alpha[digit] == max[digit])
  {
    alpha[digit] = min[digit];
    if (digit+1 < (int)alpha.size())
      return more_general_inc(alpha, min, max, digit+1);
    else
    {
      return true;
    }
  }
  return false;
}







#endif
