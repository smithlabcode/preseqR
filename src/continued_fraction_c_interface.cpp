/*    Copyright (C) 2016 University of Southern California and
 *                       Chao Deng
 *
 *    Authors: Chao Deng
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <cmath>

#include "continued_fraction.h"

using std::vector;

#include <R_ext/Arith.h>

/*
 * construct a continued fraction given a power series
 * de is the degree of the constructed continued fraction
 *
 * =============================================================================
 * di               the diagonal diagonal value of the CF
 * is_valid         an indicator to show the validness of constructed CF
 * PS_coeffs        coefficients of the derived power series from the histogram
 * PS_coeffs_l      the length of the coefficients
 * ps_coeffs        all variables below are used to store the constructed
 * ps_coeffs_l      continued fraction
 * cf_coeffs        
 * cf_coeffs_l      
 * offset_coeffs
 * diagnomal_idx
 * degree
 * =============================================================================
 *
 */
extern "C" {
  void c_PS2CF(int *di, int *de, double *PS_coeffs, int *PS_coeffs_l, 
      double *ps_coeffs, int *ps_coeffs_l, double *cf_coeffs, int *cf_coeffs_l,
      double *offset_coeffs, int *diagonal_idx, int *degree, int *is_valid)
  { 
    //ps_coeffs and ps_coeffs_l store the given power series
    const vector<double> ps(PS_coeffs, PS_coeffs + *PS_coeffs_l);

    if (*de > ps.size()) {
      ContinuedFraction empty;
      *is_valid = empty.is_valid();
      return;
    }

    ContinuedFraction CF(ps, *di, *de);  
    vector<double> estimates;
    CF.extrapolate_distinct(100, 0.05, estimates);
    // checking the function is well-defined
    for (size_t i = 1; i < estimates.size(); ++i) 
      if (!R_FINITE(estimates[i])) {
        ContinuedFraction empty;
        *is_valid = empty.is_valid();
        return;
    }


    // store informatin relevent to the constructed continued fraction
    *is_valid = CF.is_valid();

    // copy information from ContinuedFraction into variables
    // size_t converts to a int type; same things happen below
    for (size_t i = 0; i != CF.ps_coeffs.size(); i++)
      ps_coeffs[i] = CF.ps_coeffs[i];
    *ps_coeffs_l = CF.ps_coeffs.size();

    for (size_t i = 0; i != CF.cf_coeffs.size(); i++)
      cf_coeffs[i] = CF.cf_coeffs[i];
    *cf_coeffs_l = CF.cf_coeffs.size();

    for (size_t i = 0; i != CF.offset_coeffs.size(); i++)
      offset_coeffs[i] = CF.offset_coeffs[i];

    *diagonal_idx = CF.diagonal_idx;
    *degree = CF.degree;

    return;
  }
}
