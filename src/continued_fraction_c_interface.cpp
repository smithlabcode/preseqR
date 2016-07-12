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
 * the C-encode interface to calculate the value of a continued fraction given 
 * the coordinate
 *
 * =============================================================================
 * cf_coeffs        coefficients of the continued fraction (CF)
 * cf_coeffs_l      the number of coefficients  
 * offset_coeffs    the offset coefficients of the CF
 * di               the diagonal diagonal value of the CF
 * de               the degree of the CF
 * coordinate       the coordinate
 * result           store the function value given the coordinate
 * =============================================================================
 *
 */
extern "C" {
  void c_calculate_continued_fraction(double *cf_coeffs, int *cf_coeffs_l, 
                                      double *offset_coeffs, int *di, int *de,
                                      double *coordinate, double *result) 
  {
    ContinuedFraction CF;
    CF.cf_coeffs.assign(cf_coeffs, cf_coeffs + *cf_coeffs_l);
    CF.diagonal_idx = *di;
    CF.degree = *de;

    if (CF.diagonal_idx > 0)
      CF.offset_coeffs.assign(offset_coeffs, offset_coeffs + *di);
    else if (CF.diagonal_idx < 0)
      CF.offset_coeffs.assign(offset_coeffs, offset_coeffs - *di);
    *result = CF(*coordinate);

    return;
  }
}

/* 
 * the C-encoded interface to extralate from a continued fraction
 * the initial experiment is not counted
 *
 * =============================================================================
 * cf_coeffs        coefficients of the continued fraction (CF)
 * cf_coeffs_l      the number of coefficients  
 * offset_coeffs    the offset coefficients of the CF
 * di               the diagonal diagonal value of the CF
 * de               the degree of the CF
 * estimate         store values of the continued fraction
 * estimate_l       the number of stored values
 *
 * start_size       starting sample size of extrapolation
 * step_size        step size for each extrapolation
 * max_size         maximum sample size of extrapolation
 * =============================================================================
 *
 */
extern "C" 
{
  void c_extrapolate_distinct(double *cf_coeffs, int *cf_coeffs_l, 
                              double *offset_coeffs, int *di, int *de, 
                              double *start_size,
                              double *step_size, double *max_size, 
                              double *estimate, int *estimate_l) 
  { 
    double result = 0;
    vector<double> est;
    for (double t = *start_size; t <= *max_size; t += *step_size) 
    {
      c_calculate_continued_fraction(cf_coeffs, cf_coeffs_l, offset_coeffs, di, 
                                     de, &t, &result);
      est.push_back( t*result );
    }

    *estimate_l = est.size();
    for (int i = 0; i != *estimate_l; i++)
      estimate[i] = est[i];
  }
}

/*
 * a c-encoded interface to construct continued fraction given the histogram
 *
 * =============================================================================
 * hist_count       the count vector of the histogram
 * hist_count_l     the length the count vector
 * di               the diagonal diagonal value of the CF
 * mt               the maximum number of terms to try for a CF
 * is_valid         an indicator to show the validness of constructed CF
 *
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
  void c_continued_fraction_estimate(double *hist_count, int *hist_count_l, 
                                     int *di, int *mt, double *ps_coeffs, 
                                     int *ps_coeffs_l, double *cf_coeffs, 
                                     int *cf_coeffs_l, double *offset_coeffs, 
                                     int *diagonal_idx, int *degree, 
                                     int *is_valid) 
  { 
    ContinuedFractionApproximation CFA(*di, *mt);
    const vector<double> counts_hist(hist_count, hist_count + *hist_count_l);
    ContinuedFraction CF(CFA.optimal_cont_frac_distinct(counts_hist));

    // store informatin relevent to the constructed continued fraction
    *is_valid = CF.is_valid();

    // copy information from ContinuedFraction into points
    //size_t converts to a int type; same things happen below
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
