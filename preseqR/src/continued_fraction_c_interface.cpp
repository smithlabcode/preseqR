#include "continued_fraction.h"
#include <vector>
using std::vector;

/* 
 * the C-encode interface to calculate the value of a continued fraction given 
 * the coordinate
 */
extern "C" {
	void c_calculate_continued_fraction( \
	      double *cf_coeffs, // coefficients of the continued fraction (CF)
		  int *cf_coeffs_l,  // the number of coefficients
          double *offset_coeffs,  // the offset coefficients of the CF
          int *di,  // the diagonal diagonal value of the CF
		  int *de,  // the degree of the CF
          double *coordinate, // the coordinate 
          double *result) { // store the function value given the coordinate
	  ContinuedFraction CF;
	  CF.cf_coeffs.assign(cf_coeffs, cf_coeffs + *cf_coeffs_l);
	  CF.diagonal_idx = *di;
	  CF.degree = *de;
	  if (CF.degree > 0)
		  CF.offset_coeffs.assign(offset_coeffs, offset_coeffs + *di);
	  else if (CF.degree < 0)
		  CF.offset_coeffs.assign(offset_coeffs, offset_coeffs - *di);
	  *result = CF(*coordinate);
	  return;
	}
}

/* 
 * the C-encoded interface to extralate from a continued fraction
 */
extern "C" 
{
  void c_extrapolate_distinct( \
        double *cf_coeffs, // coefficients of the continued fraction (CF)
        int *cf_coeffs_l, // the number of coefficients
        double *offset_coeffs, // the offset coefficients of the CF
        int *di, // the diagonal diagonal value of the continued fraction
        int *de, // the degree of the continued fraction
        double *hist, // the count vector of the histogram
        int *hist_l, // the length the count vector
        // start_size, step_size, max_size set the number of points sampled from
		// the continued fraction
        double *start_size,
        double *step_size, 
        double *max_size, 
        double *estimate, // store values of the continued fraction
        int *estimate_l) { // the number of stored values
	//the number of distinct molecules
    double hist_sum = 0;
    for (int i = 0; i != *hist_l; i++)
      hist_sum += hist[i];

    double result = 0;
    vector<double> est;
    for (double t = *start_size; t <= *max_size; t += *step_size) {
      c_calculate_continued_fraction(cf_coeffs, cf_coeffs_l, 
				    offset_coeffs, di, de, &t, &result);
      est.push_back(hist_sum + t * result);
    }
    *estimate_l = est.size();
    for (int i = 0; i != *estimate_l; i++)
      estimate[i] = est[i];
  }
}

/*
 * a c-encoded interface to construct continued fraction given the histogram
 */
extern "C" {
  void c_continued_fraction_estimate( \
        double *hist_count, // the count vector of the histogram
        int *hist_count_l, // the length of the count vector
        int *di, // the diagonal to work with for estimates
        int *mt, // the maximum number of terms to try for a CF
        double *ss, // the step size to use when training
        double *mv, // the largest value to check when training
		// all variables below are used to store the constructed CF
        double *ps_coeffs, 
        int *ps_coeffs_l, 
        double *cf_coeffs, 
        int *cf_coeffs_l, 
        double *offset_coeffs, 
        int *diagonal_idx, 
        int *degree, 
        int *is_valid) { // an indicator to show the validness of CF
    ContinuedFractionApproximation CFA(*di, *mt, *ss, *mv);
    const vector<double> counts_hist(hist_count, hist_count + *hist_count_l);
    ContinuedFraction CF(CFA.optimal_cont_frac_distinct(counts_hist));

	// store informatin relevent to the constructed continued fraction
    *is_valid = CF.is_valid();

	// copy information from ContinuedFraction into points
	for (size_t i = 0; i != CF.ps_coeffs.size(); i++)
		ps_coeffs[i] = CF.ps_coeffs[i];
	//size_t converts to a int type; same things happen below
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
