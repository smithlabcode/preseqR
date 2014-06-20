/*    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith and Timothy Daley
 *
 *    Authors: Andrew D. Smith and Timothy Daley
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

#include "Rcontinued_fraction.hpp"
#include <cassert>
#include <numeric>

using std::vector;
using std::cerr;
using std::endl;
using std::min;

const double TOLERANCE = 1e-20;
const double DERIV_DELTA = 1e-8;


static double
GoodToulmin2xExtrap(const vector<double> &counts_hist){
    double two_fold_extrap = 0.0;
    for(size_t i = 0; i < counts_hist.size(); i++)
        two_fold_extrap += pow(-1.0, i + 1)*counts_hist[i];
    
    return two_fold_extrap;
}

static bool
check_yield_estimates_stability(const vector<double> &estimates) {
// make sure that the estimate is increasing in the time_step and
// is below the initial distinct per step_size
  for (size_t i = 1; i < estimates.size(); ++i){
    if (estimates[i] < estimates[i - 1] ){
      return false;
    }
    if(i >= 2 && (estimates[i] - estimates[i - 1] >
                  estimates[i - 1] - estimates[i - 2])){
      return false;
    }
  }
  return true;
}

//rescale value for real number
static double
get_rescale_value(const double numerator, const double denominator) {
  const double rescale_val = fabs(numerator) + fabs(denominator);
  if (rescale_val > 1.0/TOLERANCE)
    return 1.0/rescale_val;
  else if (rescale_val < TOLERANCE)
    return 1.0/rescale_val;
  return 1.0;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  QUOTIENT DIFFERENCE ALGORITHMS
////

/* quotient-difference algorithm to compute continued fraction
   coefficients
*/ 
static void
quotdiff_algorithm(const vector<double> &ps_coeffs, vector<double> &cf_coeffs) { //vector for power series coefficients & vector for continued fraction coefficients
  
  const size_t depth = ps_coeffs.size(); //degree of power series
  vector< vector<double> > q_table(depth, vector<double>(depth+1, 0.0));
  vector< vector<double> > e_table(depth, vector<double>(depth+1, 0.0));

  for (size_t j = 0; j < depth-1; j++) 
    q_table[1][j] = ps_coeffs[j + 1]/ps_coeffs[j];
  
  for (size_t j = 0; j < depth-1; j++) 
    e_table[1][j] = q_table[1][j + 1] - q_table[1][j] + e_table[0][j + 1];
    
  //using intial values of E(i)(j)'s and Q(i)(j)'s, fill rest of the q table and e table
  for (size_t i = 2; i < depth; i++) {
    for (size_t j = 0; j < depth; j++)
      q_table[i][j] = q_table[i - 1][j + 1]*e_table[i - 1][j + 1]/e_table[i - 1][j];
    
    for (size_t j = 0; j < depth; j++)
      e_table[i][j] = q_table[i][j + 1] - q_table[i][j] + e_table[i - 1][j + 1];
  }

  cf_coeffs.push_back(ps_coeffs[0]); //add first power series coefficient to end of vector for continued fraction coefficients

  //setting coefficients for continued fraction 
  for (size_t i = 1; i < depth; ++i) {
    if (i % 2 == 0) 
      cf_coeffs.push_back(-e_table[i/2][0]);
    else
      cf_coeffs.push_back(-q_table[(i + 1)/2][0]);
  }
}  


// compute CF coeffs when upper_offset > 0
//above the diagonal referring to degree of polynomial in numerator of Pade approximant is greater than degree of polynomial in the denominator 
static void
quotdiff_above_diagonal(const vector<double> &coeffs, const size_t offset,
                        vector<double> &cf_coeffs, vector<double> &offset_coeffs) {  
  //first offset coefficients set to first offset coeffs
  vector<double> holding_coeffs; 
  for (size_t i = offset; i < coeffs.size(); i++)
    holding_coeffs.push_back(coeffs[i]);
  
  // qd to determine cf_coeffs
  quotdiff_algorithm(holding_coeffs, cf_coeffs);
  for (size_t i = 0; i < offset; i++)
    offset_coeffs.push_back(coeffs[i]);
}


// calculate CF coeffs when lower_offset > 0
static void
quotdiff_below_diagonal(const vector<double> &coeffs, const size_t offset, 
                        vector<double> &cf_coeffs, vector<double> &offset_coeffs) {
  
  //need to work with reciprocal series g = 1/f, then invert
  vector<double> reciprocal_coeffs;
  reciprocal_coeffs.push_back(1.0/coeffs[0]);
  for (size_t i = 1; i < coeffs.size(); i++) {
    double holding_val = 0.0;
    for (size_t j = 0; j < i; ++j)
      holding_val += coeffs[i - j]*reciprocal_coeffs[j];
    reciprocal_coeffs.push_back(-holding_val/coeffs[0]);
  }
  
  //set offset_coeffs to 1st offset coeffs of 1/f 
  for (size_t i = 0; i < offset; i++)
    offset_coeffs.push_back(reciprocal_coeffs[i]);
  
  // qd to compute cf_coeffs using remaining coeffs
  vector<double> holding_coeffs;
  for (size_t i = offset; i < coeffs.size(); i++)
    holding_coeffs.push_back(reciprocal_coeffs[i]);
  
  quotdiff_algorithm(holding_coeffs, cf_coeffs);
}

// output new ContinuedFraction with decreased degree
// and coeffs equal to the old, but decreased in degree
ContinuedFraction
ContinuedFraction::decrease_degree(const ContinuedFraction &CF,
                                   const size_t decrement) {
  // create return ContinuedFraction
  ContinuedFraction decreasedCF;
  // properties of orig CF to decrement
  vector<double> decreased_ps_coeffs(CF.ps_coeffs);
  vector<double> decreased_cf_coeffs(CF.cf_coeffs);
  // decrease order
  for(size_t i = 0; i < decrement; i++) {
    decreased_ps_coeffs.pop_back();
    decreased_cf_coeffs.pop_back();
  }

  // just a copy
  vector<double> decreased_offset_coeffs(CF.offset_coeffs);

  // set return ContinuedFraction
  decreasedCF.ps_coeffs = decreased_ps_coeffs;
  decreasedCF.cf_coeffs = decreased_cf_coeffs;
  decreasedCF.offset_coeffs = decreased_offset_coeffs;
  decreasedCF.diagonal_idx = CF.diagonal_idx;
  decreasedCF.degree = CF.degree - decrement;

  return decreasedCF;
}

ContinuedFraction
ContinuedFraction::truncate_degree(const ContinuedFraction &CF,
				   const size_t n_terms){
  ContinuedFraction truncated_CF;
  if(CF.degree < n_terms){
    cerr << "current CF degree   = " << CF.degree << endl;
    cerr << "truncated CF degree = " << n_terms << endl; 
	try
	{
    	throw SMITHLABException("degree of truncate CF must be at least as large as current");
	}
	catch(SMITHLABException &e)
	{
		cerr << e.what() << endl;
		return truncated_CF;
	}
  }

  vector<double> truncated_ps_coeffs(CF.ps_coeffs);
  vector<double> truncated_cf_coeffs(CF.cf_coeffs);
  vector<double> truncated_offset_coeffs(CF.offset_coeffs);

  truncated_ps_coeffs.resize(n_terms);
  truncated_cf_coeffs.resize(n_terms - truncated_offset_coeffs.size());

  truncated_CF.ps_coeffs = truncated_ps_coeffs;
  truncated_CF.cf_coeffs = truncated_cf_coeffs;
  truncated_CF.offset_coeffs = truncated_offset_coeffs;
  truncated_CF.diagonal_idx = CF.diagonal_idx;
  truncated_CF.degree = n_terms;

  return truncated_CF;
}

ContinuedFraction::ContinuedFraction(const vector<double> &ps_cf, 
                                     const int di, const size_t dg) :
  ps_coeffs(ps_cf), diagonal_idx(di), degree(dg) {

  if (diagonal_idx == 0)
    quotdiff_algorithm(ps_coeffs, cf_coeffs);
  else if (diagonal_idx > 0)
    quotdiff_above_diagonal(ps_coeffs, diagonal_idx, cf_coeffs, offset_coeffs);
  else // if(cont_frac_estimate.lower_offset > 0) {
    quotdiff_below_diagonal(ps_coeffs, -diagonal_idx, cf_coeffs, offset_coeffs);
  // notice the "-" above so that -diagonal_idx > 0
}

void ContinuedFraction::copy2pointers(double * ps, int * ps_l, double * cf, int * cf_l, double * off, int * di, int * de)
{
	*ps_l = ps_coeffs.size();
	for (size_t i = 0; i != ps_coeffs.size(); i++)
		ps[i] = ps_coeffs[i];
	*cf_l = cf_coeffs.size();
	for (size_t i = 0; i != cf_coeffs.size(); i++)
		cf[i] = cf_coeffs[i];
	for (size_t i = 0; i != offset_coeffs.size(); i++)
		off[i] = offset_coeffs[i];
	*di = diagonal_idx;
	*de = degree;
	return;
}

/*
extern "C"
{
	void R_Advance_ContinueFraction(double * ps_cf, int * ps_cf_l, int * di, int * dg, double * PS_COEFFS, int * PS_COEFFS_L, double * CF_COEFFS, int * CF_COEFFS_L, double * OFFSET_COEFFS)
	{
		if (*ps_cf_length > 0 && *dg >= 0)
		{
			int last_nonzero = 0; // The number of continuous nonzero items in histogram from beginning
			for (int i = 0; i != *ps_cf_l; i++)
				if (ps_cf[i] != 0)
					last_nonzero++;
			vector<double> ps_coeffs(ps_cf, ps_cf + last_nonzero);
			ContinuedFraction::ContinuedFraction CF(ps_coeffs, *di, *dg);
			*di = CF::diagonal_idx;
			*dg = CF::degree;
			*PS_COEFFS_l = CF::ps_coeffs;
			for (int i = 0; i != *PS_COEFFS_L; i++)
				PS_COEFFS[i] = CF::ps_coeffs[i];
			*CF_COEFFS_L = CF:cf_coeffs.size();
			for (int i = 0; i != *CF_COEFFS_L; i++)
				CF_COEFFS[i] = CF::cf_coeffs[i];
			if (CP::offset_coeffs.size() > 0)
				for (int i = 0; i != CF::offset_coeffs.size(); i++)
					OFFSET_COEFFS[i] = CF::offset_coeffs[i];
		}
	}
}
*/



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  FUNCTIONS TO EVALUATE CONTINUED FRACTIONS AT A POINT
////

/* evaluate CF when upper_offset > 0 using Euler's recursion
 */
static double
evaluate_above_diagonal(const vector<double> &cf_coeffs,
                        const vector<double> &offset_coeffs,
                        const double val, const size_t depth) {
  
  double current_num = 0.0;
  double prev_num1 = cf_coeffs[0];
  double prev_num2 = 0.0;
  
  double current_denom = 0.0;
  double prev_denom1 = 1.0;
  double prev_denom2 = 1.0; 
  
  for (size_t i = 1; i < min(cf_coeffs.size(),
                             depth - offset_coeffs.size()); i++) {
    // initialize
    current_num = prev_num1 + cf_coeffs[i]*val*prev_num2;
    current_denom = prev_denom1 + cf_coeffs[i]*val*prev_denom2;
    
    prev_num2 = prev_num1;
    prev_num1 = current_num;
    
    prev_denom2= prev_denom1;
    prev_denom1 = current_denom;
    
    //rescale to avoid over- and underflow
    const double rescale_val = get_rescale_value(current_num, current_denom);
    
    current_num = current_num*rescale_val;
    current_denom = current_denom*rescale_val;
    
    prev_num1 = prev_num1*rescale_val;
    prev_num2 = prev_num2*rescale_val;

    prev_denom1 = prev_denom1*rescale_val;
    prev_denom2 = prev_denom2*rescale_val;
  }

  double offset_part = 0.0;
  for (size_t i = 0; i < offset_coeffs.size(); i++)
    offset_part += offset_coeffs[i]*pow(val, i);
  
  return offset_part + pow(val, min(depth, offset_coeffs.size()))*
    current_num/current_denom;
} 


// calculate ContinuedFraction approx when lower_offdiag > 0
static double 
evaluate_below_diagonal(const vector<double> &cf_coeffs,
                        const vector<double> &offset_coeffs,
                        const double val, const size_t depth) {
  
  //initialize
  double current_num = 0.0;
  double prev_num1 = cf_coeffs[0];
  double prev_num2 = 0.0;

  double current_denom = 0.0;
  double prev_denom1 = 1.0;
  double prev_denom2 = 1.0; 

  for (size_t i = 1; i < min(cf_coeffs.size(),
                             depth - offset_coeffs.size()); i++) {

    // recursion
    current_num = prev_num1 + cf_coeffs[i]*val*prev_num2;
    current_denom = prev_denom1 + cf_coeffs[i]*val*prev_denom2;

    prev_num2 = prev_num1;
    prev_num1 = current_num;

    prev_denom2= prev_denom1;
    prev_denom1 = current_denom;

    const double rescale_val = get_rescale_value(current_num, current_denom);

    current_num = current_num*rescale_val;
    current_denom = current_denom*rescale_val;

    prev_num1 = prev_num1*rescale_val;
    prev_num2 = prev_num2*rescale_val;
    
    prev_denom1 = prev_denom1*rescale_val;
    prev_denom2 = prev_denom2*rescale_val;
  }
  
  double offset_terms = 0.0;
  for (size_t i = 0; i < min(offset_coeffs.size(), depth); i++)
    offset_terms += offset_coeffs[i]*pow(val, i);
  
  // recall that if lower_offset > 0, we are working with 1/f, invert approx
  return 1.0/(offset_terms + pow(val, min(offset_coeffs.size(),depth))*
              current_num/current_denom);
}


// calculate ContinuedFraction approx when there is no offset
// uses euler's recursion
static double
evaluate_on_diagonal(const vector<double> &cf_coeffs, 
                     const double val, const size_t depth) {
  
  // initialize
  double current_num = 0.0;
  double prev_num1 = cf_coeffs[0];
  double prev_num2 = 0.0;

  double current_denom = 0.0;
  double prev_denom1 = 1.0;
  double prev_denom2 = 1.0; 

  for (size_t i = 1; i < min(cf_coeffs.size(), depth); i++) {
    // recursion
    current_num = prev_num1 + cf_coeffs[i]*val*prev_num2;
    current_denom = prev_denom1 + cf_coeffs[i]*val*prev_denom2;

    prev_num2 = prev_num1;
    prev_num1 = current_num;

    prev_denom2= prev_denom1;
    prev_denom1 = current_denom;

    const double rescale_val = get_rescale_value(current_num, current_denom);
    
    current_num = current_num*rescale_val;
    current_denom = current_denom*rescale_val;
    
    prev_num1 = prev_num1*rescale_val;
    prev_num2 = prev_num2*rescale_val;

    prev_denom1 = prev_denom1*rescale_val;
    prev_denom2 = prev_denom2*rescale_val;
  }
  return current_num/current_denom;
}


// calculate cont_frac approx depending on offset
double
ContinuedFraction::operator()(const double val) const {
  if (diagonal_idx > 0)
    return evaluate_above_diagonal(cf_coeffs, offset_coeffs, val, degree);
  
  if (diagonal_idx < 0)
    return evaluate_below_diagonal(cf_coeffs, offset_coeffs, val, degree);
  
  return evaluate_on_diagonal(cf_coeffs, val, degree);
}


extern "C"
{
	void R_Calculate_ContinuedFraction(double * cf_coeffs, int * cf_coeffs_l, double * offset_coeffs, int * di, int * de, double * coordinate, double * result)
	{
		vector<double> cf(cf_coeffs, cf_coeffs + *cf_coeffs_l);
		if (*di > 0)
		{
			vector<double> off(offset_coeffs, offset_coeffs + *di);
			*result = evaluate_above_diagonal(cf, off, *coordinate, *de);
		}
		else if (*di < 0)
		{
			vector<double> off(offset_coeffs, offset_coeffs - *di);
			*result = evaluate_below_diagonal(cf, off, *coordinate, *de);
		}		
		else
			*result = evaluate_on_diagonal(cf, *coordinate, *de);
	}
}


// Extrapolates the curve, for given values (step & max) and numbers
// of terms
void
ContinuedFraction::extrapolate_distinct(const vector<double> &counts_hist,
                                        const double max_value, 
                                        const double step_size,
                                        vector<double> &estimates) const {
  const double hist_sum = accumulate(counts_hist.begin(), counts_hist.end(), 0.0);
  estimates.clear();
  estimates.push_back(hist_sum);
  for (double t = step_size; t <= max_value; t += step_size)
    estimates.push_back(hist_sum + t*operator()(t));
}

extern "C"
{
	void R_Advance_extrapolate_distinct(double * cf_coeffs, int * cf_coeffs_l, double * offset_coeffs, int * di, int * de, double * hist, int * hist_l, double * start_size, double * step_size, double * max_size, double * estimate, int * estimate_l)
	{
		double hist_sum = 0;
		for (int i = 0; i != *hist_l; i++)
			hist_sum += hist[i];
		double result = 0;
		vector<double> est;
		est.push_back(hist_sum);
		for (double t = *start_size; t <= *max_size; t += *step_size)
		{
			R_Calculate_ContinuedFraction(cf_coeffs, cf_coeffs_l, offset_coeffs, di, de, &t, &result);
			est.push_back(hist_sum + t * result);
		}
		*estimate_l = est.size();
		for (int i = 0; i != *estimate_l; i++)
			estimate[i] = est[i];
	}
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////////
/////////////////
////////////////  CONTINUED FRACTION APPROXIMATION CLASS BELOW
///////////////
//////////////
/////////////
////////////

const size_t ContinuedFractionApproximation::MIN_ALLOWED_DEGREE = 4;

const double ContinuedFractionApproximation::SEARCH_MAX_VAL = 100;

const double ContinuedFractionApproximation::SEARCH_STEP_SIZE = 0.05;


// calculate cf_coeffs depending on offset
ContinuedFractionApproximation::ContinuedFractionApproximation(const int di, const size_t mt, 
                                                               const double ss, const double mv) :
  diagonal_idx(di), max_terms(mt), step_size(ss), max_value(mv) {}


/* Finds the optimal number of terms (i.e. degree, depth, etc.) of the
 * continued fraction by checking for stability of estimates at
 * specific points for yield.
 */
// New way for searching for optimal CF
ContinuedFraction
ContinuedFractionApproximation::optimal_cont_frac_distinct(const vector<double> &counts_hist) const {
  //do this outside
  // ensure that we will use an underestimate
  //  const size_t local_max_terms = max_terms - (max_terms % 2 == 1); 
 
  if(max_terms >= counts_hist.size()){
    cerr << "max terms = " << max_terms << endl;
    cerr << "hist size = " << counts_hist.size() << endl;
  } 
  assert(max_terms < counts_hist.size());
  
  // counts_sum = number of total captures
  double counts_sum  = 0.0;
  for(size_t i = 0; i < counts_hist.size(); i++)
    counts_sum += i*counts_hist[i];
  
  vector<double> full_ps_coeffs;
  for (size_t j = 1; j <= max_terms; j++)
    full_ps_coeffs.push_back(counts_hist[j]*pow(-1, j + 1));

  ContinuedFraction full_CF(full_ps_coeffs, -1, max_terms);  

  // if max terms = 4, check only that degree
  if(max_terms == 4 || max_terms == 3 || max_terms == 5 || max_terms == 6)
  {   
    vector<double> estimates;
    full_CF.extrapolate_distinct(counts_hist, SEARCH_MAX_VAL, SEARCH_STEP_SIZE, estimates);
    // return the continued fraction if it is stable
    if (check_yield_estimates_stability(estimates))
      return full_CF;
  }
  else{
    //if max terms >= 8, start at 8 and check increasing cont frac's
    size_t curr_terms = 0;
    if(max_terms % 2 == 0)
      curr_terms = 8;
    else
      curr_terms = 7;
    while (curr_terms <= max_terms) {    
      ContinuedFraction curr_cf 
	= ContinuedFraction::truncate_degree(full_CF, curr_terms);
      vector<double> estimates;
      curr_cf.extrapolate_distinct(counts_hist, SEARCH_MAX_VAL, SEARCH_STEP_SIZE, estimates);
          
    // return the continued fraction if it is stable
      if (check_yield_estimates_stability(estimates))
	return curr_cf;
    
      curr_terms += 2;
    // if not cf not acceptable, increase degree
    }
  }
   // no stable continued fraction: return null
  return ContinuedFraction();  
}

extern "C"
{
	void continuedfraction_estimate(double * hist_count, int * hist_count_l, int * di, int * mt, double * ss, double * mv, double * ps_coeffs, int * ps_coeffs_l, double * cf_coeffs, int * cf_coeffs_l, double * offset_coeffs, int * diagonal_idx, int * degree, int * is_valid)
	{				
		ContinuedFractionApproximation CFA(*di, *mt, *ss, *mv);
		const vector<double> counts_hist(hist_count, hist_count + *hist_count_l);
		ContinuedFraction CF = CFA.optimal_cont_frac_distinct(counts_hist);
		if (CF.is_valid())
			*is_valid = 1;
		else
			*is_valid = 0;
		CF.copy2pointers(ps_coeffs, ps_coeffs_l, cf_coeffs, cf_coeffs_l, offset_coeffs, diagonal_idx, degree);
		return;
	}
}


