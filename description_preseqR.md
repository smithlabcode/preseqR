#Contents, Interface Design, Principles, Implementation and  Difficulties

##Contents
Three files are included: Rcontinuedfraction.hpp, Rcontinuedfraction.cpp and 
Rcontinuedfraction.R. The Rcontinuedfraction.cpp is a source code containing 
all required functions written in c++. The Rcontinuedfraction.hpp is its 
header. The Rcontinuedfraction.R includes all functions written in R, which 
can call c++ functions, providing tools for capture and re-capture model.

##Interface of the preseqR
Four functions are supported by preseqR for common users:

  1.  preseqR.continued.fraction.estimate(): given a histogram, preseqR produces 
	  a continued rational function (CRF);
  2.  preseqR.bootstrap.complexity.curve(): given a histogram, preseqR produces
	  a complexity curve of the library;
  3.  preseqR.print(): write down results into a file after running functions
	  in preseqR
  4.  preseqR.calculate.continued.fraction(): given a CRF and its coordinates, 
      preseqR calculates the function value.

Two functions are supported by preseqR for advanced users:

  1.  preseqR.extrapolate.distinct(): extrapolate distinct molecules given a CRF;
  2.  preseqR.interpolation(): interpolate distinct molecules given a histogram;

All the functions are in the Rcontinuedfraction.R. They are written in R, but they 
may call c++ function during implementation.

Other functions in the Rcontinuedfraction.R are helpers, which supposed to be 
used by no one.

  1.  read.hist(): read a histogram file and output the count vector of 
	  the histogram
  2.  replace.sampling(): re-sampling histograms givan a histogram. It is based
	  on replacement sampling. 
  3.  nonreplace.sampling(): do with replacement sampling given a histogram
  4.  count.distinct(): count the number of distinct spicies/reads given
	  a sample data
  5.  goodtoulmin.2x.extrap() : check the goodness of the sample through its 
	  histogram based on Good & Toulmin's model

##Principles and implementation

1.	In order to pack all functions and make an available R package, all 
	platform-dependent codes should be removed or rewritten. For example codes call 
	*gsl library* or *unistd.h*. I am not sure whether smithlab_cpp requires any extra 
	libraries. So I just grabbed the piece of codes required by the
	Rcontinuedfraction.cpp to make sure that the Rcontinuedfraction.cpp can be 
	compiled properly and no other plat-specific libraries are needed. The only reason
	I do it this way is that it is simple and easy. We can discuss about how to 
	design. The sampling function based on *gsl library* has been rewritten in R 
	language. Thus gsl is not a required library for continuedfraction class.

2.	R functions could only call functions from the c code. To explain it from 
	what I learn, R could only call the function written in c with the format
	```
		void * function(type1 *pointer, type2 *pointer, type3 *pointer, ...)
	```
	and each type should be valid in c. In order to call a c++ function through R, 
	extra functions encoded in c language are required to achieve a bridge between
	c++ and R. The basic strategy I followed is like this: first, the R function calls
	the c function and pass R's parameters to the c function; then the c function 
	constructs all parameters required by the c++ function and call the function. 
	Results of the c++ function are recorded in all pointers of the c function. 
	Finally the R function gets all results through the pointers. 

Two tools are provided by R to imply this process. ```Extern "C"{ }``` is to
tell the compile the following code is a c function. ```.C( )``` is used by R
to call a function encoded by c language. All c functions should title with 
```Extern "C"```. See the code for the usage. 

It is trivial to call a c function through R. But one should notice that it is 
required to explicitly convert parameters in R to data types which can be 
recognized by the c code. In R, I use ```as.datatype``` to convert. The reason
for this is that, unlike R, each data type of the parameter should be 
predetermined in a c-encoded function. Since I do not find a way to convert 
the data type `struct` in c, all c functions I wrote can only take primitive 
data types.

Calling a c++ function through a c function involves two problems. One is how to
call c++ function, which contains non-trivial parameter data types. The other is
how to store results after calling a c++ function. In the continuedfraction.cpp
file, the only parameter type that can not be directly converted into a c-encoded 
type is the class ```vector``` derived from standard c++ library. My strategy is 
to set two parameters in the c function, one parameter is the pointer to an 
array, the other parameter is the length of the array. Then I used these two 
parameters to construct the vector and call the c++ function. The drawback of 
this strategy is that the number of parameters in the c function could be 
extended by factor two. Using ```struct pointer``` should solve this problem, but
I do not know how to do it yet. 

Extra parameters are also needed in a c-encoded function in order to store 
returning results from a c++ function. The same as above, they will further 
extend the number of parameters of the c-encoded function . A more critical 
issue on storing results from the c++ function is how to allocate memory. One
option is to let the c function allocate memory, use pointers from parameters to
point the address of the memory, copy all results into assigned memory and pass
them to the R function. However, I feel this may lead to leak memory. Because R
should not be capable to free allocated memory as it was created by the c 
function. To avoid this, I use the R function to allocate enough memory before 
the R function passes its parameters to the c function. The drawback of this 
method is that it wastes memory.

##Difficulties
I am not clear which types of pointers in a c-encoded function can pass to R 
except some primary types. For example, a data type `size_t` in c++, I do not 
know which the type of pointer corresponds in c language. Currently I just use
`int*` type to point the size_t type. Memory allocation is also an issue. 
Whether c function or R function should manager memory still bothers me. The last
thing is to make interfaces written in the c code neat and clean. How to pass 
struct pointers to a R function is the key. But I would leave it after 
successfully building the R package. 

##Suggestion
It would be great if preseqR can be designed as a sub module of preseq. I guess 
it is possible. All c-encoded functions can be extracted into a separated file. 
A few changes can keep continuedfraction class away from smithlab.cpp. The main 
problem is that compiling continuedfraction class requires the gsl library. I
have already written a R function to repalce the original sampling funtion, 
which relies on gsl library. So preseqR should not need these functions any more.
I would suggest that create a different Makefile file corresponding to compile 
preseqR, coupling with adding some `macros` in the continuedfraction.cpp file, 
which can skip all functions that require gsl for compiling. It could make the 
structure dense and less duplicated codes. But it is more like a pesudo 
sub-module. PreseqR is still an independent sub module in preseq. 
Downsample in preseq is written through gsl, while downsample in preseqR is 
written through R libraries. 

##Documents of global.variables and functions in preseqR

#### global.variables

  1.  MAXLENGTH = 10000000: the size of each pre-allocated vector by R to store
      results from calling C++ function
  2.  MULTINOMIAL.SAMPLE.TIMES = 19: a number to define the number of random
      vectors to draw from a multinomial distribution. It also defines the number
	  of sampling times in each iteration in bootstrapping. 
  3.  MINOR.correction = 1e-1: a small positive number to avoid randomness when
      two double type numbers with almost identical values are compared.
  4.  BOOTSTRAP.factor = 0.1: a number to deine the efficiency of bootstrap. See
      `preseqR.bootstrap.complexity.curve` for details.

#### R-functions

  1.  read.hist(hist.file): a function to read a histogram file and output a 
	  count vector of the histogram. The hisogram must have two columns and no 
	  headers are allowed. In the first column, each number x represents a 
	  sequencing read/species occurs x times in the experiment. The corresponding
	  number y in the second column represents there are y reads/species, each 
	  of which occurs x times.
  2.  preseqR.calculate.continued.fraction(CF, x): a function to calculate the
	  value of a continued fraction given its coordinates. CF is a continued
	  fraction with `CF` attribute. `x` is a non-negtive real number. It calls 
	  a C-function ```c_calculate_continued_fraction``` to calculate the
	  value.
  3.  preseqR.extrapolate.distinct(hist.count, CF, start.size = NULL, step.size = NULL, 
	  max.size = NULL): a function to extrapolate through a continued fraction.
	  It returns a list containing various sample sizes and their extrapolation 
	  values. `hist.count` is a count vector of a histogram. The ith number x in
      the vector represents there are x reads/species occurring i times in the 
	  experiment. `CF` is a continued fraction. `start.size` is a starting 
	  sample size for extrapolation. `step.size` is a increasement size in each 
	  step for extrapolation. `max.size` is a upper bound for extrapolation. It
	  calls a C-function `c_extrapolate_distinct` to implement.
  4.  replace.sampling(n, hist.count): a function to re-sample histograms given 
	  a histogram. `n` is the number of histograms to draw. `hist.count` is a 
	  count vector of a histogram. It calls a built-in R-function `rmultinomial()`
	  to implement.
  5.  nonreplace.sampling(size, hist.count): a function to do without replacement 
	  sampling. It returns a sample vector composed by positive number. Each 
	  positive number is considered as a id for a distinct item. `size` is a 
	  non-negative integer giving the number of items to choose. `hist.count` is
      a count vector of a histogram. See above for detail. It calls a built-in 
	  R-function `sample()` to implement.
  6.  count.distinct(sample): a function to count the number of distinct items
	  given a sample result. `sample` is a positive number vector. Each distinct 
	  item is uniquely represented by a number. It builds a vector and set the
	  value to one whenever its index appears in the sample. Then `sum()` is
	  used to get a distinct number of items.
  7.  preseqR.interpolate.distinct(hist.count, ss): a function to interpolate
	  given a histogram. It tells how many distinct items if the sample size is
	  less than the initial experiment. The returning value is a list that 
	  contains various sample sizes and their interpolation values. It also 
	  contains a starting sample size, which is larger than the initial
	  experiment and cannot be determined by interpolation. `hist.count` is a 
	  count vector of a histogram. `ss` is a step size. It is a increasement 
	  size in each step for interpolation. For each givan sample size, it do
	  sampling without replacement. The histogram is used as an underlining
	  distribution. Then it calls `count.distinct()` on the sample result to 
	  to count the number of distinct items.
  8.  goodtoulmin.2x.extrap(hist.count): a function to check quality of a
	  given histogram. See Good, I. J., and G. H. Toulmin 1956. `hist.count` is
	  a count vector of a histogram.
  9.  preseqR.continued.fraction.estimate(hist, di = 0, mt = 100, ss = 1e6, 
	  mv = 1e10,  max.extrapolation = 1e10, step.adjust=TRUE): a function to
	  construct a continued fraction given a histogram. The continued fraction,
	  which is based on a capture-recapture model, is used to predict the number
	  of distinct items if additional experiment conducts. The function returns
	  a constructed continued fraction and estimated number of distinct items
	  given various sample sizes. `hist` could be either a histogram file or a 
	  count vector of a histogram. `di` is the diagonal value for a continued
	  fraction. `mt` is the maximum number of parameters allowed in a constructed
	  continued fraction. `ss` is the step size, which defines the distance
	  between two neighbors in a sample size vector. `mv` is the maximum value 
	  for testing a continued fraction. `max.extrapolation` is the upper bound
	  for extrapolation. `step.adjust` is a logic value. When it set TRUE, the 
	  step size of sampling points could adapt in order to correspond to the size
	  of an inital experiment. The function  calls a C-function 
	  `c_continued_fraction_estimate` to approximaate a continued fraction. 
	  `preseqR.interpolate.distinct` is used to count the number of distinct 
	  items when the sample size is less than the inital experiment. 
	  `preseqR.extrapolate.distinct` is used to predict the number of distinct
	  items when the sample size is larger than the size of the experiment. 
	  MINOR.correction is used to make sure for same step size and maximum 
	  extrapolation value, the number of extrapolation points should be same.
  10.  preseqR.bootstrap.complexity.curve(hist, bootstrap.times = 100, di = 0, 
	   mt = 100, ss = 1e6, mv = 1e10, max.extrapolation = 1e10, step.adjust=TRUE):
	   a function to estimate the number of distinct items given various sample
	   sizes. Bootstrap is used to improve estimation and builds a confident
	   interval. It returns approximated number of distinct items given various
	   sample size plus a confident interval for each estimation. All input 
	   parameters are same as parameters in `preseqR.continued.fraction.estimate()` 
	   except bootstrap.times, which defines the minimal successful estimation
	   times for bootstrapping. For each iteration, the function generates
	   n = MULTINOMIAL.SAMPLE.TIMES histograms. Then 
	   `preseqR.continued.fraction.estimate()` is called by each histogram to
	   make an estimateion. Because the implementation is in a vector way, setting
	   a proper value for MULTINOMIAL.SAMPLE.TIMES could sppedup the function. 
	   The function will stop under two situations. One is
	   that times of successful estimation reach a number defined by bootstrap.times.
	   The other is that total resampling times beyond a number, defined as 
	   ```bootstrap.times / BOOTSTRAP.factor```. 
  11.  preseqR.printout(file.prefix = NULL, ...): a function to write down 
	   results from other functions to a file. `file.prefix` defined the prefix 
	   name of created files. Results include returning results from 
	   ```preseqR.interpolation(), preseqR.extrapolate.distinct(), 
	      preseqR.continued.fraction.estimate(), 
		  preseqR.bootstrap.complexity.curve()```. It can also print out a 
	   continued fraction.

#### C-functions

  1.  void c_calculate_continued_fraction( double *cf_coeffs, int *cf_coeffs_l,
	  double *offset_coeffs, int *di, int *de, double *coordinate, double *result): 
	  a function to calculate the value of a continued fraction given its 
	  coordinates. It is a C-encoded interface for C++ function `operator()` in
	  `continued_fraction.cpp` and `preseqR.calculate.continued.fraction` calls 
	  it to calculate the value. `cf_coeffs` stores coefficients of a continued
	  fraction. `cf_coeffs_l` is the number of coefficients. `offset_coeffs` 
	  stores offset coefficients of a continued fraction if the diagonal value
	  of the continued fraction is nonzero. `di` is the diagonal value of the
	  continued fraction. `de` is the degree of the continued fraction.
	  `coordinate` is a positive real value. `result` stores calculated result.
	  First, the function construct a continued fraction through information
	  stored in ```cf_coeffs, cf_coeffs_l, offset_coeffs, di, de```. Then it
	  calls a C++ function `operator()` in `continued_fraction.cpp` to calculate
	  the value.

  2.  void c_extrapolate_distinct( double *cf_coeffs, int *cf_coeffs_l,
	  double *offset_coeffs, int *di, int *de, double *hist, int *hist_l,
	  double *start_size, double *step_size, double *max_size, double *estimate,
	  int *estimate_l): a function to extrapolate through a continued fraction.
	  It calculates extrapolation values given various sample size and returns
	  results to `estimate` and `estimate_l`. `start.size` is a starting sample 
	  size for extrapolation. `step.size` is a increasement size in each step 
	  for extrapolation. `max.size` is a upper bound for extrapolation. `estimate`
	  is used to store extrapolation values. `estimate_l` stores the number of
	  extrapolation points.	For each sample size, the function calls 
	  `c_calculate_continued_fraction` to calculate the extrapolation value.
	  Then it save results into `estimate` and `estimate_l`.
	  `preseqR.extrapolate.distinct` calls the function to calculate the
	  extrapolation values.

  3.  void c_continued_fraction_estimate(double *hist_count, int *hist_count_l,
	  int *di, int *mt, double *ss, double *mv, double *ps_coeffs, 
	  int *ps_coeffs_l, double *cf_coeffs, int *cf_coeffs_l, 
	  double *offset_coeffs, int *diagonal_idx, int *degree, int *is_valid)(): 
	  a function to estimate a continued fraction given a histogram. `hist_count`
	  stores a count vector of a histogram. `hist_count_l` is the number of
	  items in `hist_count`. `is_valid` is a logic value to indicate the 
	  validness of the constructed continued fraction. All the information of 
	  the constructed continued fraction is stores in `ps_coeffs, ps_coeffs_l, 
	  cf_coeffs, cf_coeffs_l, offset_coeffs, diagonal_idx, degree`. The function
	  calls a c++ function `optimal_cont_frac_distinct` in `continued_fraction.cpp`
	  for constructing a continued fraction. It is the C-interface for 
	  constructing a continued fraction.
