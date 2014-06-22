#Contents, Interface Design, Principles, Implementation,  Difficulties

##Contents
Three files are included: Rcontinuedfraction.hpp, Rcontinuedfraction.cpp, 
Rcontinuedfraction.R. Rcontinuedfraction.cpp is the source code containing 
all required functions written in c++. Rcontinuedfraction.hpp is its 
header. Rcontinuedfraction.R includes all functions written in R, which 
could call c++ functions to archieve capture and re-capture.

##Interface of the preseqR
Four functions are supported by preseqR for common users:
1. preseqR_continuedfraction_estimate(): given a histogram, preseqR produces a 
continued rational function (CRF); 
2. preseqR_complex_curve(): given a histogram, preseqR produces a complexity 
curve of the library;
3. preseqR_print_continuedfraction(): print out the continued rational function 
in a friendly way such that user knows what it is; 
4. preseqR_calculate_ContinueFraction(): given a CRF and the coordinates, 
preseqR calculates the function values.

Two functions are supported by preseqR for advanced users:
1. preseqR_extrapolate_distinct(): extrapolate distinct molecules given a CRF;
2. preseqR_interpolation(): interpolate distinct molecules given a histogram;

All the functions are in Rcontinuedfraction.R. They are written in R, but they 
may call c++ function during implementation.

Other functions in Rcontinuedfraction.R are helpers, which should be used by no 
one except other functions.
1. R_read_hist(): read a histogram file and output a histogram count vector;
2. preseqR_sample(): do with/without replacement sampling based on a histogram;
3. sample2hist_count(): convert a sample vector into a histogram count vector;

##Principles and implementation
1. In order to pack all functions and make an available R package, all 
platform-dependent codes should be removed or rewritten. For example codes call 
gsl library or unistd.h. I am not sure smithlab_cpp requires any extra 
libraries. So I just grabbed the piece of codes required by 
Rcontinuedfraction.cpp to make sure that Rcontinuedfraction.cpp can be compiled 
properly and no other plat-specific libraries are needed. The only reason I do 
it this way is because it is simple and easy. We can discuss more on how to 
design. The sampling function based on gsl library is also rewritten in R 
language. Thus gsl is not a required library for continuedfraction class.

2. R functions could only call functions in c code. To explain it from what I
learn, R could only call function written in c with the format 
```void * function(type1 *pointer, type2 *pointer, type3 *pointer, ...)```. 
and each type should be valid in c. In order to call c++ function through R, 
extra c code functions are required to achieve a bridge between c++ and R. The 
basic idea I am doing is like this: first, R function calls c function and pass
R's parameters to c function; then c function constructs all parameters required
by c++ function and call c++ function. The results of c++ function recorded in 
all pointers of c function. Finally R function gets all results through the 
pointers. Two tools are provided by R to imply this process. *Extern "C"* is to 
tell the compile the following code is a c function. *.C* is used by R to call
c function. All c function should title with *Extern "C"*. See the code for the
usage. It is trivial to call c function through R. But one thing should be 
emphasized that you still need to explicitly convert the data types in R to the
data type which can be recognized by c code. In R, I use *as.datatype* to 
convert. Unlike R, for c function, each data type of the parameter should be
predetermined. Since I do not find a way to convert the data type to struct in 
c, all c functions I wrote can only take primitive data types.  Calling c++ 
function through c function contains two problems. One is how to call c++ 
function, which has non-trivial parameter data type. The other is  how to store 
results from c++ function. In continuedfraction.cpp file, the only parameter
type that can not be directly converted into c code is the class vector from
standard c++ library. My strategy is to set two parameters in c function, one
parameter is the pointer to an array, the other parameter is the length of the 
array. Then I used these two parameters to construct the vector and call the c++
function. The drawback of the strategy is that parameters in c function could be
unnecessarily extended. Use struct pointer should solve this problem, but I do 
not know how to do it yet. Extra parameters are also needed in c code in order
to store returning results from c++ function. The same as above, they will 
further extend the parameters of c function. More critical issue on storing
results from c++ function is how to allocate memory. One option is to let c
function allocate memory, use pointers in parameters to point the address of the
memory, copy all results into assigned memory and pass it to R function. 
However, I feel this may lead to leak memory. Because R should not be capable 
free allocated memory as they were created by c function. To avoid this, I use R
function to allocate enough memory before it passes the parameters to c 
function. The drawback of this method is that it wastes memory.

##Difficulties
I am not clear which types of pointers in c function could pass to R except some
primary types. For example, data type size_t c++, I do not know which type of 
pointer I should use in c. Currently I just use int * for pointing size_t. 
Memory allocation is also an issue. Whether c function or R function should 
manager memory still bothers me. The last thing is to make the interfaces 
written in c code neat and clean. How to pass struct pointers to R function is 
the key. But I would leave it after successfully building the R package. 

##Suggestion
It would be great if preseqR can be designed as a sub module of preseq. I guess 
it is possible. All c functions can be extracted into a separated file. A few 
changes can make continuedfraction class independent of smithlab.cpp. The main 
problem is that compiling continuedfraction class requires gsl library. I have 
already written R function to repalce sampling compiled by gsl. So preseqR do 
not need these functions any more. I would suggest create a different Makefile 
file corresponding to compile preseqR, coupling by adding some macros in 
continuedfraction.cpp, which can skip all functions require gsl. It could make 
the structure dense and less duplicated code produced. But in this way preseqR 
is still an independent sub module in preseq. Downsample in preseq is written 
through gsl, while downsample in preseqR is written through R libraries. 

