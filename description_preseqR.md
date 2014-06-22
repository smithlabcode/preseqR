#Contents, Interface Design, Principles, Implementation,  Difficulties

##Contents
Three files are included: Rcontinuedfraction.hpp, Rcontinuedfraction.cpp, 
Rcontinuedfraction.R. Rcontinuedfraction.cpp is a source code containing 
all required functions written in c++. Rcontinuedfraction.hpp is its 
header. Rcontinuedfraction.R includes all functions written in R, which 
can call c++ functions, providing tools for capture and re-capture.

##Interface of the preseqR
Four functions are supported by preseqR for common users:

  1.  preseqR_continuedfraction_estimate(): given a histogram, preseqR produces 
	  a continued rational function (CRF);
  2.  preseqR_complex_curve(): given a histogram, preseqR produces a complexity 
      curve of the library;
  3.  preseqR_print_continuedfraction(): print out the construct continued 
	  rational function in a friendly way;
  4.  preseqR_calculate_ContinueFraction(): given a CRF and its coordinates, 
      preseqR calculates the function value.

Two functions are supported by preseqR for advanced users:

  1.  preseqR_extrapolate_distinct(): extrapolate distinct molecules given a CRF;
  2.  preseqR_interpolation(): interpolate distinct molecules given a histogram;

All the functions are in Rcontinuedfraction.R. They are written in R, but they 
may call c++ function during implementation.

Other functions in Rcontinuedfraction.R are helpers, which should be used by no 
one except other functions.

  1.  R_read_hist(): read a histogram file and output a histogram count vector;
  2.  preseqR_sample(): do with/without replacement sampling based on a histogram;
  3.  sample2hist_count(): convert a sample vector into a histogram count vector;

##Principles and implementation
1. In order to pack all functions and make an available R package, all 
platform-dependent codes should be removed or rewritten. For example codes call 
*gsl library* or *unistd.h*. I am not sure smithlab_cpp requires any extra 
libraries. So I just grabbed the piece of codes required by 
Rcontinuedfraction.cpp to make sure that Rcontinuedfraction.cpp can be compiled 
properly and no other plat-specific libraries are needed. The only reason I do 
it in this way is because it is simple and easy. We can discuss more on how to 
design. The sampling function based on *gsl library* has been rewritten in R 
language. Thus gsl is not a required library for continuedfraction class.

2. R functions could only call functions from the c code. To explain it from 
what I learn, R could only call the function written in c with the format 
```void * function(type1 *pointer, type2 *pointer, type3 *pointer, ...)```. 
and each type should be valid in c. In order to call a c++ function through R, 
extra functions encoded in c language are required to achieve a bridge between
c++ and R. The basic idea I am doing is like this: first, the R function calls
the c function and pass R's parameters to the c function; then the c function 
constructs all parameters required by the c++ function and call the function. 
Results of the c++ function are recorded in all pointers of the c function. 
Finally the R function gets all results through the pointers. 

Two tools are provided by R to imply this process. ```Extern "C"{ }``` is to
tell the compile the following code is a c function. ```.C( )``` is used by R
to call a function encoded by c language. All c functions should title with 
```Extern "C"```. See the code for the usage. 

It is trivial to call a c function through R. But one should notice that you 
still need to explicitly convert parameters in R to data types which can be 
recognized by the c code. In R, I use ```as.datatype``` to convert. The reason
for this is that, unlike R, each data type of the parameter should be 
predetermined in a c-encoded function. Since I do not find a way to convert 
the data type ```struct``` in c, all c functions I wrote can only take 
primitive data types.

Calling a c++ function through a c function involves two problems. One is how to
call c++ function, which contains non-trivial parameter data types. The other is
how to store results after calling a c++ function. In continuedfraction.cpp file
, the only parameter type that can not be directly converted into a c-encoded 
type is the class ```vector``` derived from standard c++ library. My strategy is 
to set two parameters in the c function, one parameter is the pointer to an 
array, the other parameter is the length of the array. Then I used these two 
parameters to construct the vector and call the c++ function. The drawback of 
this strategy is that the number of parameters in the c function could be 
unnecessarily extended. Using ```struct pointer``` should solve this problem, but
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
except some primary types. For example, a data type size_t in c++, I do not know
which the type of pointer corresponds in c language. Currently I just use int *
type to point the size_t type. Memory allocation is also an issue. Whether c 
function or R function should manager memory still bothers me. The last thing is
to make interfaces written in the c code neat and clean. How to pass struct 
pointers to a R function is the key. But I would leave it after successfully
building the R package. 

##Suggestion
It would be great if preseqR can be designed as a sub module of preseq. I guess 
it is possible. All c-encoded functions can be extracted into a separated file. 
A few changes can keep continuedfraction class from smithlab.cpp. The main 
problem is that compiling continuedfraction class requires the gsl library. I
have already written a R function to repalce the sampling funtion, which 
requires gsl library. So preseqR should not need these functions any more. I 
would suggest that create a different Makefile file corresponding to compile 
preseqR, coupling with adding some macros in the continuedfraction.cpp file, 
which can skip all functions that require gsl for compiling. It could make the 
structure dense and less duplicated codes. But it is more like a pesudo 
sub-module. ```preseqR``` is still an independent sub module in preseq. 
Downsample in preseq is written through gsl, while downsample in preseqR is 
written through R libraries. 
