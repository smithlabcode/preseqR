## MAXLENGTH the allocate size for storing complexity curve
## MULTINOMIAL.SAMPLE.TIMES number of random vectors to draw each time
## MINOR.correction a very small number to correct comparison result between
## two double type numbers when precisions can bias the result
## BOOTSTRAP.factor the cut off ratio of success times / total bootstrap times
## BOOTSTRAP.times default number of times for bootstrap
MAXLENGTH = 10000000
MULTINOMIAL.SAMPLE.TIMES = 10
MINOR.correction = 1e-5
BOOTSTRAP.factor = 0.1
BOOTSTRAP.times = 100

## read a histogram file; return the histogram count vector
## count vector represent frequencies of indexes. For those indexes not showing
## in the histogram file, use zeros to represent their values
preseqR.read.hist <- function(hist.file)
{
	#the first column is the index of the histogram
	#the second column is the frequency of indexes
	hist.table = read.table(hist.file);
	hist.count = vector(mode = 'numeric', length = 0);
	pre.pos = 0;
	for (i in 1:length(hist.table[, 1]))
	{
		#Indexes of the histogram should be increasing order
		if (i > 1 && hist.table[i, 1] <= hist.table[i - 1, 1])
		{
			stop("The input histogram is not sorted in increasing order");
		}
		#set the frequency of each index; 
		#fill the gaps between indexes the zeros
		hist.count = c(hist.count, rep(0, hist.table[i, 1] - pre.pos));
		hist.count[hist.table[i, 1]] = hist.table[i, 2];
		pre.pos = hist.table[i, 1];
	}
	return(hist.count);
}


## calculate the value of the continued fraction CF given the coordinate x 
## call c-encoded function "c.calculate.continued.fraction()" through 
## preseqR.calculate.continue.fraction
preseqR.calculate.continued.fraction <- function(CF, x)
{
	out <- .C("c_calculate_continued_fraction", 
			  cf = as.double(CF$cf.coeffs), 
	   		  cf.l = as.integer(length(CF$cf.coeffs)),
			  off = as.double(CF$offset.coeffs), 
			  di = as.integer(CF$diagonal.idx), 
			  de = as.integer(CF$degree), 
			  coordinate = as.double(x), 
		  	  result = as.double(0));
	return(out$result);
}

## extrapolate given a histogram and a continued fraction
preseqR.extrapolate.distinct <- function(hist.count, CF, start.size = NULL,
		 step.size = NULL, max.size = 100)
{
	cf.coeffs = as.double(CF$cf.coeffs);
	cf.coeffs.l = as.integer(length(CF$cf.coeffs));
	offset.coeffs = as.double(CF$offset.coeffs);
	di = as.integer(CF$diagonal.idx);
	de = as.integer(CF$degree);
	hist.count = as.double(hist.count);
	# record the size of the sample based on the histogram count
	total.reads = 0.0;
	for (i in 1:length(hist.count))
		total.reads <- total.reads + i * as.integer(hist.count[i]);

	# the styles of the histogram count vector are different between R code 
	# and c++ code; The first line of the histogram is always [0  0] in c++ 
	# but the line is removed in R-encoded function
	hist.count = c(0, hist.count);
	hist.count.l = as.integer(length(hist.count));

	# set start.size and step.size if they are not defined by user
	if (is.null(start.size))
		start.size = total.reads;
	if (start.size > max.size)
	{
		write("start position has already beyond the maximum prediction",
			  stderr());
		return(NULL);
	}
	if (is.null(step.size))
		step.size = total.reads;

	# allocate memory to store extrapolation results
	# first "c.extrapolate.distinct" stores the observed number of distinct 
	# molecules into estimate, then it stores the extrapolation values
	# thus the allocated memory size is 1 plus the size of extrapolation values,
	# which is (max.size - start.size) / step.size) + 1
	extrap.size = as.integer((max.size - start.size) / step.size) + 1;
	out <- .C("c_extrapolate_distinct", cf.coeffs, cf.coeffs.l, offset.coeffs,
	   di, de, hist.count, hist.count.l, as.double(start.size), 
	   as.double(step.size), as.double(max.size), 
	   estimate = as.double(vector(mode = 'numeric', extrap.size + 1)), 
	   estimate.l = as.integer(0));
	return(out$estimate[ 1: out$estimate.l ]);
}
	
#do with/without replacement of random sampling given a histogram count vector
#size is a user defined sample size
preseqR.hist.sample <- function(hist.count, size, replace = NULL)
{
	if (replace == FALSE)
	{
		total.sample = 0;
		i = 1;
		# calculate total number of sample
		freq = 1:length(hist.count);
		total.sample = freq %*% hist.count;
		#construct a sample space X 
		distinct.sample = sum(hist.count);
		# identities for each distinct read
		ind = 1:as.integer(distinct.sample);
		# the size of each read in the library
		n = rep(freq, as.integer(hist.count))
		# the whole library represents by its indexes. If a read presents t
		# times in the library, its indexes presents t times in X
		X = rep(ind, n); 
		return(sample(X, size, replace = FALSE));
	}
	else if (replace == TRUE) {
		return(rmultinom(MULTINOMIAL.SAMPLE.TIMES, size, hist.count));
	} else {
		write("Specify the sampling methods(wit/without replacement)", stderr())
		return();
	}
}

## convert sample points into a count vector
## for different sampling methods, the form of sampled points are different
## thus the preprocesses are different between two methods.
preseqR.sample2hist.count <- function(sample.points, replace = NULL)
{
	# V is the correponding histogram
	if (replace == TRUE) {
		V = table(sample.points);
	    # index "0" corresponds to the number of molecules not being sampled. 
	    # remove this term if it exists
	    if (names(V)[1] == '0')
			V = V[-1]
	}
	else if (replace == FALSE) {
		V = table(table(sample.points));
	}
	# convert the histogram to its count vector
	freq = as.integer(names(V));
	hist.count = vector(mode = 'numeric', length = max(freq));
	hist.count[freq] = as.integer(V);
	return(hist.count);
}

# given a sample vector, the function counts the number of distinct molecules
count.distinct <- function(sample)
{
	max.value = max(sample);
	sample.table = vector(mode = "numeric", length = max.value);
	sample.table[sample] = 1;
	return(sum(sample.table));
}
# interpolate when the sample size is no more than the size of 
# the initial experiment
# return the interpolated estimations and the sample size beyond 
# the initial experiment
preseqR.interpolate.distinct <- function(hist.count, ss)
{
	# calculate total number of sample
	freq = 1:length(hist.count);
	total.sample = freq %*% hist.count;

	inital.distint = sum(hist.count);
	upper.limit = as.integer(total.sample)
	step = ss;
	sample = step;
	# l is the number of interpolation points
	l = as.integer(upper.limit / sample);
	# if the sample size is larger than the size of experiment, return NULL
	if (l == 0)
		return();
	# sample size vector
	x = sample * ( 1:l );
	# dimesion must be defined in order to use R apply
	dim(x) = length(x);
	# do sampling without replacement 
	s = apply(x, 1, function(x) preseqR.hist.sample(hist.count, x, replace=FALSE));
	# calculate the number of distinct reads based on each sample
	dim(s) = length(s);
	yield.estimates = sapply(s, count.distinct);
	# sample stores the starting sample size for extrapolation
	sample <- sample + sample * l
	out = list(yield.estimates, sample);
	names(out) = c("yield.estimates", "sample.size")
	return(out);
}

#check the goodness of the sample based on good Good & Toulmin's model
goodtoulmin.2x.extrap <- function(hist.count)
{
	two_fold_extrap = 0.0;
	for ( i in 1:length(hist.count) )
		two_fold_extrap <- two_fold_extrap + (-1.0)^(i + 1) * hist.count[i];
    return(two_fold_extrap);
}

## estimate a continued fraction given a the count vector of the histogram
## di = diagonal, mt = max_terms, ss = step_size, 
## mv = max_value for training
## step.adjust is an indicator for whether or not to adjust step.size
preseqR.continued.fraction.estimate <- function(hist, di = 0, mt = 100,
	   	ss = 1e6, mv = 1e10,  max.extrapolation = 1e10, step.adjust=TRUE)
{
	# input could be either histogram file or count vector of the histogram
	if (mode(hist) == "character") {
		hist.count = preseqR.read.hist(hist);
	}
	else {
		hist.count = hist;
	}

	# minimum required number of terms of power series in order to construct
	# continued fraction
	MIN_REQUIRED_TERMS = 4
	# calculate total number of sample
	freq = 1:length(hist.count);
	total.sample = freq %*% hist.count;
	step.size = ss;
	# no interpolation if step.size is larger than the size of experiment
	# set the starting sample size as the step.size
	if (step.size > total.sample)
	{
		yield.estimates = vector(mode = 'numeric', length = 0);
		starting.size = step.size;
	}
	# interpolation when sample size is no more than total sample size
	else 
	{
		# adjust step.size when it is too small
		if (step.adjust == TRUE && step.size < (total.sample / 20))
		{
			step.size = max(step.size, step.size*round(total.sample/(20*step.size)));
			# output the adjusted step size to stderr
			m = paste("adjust step size to", toString(step.size), '\n', sep = ' ');
			write(m, stderr());
		}
		# interpolate and set the size of sample for initial extrapolation
		out = preseqR.interpolate.distinct(hist.count, step.size)
		yield.estimates = out$yield.estimates;
		starting.size = out$sample.size
	}

	counts.before.first.zero = 1;
	while (as.integer(counts.before.first.zero) <= length(hist.count) && 
		   			  hist.count[counts.before.first.zero] != 0)
		counts.before.first.zero <- counts.before.first.zero + 1;

	# starting sample size for extrapolation

	# continued fraction with even degree conservatively estimates
	mt = min(mt, counts.before.first.zero - 1);
	mt = mt - (mt %% 2);

	# pre-check to make sure the sample is good for prediction
	if (mt < MIN_REQUIRED_TERMS)
	{
		m = paste("max count before zero is les than min required count (4), ",
				  "sample not sufficiently deep or duplicates removed",
				  sep = '')
		write(m, stderr());
		return();
	}
	if(goodtoulmin.2x.extrap(hist.count) < 0.0)
	{
		m = paste("Library expected to saturate in doubling of size, ",
				  "unable to extrapolate", sep = '');
		write(m, stderr());
		return();
	}

	# adjust the format of count vector of the histogram in order to
	# call c-encoded function
	hist.count = c(0, hist.count);
	# allocate spaces to store constructed continued fraction
	# construct a continued fraction with minimum degree
	out <- .C('c_continued_fraction_estimate', as.double(hist.count), 
			  as.integer(length(hist.count)), as.integer(di), as.integer(mt), 
			  as.double(step.size), as.double(mv), 
			  ps.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
			  ps.coeffs.l = as.integer(0), 
			  cf.coeffs = as.double(vector(mode = 'numeric', length=MAXLENGTH)),
			  cf.coeffs.l = as.integer(0), 
			  offset.coeffs =as.double(vector(mode='numeric',length=MAXLENGTH)),
			  diagonal.idx = as.integer(0), 
			  degree = as.integer(0),
			  is.valid = as.integer(0));
	if (!out$is.valid)
	{
		write("Fail to construct continued fraction", stderr());
		return();
	}
	length(out$ps.coeffs) = out$ps.coeffs.l;
	length(out$cf.coeffs) = out$cf.coeffs.l;
	length(out$offset.coeffs) = as.integer(abs(out$diagonal.idx));
	CF = list(out$ps.coeffs, out$cf.coeffs, out$offset.coeffs, out$diagonal.idx,
			  out$degree);
	names(CF) = c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx', 
				  'degree');
	# if the sample size is larger than max.extrapolation
	# stop extrapolation
	# prevent machinary precision from biasing comparison result
	if (starting.size > (max.extrapolation + MINOR.correction))
	{
		index = as.double(step.size) * (1:length(yield.estimates));
		yield.estimates = list(sample.size = index, yields = yield.estimates);
		result = list(continued.fraction = CF, yield.estimates = yield.estimates);
		return(result);
	}
	est <- preseqR.extrapolate.distinct( hist.count, CF, 
		       (starting.size - total.sample) / total.sample, step.size / total.sample, 
		       (max.extrapolation + MINOR.correction- total.sample) / total.sample);
	# est[1] is the number of the distinct molecules from experiments
	# est[-1] are extrapolation results
	est = est[-1]
	yield.estimates = c(yield.estimates, est);
	index = as.double(step.size) * (1: length(yield.estimates));
	yield.estimates = list(sample.size = index, yields = yield.estimates);
	result = list(continued.fraction = CF, yield.estimates = yield.estimates)
	return(result);
}

print.continuedfraction <- function(CF)
{
	print(CF$cf.coeffs);
}

## generate complexity curve through bootstrap the histogram
bootstrap.complex.curve <- function(hist, times = 100, di = 0, mt = 100,
									ss = 1e6, mv = 1e10, 
									max.extrapolation = 1e10, step.adjust=TRUE)
{
	# input could be either histogram file or count vector of the histogram
	if (mode(hist) == "character") {
		hist.count = preseqR.read.hist(hist);
	}
	else {
		hist.count = hist;
	}
	# calculate total number of sample
	freq = 1:length(hist.count);
	total.sample = freq %*% hist.count;
	# adjust step.size for sampling complexity curve
	step.size = ss;
	if (step.adjust == TRUE && step.size < (total.sample / 20))
	{
		step.size = max(step.size, step.size*round(total.sample/(20*step.size)));
		# output the adjusted step size to stderr
		m = paste("adjust step size to", toString(step.size), '\n', sep = ' ');
		write(m, stderr());
	}
	# calculate the distinct number of sample
	distinct.sample = sum(hist.count)
	# resampled vector count of a histogram
	re.hist.count = vector(mode = "numeric", length = length(hist.count));
	# index of nonzero items in hist.count
	index.hist.count = which(hist.count != 0);
	if (times == 1) {
		out <- preseqR.continued.fraction.estimate(hist.count, di, 
				                          mt, ss, mv, max.extrapolation);
		if (!is.null(out)) {
			return(out$yield.estimates);
		}
		else {
			return();
		}
	} else if (times > 1) {
		# the actually step.size preseqR.continued.fraction.estimate uses
		step.size = 0

		yield.estimates = vector(mode = "numeric", length = 0);
		for (i in 1:as.integer(times))
		{
			# do sampling with replacement 
			sample = preseqR.hist.sample(hist.count, as.integer(distinct.sample), 
								replace = TRUE);
			re.hist.count[index.hist.count] = sample;
			# build count vector of the histogram based on sampling results
			out <- preseqR.continued.fraction.estimate(re.hist.count, di, mt, ss,
				                                    	mv, max.extrapolation);
			if (!is.null(out))
			{
				step.size = out$step.size
				yield.estimates = cbind(yield.estimates, out$yield.estimates$yields);
			}
		}
		# check the number times of success beyond a threshond
		if (is.null(dim(yield.estimates)) || dim(yield.estimates)[2] < BOOTSTRAP.factor * times)
		{
			write("fail to bootstrap since the histogram is poor", stderr());
			return();
		}
		# the number of sampled points for complexity curve
		n = dim(yield.estimates)[1];
		# sample sizes
		index = as.double(step.size) * ( 1:n )
		# mean values are used as complexity curve
		mean = apply(yield.estimates, 1, mean)
		variance = apply(yield.estimates, 1, var)
		# 95% confident interval based on normal distribution
		left.interval = mean - qnorm(0.975) * sqrt(variance / n);
		right.interval = mean + qnorm(0.975) * sqrt(variance / n);
		yield.estimates = list(sample.size = index, yields = mean)
		result = list(yield.estimates, left.interval, right.interval);
		names(result) = c("yield.estimates", "LOWER_0.95CI", 
						  "UPPER_0.95CI")
		return(result);
	} else {
		write("the paramter times should be at least one", stderr());
		return();
	}
}

bootstrap.complexity.curve <- function(hist, bootstrap.times = 100, di = 0, 
									   mt = 100, ss = 1e6, mv = 1e10,
									   max.extrapolation = 1e10, step.adjust=TRUE)
{
	if (mode(hist) == 'character') {
		hist.count = preseqR.read.hist(hist);
	} else {
		hist.count = hist;
	}
	# calculate the total number of sample
	freq = 1:length(hist.count);
	total.sample = freq %*% hist.count;
	# calculate the distinct number of sample
	distinct.sample = sum(hist.count)
	# adjust step.size for sampling complexity curve
	step.size = ss;
	if (step.adjust == TRUE && step.size < (total.sample / 20))
	{
		step.size = max(step.size, step.size*round(total.sample/(20*step.size)));
		# output the adjusted step size to stderr
		m = paste("adjust step size to", toString(step.size), '\n', sep = ' ');
		write(m, stderr());
	}
	# record count vector of a resampled histogram
	re.hist.count = matrix(data = 0, nrow = length(hist.count), 
						   ncol = MULTINOMIAL.SAMPLE.TIMES)
	# nonzero indexes in the hist.count
	nonzero.index = which(hist.count != 0);
	# nonzero values in the hist.count
	nonzero.hist.count = hist.count[nonzero.index];
	# the number of resampling times
	counter = 0;
	yield.estimates = vector(mode = "numeric", length = 0);
	# upperbound of times of iterations for bootstrapping
	upper.limit = bootstrap.times / BOOTSTRAP.factor
	while (bootstrap.times > 0) {
		# do sampling with replacement
		resample = preseqR.hist.sample(nonzero.hist.count, as.integer(distinct.sample), 
									   replace = TRUE);
		# reconstruct count vectors of histograms
		re.hist.count[ nonzero.index, 1:MULTINOMIAL.SAMPLE.TIMES ] = resample;
		# make estimation for each histogram
		out = apply(re.hist.count, 2, function(x) preseqR.continued.fraction.estimate(
					  x, di, mt, step.size, mv, max.extrapolation + MINOR.correction, 
					  step.adjust=FALSE))
		# eliminate NULL items in results
		out[sapply(out, is.null)] <- NULL
		# extract yields estimation from each estimation result. 
		yields = sapply(out, function(x) x$yield.estimates$yields)
		# the number of successful estimation time
		if ( !is.null( dim(yields) ) )
		{
			success.times = dim(yields)[2];
			bootstrap.times <- bootstrap.times - success.times;
			yield.estimates = cbind(yield.estimates, yields);
		}
		counter <- counter + MULTINOMIAL.SAMPLE.TIMES;
		if (counter > upper.limit)
			break;
	}
	# boostrap succeeds 
	if (bootstrap.times <= 0) {
		# the number of sampled points for complexity curve
		n = dim(yield.estimates)[1];
		# sample sizes
		index = as.double(step.size) * ( 1:n )
		# mean values are used as complexity curve
		mean = apply(yield.estimates, 1, mean)
		variance = apply(yield.estimates, 1, var)
		# 95% confident interval based on normal distribution
		left.interval = mean - qnorm(0.975) * sqrt(variance / n);
		right.interval = mean + qnorm(0.975) * sqrt(variance / n);
		yield.estimates = list(sample.size = index, yields = mean)
		result = list(yield.estimates, left.interval, right.interval);
		names(result) = c("yield.estimates", "LOWER_0.95CI", 
						  "UPPER_0.95CI")
		return(result);
	} else {
		write("fail to bootstrap!", stderr());
		return();
	}
}


