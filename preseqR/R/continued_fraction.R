## MAXLENGTH the allocate size for storing complexity curve
## MULTINOMIAL.SAMPLE.TIMES number of random vectors to draw each time
## MINOR.correction a very small number to correct comparison result between
## two double type numbers when precisions can bias the result
## BOOTSTRAP.factor the cut off ratio of success times / total bootstrap times
## BOOTSTRAP.times default number of times for bootstrap
MAXLENGTH = 10000000
MULTINOMIAL.SAMPLE.TIMES = 19
MINOR.correction = 1e-1
BOOTSTRAP.factor = 0.1

## read a histogram file, a two-column table, or a count vector of the histogram
## return the histogram count vector
## count vector represent frequencies of indexes. For those indexes not showing
## in the histogram file, use zeros to represent their values
read.hist <- function(hist.file, header = FALSE)
{
	if (class(hist.file) == "character") {
		hist.table = read.table(hist.file, header = header);
	} else if (!is.null(ncol(hist.file)) && ncol(hist.file) == 2) {
		if (is.numeric(hist.file[, 1]) && is.numeric(hist.file[, 2])) {
			hist.table = hist.file;
		} else {
			stop("All items in the input should be numeric values");
		}
	} else {
		stop("input should be a variable or a file of a two-column table")
	}
	# the first column is the frequencies of observed items
	freq = hist.table[, 1];
	# the second column is the number of observed items for each frequency
	number.items = hist.table[, 2];
	# check whether frequencies are at least one and the histogram is sorted
	for (i in 1:length(freq))
		if (freq[i] <= 0 || freq[i] != floor(freq[i]) || number.items[i] < 0) {
			stop("frequencies should not be positive integers!")
		} else {
			if (i > 1 && freq[i - 1] >= freq[i])
				stop("The input histogram is not sorted in increasing order");
		}
	# hist.count is the count vector of the histogram
	hist.count = vector(mode = 'numeric', length = max(freq));
	hist.count[freq] = number.items;
	return(hist.count);
}

## calculate the value of the continued fraction CF given the coordinate x 
## call c-encoded function "c.calculate.continued.fraction()" through 
## preseqR.calculate.continue.fraction
preseqR.calculate.continued.fraction <- function(CF, x)
{
	if (class(CF) != "CF")
		return();
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
		 step.size = NULL, max.size = NULL)
{
	if (class(CF) != "CF")
		return();
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

	# set start.size, step.size, max.size if they are not defined by user
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
	if (is.null(max.size))
	# 1000 is a magic number; 
		max.size = 1000 * total.reads;

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

	extrapolation = out$estimate[ 1: out$estimate.l ]
	# sample size vector for extrapolation
	sample.size = start.size + step.size * ( 1:length(extrapolation) - 1 );
	# put sample.size and extrapolation results together into a matrix
	result = matrix(c(sample.size, extrapolation), ncol = 2, byrow = FALSE);
	colnames(result) = c('sample.size', 'extrapolation');
	return(result);
}

## do withouat replacement of random sampling given a count vector of a histogram
## size is a user defined sample size
nonreplace.sampling <- function(size, hist.count)
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

## the function samples n histograms given a histogram
## it is based on sampling with replacement (the multinomial distribution)
replace.sampling <- function(n, hist.count)
{
	# record count vector of a resampled histogram
	re.hist.count = matrix(data = 0, nrow = length(hist.count), ncol = n) 
	# nonzero indexes in the hist.count
	nonzero.index = which(hist.count != 0);
	# nonzero values in the hist.count
	nonzero.hist.count = hist.count[nonzero.index];
	# calculate the distinct number of sample
	distinct.sample = sum(hist.count);
	# do sampling with replacement
	resample = rmultinom(n, as.integer(distinct.sample), hist.count);
	# reconstruct count vectors of histograms
	re.hist.count[ nonzero.index, 1:n ] = resample;
	return(re.hist.count);
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
	s = lapply(x, function(x) nonreplace.sampling(x, hist.count));
	# calculate the number of distinct reads based on each sample
	dim(s) = length(s)
	yield.estimates = sapply(s, function(x) count.distinct(x));
    # put sample.size and yield.estimates together into a matrix
	result = matrix(c(x, yield.estimates), ncol = 2, byrow = FALSE);
	colnames(result) = c('sample.size', 'interpolation')
	return(result);
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
## di = diagonal, mt = max_terms, 
## step.adjust is an indicator for whether or not to adjust step.size
preseqR.continued.fraction.estimate <- function(hist, di = 0, mt = 100, 
		ss = NULL,  max.extrapolation = NULL, step.adjust=TRUE, header = FALSE)
{
	hist.count = read.hist(hist, header);
	# minimum required number of terms of power series in order to construct
	# continued fraction
	MIN_REQUIRED_TERMS = 4
	# calculate total number of sample
	freq = 1:length(hist.count);
	total.sample = freq %*% hist.count;
	if (is.null(ss)) {
		ss = floor(total.sample);
		step.size = ss;
	} else {
		step.size = ss;
	}
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
		out = preseqR.interpolate.distinct(hist.count, step.size);
		yield.estimates = out[, 2];
		starting.size = (as.integer(total.sample / step.size) + 1) * step.size;
	}
	if (is.null(max.extrapolation)) {
		# extrapolation 100 times
		max.extrapolation = 100 * step.size;
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
	class(CF) = 'CF'
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
		       (max.extrapolation+MINOR.correction-total.sample) / total.sample);

	yield.estimates = c(yield.estimates, est[, 2]);
	index = as.double(step.size) * (1: length(yield.estimates));
	# put index and estimated yields together into a two-colunm matrix
	yield.estimates = matrix(c(index, yield.estimates), ncol = 2, byrow = FALSE);
	colnames(yield.estimates) = c('sample.size', 'yield.estimates');
	result = list(continued.fraction = CF, yield.estimates = yield.estimates)
	return(result);
}

## generate complexity curve through bootstrapping the histogram
preseqR.bootstrap.complexity.curve <- function(hist, bootstrap.times = 100, 
		di = 0, mt = 100, ss = NULL, max.extrapolation = NULL, step.adjust=TRUE,
	   	header = FALSE)
{
	hist.count = read.hist(hist, header);
	# calculate the total number of sample
	freq = 1:length(hist.count);
	total.sample = freq %*% hist.count;
	# calculate the distinct number of sample
	distinct.sample = sum(hist.count)
	if (is.null(ss)) {
		ss = floor(total.sample);
		step.size = ss;
	} else {
		step.size = ss;
	}
	# adjust step.size for sampling complexity curve
	if (step.adjust == TRUE && step.size < (total.sample / 20))
	{
		step.size = max(step.size, step.size*round(total.sample/(20*step.size)));
		# output the adjusted step size to stderr
		m = paste("adjust step size to", toString(step.size), '\n', sep = ' ');
		write(m, stderr());
	}
	if (is.null(max.extrapolation)) {
		# extrapolation 100 times
		max.extrapolation = 100 * step.size;
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
	f <- function(x)
	{
		# convert a count vector of a histogram into a histogram table
		freq = which(x != 0);
		number.items = x[freq];
		hist.table = matrix(c(freq, number.items), ncol = 2, byrow = FALSE);
		preseqR.continued.fraction.estimate(hist.table, di, mt, step.size, 
				max.extrapolation, step.adjust = FALSE);
	}

	while (bootstrap.times > 0) {
		# do sampling with replacement
		# the piece of code achives the same function as nonreplace.sampling()
		# I extend the code of the function for speedup
		# the cleaness of code is sacrificed to trade off efficiency
		resample = rmultinom(MULTINOMIAL.SAMPLE.TIMES, 
							 as.integer(distinct.sample), nonzero.hist.count);
		# reconstruct count vectors of histograms
		re.hist.count[ nonzero.index, 1:MULTINOMIAL.SAMPLE.TIMES ] = resample;
		# make estimation for each histogram
		out = apply(re.hist.count, 2, f)
		# eliminate NULL items in results
		out[sapply(out, is.null)] <- NULL
		# extract yields estimation from each estimation result. 
		yields = sapply(out, function(x) x$yield.estimates[, 2])
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
		index = as.double(step.size) * ( 1:n );
		# mean values are used as complexity curve
		median.estimate = apply(yield.estimates, 1, median);
		variance = apply(yield.estimates, 1, var);

		# 95% confident interval based on lognormal distribution
		C = exp(qnorm(0.975) * sqrt(log(1.0 + variance / (median.estimate^2))))
		left.interval = median.estimate / C;
		right.interval = median.estimate * C;

		result = matrix(c(index, median.estimate, left.interval, right.interval), ncol = 4,
				byrow = FALSE)
		colnames(result) = c('sample.size', 'yield.estimates','lower.0.95CI','upper.0.95CI')
		return(result);
	} else {
		write("fail to bootstrap!", stderr());
		return();
	}
}

print.continued.fraction <- function(X, filename)
{
	# use the variable name as the name of the continued fraction
	s = paste("CONTINUED FRACTION ", deparse(substitute(X)), ':\n', sep = '');
	write(s, filename);
	# print the degree of the continued fraction
	s = paste("DEGREE", toString(X$degree), "\n", sep = '\t');
	write(s, filename, append = TRUE);
	# print the diagonal value
	s = paste("DIAGONAL VALUE", toString(X$diagonal.idx), "\n", sep = '\t');
	write(s, filename, append = TRUE);
	write("COEFFICIENTS:", filename, append = TRUE)
	di = abs(X$diagonal.idx); 
	# print offset values if any
	if (di > 0)
	{
		index = 1:di;
		dim(index) = di;
		s = apply(index, 1, function(x)  paste('a_', toString(x - 1), ' = ', 
				               toString(X$offset.coeffs[x]), sep = ''));
		write(s, filename, append = TRUE);
	}
	# print coeffients if any
	if (length(X$cf.coeffs) > 0)
	{
		index = 1:length(X$cf.coeffs);
		dim(index) = length(X$cf.coeffs);
		s = apply(index, 1,function(x) paste('a_', toString(x - 1 + di),
					              ' = ', toString(X$cf.coeffs[x]), sep = ''));
		write(s, filename, append = TRUE);
	}
}

print.yield.estimates <- function(X, filename, digit = 0)
{
	if (!is.null(colnames(X))) {
		s = formatC(toupper( colnames(X) ), width = 15, format = 's', flag = '-')
		write(paste(s, collapse = ''), filename);
	}
	if (!is.null( dim(X) )) {
		f <- function(x) {
			s = formatC(x, digit, width = 15, format = 'f', flag = '-')
			write(paste(s, collapse = ''), filename, append = TRUE);
		}
		apply(X, 1, function(x) f(x));
	}
}

preseqR.print2file <- function(X, prefix = '', digit = 0)
{
	if ( !is.null(class(X)) ) {
		# check if X is a continued fraction
		if ( class(X) == "CF" ) {
			filename = paste(prefix, "_continued_fraction.txt", sep = '');
			print.continued.fraction(X, filename);
			return(0)
		  # check if X is a matrix
		} else if (class(X) == "matrix") {
			filename.YE = paste(prefix, "_yield_estimates.txt", sep = '');
			print.yield.estimates(X, filename.YE, digit);
			return(0)
		} else if (class(X) == "list") {
			# check if X is a result from preseqR.continued.fraction.estimate
			if (!is.null( names(X) ) && length(names(X)) == 2 && 
					all(names(X) == c("continued.fraction", "yield.estimates")))
			{	
				filename.CF = paste(prefix, "_continued_fraction.txt", sep = '');
				filename.YE = paste(prefix, "_yield_estimates.txt", sep = '');
				print.continued.fraction(X$continued.fraction, filename.CF);
				print.yield.estimates(X$yield.estimates, filename.YE, digit);
				return(0)
			}
		}
	}

	# invalid parameter X
	write("unknown input variables!", stderr());
	return(1);
}
