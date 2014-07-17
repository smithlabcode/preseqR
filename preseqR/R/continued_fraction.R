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

## read a histogram file; return the histogram count vector
## count vector represent frequencies of indexes. For those indexes not showing
## in the histogram file, use zeros to represent their values
read.hist <- function(hist.file)
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
	return(list(sample.size = sample.size, yield.estimates = extrapolation))
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
	s = apply(x, 1, function(x) nonreplace.sampling(x, hist.count));
	# calculate the number of distinct reads based on each sample
	dim(s) = length(s);
	yield.estimates = sapply(s, count.distinct);
    # yield.estimates
	yield.estimates = list(sample.size = x, yields = yield.estimates);
	return(yield.estimates);
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
		hist.count = read.hist(hist);
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
		out = preseqR.interpolate.distinct(hist.count, step.size);
		yield.estimates = out$yields;
		starting.size = (as.integer(total.sample / step.size) + 1) * step.size;
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

	yield.estimates = c(yield.estimates, est$yield.estimates);
	index = as.double(step.size) * (1: length(yield.estimates));
	yield.estimates = list(sample.size = index, yields = yield.estimates);
	result = list(continued.fraction = CF, yield.estimates = yield.estimates)
	return(result);
}

## generate complexity curve through bootstrapping the histogram
preseqR.bootstrap.complexity.curve <- function(hist, bootstrap.times = 100, di = 0, 
									   mt = 100, ss = 1e6, mv = 1e10,
									   max.extrapolation = 1e10, step.adjust=TRUE)
{
	if (mode(hist) == 'character') {
		hist.count = read.hist(hist);
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
		# the piece of code achives the same function as nonreplace.sampling()
		# I extend the code of the function for speedup
		# the cleaness of code is sacrificed to trade off efficiency
		resample = rmultinom(MULTINOMIAL.SAMPLE.TIMES, 
							 as.integer(distinct.sample), nonzero.hist.count);
		# reconstruct count vectors of histograms
		re.hist.count[ nonzero.index, 1:MULTINOMIAL.SAMPLE.TIMES ] = resample;
		# make estimation for each histogram
		out = apply(re.hist.count, 2, function(x) preseqR.continued.fraction.estimate(
					  x, di, mt, step.size, mv, max.extrapolation, step.adjust=FALSE))
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
		# successful resampling times
		resampling.n = dim(yield.estimates)[2];
		# sample sizes
		index = as.double(step.size) * ( 1:n );
		# mean values are used as complexity curve
		mean = apply(yield.estimates, 1, mean);
		variance = apply(yield.estimates, 1, var);
		# 95% confident interval based on normal distribution
		left.interval = mean - qnorm(0.975) * sqrt(variance / resampling.n);
		right.interval = mean + qnorm(0.975) * sqrt(variance / resampling.n);
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

## density function of truncated zero negative binomial distribution
preseqR.zerotruncated.dnbinom <- function(x, size, mu, log = FALSE)
{
	# the density of x in negative binomial
	p = dnbinom(x, size = size, mu = mu, log = log);
	# set zeros in x with zero probability
	p[ which(x == 0) ] = 0;
	# the density of non-zero in negative binomial
	q = 1 - dnbinom(0, size = size, mu = mu, log = log);
	# normalize all non-zero values in negrative binomial to generate ZTNB
	if (log == FALSE) {
		return(p / q);
	} else {
		return(log(p / q));
	}
}

## negative loglikelyhood 
zerotruncated.minus.log.likelyhood <- function(x, size, mu)
{
	prob = zerotruncated.dnbinom(1:length(x), size, mu, log = TRUE);
	# minus log likelyhood
	prob = -prob;
	# negative loglikelyhood
	return( x %*% prob)
}

## MLE
preseqR.zerotruncated.mle <- function(hist, size.init = NULL, 
									  mu.init = NULL)
{
	if (mode(hist) == 'character') {
		hist.count = read.hist(hist);
	} else {
		hist.count = hist;
	}
	if (mu.init == NULL)
	{
		total.sample = (1:length(hist.count) %*% hist.count);
		distinct.sample = sum(hist.count);
		mu.init = total.sample / distinct.sample;
	}
	if (size.init == NULL)
		size.init = total.sample;
	f <- function(size, mu) zerotruncated.minus.log.likelyhood(hist.count,size,mu);
	return(optim(c(size.init, mu.init), f, NULL, method = "L-BFGS-B", 
				lower = c(0, 0)), upper = c(1e10, 1e10))
}

## predict the number of distinct items using zero truncated negative binomial 
## distribution
## n is the size of experiment
preseqR.zerotruncated.estimate <- function(hist.count, n)
{
	total.sample = (1:length(hist.count) %*% hist.count);
	distinct.sample = sum(hist.count);

	opt <- preseqR.zerotruncated.mle(hist.count);
	size = opt$par[1];
	mu = opt$par[2];
	# the probability of being sampled in the initial experiment
	p = 1 - dnbinom(0, size = size, mu = mu);
	# L is the estimated total number of distinct items
	L = distinct / p;
	# the parameters of negative binomial in the experiment with size n
	size.new = size;
	mu.new = mu * as.double(n) / total.sample;
	# the probability of being sampled under the new experiment
	P = 1 - dnbinom(0, size = size.new, mu = mu.new, log);
	# the expected number of distinct items under the new experiment
	return(L * P);
}

preseqR.zerotruncated.complexity.curve <- function(hist, ss = 1e6, 
												max.extrapolation = 1e10)
{
	if (mode(hist) == 'character') {
		hist.count = read.hist(hist);
	} else {
		hist.count = hist;
	}
	# n is the number of experiments
	n = as.integer(max.extrapolation / ss);
	sample.size = as.double(ss) * (1: n);
	dim(sample.size) = n;
	return(apply(sample.size, 1, function(x) preseqR.zerotruncated.estimate(hist.count, x)));
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

print.yield.estimates <- function(X, filename)
{
	if (!is.null(names(X)))
		write(toupper(paste(names(X), collapse = '\t')), filename);
	l = length(X)
	if (l > 0)
	{
		X = matrix(unlist(X), ncol = l, byrow = FALSE);
		apply(X, 1, function(x) write(paste(x, collapse = '\t'), 
								filename, append = TRUE));
	}
}

preseqR.print2file <- function(X, prefix = '')
{
	# check if X is a continued fraction
	if (class(X) == "CF") {
		filename = paste(prefix, "_continued_fraction.txt", sep = '');
		print.continued.fraction(X, filename);
	} else if (!is.null(names(X))) {
		# check if X is the result from preseqR.continued.fraction.estimate
		if ( length(names(X)) == 2 && 
				all( names(X) == c("continued.fraction", "yield.estimates") ) )
		{
			filename.CF = paste(prefix, "_continued_fraction.txt", sep = '');
			filename.YE = paste(prefix, "_yield.estimates.txt", sep = '');
			print.continued.fraction(X$continued.fraction, filename.CF);
			print.yield.estimates(X$yield.estimates, filename.YE);
		  # check if X is the result from preseqR.bootstrap.complexity.curve
		} else if( length(names(X)) == 3 && 
					all( names(X) == c("yield.estimates", "LOWER_0.95CI", "UPPER_0.95CI") ) )
		{
			X = list(SAMPLE.SIZE = X$yield.estimates$sample.size,
					 YIELDS = X$yield.estimates$yields,
					 LOWER_0.95CI = X$LOWER_0.95CI,
					 UPPER_0.95CI = X$UPPER_0.95CI)
			filename.YE = paste(prefix, "_yield.estimates.txt", sep = '');
			print.yield.estimates(X, filename.YE);
		}
	} else
	{
		# invalid parameter X
		write("unknown input variables!", stderr());
	}
}
