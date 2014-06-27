## MAXLENGTH the allocate size for storing complexity curve
## BOOTSTRAP.factor the cut off ratio of success times / total bootstrap times
## default number of times for bootstrap
MAXLENGTH = 10000000
BOOTSTRAP.factor = 0.7
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
preseqR.hist.sample <- function(hist.count, size, replace)
{
	if (replace == FALSE)
	{
		total.reads = 0;
		i = 1;
		while (i <= length(hist.count))
		{
			total.reads = i * as.integer(hist.count[i]);
			i <- i + 1;
		}
		#construct a sample space X 
		X = vector(mode = 'numeric', length = as.integer(total.reads));
		pos = 1;
		p = 1;
		value = 1;
		for (l in hist.count)
		{
			if (as.integer(l) > 0)
			{
				X[pos: (pos + as.integer(l) * p - 1)] <- 
					rep(value: (value + as.integer(l) - 1), p);
				value <- value + as.integer(l);
				pos <- pos + as.integer(l) * p;
			}
			p <- p + 1;
		}
		return(sample(X, size, replace = FALSE));
	}
	else if (replace == TRUE) {
		distinct.reads = as.integer(sum(hist.count));
		#construct the pdf of the multinomial distribution
		prob = vector(mode = 'numeric', length = distinct.reads);
		p = 1;
		pos = 1;
		for (l in hist.count)
		{
			if (as.integer(l) > 0)
			{
				prob[pos: (pos + as.integer(l) - 1)] = p;
				pos <- pos + as.integer(l);
			}
			p <- p + 1;
		}
		return(rmultinom(1, size, prob));
	} else {
		write("Specify the sampling methods(wit/without replacement)", stderr())
		return();
	}
}

## convert points from without replacement sampling into a count vector of the 
## histogram 
preseqR.nonreplace.sample2hist.count <- function(sample.points)
{
	V = table(sample.points);
	value = as.vector(V);
	hist.count = vector(mode = 'numeric', length = max(value));
	for (i in value)
	{
		hist.count[i] <- hist.count[i] + 1;
	}
	return(hist.count)
}

## convert points from replacement sampling into a count vector of the histogram
preseqR.replace.sample2hist.count <- function(sample.points)
{
	hist.count = vector(mode = 'numeric', length = max(sample.points));
	for (v in sample.points)
	{
		if (v != 0) hist.count[v] <- hist.count[v] + 1
	}
	return(hist.count)
}

## combine two sample2hist methods into one; repalce represents points are 
## sampled through replacement sampling or non replacement sampling
preseqR.sample2hist.count <- function(sample.points, replace)
{
	if (replace == TRUE) {
		preseqR.replace.sample2hist.count(sample.points);
	} else {
		preseqR.nonreplace.sample2hist.count(sample.points);
	}
}
# interpolate when the sample size is no more than the size of 
# the initial experiment
# return the interpolated estimations and the sample size beyond 
# the initial experiment
preseqR.interpolate.distinct <- function(hist.count, ss)
{
	total.sample = 0.0;
	for (i in 1:length(hist.count))
		total.sample <- total.sample + i * hist.count[i];
	inital.distint = sum(hist.count);
	upper.limit = as.integer(total.sample)
	step = ss;
	sample = step;
	# l is the number of interpolation points
	l = as.integer(upper.limit / sample);
	# if the sample size is larger than the size of experiment, return NULL
	if (l == 0)
		return();

	yield.estimates = as.double(vector(mode = 'numeric', length = l));
	for (i in 1:l)
	{
		s = preseqR.hist.sample(hist.count, sample, replace = FALSE);
		yield = sum(preseqR.sample2hist.count(s, replace = FALSE));
		yield.estimates[i] = yield;
		sample <- sample + step;
	}
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
preseqR.continued.fraction.estimate <- function(hist, di = 0, mt = 100,
	   	ss = 1e6, mv = 1e10,  max.extrapolation = 1e10)
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
	total.sample = 0.0;
	for (i in 1:length(hist.count))
		total.sample <- total.sample + i * hist.count[i];
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
		if (step.size < (total.sample / 20))
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
			  step.size = as.double(step.size), as.double(mv), 
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
		write("Fail to construct and need to bootstrap to obtain estimates", 
			  stderr());
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
	if (starting.size > max.extrapolation)
	{
		index = as.double(ss) * (1:length(yield.estimates));
		yield.estimates = list(sample.size = index, yields = yield.estimates);
		result = list(CF, yield.estimates, out$step.size);
		names(result) = c("continued.fraction","yield.estimates", "step.size");
		return(result);
	}
	est <- preseqR.extrapolate.distinct(
				hist.count, CF, (starting.size - total.sample) / total.sample,
				out$step.size / total.sample, 
				(max.extrapolation - total.sample) / total.sample);
	# est[1] is the number of the distinct molecules from experiments
	# est[-1] are extrapolation results
	est = est[-1]
	yield.estimates = c(yield.estimates, est);
	index = as.double(out$step.size) * (1: length(yield.estimates));
	yield.estimates = list(sample.size = index, yields = yield.estimates);
	result = list(CF, yield.estimates, out$step.size)
	names(result) = c("continued.fraction",	"yield.estimates", "step.size");
	return(result);
}

print.continuedfraction <- function(CF)
{
	print(CF$cf.coeffs);
}

## generate complexity curve through bootstrap the histogram
bootstrap.complex.curve <- function(hist, times = 100, di = 0, mt = 100,
									ss = 1e6, mv = 1e10, 
									max.extrapolation = 1e10)
{
	# input could be either histogram file or count vector of the histogram
	if (mode(hist) == "character") {
		hist.count = preseqR.read.hist(hist);
	}
	else {
		hist.count = hist;
	}
	total.sample = 0.0;
	for (i in 1:length(hist.count))
		total.sample <- total.sample + i * hist.count[i];
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
		# the number of sampled points for complexity curve
		N = 0
		# the number of times bootstrap success
		count = 0
		# the actually step.size preseqR.continued.fraction.estimate uses
		step.size = 0
#		estimates = list();
		estimates = matrix(data = NA, nrow = max.extrapolation / ss, 
				           ncol = times, byrow = FALSE);
		for (i in 1:as.integer(times))
		{
			# do sampling with replacement 
			sample = preseqR.hist.sample(hist.count, as.integer(total.sample), 
								replace = TRUE);
			# build count vector of the histogram based on sampling results
			hist = preseqR.sample2hist.count(sample, replace = TRUE);
			out <- preseqR.continued.fraction.estimate(hist, di, mt, ss,
				                                    	mv, max.extrapolation);
			if (!is.null(out))
			{
				count <- count + 1;
				N = length(out$yield.estimates$yields);
				step.size = out$step.size
				estimates[, i][1: N] = out$yield.estimates$yields;
			}
		}
		# return NULL if all bootstrap failed to construct a continued fraction
		if (N == 0)
		{
			write("can not make prediction based on the given histogram", stderr());
			return()
		}
		if (count < BOOTSTRAP.factor * times)
		{
			write("fail to bootstrap since the histogram is poor", stderr());
			return();
		}
		index = as.double(step.size) * ( 1:N )
		# mean values are used as complexity curve
		mean = apply(estimates[1:N, ], 1, mean, na.rm = TRUE)
		variance = apply(estimates[1:N, ], 1, var, na.rm = TRUE)
		# count the number of values not zero
		n = as.vector(apply(estimates, 1, function(x) length(which(!is.na(x)))))
		n = n[1: N]
		# sample size list
		# 95% confident interval based on normal distribution
		left.interval = mean - qnorm(0.975) * sqrt(variance / n);
		right.interval = mean + qnorm(0.975) * sqrt(variance / n);
		yield.estimates = list(sample.size = index, yields = yield.estimates)
		result = list(yield.estimates, left.interval, right.interval);
		names(result) = c("yield.estimates", "LOWER_0.95CI", 
						  "UPPER_0.95CI")
		return(result);
	} else {
		write("the paramter times should be at least one", stderr());
		return();
	}
}
