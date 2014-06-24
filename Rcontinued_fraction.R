MAXLENGTH = 10000000

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


#R.ADVANCE.ContinueFraction <- function(hist.count, diagonal.idx, degree)
#{
#	hist.count = as.double(hist.count);
## For hist.count as parameter in C function, the first row should be (0, 0)
## For hist.count as parameter in R function, the first row should be (1, x)
#	hist.count = c(0, hist.count);
#	l = as.integer(length(hist.count));
#	if (l == 0) stop('empty hist.count');
#	di = as.integer(diagonal);
#	de = as.integer(degree);
#	PS.COEFFS = as.double(vector(mode = 'numeric', length = l));
#	PS.COEFFS.L = as.integer(1);
#	CF.COEFFS = as.double(vector(mode = 'numeric', length = l));
#	CF.COEFFS = as.integer(1);
#	OFFSET.COEFFS = as.double(vector(mode = 'numeric', length = l));

#	out <- .C("R.Advance.ContinueFraction", hist.count = c(0, hist.count), l, \
#   di, de, PS.COEFFS, PS.COEFFS.L, CF.COEFFS, CF.COEFFS.L, OFFSET.COEFFS);
#
#	length(PS.COEFFS) = PS.COEFFS.L;
#	length(CF.COEFFS) = CF.COEFFS.L;
#	length(OFFSET.COEFFS) = as.integer(abs(di));
#	CF = list(PS.COEFFS, CF.COEFFS, OFFSET.COEFFS, di, de);
#	names(CF) = c(ps.coeffs, cf.coeffs, offset.coeffs, diagonal.idx, degree);
#	return(CF);
#}

## calculate the value of the continued fraction CF given the coordinate x 
## call c-encoded function "c.calculate.continued.fraction()" through 
## preseqR.calculate.continue.fraction
preseqR.calculate.continue.fraction <- function(CF, x)
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
	   estimate.l = as.integer(0), DUP = FALSE);
	return(out$estimate[ 1: out$estimate.l ]);
}
	
#do with/without replacement of random sampling given a histogram count vector
#size is a user defined sample size
preseqR.hist.sample <- function(hist.count, size, replace = FALSE)
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
	else
	{
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
	}
}

#convert sampled points into a count vector of the histogram 
preseqR.sample2hist.count <- function(sample.points)
{
	V = table(sample.points);
	if (names(V)[1] == '0')
	{
		V = V[-1];
	}
	value = as.vector(V);
	hist.count = vector(mode = 'numeric', length = max(value));
	for (v in value)
	{
		hist.count[v] <- hist.count[v] + 1;
	}
	return(hist.count)
}


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
	yield.estimates = as.double(vector(mode = 'numeric', length = l));
	for (i in 1:l)
	{
		s = preseqR.hist.sample(hist.count, sample);
		yield = sum(preseqR.sample2hist.count(s));
		yield.estimates[i] = yield;
		sample <- sample + step;
	}
	return(yield.estimates);
}


## estimate a continued fraction given a the count vector of the histogram
preseqR.continued.fraction.estimate <- function(hist.count, di, mt, ss, mv, 
										 max.extrapolation)
{
	total.sample = 0.0;
	for (i in 1:length(hist.count))
		total.sample <- total.sample + i * hist.count[i];
	#interpolation when sample size is no more than total sample size
	yield.estimates = preseqR.interpolate.distinct(hist.count, ss)

	counts.before.first.zero = 1;
	while (as.integer(counts.before.first.zero) <= length(hist.count) && 
		   			  hist.count[counts.before.first.zero] != 0)
		counts.before.first.zero <- counts.before.first.zero + 1;

	# starting sample size for extrapolation
	sample = ss * as.integer(total.sample / ss) + ss;
	# continued fraction with even degree conservatively estimates
	mt = min(mt, counts.before.first.zero - 1);
	mt = mt - (mt %% 2);

	hist.count = c(0, hist.count);
	# allocate spaces to store constructed continued fraction
	ps.coeffs = as.double(vector(mode = 'numeric', length = MAXLENGTH));
	cf.coeffs = as.double(vector(mode = 'numeric', length = MAXLENGTH));
	offset.coeffs = as.double(vector(mode = 'numeric', length = MAXLENGTH));
	diagonal.idx = as.integer(0);
	degree = as.integer(0);
	is.valid = as.integer(0);
	# construct a continued fraction with minimum degree
	out <- .C('c_continued_fraction_estimate', as.double(hist.count), 
			  as.integer(length(hist.count)), as.integer(di), as.integer(mt), 
			  as.double(ss), as.double(mv), ps.coeffs, 
			  ps.coeffs.l = as.integer(0), cf.coeffs, 
			  cf.coeffs.l = as.integer(0), offset.coeffs, 
			  diagonal.idx, degree, is.valid, DUP = FALSE);
	if (!is.valid)
	{
		write("Fail to construct and need to bootstrap to obtain estimates", 
			  stderr());
		return;
	}
	length(ps.coeffs) = out$ps.coeffs.l;
	length(cf.coeffs) = out$cf.coeffs.l;
	length(offset.coeffs) = as.integer(abs(diagonal.idx));
	CF = list(ps.coeffs, cf.coeffs, offset.coeffs, diagonal.idx, degree);
	names(CF) = c('ps.coeffs', 'cf.coeffs', 'offset.coeffs', 'diagonal.idx', 
				  'degree');
	# the value of the step.size is equal to the size of the sample from the 
    # histogram. thus set step.size equal to one
	# extrapolation when sample size is larger than the inital experiment
	est <- preseqR.extrapolate.distinct(hist.count, CF, sample / total.sample,
		 step.size = 1, max.extrapolation / total.sample);
	yield.estimates = c(yield.estimates, est);
	result = list(CF, yield.estimates)
	names(result) = c("continued.fraction", "yield.estimates");
	return(result);
}

print.continuedfraction <- function(CF)
{}

bootstrap.complex.curve <- function(hist)
{}

