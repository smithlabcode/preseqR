## density function of truncated zero negative binomial distribution
## size and mu are two parameters for negative binomial
zerotruncated.dnbinom <- function(x, size, mu, log = FALSE)
{
	# the density of x in negative binomial
	p = dnbinom(x, size = size, mu = mu, log = log);
	# set zeros in x with zero probability
	if (log == FALSE) {
		p[ which(x == 0) ] = 0;
	} else {
		p[ which(x == 0) ] = -Inf;
	}
	# the density of non-zero in negative binomial
	q = 1 - dnbinom(0, size = size, mu = mu);
	# normalize all non-zero values in negrative binomial to generate ZTNB
	if (log == FALSE) {
		return(p / q);
	} else {
		return(p - log(q));
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
preseqR.zerotruncated.mle <- function(hist, size = SIZE.INIT, 
									  mu = MU.INIT)
{
	if (mode(hist) == 'character') {
		hist.count = read.hist(hist);
	} else {
		hist.count = hist;
	}
	total.sample = (1:length(hist.count) %*% hist.count);
	distinct.sample = sum(hist.count);
	f <- function(x) zerotruncated.minus.log.likelyhood(hist.count, 
														size = x[1], mu = x[2]);
	return(optim(c(size, mu), f, NULL, method = "L-BFGS-B", 
				lower = c(0.0001, 0.0001), upper = c(10000, 10000)))
}

## predict the number of distinct items using zero truncated negative binomial 
## distribution
## n is the size of experiment
preseqR.zerotruncated.estimate <- function(hist, n)
{
	if (mode(hist) == 'character') {
		hist.count = read.hist(hist);
    } else {
		hist.count = hist;
	}   

	total.sample = (1:length(hist.count) %*% hist.count);
	distinct.sample = sum(hist.count);

	opt <- preseqR.zerotruncated.mle(hist.count);
	size = opt$par[1];
	mu = opt$par[2];
	# the probability of being sampled in the initial experiment
	p = 1 - dnbinom(0, size = size, mu = mu);
	# L is the estimated total number of distinct items
	L = distinct.sample / p;
	# update parameters of negative binomial in the experiment with size n
	mu = mu * n / as.double(total.sample);
	# the probability of being sampled under the new experiment
	P = 1 - dnbinom(0, size = size, mu = mu);
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
	total.sample = (1:length(hist.count) %*% hist.count);
	distinct.sample = sum(hist.count);
	# n is the number of experiments
	n = as.integer(max.extrapolation / ss);
	sample.size = as.double(ss) * (1: n);
	dim(sample.size) = n;

	# estimate parameters
	opt <- preseqR.zerotruncated.mle(hist.count)
	size = opt$par[1];
	mu = opt$par[2];
	# the probability of being sampled in the initial experiment
	p = 1 - dnbinom(0, size = size, mu = mu);
	# L is the estimated total number of distinct items
	L = distinct.sample / p;
	# estimate the item being sampled under new experiments with different size
	t = sample.size / as.double(total.sample);
	dim(t) = length(t)
	P = apply(t, 1, function(x) 1 - dnbinom(0, size, mu = x * mu))
	yield.estimates = L * P;
	yield.estimates = list(sample.size = sample.size, 
						   yield.estimates = yield.estimates);
	return(yield.estimates)
}


## calculate the negative binomial loglikelyhood
## hist.count is a count vector of a histogram of observed items
## zero.items is number of items unobserved
## size and mu are parameters in a negative binomial distribution
nb.loglikelyhood <- function(hist.count, zero.items, size, mu)
{
	# likelyhood of nonzero terms
	log.prob = dnbinom(1:length(hist.count), size = size, mu = mu, log = TRUE);
	loglikelyhood = log.prob %*% hist.count
	# add items with zero count
	log.zero.prob = dnbinom(0, size = size, mu = mu, log = TRUE)
	loglikelyhood <- loglikelyhood + zero.items * log.zero.prob
	return(loglikelyhood)
}


## EM algorithm to fit a 
## hist only includes information for observation
## the number of unobserved items is missing data
preseqR.nbinom.em <- function(hist, size = SIZE.INIT, mu = MU.INIT)
{
	if (mode(hist) == 'character') {
		hist.count = read.hist(hist);
	} else {
		hist.count = hist;
	}
	# setting the number of unobserved items as 0
	zero.prob = exp(dnbinom(0, size = size, mu = mu, log = TRUE))
	# estimate the total number of distinct items
	observed.items = sum(hist.count);
	L = observed.items / (1 - zero.prob);
	# expected the number of unobservations
	zero.items = L * zero.prob
	# convert zero.items into an integer
	zero.items = floor(zero.items)
	# estimated mean 
	m = (1:length(hist.count) %*% hist.count) / L
	# estimated variance
	v = (((1:length(hist.count) - m)^2) %*% hist.count + m^2 * zero.items) / (L - 1)
	# make sure each item in histogram is an integer
	hist.count = floor(hist.count)
	loglikelyhood.pre = Inf;
	f <- function(x) -nb.loglikelyhood(hist.count, zero.items, size = x, mu = m) / L
	# derivative of f
	gr <- function(x) 
	{
		freq = 1:length(hist.count)
		first.term = (digamma(x) * zero.items + (digamma(freq + size) %*% hist.count)) / L
		second.term = digamma(x)
		third.term = log(x) - log(x + m)
		result = first.term - second.term + third.term
		# f is negative loglikelyhood
		return(-result)
	}
	# estimate size and mu based on first and second moments
	if (v > m) {
		res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B", 
					 lower = 0.0001, upper = 10000)
	} else {
		res <- optim(size, f, gr, method = "L-BFGS-B", 
					 lower = 0.0001, upper = 10000)
	}
	# count the times of iteration
	iter = as.double(1)
	# zerotruncated loglikelyhood 
	loglikelyhood = zerotruncated.minus.log.likelyhood(hist.count, res$par, m)
	# make sure EM algorithm could terminate
	while (abs(loglikelyhood.pre - loglikelyhood) / observed.items > TOLERANCE 
			&& iter < ITER.TOLERANCE)
	{
		# update minus loglikelyhood
		loglikelyhood.pre = loglikelyhood
		# update parameters
		size = res$par;
		mu = m;
		# update the probility an item unobserved
		zero.prob = exp(dnbinom(0, size = size, mu = mu, log = TRUE))
		# estimate the total number of distinct items
		L = observed.items / (1 - zero.prob)
		# update expected number of unobserved items
		zero.items = L * zero.prob
		# convert zero.items into an integer
		zero.items = floor(zero.items)
		# estimated mean 
		m = (1:length(hist.count) %*% hist.count) / L
		# estimated variance
		v = (((1:length(hist.count) - m)^2) %*% hist.count + m^2 * zero.items) / (L - 1)
		# M step: estimate the parameters size and mu
		if (v > m) {
			res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B", 
						 lower = 0.0001, upper = 10000)
		} else {
			res <- optim(size, f, gr, method = "L-BFGS-B", 
						 lower = 0.0001, upper = 10000)
		}
		iter <- iter + 1
		# zerotruncated loglikelyhood 
		loglikelyhood = zerotruncated.minus.log.likelyhood(hist.count, res$par, m)
	}
	return(list(size = size, mu = m, loglik = -loglikelyhood.pre))
}
