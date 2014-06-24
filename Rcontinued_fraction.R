MAXLENGTH = 10000000

## read a histogram file; return the histogram count vector
## count vector represent frequencies of indexes. For those indexes not showing
## in the histogram file, use zeros to represent their values
preseqR_read_hist <- function(hist_file)
{
	#the first column is the index of the histogram
	#the second column is the frequency of indexes
	hist_table = read.table(hist_file);
	hist_count = vector(mode = 'numeric', length = 0);
	pre_pos = 0;
	for (i in 1:length(hist_table[, 1]))
	{
		#Indexes of the histogram should be increasing order
		if (i > 1 && hist_table[i, 1] <= hist_table[i - 1, 1])
		{
			stop("The input histogram is not sorted in increasing order");
		}
		#set the frequency of each index; 
		#fill the gaps between indexes the zeros
		hist_count = c(hist_count, rep(0, hist_table[i, 1] - pre_pos));
		hist_count[hist_table[i, 1]] = hist_table[i, 2];
		pre_pos = hist_table[i, 1];
	}
	return(hist_count);
}


#R_ADVANCE_ContinueFraction <- function(hist_count, diagonal_idx, degree)
#{
#	hist_count = as.double(hist_count);
## For hist_count as parameter in C function, the first row should be (0, 0)
## For hist_count as parameter in R function, the first row should be (1, x)
#	hist_count = c(0, hist_count);
#	l = as.integer(length(hist_count));
#	if (l == 0) stop('empty hist_count');
#	di = as.integer(diagonal);
#	de = as.integer(degree);
#	PS_COEFFS = as.double(vector(mode = 'numeric', length = l));
#	PS_COEFFS_L = as.integer(1);
#	CF_COEFFS = as.double(vector(mode = 'numeric', length = l));
#	CF_COEFFS = as.integer(1);
#	OFFSET_COEFFS = as.double(vector(mode = 'numeric', length = l));

#	out <- .C("R_Advance_ContinueFraction", hist_count = c(0, hist_count), l, \
#   di, de, PS_COEFFS, PS_COEFFS_L, CF_COEFFS, CF_COEFFS_L, OFFSET_COEFFS);
#
#	length(PS_COEFFS) = PS_COEFFS_L;
#	length(CF_COEFFS) = CF_COEFFS_L;
#	length(OFFSET_COEFFS) = as.integer(abs(di));
#	CF = list(PS_COEFFS, CF_COEFFS, OFFSET_COEFFS, di, de);
#	names(CF) = c(ps_coeffs, cf_coeffs, offset_coeffs, diagonal_idx, degree);
#	return(CF);
#}

## calculate the value of the continued fraction CF given the coordinate x 
## call c-encoded function "c_calculate_continued_fraction()" through 
## preseqR_calculate_continue_fraction
preseqR_calculate_continue_fraction <- function(CF, x)
{
	out <- .C("c_calculate_continued_fraction", 
			  cf = as.double(CF$cf_coeffs), 
	   		  cf_l = as.integer(length(CF$cf_coeffs)),
			  off = as.double(CF$offset_coeffs), 
			  di = as.integer(CF$diagonal_idx), 
			  de = as.integer(CF$degree), 
			  coordinate = as.double(x), 
		  	  result = as.double(0));
	return(out$result);
}

## extrapolate given a histogram and a continued fraction
preseqR_extrapolate_distinct <- function(hist_count, CF, start_size = NULL,
		 step_size = NULL, max_size = 100)
{
	cf_coeffs = as.double(CF$cf_coeffs);
	cf_coeffs_l = as.integer(length(CF$cf_coeffs));
	offset_coeffs = as.double(CF$offset_coeffs);
	di = as.integer(CF$diagonal_idx);
	de = as.integer(CF$degree);
	hist_count = as.double(hist_count);
	# record the size of the sample based on the histogram count
	total_reads = 0.0;
	for (i in 1:length(hist_count))
		total_reads <- total_reads + i * as.integer(hist_count[i]);

	# the styles of the histogram count vector are different between R code 
	# and c++ code; The first line of the histogram is always [0  0] in c++ 
	# but the line is removed in R-encoded function
	hist_count = c(0, hist_count);
	hist_count_l = as.integer(length(hist_count));

	# set start_size and step_size if they are not defined by user
	if (is.null(start_size))
		start_size = total_reads;
	if (start_size > max_size)
	{
		write("start position has already beyond the maximum prediction",
			  stderr());
		return(NULL);
	}
	if (is.null(step_size))
		step_size = total_reads;

	# allocate memory to store extrapolation results
	# first "c_extrapolate_distinct" stores the observed number of distinct 
	# molecules into estimate, then it stores the extrapolation values
	# thus the allocated memory size is 1 plus the size of extrapolation values,
	# which is (max_size - start_size) / step_size) + 1
	extrap_size = as.integer((max_size - start_size) / step_size) + 1;
	out <- .C("c_extrapolate_distinct", cf_coeffs, cf_coeffs_l, offset_coeffs,
	   di, de, hist_count, hist_count_l, as.double(start_size), 
	   as.double(step_size), as.double(max_size), 
	   estimate = as.double(vector(mode = 'numeric', extrap_size + 1)), 
	   estimate_l = as.integer(0), DUP = FALSE);
	return(out$estimate[ 1: out$estimate_l ]);
}
	
#do with/without replacement of random sampling given a histogram count vector
#size is a user defined sample size
preseqR_hist_sample <- function(hist_count, size, replace = FALSE)
{
	if (replace == FALSE)
	{
		total_reads = 0;
		i = 1;
		while (i <= length(hist_count))
		{
			total_reads = i * as.integer(hist_count[i]);
			i <- i + 1;
		}
		#construct a sample space X 
		X = vector(mode = 'numeric', length = as.integer(total_reads));
		pos = 1;
		p = 1;
		value = 1;
		for (l in hist_count)
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
		distinct_reads = as.integer(sum(hist_count));
		#construct the pdf of the multinomial distribution
		prob = vector(mode = 'numeric', length = distinct_reads);
		p = 1;
		pos = 1;
		for (l in hist_count)
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
preseqR_sample2hist_count <- function(sample_points)
{
	V = table(sample_points);
	if (names(V)[1] == '0')
	{
		V = V[-1];
	}
	value = as.vector(V);
	hist_count = vector(mode = 'numeric', length = max(value));
	for (v in value)
	{
		hist_count[v] <- hist_count[v] + 1;
	}
	return(hist_count)
}

#estimate a continued fraction given a the count vector of the histogram
preseqR_continued_fraction_estimate <- function(hist_count, di, mt, ss, mv, 
										 max_extrapolation)
{
	total_sample = 0.0;
	for (i in 1:length(hist_count))
		total_sample <- total_sample + i * hist_count[i];
	inital_distint = sum(hist_count);
	upper_limit = as.integer(total_sample)
	step = ss;
	sample = step;
	yield_estimate = as.double(vector(mode = 'numeric', length = 0));
	while(sample <= upper_limit)
	{
		s = preseqR_hist_sample(hist_count, sample);
		yield = sum(preseqR_sample2hist_count(s));
		yield_estimate = c(yield_estimate, yield);
		sample <- sample + step;
	}
	counts_before_first_zero = 1;
	while (as.integer(counts_before_first_zero) <= length(hist_count) && 
		   			  hist_count[counts_before_first_zero] != 0)
		counts_before_first_zero <- counts_before_first_zero + 1;

	## continued fraction with even degree conservatively estimates
	mt = min(mt, counts_before_first_zero - 1);
	mt = mt - (mt %% 2);

	hist_count = c(0, hist_count);
	## allocate spaces to store constructed continued fraction
	ps_coeffs = as.double(vector(mode = 'numeric', length = MAXLENGTH));
	cf_coeffs = as.double(vector(mode = 'numeric', length = MAXLENGTH));
	offset_coeffs = as.double(vector(mode = 'numeric', length = MAXLENGTH));
	diagonal_idx = as.integer(0);
	degree = as.integer(0);
	is_valid = as.integer(0);
	out <- .C('c_continued_fraction_estimate', as.double(hist_count), 
			  as.integer(length(hist_count)), as.integer(di), as.integer(mt), 
			  as.double(ss), as.double(mv), ps_coeffs, 
			  ps_coeffs_l = as.integer(0), cf_coeffs, 
			  cf_coeffs_l = as.integer(0), offset_coeffs, 
			  diagonal_idx, degree, is_valid, DUP = FALSE);
	if (!is_valid)
	{
		write("Fail to construct and need to bootstrap to obtain estimates", 
			  stderr());
		return;
	}
	length(ps_coeffs) = out$ps_coeffs_l;
	length(cf_coeffs) = out$cf_coeffs_l;
	length(offset_coeffs) = as.integer(abs(diagonal_idx));
	CF = list(ps_coeffs, cf_coeffs, offset_coeffs, diagonal_idx, degree);
	names(CF) = c('ps_coeffs', 'cf_coeffs', 'offset_coeffs', 'diagonal_idx', 
				  'degree');
	# the value of the step_size is equal to the size of the sample from the 
    # histogram
	est <- preseqR_extrapolate_distinct(hist_count, CF, sample / total_sample,
		 step_size = 1, max_extrapolation / total_sample);
	yield_estimate = c(yield_estimate, est);
	result = list(CF, yield_estimate)
	names(result) = c("continued_fraction", "yield_estimates");
	return(result);
}

print_continuedfraction <- function(CF)
{}

bootstrap_complex_curve <- function(hist)
{}

