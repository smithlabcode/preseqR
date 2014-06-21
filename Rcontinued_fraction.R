MAXLENGTH = 10000000

R_read_hist <- function(hist_file)
{
	hist_table = read.table(hist_file);
	hist_count = vector(mode = 'numeric', length = 0);
	pre_pos = 0;
	for (i in 1:length(hist_table[, 1]))
	{
		if (i > 1 && hist_table[i, 1] <= hist_table[i - 1, 1])
		{
			stop("The input histogram is not sorted in increasing order");
		}
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

#	out <- .C("R_Advance_ContinueFraction", hist_count = c(0, hist_count), l, di, de, PS_COEFFS, PS_COEFFS_L, CF_COEFFS, CF_COEFFS_L, OFFSET_COEFFS);
#
#	length(PS_COEFFS) = PS_COEFFS_L;
#	length(CF_COEFFS) = CF_COEFFS_L;
#	length(OFFSET_COEFFS) = as.integer(abs(di));
#	CF = list(PS_COEFFS, CF_COEFFS, OFFSET_COEFFS, di, de);
#	names(CF) = c(ps_coeffs, cf_coeffs, offset_coeffs, diagonal_idx, degree);
#	return(CF);
#}


R_CAL_ContinueFraction <- function(CF, value)
{
	cf = as.double(CF$cf_coeffs);
	cf_l = as.integer(length(cf));
	off = as.double(CF$offset_coeffs);
	di = as.integer(CF$diagonal_idx);
	de = as.integer(CF$degree);
	coordinate = as.double(value);
	out <- .C("R_Calculate_ContinuedFraction", cf, cf_l, off, di, de, coordinate, r = as.double(1))
	return(out$r);
}

R_ADVANCE_extrapolate_distinct <- function(hist_count, CF, start_size, step_size, max_size)
{
	cf_coeffs = as.double(CF$cf_coeffs);
	cf_coeffs_l = as.integer(length(CF$cf_coeffs));
	offset_coeffs = as.double(CF$offset_coeffs);
	di = as.integer(CF$diagonal_idx);
	de = as.integer(CF$degree);
	hist_count = as.double(hist_count);
	hist_count = c(0, hist_count);
	hist_count_l = as.integer(length(hist_count));
	estimate = as.double(vector(mode = 'numeric', length = max_size / step_size + 1));
	estimate_l = as.integer(1);
	.C("R_Advance_extrapolate_distinct", cf_coeffs, cf_coeffs_l, offset_coeffs, di, de, hist_count, hist_count_l, as.double(start_size), as.double(step_size), as.double(max_size), estimate, estimate_l);
	length(estimate) = estimate_l;
	return(estimate);
}
	
R_ADVANCE_sample <- function(hist_count, size, replace = FALSE)
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
		x = vector(mode = 'numeric', length = as.integer(total_reads));
		pos = 1;
		p = 1;
		value = 1;
		for (l in hist_count)
		{
			if (as.integer(l) > 0)
			{
				x[pos: (pos + as.integer(l) * p - 1)] = rep(value: (value + as.integer(l) - 1), p);
				value <- value + as.integer(l);
				pos <- pos + as.integer(l) * p;
			}
			p <- p + 1;
		}
		return(sample(x, size, replace = FALSE));
	}
	else
	{
		distinct_reads = as.integer(sum(hist_count));
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

sample2hist_count <- function(vec)
{
	V = table(vec);
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

R_continuedfraction_estimate <- function(hist_count, di, mt, ss, mv)
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
		s = R_ADVANCE_sample(hist_count, sample);
		yield = sum(sample2hist_count(s));
		yield_estimate = c(yield_estimate, yield);
		sample <- sample + step;
	}
	counts_before_first_zero = 1;
	while (as.integer(counts_before_first_zero) <= length(hist_count) && hist_count[counts_before_first_zero] != 0)
		counts_before_first_zero <- counts_before_first_zero + 1;

## continued fraction with even degree conservatively estimates
	mt = min(mt, counts_before_first_zero - 1);
	mt = mt - (mt %% 2);

	hist_count = c(0, hist_count);
	ps_coeffs = as.double(vector(mode = 'numeric', length = MAXLENGTH));
	cf_coeffs = as.double(vector(mode = 'numeric', length = MAXLENGTH));
	offset_coeffs = as.double(vector(mode = 'numeric', length = MAXLENGTH));
	diagonal_idx = as.integer(0);
	degree = as.integer(0);
	is_valid = as.integer(0);
	out <- .C('continuefraction_estimate', as.double(hist_count), as.integer(length(hist_count)), as.integer(di), as.integer(mt), as.double(ss), as.double(mv), ps_coeffs, ps_coeffs_l = as.integer(length(ps_coeffs)), cf_coeffs, cf_coeffs_l = as.integer(length(cf_coeffs)), offset_coeffs, diagonal_idx, degree, is_valid);
	if (!is_valid)
	{
		write("Fail to construct and need to bootstrap to obtain estimates", stderr);
		return;
	}
	est = R_ADVANCE_extrapolate_distinct(hist_count, CF, sample / total_sample, step_size / total_sample, max_size / total_sample);
	yield_estimate = c(yield_estimate, est);
	length(ps_coeffs) = ps_coeffs_l;
	length(cf_coeffs) = cf_coeffs_l;
	length(offset_coeffs) = as.integer(abs(diagonal_idx));
	CF = list(ps_coeffs, cf_coeffs, offset_coeffs, diagonal_idx, degree);
	names(CF) = c('ps_coeffs', 'cf_coeffs', 'offset_coeffs', 'diagonal_idx', 'degree');
	result = c(CF, yield_estimate)
	return(result);
}

print_continuedfraction <- function(CF)
{}

bootstrap_complex_curve <- function(hist)
{}

