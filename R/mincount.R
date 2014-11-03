## calculate expected number of species appearing at least k times
## (k >= 2) through non-parameteric emprical estimator

## map the histogram into a power series;s
## read the calculated histogram through the power series;
## the parameter hist is the original histogram; k is the minimum required freq


get.hist.freq <- function(hist, j) {
	if (j %in% hist[, 1]) {
		return(hist[j, 2])
	} else {
		return(0)
	}
}

mincount.hist <- function(hist, k) {
	max.freq <- max(hist[, 1])
	res <- vector(length = max.freq, mode = 'numeric')
	for (i in 1:max.freq) {
		co.eff <- 0
		for ( l in 0:(k-1) ) {
			for (j in 0:i) {
				index <- j + l
				if (i <= index) {
					co.eff <- co.eff + (-1)^(j + 1) * get.hist.freq(hist, index) *
						       exp( lgamma(j + l + 1) - lgamma(j + 1) - lgamma(l + 1) ) *
									 exp( lgamma(l+1) - lgamma(i-j+1) - lgamma(l - i + j + 1) )
				}
			}
		}
		res[i] <- co.eff
	}
	res
}

