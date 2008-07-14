## ========================================================
## Creates a vector of 0 and 1 used for plotting a sequence
## ========================================================

seqgbar <- function(seq, seql, statl, nbstat) {
	gbar <- vector("integer", seql*(nbstat+1))

	for (j in 1:seql) {
		if (!is.na(seq[j])) 
			gbar[((j-1)*(nbstat+1))+which(statl==seq[j])] <- 1
		else 
			gbar[((j-1)*(nbstat+1))+nbstat+1] <- 1
	}

	return(gbar)
}
	
