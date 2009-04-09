## ===========================
## Treatment of missing values
## ===========================

seqprep <- function(seqdata, left=NA, right="DEL", gaps=NA, 
	neutral="#", missing=NA, void="%", nr="*") {

	nbseq <- nrow(seqdata)
	sl <- ncol(seqdata)

	message(" [>] preparing ",nbseq, " sequences")
	message(" [>] using ", void, " for void elements and ", nr, " for missing values")

	if (is.na(missing)) mstate <- is.na(seqdata)
	else mstate <- seqdata==missing

	for (i in 1:nbseq) 
		seqdata[i,] <- TraMineR.trunc(seqdata[i,], mstate[i,], sl,
			left=left, right=right, gaps=gaps, 
			neutral=neutral, void=void)

	## Setting a new code for missing statuses
	if (is.na(missing)) seqdata[is.na(seqdata)] <- nr
	else seqdata[seqdata==missing] <- nr

	return(seqdata)
} 


