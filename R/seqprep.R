## ===========================
## Treatment of missing values
## ===========================

seqprep <- function(seqdata, left=NA, right="DEL", gaps=NA, 
	neutral="#", missing=NA, void="%", nr="*") {

	nbseq <- seqdim(seqdata)[1]

	message(" [>] preparing ",nbseq, " sequences")
	message(" [>] using ", void, " for void elements and ", nr, " for missing values")

	for (i in 1:nbseq)
		seqdata[i,] <- TraMineR.trunc(seqdata[i,], 
			left=left, right=right, gaps=gaps, 
			neutral="#", missing=missing, void=void)

	## Setting a new code for missing statuses
	if (is.na(missing)) seqdata[is.na(seqdata)] <- nr
	else seqdata[seqdata==missing] <- nr

	return(seqdata)
} 


