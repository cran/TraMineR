## ============================
## Number of matching positions
## ============================

seqmpos <- function(iseq,jseq) {

	clength <- min(length(iseq), length(jseq))

	mpos <- 0

	for (isym in 1:clength) {
		c <- iseq[isym]==jseq[isym]
		
		if (!is.na(c)) {
			if (c==TRUE) mpos <- mpos+1
		}
	}

	return(mpos)
}
