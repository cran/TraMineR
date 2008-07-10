## ====================================
## Calculates the Length of the Longest
## Common Prefix
## ====================================

seqLCP <- function(iseq,jseq) {

	clength <- min(length(iseq), length(jseq))

	dlength <- 0

	for (isym in 1:clength) {
		c <- iseq[isym]==jseq[isym]
		
		if (is.na(c) | c==FALSE) return(dlength)
			else dlength <- isym
			}

	return(dlength)
}







