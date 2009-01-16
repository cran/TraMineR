## =======================================
## Extracts distinct states from sequences
## =======================================

seqdss <- function(seqdata) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	nbseq <- seqdim(seqdata)[1]
	maxsl <- max(seqlength(seqdata))
	trans <- matrix(nrow=nbseq, ncol=maxsl)
	statl <- attr(seqdata,"alphabet")

	for (i in 1:nbseq) {
		idx <- 1
		j <- 1

		tmpseq <- seqdata[i,]
		sl <- seqlength(tmpseq)
		tmpseq <- as.integer(seqdata[i,])
		
		while (idx <= sl) {
			iseq <- tmpseq[idx]
			dur <- 1

			while (idx < sl & tmpseq[idx+1]==iseq) { 
					idx <- idx+1
			}

			trans[i,j] <- statl[iseq]

			j <- j+1
			idx <- idx+1
		}
	}
	
	trans <- suppressMessages(seqdef(trans, alphabet=statl, cnames=paste("ST",seq(1:maxsl),sep="")))

	return(trans)
}

