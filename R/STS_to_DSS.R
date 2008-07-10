## ==============================
## Convert from STS to DSS format
## ==============================

STS_to_DSS <- function(seqdata) {

	nbseq <- seqdim(seqdata)[1]
	out <- matrix("", nrow=nbseq, ncol=1)
	rownames(out) <- paste("[",seq(1:nbseq),"]",sep="")
	colnames(out) <- "DSS sequence"

	for (i in 1:nbseq) {
		idx <- 1
		tmpseq <- strsplit(seqdata[i],split="-")[[1]]

		while (idx <= length(tmpseq)) {
			iseq <- tmpseq[idx]
			
			if (idx==1) out[i] <- paste(out[i],iseq,sep="")
			else out[i] <- paste(out[i],"-",iseq,sep="")

			while (idx < length(tmpseq) & tmpseq[idx+1]==iseq) { 
				idx <- idx+1
			}
	
			idx <- idx+1
		}
	}
	return(out)
}

