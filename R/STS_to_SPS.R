## ==============================
## Convert from STS to SPS format
## ==============================

STS_to_SPS <- function(seqdata,type) {

	nbseq <- seqdim(seqdata)[1]
	out <- matrix("", nrow=nbseq, ncol=1)

	rownames(out) <- paste("[",seq(1:nbseq),"]",sep="")
	colnames(out) <- "SPS sequence"

	if (type==1) {
		prefix <- "("
		stdursep <- ","
		suffix <- ")"
	}
	else {
		prefix <- ""
		stdursep <- "/"
		suffix <- ""
	}		

	for (i in 1:nbseq) {
		idx <- 1
		tmpseq <- strsplit(seqdata[i],split="-")[[1]]
		
		while (idx <= length(tmpseq)) {
			iseq <- tmpseq[idx]

			if (idx==1) out[i] <- paste(out[i],prefix,iseq,stdursep,sep="")
			else out[i] <- paste(out[i],"-",prefix,iseq,stdursep,sep="")

			dur <- 1
			while (idx < length(tmpseq) & tmpseq[idx+1]==iseq) { 
				dur <- dur+1
				idx <- idx+1
			}

			## adding suffix
			out[i] <- paste(out[i],dur,suffix,sep="")
			idx <- idx+1
		}

	}

	return(out)
}

