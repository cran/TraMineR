## ==============================
## Convert from SPS to STS format
## ==============================

SPS_to_STS <- function(seqdata,type,stsep) {

	nbseq <- seqdim(seqdata)[1]
	trans <- matrix("", nrow=nbseq, ncol=1)
	
	for (i in 1:nbseq) {
		tmpseq <- strsplit(seqdata[i],split=stsep)[[1]]
		for (s in 1:length(tmpseq)) {
			if (type==1) sps <- strsplit(gsub("[()]","",tmpseq[s]),split=",")[[1]]
			else sps <- strsplit(tmpseq[s],split="/")[[1]]

			seq <- sps[[1]]
			dur <- as.integer(sps[[2]])

			if (s==1) trans[i] <- paste(trans[i],seq,sep="") 
			else trans[i] <- paste(trans[i],seq,sep="-")

			if (dur>1)
				for (r in 2:dur) trans[i] <- paste(trans[i],"-",seq,sep="")
		}
	}

	return(trans)
}
