## ================================
## Convert from SPELL to STS format
## ================================

TSE_to_STS <- function(seqdata, stdur=TRUE, start=NULL,istate=NULL, end, eqtime=0) {

	lid <- unique(seqdata[,1])
	nbseq <- length(lid)

	trans <- matrix("", nrow=nbseq, ncol=1)
		
	for (i in 1:nbseq) {
		## OBLIGE D'UTILISER L'INDEX POUR LES COLONNES SINON NE MARCHE PAS
		spell <- seqdata[seqdata[,1]==lid[i],]
		if (!is.null(istate)) {
			s0 <- cbind(lid[i], start, istate)
			spell <- rbind(s0, spell)
		}
		idxmax <- nrow(spell)
		
		if (idxmax>0) {
			for (j in 1:idxmax) {
				if (stdur==TRUE) {
					if (j<idxmax) dur <- spell[j+1,2]-spell[j,2]
					else dur <- end-spell[j,2]
				} else dur <- 1

				if (dur==0) {
					if (eqtime==0) seq <- "NA"
					if (eqtime==1) seq <- spell[j,3]
					if (eqtime==2) seq <- spell[j+1,3]
					if (eqtime==3) seq <- ""
				} else {
					if (dur<0 | is.na(dur)) {
						warning("negative or invalid duration for id ",lid[i],", spell ",j)
						seq <- "NA"
					}
					else seq <- paste(rep(spell[j,3],dur),collapse="-")
				}

				if (j==1) trans[i] <- paste(trans[i],seq,sep="") 
				else {
					if (seq!="") trans[i] <- paste(trans[i],seq,sep="-")
				}
			}
		}
	}

	return(trans)
}
