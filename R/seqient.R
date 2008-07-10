## =======================
## Within Sequence Entropy
## =======================

seqient <- function(seqdata, norm=TRUE) {

	if (!inherits(seqdata,"stslist")) {
		cat(" [!] data is NOT a sequence object, see seqdef function to create one\n")
		return()
		}

	statl <- attr(seqdata,"alphabet")
	nbstat <- length(statl)
	iseqtab <- matrix(nrow=seqdim(seqdata)[1],ncol=nbstat)

	cat(" [>] computing within sequence entropy for",seqdim(seqdata)[1],"sequences... \n")

	for (i in 1:nbstat) {
		iseqtab[,i] <- apply(seqdata,1,function(x) sum(x==statl[i],na.rm=TRUE))
	}			
	
	ient <- apply(iseqtab,1,entropy)
	ient <- as.matrix(ient)
	if (norm==TRUE) {
		emax <- entropy(rep(1/nbstat,nbstat)) 
		ient <- ient/emax
		}

	colnames(ient) <- "Entropy"
	rownames(ient) <- paste("[",seq(1:length(ient)),"]",sep="")

	return(ient)

	}	
