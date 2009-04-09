## =======================
## Within Sequence Entropy
## =======================

seqient <- function(seqdata, norm=TRUE, with.miss=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	statl <- attr(seqdata,"alphabet")
	nbstat <- length(statl)

	message(" [>] computing entropy for ",nrow(seqdata)," sequences, please wait...")

	iseqtab <- seqistatd(seqdata, with.miss)
	
	ient <- apply(iseqtab,1,entropy)
	ient <- as.matrix(ient)
	if (norm==TRUE) {
		emax <- entropy(rep(1/nbstat,nbstat)) 
		ient <- ient/emax
		}

	colnames(ient) <- "Entropy"
	rownames(ient) <- rownames(seqdata)

	return(ient)

	}	
