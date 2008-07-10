## ========================================
## Change the alphabet of a sequence object
## ========================================

seqnum <- function(seqdata) {
	if (!inherits(seqdata,"stslist")) {
		stop(" [!] data is not a sequence object, see seqdef function to create one\n")
	return()
	}

	old <- attr(seqdata,"alphabet")
	nbstat <- length(old)

	alphabet <- 0:(nbstat-1)

	if (length(alphabet)!=nbstat) {
		stop("lengths of old and new alphabet are different")
		} 

	for (i in 1:seqdim(seqdata)[2])
		seqdata[,i] <- factor(seqdata[,i], levels=old, labels=alphabet)

	attr(seqdata,"alphabet") <- alphabet

	return(seqdata)

	}
