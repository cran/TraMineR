## States frequency by time unit

seqstatd <- function(seqdata, digits=2, norm=TRUE) {

	if (!inherits(seqdata,"stslist")) {
		cat(" => data is not a sequence object, see function seqdef to create one\n")
		return()
		}
	
	## Retrieving the alphabet
	statl <- attr(seqdata,"alphabet")
	nbstat <- length(statl)

	seql <- seqdim(seqdata)[2]

	sd <- matrix(nrow=nbstat,ncol=seql)
	row.names(sd) <- statl
	colnames(sd) <- colnames(seqdata)

	for (i in 1:nbstat)
		for (j in 1:seql)
			sd[i,j] <- sum(seqdata[,j]==statl[i],na.rm=TRUE)

	## sd <-	apply(seqdata,2,table)
	N <- apply(sd,2,sum)
	for (i in 1:seql) sd[,i] <- sd[,i]/N[i]

	E <- apply(sd,2,entropy)
	## Maximum entropy is the entropy of the alphabet
	if (norm==TRUE) {
		E.max <- entropy(rep(1/nbstat,nbstat)) 
		E <- E/E.max
	}

	if (!is.null(digits)) {
		sd <- round(sd,digits)
		E <- round(E,digits)
	}

	sd <- list(sd,N,E)
	names(sd) <- c("Frequencies","ValidStates","Entropy")
	return(sd)
	}
