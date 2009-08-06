## =========================
## Computes transition rates
## =========================

seqtrate <- function(seqdata, statl=NULL, time.varying=FALSE, weighted=TRUE) {

	if (!inherits(seqdata,"stslist")) 
		stop(" [!] seqdata is NOT a sequence object, see seqdef function to create one")

	## State list if not specified
	if (is.null(statl)) statl <- attr(seqdata,"alphabet")
	nr <- attr(seqdata,"nr")
	void <- attr(seqdata,"void")

	## Weights
	weights <- attr(seqdata, "weights")
	
	if (!weighted || is.null(weights)) {
		weights <- rep(1, nrow(seqdata))
	}
	nbetat <- length(statl)

	sdur <- ncol(seqdata)

	seqdata <- as.matrix(seqdata)
	
	numtransition <- sdur-1
	## =================================
	## time varying transition rates
	## =================================
	if (time.varying) {
		message(" [>] computing time varying transition rates for states ",paste(statl,collapse="/")," ...",sep="")
		## Dimension names
		dnames <- list()
		dnames[[1]] <- paste("[",statl," ->]",sep="")
		dnames[[2]] <- paste("[-> ",statl,"]",sep="")
		dnames[[3]] <- colnames(seqdata)[-sdur]
		
		tmat <- array(0, dim=c(nbetat, nbetat, numtransition))
		dimnames(tmat) <- dnames
		for (sl in 1:numtransition) {
			for (x in 1:nbetat) {
				PA <- sum(weights[seqdata[,sl]==statl[x] 
					& seqdata[,sl+1]!=nr & seqdata[,sl+1]!=void])
					
				if (PA == 0) {
					tmat[x,,sl] <- 0
				}
				else {
					for (y in 1:nbetat) {
						PAB <- sum(weights[seqdata[,sl]==statl[x] & seqdata[,sl+1]==statl[y]])
						tmat[x,y,sl] <- PAB/PA
					}
				}
			}
		}
	}
	## =================================
	## Non time varying transition rates
	## =================================
	else {
		message(" [>] computing transition rates for states ",paste(statl,collapse="/")," ...",sep="")
		tmat <- matrix(nrow=nbetat, ncol=nbetat)
		row.names(tmat) <- paste("[",statl," ->]",sep="")
		colnames(tmat) <- paste("[-> ",statl,"]",sep="")

		for (x in 1:nbetat) {
			## Count
			PA <- 0
			for (sl in 1:numtransition) {
				PA <- PA + sum(weights[seqdata[,sl]==statl[x] 
					& seqdata[,sl+1]!=nr & seqdata[,sl+1]!=void])
			}
			for (y in 1:nbetat) {
				PAB <- 0
				for (i in 1:numtransition) {
					PAB <- PAB + sum(weights[seqdata[,i]==statl[x] & seqdata[,i+1]==statl[y]])
					}
				if (PA==0) tmat[x,y] <- 0 else tmat[x,y] <- PAB/PA
			}
		}
	}
	return(tmat)
}


