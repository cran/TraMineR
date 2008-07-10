## =============================================
## Number of distinct subsequences in a sequence
## =============================================

nsubs <- function (x, nbstat, statlist) {
		l <- vector(mode="integer", nbstat)
		slength <- sum(!is.na(x))
		N <- vector(mode="integer",(slength+1))
		N[1] <- 1

		for (i in 2:(slength+1)) {
			N[i] <- 2*N[i-1]
			cidx <- which(statlist==x[i-1])

			if (l[cidx] > 0) N[i] <- N[i] - N[l[cidx]]
			l[cidx] <- i-1
		}
		return(N[(slength+1)])
	}


seqsubsn <- function(seqdata, DSS=TRUE) {

	if (!inherits(seqdata,"stslist")) {
		stop("data is not a sequence object, see seqdef function to create one")
	}

	nbseq <- seqdim(seqdata)[1]

	if (DSS==TRUE) seqdata <-  suppressMessages(seqdss(seqdata))

	## alphabet
	sl <- attr(seqdata,"alphabet")
	ns <- length(sl)
		
	result <- apply(seqdata,1,nsubs,nbstat=ns,statlist=sl)
	
	names(result) <- NULL
	return(result)
}

