## ========================================
## Computes the sequence turbulence measure
## proposed by Elzinga
## ========================================

turb <- function(x) {
		phi <- x[1] 
		n <- x[2]
		## sum.tx <- x[3]
		mean.tx <- x[3]
		s2.tx <- x[4] 

		s2max <- (n-1) * (1-mean.tx)^2
		
		Tux <- log2(phi* ((s2max+1)/(s2.tx+1)))
		return(Tux)
}

## the var function in R gives the unbiased variance with the n-1 denominator
## but we need the real variance here
realvar <- function(x) {
	n <- sum(!is.na(x))
	var <- 1/n*sum((x - mean(x,na.rm=TRUE))^2,na.rm=TRUE)
	return(var)
	}
	

seqST <- function(seqdata) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see 'seqdef' function to create one")

	message(" [>] extracting symbols and durations...")
	states <- seqdss(seqdata)
	dur <- seqdur(seqdata)

	message(" [>] Computing turbulence for ",seqdim(seqdata)[1],"sequences, please wait...")
	phi <- seqsubsn(states, DSS=FALSE)
	s2.tx <- apply(dur, 1, realvar)
	mean.tx <- apply(dur, 1, mean, na.rm=TRUE)
	sum.tx <- apply(dur, 1, sum, na.rm=TRUE)
	n <- apply(states, 1, function(x) sum(x!=attr(seqdata,"void")))

	tmp <- data.frame(phi, n, mean.tx, s2.tx)
	Tx <- apply(tmp, 1, turb)
	Tx <- as.matrix(Tx)

	rownames(Tx) <- paste("[",seq(1:length(Tx)),"]",sep="")
	colnames(Tx) <- "Turbulence"

	return(Tx)
}

	

 
