## =========================
## Computes transition rates
## =========================

seqtrate <- function(seqdata, statl=NULL) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is NOT a sequence object, see seqdef function to create one")

	## State list if not specified
	if (is.null(statl)) statl <- attr(seqdata,"alphabet")
	nr <- attr(seqdata,"nr")

	nbetat <- length(statl)
	tmat <- matrix(nrow=nbetat, ncol=nbetat)
	row.names(tmat) <- paste("[",statl," ->]",sep="")
	colnames(tmat) <- paste("[-> ",statl,"]",sep="")

	sdur <- seqdim(seqdata)[2]

	message(" [>] computing transition rates for states ",paste(statl,collapse="/")," ...\n",sep="")

	seqdata <- as.matrix(seqdata)

	for (x in 1:nbetat) {
		## Count
		PA <- 0
		for (sl in 1:(sdur-1))
			PA <- PA + sum(seqdata[,sl]==statl[x] & seqdata[,sl+1]!=nr,na.rm=TRUE)
		for (y in 1:nbetat) {
			PAB <- 0
			for (i in 1:sdur-1) {
				PAB <- PAB + sum(seqdata[,i]==statl[x] & seqdata[,i+1]==statl[y],na.rm=TRUE)
				}
			if (PA==0) tmat[x,y] <- 0 else tmat[x,y] <- PAB/PA
		}
	}
	return(tmat)
}


