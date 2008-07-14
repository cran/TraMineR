## ============================
## Matrix of substitution costs
## ============================

seqsubm <- function(seqdata, method, cval=NA) {

	if (!inherits(seqdata,"stslist")) {
		stop("data is NOT a sequence object, see seqdef function to create one")
		}

	metlist <- c("CONSTANT","TRATE")
	if (!method %in% metlist) {
			stop("method must be one of: ", paste(metlist,collapse=" "))
		}

	alphabet <- attr(seqdata,"alphabet")
	alphsize <- length(alphabet)

	if (method=="CONSTANT") {
		message("creating ",alphsize,"x",alphsize," substitution-cost matrix using ",cval," as constant value")
		if (is.na(cval)) {
				stop("no value for the constant substitution-cost")
			}
		costs <- matrix(cval,nrow=alphsize,ncol=alphsize,dimnames=list(alphabet,alphabet))
		diag(costs) <- 0
	}

	if (method=="TRATE") {
		message(" [>] creating substitution-cost matrix using transition rates ...")
		tr <- seqtrate(seqdata)
		tmat <- nrow(tr)
		costs <- matrix(nrow=tmat,ncol=tmat,dimnames=list(alphabet,alphabet))
		for (i in 1:tmat) {
			for (j in 1:tmat) {
				if (i==j) costs[i,j] <- 0
				else costs[i,j] <- 2-tr[i,j]-tr[j,i]
			}
		}
	}
	return(costs)
}

