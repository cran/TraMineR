## ============================
## Matrix of substitution costs
## ============================

seqsubm <- function(seqdata, method, cval=2, with.miss=FALSE, miss.cost=2) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	metlist <- c("CONSTANT","TRATE")
	if (!method %in% metlist)
		stop("method must be one of: ", paste(metlist,collapse=" "))

	alphabet <- attr(seqdata,"alphabet")
	if (with.miss) {
		message(" [>] setting ",miss.cost," as substitution cost for missing values")
		alphabet <- c(alphabet,attr(seqdata,"nr"))
	}

	alphsize <- length(alphabet)

	if (method=="CONSTANT") {
		message(" [>] creating ",alphsize,"x",alphsize," substitution-cost matrix using ",cval," as constant value")
		if (is.na(cval))
				stop("no value for the constant substitution-cost")

		costs <- matrix(cval,nrow=alphsize,ncol=alphsize)
		diag(costs) <- 0
	}

	if (method=="TRATE") {
		message(" [>] creating substitution-cost matrix using transition rates ...")
		tr <- seqtrate(seqdata)
		tmat <- nrow(tr)
		costs <- matrix(nrow=alphsize,ncol=alphsize)
		for (i in 1:tmat) {
			for (j in 1:tmat) {
				if (i==j) costs[i,j] <- 0
				else costs[i,j] <- 2-tr[i,j]-tr[j,i]
			}
		}
	}
	
	if (with.miss) {
			costs[alphsize,] <- miss.cost
			costs[,alphsize] <- miss.cost
	}

	## Setting rows and columns labels
	rclab <- paste(alphabet,"->",sep="")
	dimnames(costs) <- list(rclab,rclab)

	return(costs)
}

