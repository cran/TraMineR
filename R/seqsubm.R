## ============================
## Matrix of substitution costs
## ============================

seqsubm <- function(seqdata, method, cval=2, with.missing=FALSE, 
	miss.cost=2, time.varying=FALSE, weighted=TRUE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	metlist <- c("CONSTANT","TRATE")
	if (missing(method) || !method %in% metlist)
		stop("method must be one of: ", paste(metlist,collapse=" "))

	alphabet <- attr(seqdata,"alphabet")

	## Adding an entry for for missing state
	if (with.missing) {
		message(" [>] setting ",miss.cost," as substitution cost for missing values")
		alphabet <- c(alphabet,attr(seqdata,"nr"))
	}

	alphsize <- length(alphabet)

	if (method=="CONSTANT") {
		if (is.na(cval)){
			stop("no value for the constant substitution-cost")
		}
		if (time.varying) {
			time <- ncol(seqdata)
			
			message(" [>] creating ",alphsize,"x",alphsize,"x",time,
				" time varying substitution-cost matrix using ",
				cval," as constant value")
			costs <- array(cval, dim=c(alphsize, alphsize, time))
			for (i in 1:time) {
				diag(costs[,,i]) <- 0
			}
		}
		else {
			message(" [>] creating ",alphsize,"x",alphsize,
				" substitution-cost matrix using ",cval," as constant value")

			costs <- matrix(cval,nrow=alphsize,ncol=alphsize)
			diag(costs) <- 0
		}
	}

	if (method=="TRATE") {
		if (time.varying) {
			message(" [>] creating time varying substitution-cost matrix using transition rates ...")
			tr <- seqtrate(seqdata, time.varying=TRUE, weighted=weighted)
			tmat <- nrow(tr)
			time <- ncol(seqdata)
			costs <- array(0, dim=c(alphsize, alphsize, time))
			## Function to compute the cost according to transition rates
			tratecost <- function(trate, time, state1, state2, debut, fin){
				cost <- 0
				if (!debut) { ## Premier état
					cost <- cost - trate[state1,state2, time-1] - trate[state2, state1, time-1]
				}
				if (!fin) { ##Dernier Etat
					cost <- cost - trate[state1,state2, time] - trate[state2, state1, time]
				}
				if (!debut && !fin) {
					return(cost + 4)
				}
				else{
					return(4 + 2*cost)
				}
			}
			for (t in 1:time) {
				for (i in 1:(tmat-1)) {
					for (j in (i+1):tmat) {
						cost <- tratecost(tr, t, i, j, t==1, t==time)
						costs[i, j, t] <- cost
						costs[j, i, t] <- cost
					}
				}
			}
		}
		else {
			message(" [>] creating substitution-cost matrix using transition rates ...")
			tr <- seqtrate(seqdata, time.varying=FALSE, weighted=weighted)
			tmat <- nrow(tr)
			costs <- matrix(nrow=alphsize,ncol=alphsize)
			diag(costs) <- 0
			for (i in 1:tmat) {
				for (j in 1:tmat) {
					if (i!=j)
						costs[i,j] <- 2-tr[i,j]-tr[j,i]
				}
			}
		}
	}

	## 
	if (with.missing) {
		if (time.varying) {
			costs[alphsize,1:(alphsize-1),] <- miss.cost
			costs[1:(alphsize-1),alphsize,] <- miss.cost
		}
		else {
			costs[alphsize,1:(alphsize-1)] <- miss.cost
			costs[1:(alphsize-1),alphsize] <- miss.cost
		}
	}

	## Setting rows and columns labels
	rclab <- paste(alphabet,"->",sep="")
	if (time.varying) {
		dimnames(costs) <- list(rclab, rclab, colnames(seqdata))
	}
	else {
		dimnames(costs) <- list(rclab,rclab)
	}

	return(costs)
}

