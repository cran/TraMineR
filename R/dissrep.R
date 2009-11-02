## =============================================
## Representative sequence of a set of sequences
## =============================================

dissrep <- function(dist.matrix, criterion="density", score=NULL, decreasing=TRUE, trep=0.25, nrep=NULL, tsim=0.10, 
	dmax=NULL) {

	nbobj <- nrow(dist.matrix)

	## Max theoretical distance
	if (is.null(dmax)) dmax <- max(dist.matrix)

	if (tsim<0 || tsim>1)
		stop("tsim must be between 0 and 1", call.=FALSE)
	tsim <- dmax*tsim

	message(" [>] max. distance: ", round(dmax,2))
	message(" [>] neighborhood diameter: ", round(tsim,2))

	## =====================
	## Neighbourhood density
	## =====================
	if (criterion=="density") {
 
		neighbours <- dist.matrix<tsim
		score <- rowSums(neighbours)
		decreasing <- TRUE
	} 
	## =====================
	## Neighbourhood density
	## =====================
	else if (criterion=="freq") {
 
		neighbours <- dist.matrix==0
		score <- rowSums(neighbours)

		decreasing <- TRUE
	} 
	## ============
	## Min distance
	## ============
	else if (criterion=="dist") {
		## Sum of distances for all distinct sequences	
		score <- rowSums(dist.matrix)
		decreasing <- FALSE
	}
	else if (is.null(score))
		stop("Unknown criterion / no score provided")

	## ===========================
	## Sorting candidates by score
	## ===========================
	if (length(score)!=nbobj)
		stop("Score must be a vector of length equal to",nbobj)

	score.sort <- order(score, decreasing=decreasing)
	## score <- score[score.sort]

	rep.dist <- dist.matrix[score.sort, score.sort]
	## rep.dist <- rep.dist[score.sort, score.sort]
	
	## ==========================
	## Selecting representatives
	## ==========================
	idx <- 1
	idxrep <- 1
	
	## Coverage fixed	
	if (is.null(nrep) && trep>0) {
		pctrep <- 0

		while (pctrep<trep && idx < nbobj) {
			tempm <- as.matrix(rep.dist[, idxrep])

			nbnear <- sum(rowSums(tempm<tsim)>0)
			pctrep <- nbnear/nbobj

			## Searching for next non-redundant sequence in candidate list
			idx <- idx+1
			while (any(rep.dist[idx, idxrep]<tsim) && idx < nbobj) {
				## cat("distance to next:", round(rep.dist[idxrep, idx],2), "\n")
				idx <- idx+1 
			}
			## message(" [>] Adding object with index", score.sort[idx], " to representative set")
			idxrep <- c(idxrep, idx)
		}

		## If the last object is non-redundant, we add it
		if (idx==nbobj && all(rep.dist[idx, idxrep]>tsim))
			idxrep <- c(idxrep, idx)		

		nbkeep <- length(idxrep)
		message(" [>] ", nbkeep, " representative(s) selected, coverage=",
			round(pctrep,2)*100,"% (threshold=",round(trep,2)*100,"%)")
	}
	## Number of desired representative fixed
	else {
		repcount <- 1

		while (repcount<nrep && idx<=nbobj) {
			## Searching for next non-redundant sequence in candidate list
			idx <- idx+1
			while (any(rep.dist[idx, idxrep]<tsim)) {
				## cat("distance to next:", round(rep.dist[idxrep, idx],2), "\n")
				idx <- idx+1
			}
			## message(" [>] Adding object with index", score.sort[idx], " to representative set")
			idxrep <- c(idxrep, idx)
			repcount <- repcount+1
		}

		nbkeep <- length(idxrep)
		message(" [>] ", nbkeep, " representative(s) selected")
	}

	## On force avec as.matrix car sinon il y a une erreur
	## si 1 seule colonne
	dist.repseq <- as.matrix(dist.matrix[,score.sort[idxrep]])

	## ================
	## Quality measures
	## ================

	## Keeping distance to the nearest representative sequence only
	dc.tot <- disscenter(dist.matrix)
	if (nbkeep>1) {
		tied <- 0
		minidx <- apply(dist.repseq,1, which.min)
		
		for (i in 1:nbobj) {
			dist.repseq[i,-minidx[i]] <- NA
		}
	
		## message(" [>] ",tied," observations equaly distant from two or more representatives")

		na <- colSums(!is.na(dist.repseq))
		MD <- colMeans(dist.repseq, na.rm=TRUE)
		SD <- colSums(dist.repseq, na.rm=TRUE)

		## Number of similar sequences dist.repsequred by each representative
		nb <- colSums(dist.repseq < (dmax*0.10), na.rm=TRUE)
		
		## Sum of distances to global center
		DC <- matrix(nrow=nbkeep,ncol=1)
		for (i in 1:nbkeep)
			DC[i] <- sum(dc.tot[!is.na(dist.repseq[,i])])
		
		## Inertia
		V <- matrix(nrow=nbkeep,ncol=1)
		for (i in 1:nbkeep) {
			tmp <- dist.matrix[!is.na(dist.repseq[,i]),!is.na(dist.repseq[,i])]
			tmp <- as.matrix(tmp)	
			V[i] <- mean(disscenter(tmp))
		}

	} else {
		na <- sum(!is.na(dist.repseq))
		MD <- mean(dist.repseq, na.rm=TRUE)
		SD <- sum(dist.repseq, na.rm=TRUE)
		nb <- sum(dist.repseq < tsim, na.rm=TRUE)
		DC <- sum(dc.tot)
		V <- mean(dc.tot)
	}

	quality <- (sum(dc.tot)-sum(dist.repseq, na.rm=TRUE))/sum(dc.tot)
	
	## Overall stats
	na <- c(na,sum(na))
	nb <- c(nb,sum(nb))
	MD <- c(MD,mean(dist.repseq, na.rm=TRUE))
	SD <- c(SD, sum(dist.repseq, na.rm=TRUE))
	DC <- c(DC, sum(dc.tot))
	V <- c(V, mean(dc.tot))

	pcta <- (na/nbobj)*100
	pctb <- (nb/nbobj)*100
	Q <- (DC-SD)/DC*100	

	stats <- data.frame(na, pcta, nb, pctb, SD, MD, DC, V, Q)
	colnames(stats) <- c("na", "na(%)", "nb", "nb(%)", "SD", "MD", "DC", "V", "Q")
	rownames(stats) <- c(paste("r",1:nbkeep,sep=""), "Total")

	## ============
	## Final object
	## ============
	res <- score.sort[idxrep]
	class(res) <- c("diss.rep", class(res))

	attr(res, "n") <- nbobj
	attr(res, "criterion") <- criterion
	attr(res, "dmax") <- dmax
	attr(res, "Scores") <- score
	attr(res, "Distances") <- dist.repseq
	attr(res, "Statistics") <- stats
	attr(res, "Quality") <- quality 

	return(res)
}	
	
