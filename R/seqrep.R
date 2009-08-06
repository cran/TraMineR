## =============================================
## Representative sequence of a set of sequences
## =============================================

seqrep <- function(seqdata, criterion="freq", trep=0.25, nrep=NULL, tsim=0.10, 
	dmax=NULL, dist.matrix=NULL, ...) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one", call.=FALSE)

	slength <- ncol(seqdata)
	statelist <- alphabet(seqdata)
	nbseq <- nrow(seqdata)

	## State distribution
	freq <- seqstatd(seqdata)$Frequencies

	## Distance matrix
	if (missing(dist.matrix) || is.null(dist.matrix))
		dist.matrix <- seqdist(seqdata, ...)

	## ==============
	## Modal sequence
	## ==============
	if (criterion=="freq") {
		useq <- unique(seqdata)
		mcorr <- match(seqconc(useq), seqconc(seqdata))
		
		## Finding the most frequent sequence	
		ctype <- seqtab(seqdata, tlim=0)

		tab <- attr(ctype,"freq")
		score <- tab$Percent
		names(score) <- rownames(seqdata)[mcorr]

		## Converting to a sequence object
		## ctype <- suppressMessages(seqdecomp(row.names(tab)))
		## ctype <- suppressMessages(seqdef(ctype, alphabet=alphabet(seqdata), labels=attr(seqdata,"labels")))

		## Tri des séquences uniques selon leur pourcentage
		score.sort <- match(seqconc(ctype), seqconc(useq))
		ctype <- useq[score.sort,]
	}  
	## ====================
	## Max representativity
	## ====================
	else if (criterion=="mscore") {
		message(" [!] criterion still under development")

		## Computing the scores
		useq <- unique(seqdata)
		mcorr <- match(seqconc(useq), seqconc(seqdata))

		score <- apply(useq,1, TraMineR:::TraMineR.mscore, slength, statelist, freq)
		score.sort <- order(score, decreasing=TRUE)
		
		## Ordering by score
		score <- score[score.sort]
		ctype <- useq[score.sort,]
	} 
	## ============
	## Min distance
	## ============
	else if (criterion=="dist") {

		## Selecting distinct sequences only and saving the indexes 
		useq <- unique(seqdata)
		mcorr <- match(seqconc(useq), seqconc(seqdata))

		## Sum of distances for all distinct sequences	
		score <- apply(dist.matrix[mcorr,], 1, sum)

		score.sort <- order(score, decreasing=FALSE)

		score <- score[score.sort]
		ctype <- useq[score.sort,]
	}
	## ===============
	## Max probability
	## ===============
	else if (criterion=="prob") {

		logprob <- TraMineR:::seqlogp(seqdata)

		## Selecting distinct sequences only and saving the indexes 
		useq <- unique(seqdata)
		mcorr <- match(seqconc(useq), seqconc(seqdata))

		score <- logprob[mcorr]
		score.sort <- order(score, decreasing=FALSE)
		
		## Ordering by score
		score <- score[score.sort]
		ctype <- useq[score.sort,]

	} 
	## =====================
	## Neighbourhood density
	## =====================
	else if (criterion=="density") {
		useq <- unique(seqdata)
		mcorr <- match(seqconc(useq), seqconc(seqdata))

		neighbours <- dist.matrix[mcorr,]<tsim
		score <- apply(neighbours, 1, sum)

		score.sort <- order(score, decreasing=TRUE)

		score <- score[score.sort]
		ctype <- useq[score.sort,]
	}
	else 
		stop("Unknown criterion")
	
	## ============
	## Final object
	## ============
	
	## Occurence of the representative sequence
	nds <- nrow(ctype)
	message(" [>] ", nds, " distinct sequence(s)")

	## ==========================
	## Size of the candidate list
	## ==========================
	if (is.null(nrep) && trep>0) {
		s <- round(nbseq*trep,0)
		cum <- 0
		idx <- 1
		while (cum<s) {
			cum <- cum+length(seqfind(ctype[idx,], seqdata))
			idx <- idx+1
		}

		message(" [>] Selecting ", idx, " distinct representative sequence(s) (threshold=",trep,")")
	}
	else {
		idx <- nrep
		message(" [>] Selecting ", idx, " distinct representative sequence(s) (nrep=",nrep,")")
	}

	ctype.dist <- dist.matrix[mcorr, mcorr]
	ctype.dist <- ctype.dist[score.sort, score.sort]
	ctype.dist <- ctype.dist[1:idx,1:idx]

	ctype <- ctype[1:idx,]
	if (is.null(dmax)) dmax <- max(dist.matrix)
	keep <- 1

	## =============================	
	## Excluding redundant sequences
	## =============================
	if (idx>1) {
		for (i in 2:idx) {
			if (all(ctype.dist[i, keep] > (dmax*tsim)))
				keep <- c(keep, i)
		}
	}

	nbkeep <- length(keep)

	message(" [>] ", nbkeep, " distinct representative sequence(s) remaining")
	ctype <- ctype[keep,]		
	
	## Index(es) of the representative sequences 	
	repidx <- seqfind(ctype, seqdata)
	repmatch <- match(seqconc(ctype), seqconc(seqdata))

	## message(" [>] ", length(ctype), " representative sequence(s) present in the data")

	## On force avec as.matrix car sinon il y a une erreur
	## si 1 seule colonne
	dist.repseq <- as.matrix(dist.matrix[,repmatch])

	## ================
	## Quality measures
	## ================

	## Keeping distance to the nearest representative sequence only
	dc <- disscenter(dist.matrix)
	if (nbkeep>1) {
		tied <- 0
		for (i in 1:nbseq) {
			minidx <- which(dist.repseq[i,]==min(dist.repseq[i,]))
			if (length(minidx)>1)
				tied <- tied+1
			dist.repseq[i,-minidx[1]] <- NA
		}
	
		message(" [>] ",tied," observations equaly distant from two or more representatives")

		nbcapt <- colSums(!is.na(dist.repseq))
		dmean <- colMeans(dist.repseq, na.rm=TRUE)
		dsum <- colSums(dist.repseq, na.rm=TRUE)

		## Number of similar sequences dist.repsequred by each representative
		nbprox <- colSums(dist.repseq < (dmax*0.10), na.rm=TRUE)
		
		## Sum of distances to global center
		dcenter <- matrix(nrow=nbkeep,ncol=1)
		for (i in 1:nbkeep)
			dcenter[i] <- sum(dc[!is.na(dist.repseq[,i])])
		
		## Inertia
		inertia <- matrix(nrow=nbkeep,ncol=1)
		for (i in 1:nbkeep) {
			tmp <- dist.matrix[!is.na(dist.repseq[,i]),!is.na(dist.repseq[,i])]
			tmp <- as.matrix(tmp)	
			inertia[i] <- mean(disscenter(tmp))
		}

	} else {
		nbcapt <- sum(!is.na(dist.repseq))
		dmean <- mean(dist.repseq, na.rm=TRUE)
		dsum <- sum(dist.repseq, na.rm=TRUE)
		nbprox <- sum(dist.repseq < (dmax*0.10), na.rm=TRUE)
		dcenter <- sum(dc)
		inertia <- mean(dc)
	}

	## Overall stats
	nbcapt <- c(nbcapt,sum(nbcapt))
	nbprox <- c(nbprox,sum(nbprox))
	dmean <- c(dmean,mean(dist.repseq, na.rm=TRUE))
	dsum <- c(dsum, sum(dist.repseq, na.rm=TRUE))
	dcenter <- c(dcenter, sum(dc))
	inertia <- c(inertia, mean(dc))

	pctcapt <- (nbcapt/nbseq)*100
	pctprox <- (nbprox/nbseq)*100
	rindex <- (dcenter-dsum)/dcenter*100	

	quality <- data.frame(nbcapt, pctcapt, nbprox, pctprox, dmean, dsum, dcenter, inertia, rindex)
	rownames(quality) <- c(1:nbkeep, "Total")

	rindex <- (sum(dc)-sum(quality$dsum))/sum(dc)
	
	## ============
	## Final object
	## ============
	res <- ctype
	rownames(res) <- paste("[",1:nrow(ctype),"]", sep="")
	class(res) <- c("stslist.rep", class(res))

	attr(res, "nbseq") <- nrow(seqdata)
	attr(res, "criterion") <- criterion
	attr(res, "dmax") <- dmax
	attr(res, "Index") <- repidx
	attr(res, "Scores") <- score
	attr(res, "Distances") <- dist.repseq
	attr(res, "Quality") <- quality
	attr(res, "rindex") <- rindex 

	return(res)
}	
	
