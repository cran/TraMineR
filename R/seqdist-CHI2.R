# Should only be used through seqdist()

CHI2 <- function(seqdata, breaks=NULL, step=1, with.missing=FALSE, norm=TRUE, weighted=TRUE, overlap=FALSE, notC=FALSE, euclid=FALSE){
  if(euclid){
		weighted <- FALSE
	}
	if(is.null(breaks)){
		breaks <- list()
		if(step==1){
			for(i in 1:ncol(seqdata)){
				breaks[[i]] <- c(i, i)
			}
		}else if (step==ncol(seqdata)){
			breaks[[1]] <- c(1, ncol(seqdata))
		}else{
			bb <- seq(from=1, to=ncol(seqdata), by=step)
			if(bb[length(bb)]!=(1+ncol(seqdata)-step)){
				warning("[!] last episode is shorter than the other")
			}
			bb <- c(bb, ncol(seqdata)+1)

			bi <- 1
			if(overlap){
				breaks[[bi]] <- c(1, 1+step/2)
				bi <- 2
			}
			for(i in 2:length(bb)){
				breaks[[bi]] <- c(bb[(i-1)], bb[i]-1)
				bi <- bi +1
				if(overlap){
					breaks[[bi]] <- pmin(breaks[[bi-1]]+step/2, ncol(seqdata))
					bi <- bi +1
				}
			}

		}
	}
	if(!is.list(breaks)){
		stop(" [!] breaks should be a list of couples: breaks=list(c(start1, end1), c(start2, end2), ...))")
	}

	labels <- attr(seqdata, "labels")
	cpal <- cpal(seqdata)
	alph <- alphabet(seqdata)
	if(with.missing){
		labels <- c(labels, "missing")
		cpal <- c(cpal, attr(seqdata, "missing.color"))
		alph <- c(alph, attr(seqdata, "nr"))
	}
	nalph <- length(alph)
	weights <- attr(seqdata, "weights")
	if(is.null(weights)|| !weighted){
		weights <- rep(1, nrow(seqdata))
	}
	seqdata <- as.matrix(seqdata)

	dummies <- function(b){
		lastbi <- breaks[[b]][1]
		bi <- breaks[[b]][2]
		bindice <- lastbi:bi
		blength <- ifelse(norm, length(bindice), 1)
		mat <- matrix(0, nrow=nrow(seqdata), ncol=nalph)
		bseq <- seqdata[, bindice]
		myrowSums <- function(x){
			if(!is.null(ncol(x))){
				return(weights*rowSums(x))
			}else{
				return(weights*x)
			}
		}
		for(i in 1:nalph){
			mat[, i] <- myrowSums(bseq==alph[i])/blength
		}
		return(mat)
	}
	allmat <- list()
	for(b in 1:length(breaks)){
		allmat[[b]] <- dummies(b)
	}
	allmat <- do.call(cbind, allmat)
	ndotj <- colSums(allmat)
	cond <- ndotj>0
	allmat <- allmat[, cond]
	if(euclid){
		ndotj <- rep(1.0, ncol(allmat))
	} else{
		ndotj <- ndotj[cond]
	}
	chdist <- function(i, j){
		return(sqrt(sum((allmat[i, ]-allmat[j, ])^2/ndotj)))
	}
	n <- nrow(seqdata)
	if(notC){ ##Unused
		dd <- numeric(n*(n-1)/2)
		for(i in 1:(n-1)){
			for(j in (i+1):n){
				dd[n*(i-1) - i*(i-1)/2 + j-i] <- chdist(i, j)
			}
		}
	}else{
		## SEXP tmrChisq(SEXP ChiTableS, SEXP tdimS, SEXP margeS)
		dd <- .Call(C_tmrChisq, as.double(allmat), as.integer(dim(allmat)), as.double(ndotj))
	}
	attributes(dd) <- list(Size = n, Labels = rownames(seqdata), Diag = FALSE,
        Upper = FALSE, method = "Chi square sequence", call = match.call(),
        class = "dist")
	return(dd)
}
