############################
## Compute distance to center for a group
############################

disscenter <- function(diss, group=NULL, medoids.index=NULL, allcenter=FALSE) {
	if(is.logical(medoids.index)){
		if(medoids.index){
			medoids.index <- "First"
		}else{
			medoids.index <- NULL
		}
	}
	retmedoids <- !is.null(medoids.index)
	if (retmedoids) {
		allcenter <- FALSE
	}
	allmedoids <- FALSE
	if(!is.null(medoids.index)) {
		if(medoids.index=="all"){
			allmedoids <- TRUE
		} else if (medoids.index!="first") {
			stop('medoids.index argument should be one of "first", "all" or NULL')
		}
	}
	
	trim <- 0
	## max.iter <- 20
	if (inherits(diss, "dist")) {
		diss <- as.matrix(diss)
	}
	if (is.null(group)) {
		group <- integer(nrow(diss))
		group[] <- 1
	}
	ind <- 1:nrow(diss)
	grp <- factor(group)
	lgrp <- levels(grp)
	if(allcenter){
		ret <- data.frame(numeric(nrow(diss)))
	}else{
		ret <- numeric(nrow(diss))
	}
	if (retmedoids) {
		if (allmedoids) {
			medoids <- list()
		}
		else {
			medoids <- numeric(length(lgrp))
		}
	}
	keep=1-trim
	## pour chaque valeur du groupe
	for (i in 1:length(lgrp)) {
		## on crée le groupe en question
		cond <- grp==lgrp[i]
		grpindiv <- sort(ind[cond])
		## on calcul la contribution a l'inertie intraclasse
		if (allcenter) {
			ret[, i] <- 0
			others <- sort(ind[!cond])
			dT <- .Call("tmrinertiacontribext", diss, grpindiv, others, PACKAGE="TraMineR")
			ret[grpindiv, i] <- dT[1:sum(cond)]
			ret[others, i] <- dT[-(1:sum(cond))]
		}
		else {
			dc <- .Call("tmrinertiacontrib", diss, as.integer(grpindiv), PACKAGE="TraMineR")
			dc <- dc-mean(dc)/2
			
			ret[grpindiv] <- dc
			if (retmedoids) {
				if (trim>0) {
					maxdist <- quantile(dc, probs=keep)
					trimmedcond <- dc<=maxdist
					dT <- .Call("tmrinertiacontribext", diss, grpindiv[cond], grpindiv[!cond], PACKAGE="TraMineR")
					ntrimmed <- sum(trimmedcond)
					dT <- dT-mean(dT[1:ntrimmed])/2
					mindc <- min(dT)
					ret[grpindiv[cond]] <- dT[1:ntrimmed]
					ret[grpindiv[!cond]] <- dT[(ntrimmed+1):length(grpindiv)]
				}
				mindc <- min(dc)
				if (allmedoids) {
					medoids[[i]] <- which(ret==mindc & cond)
				}
				else {
					medoids[i] <- sort(which(ret==mindc & cond))[1]
				}
			}
		}
	}
	if (retmedoids) {
		## No group
		if (length(lgrp)==1) {
			return(medoids[[1]])
		}
		names(medoids) <- lgrp
		return(medoids)
	}
	if (allcenter) {
		names(ret) <- lgrp
	}
	return(ret)
}