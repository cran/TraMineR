############################
## Compute distance to center for a group
############################

disscenter <- function(diss, group=NULL, medoids.index=FALSE) {
	trim <- 0
	## max.iter <- 20
	if (inherits(diss, "dist")) {
		diss <- as.matrix(diss)
	}
	if (is.null(group)) {
		group <- integer(nrow(diss))
		group[] <- 1
	}
	ret <- numeric(nrow(diss))
	ind <- 1:nrow(diss)
	grp <- factor(group)
	lgrp <- levels(grp)
	medoids <- numeric(length(lgrp))
	keep=1-trim
	## pour chaque valeur du groupe
	for (i in 1:length(lgrp)) {
		## on crée le groupe en question
		cond <- grp==lgrp[i]
		grpindiv <- sort(ind[cond])
		## on calcul la contribution a l'inertie intraclasse
		dc <- .Call("tmrinertiacontrib", diss, as.integer(grpindiv), PACKAGE="TraMineR")
		dc <- dc-mean(dc)/2
		ret[grpindiv] <- dc
		mindc <- min(dc)
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
		medoids[i] <- sort(which(ret==mindc&cond))[1]
	}
	names(medoids) <- lgrp
	if (medoids.index) return(medoids)
	return(ret)
}