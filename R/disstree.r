#################
##DisTreeNode
#################
# Donnée accessible
# predictor : liste des prédicteurs
# dissmatrix : matrices des dissimilarités internes
#

#Donnée interne
# split : predicteur choisi(NULL pour noeux terminaux)
# vardis : variabilité interne
# children : noeud enfant (NULL pour noeux terminaux)
# ind: liste des index des individus du noeuds.
# depth: profondeur du noeud
# label: label du noeud (valeurs du predicteur
# R2: R2 du split, (NULL pour noeux terminaux)
#
DTNInit <- function(ind, vardis, depth, label) {
	node <- list()
	node$label <- label
	node$children <- NULL
	node$split <- NULL
	node$ind <- ind
	node$depth <- depth
	node$vardis <- vardis
	class(node) <- c("DissTreeNode", "list")
	return(node)
}

###########################
## Retrieve leaf belonging
###########################

disstreeleaf <- function(tree) {
	if (!inherits(tree, "disstree"))stop("tree should be a disstree object")
	categorie <- numeric(length(tree$root$ind))
	categorie[] <- -1
	counter <- 0
	catobject <- list(categorie=categorie, counter=counter)
	return(factor(DTNdisstreeleaf(tree$root, catobject)$categorie))
}

###########################
## Internal recursion
###########################
DTNdisstreeleaf <- function(node, co) {
	if (is.null(node$children)) {
		co$categorie[node$ind] <- (co$counter+1)
		catobject <- list(categorie=co$categorie, counter=co$counter+1)
		return(catobject)
	} else {
		co1 <- DTNdisstreeleaf(node$children$left, co)
		return(DTNdisstreeleaf(node$children$right, co1))
	}
}

###########################
## Internal test for computing association
###########################
DTNdissassoc <- function(dmat, grp, indiv, R) {
	## Significativity can be computed comparing only SCres, sort only starting at 750 indiv
	internalDTNdissassoc <- function(grp, ind, indiv, dmat, grp2, use.sort) {
		 if (use.sort) {
				groupe1 <- sort.int(indiv[ind[grp]], method="quick")
				groupe2 <- sort.int(indiv[ind[grp2]], method="quick")
		 } else {
			 groupe1 <- indiv[ind[grp]]
			 groupe2 <- indiv[ind[grp2]]
		 }
		 r1 <- .Call("tmrsubmatrixinertia", dmat, as.integer(groupe1), PACKAGE="TraMineR")
		 r2 <- .Call("tmrsubmatrixinertia", dmat, as.integer(groupe2), PACKAGE="TraMineR")
		 return(-(r1+r2))
	}
	bts <- boot(grp, internalDTNdissassoc, R=R, sim="permutation", stype="i",
			dmat=dmat, indiv=indiv, grp2=(!grp), use.sort=(length(grp)>750))
	return(sum(bts$t[, 1]>bts$t0[1])/bts$R)
}

###########################
## disstree main function
###########################
disstree <- function(formula, data=NULL, minSize=0.05, maxdepth=5, R=1000, pval=0.01) {
	formula.call <- formula
	dissmatrix <- eval(formula[[2]], data, parent.frame()) # to force evaluation
	if (inherits(dissmatrix, "dist")) {
		dissmatrix <- as.matrix(dissmatrix)
 	}
	formula[[2]] <- NULL
	## Model matrix from forumla
	predictor <- model.frame(formula, data, drop.unused.levels = TRUE, na.action=NULL)

	pop <- nrow(dissmatrix)
	if (minSize<1) {
		minSize <- round(pop*minSize)
	}
	if (pop!=nrow(predictor)) {
		stop("dissimilarity matrix and data should be of the same size")
	}
	vardis <- sum(dissmatrix)/(2*pop*pop)
	tree <- list(formula=formula.call)
	tree$root <- DTNBuildNode(dissmatrix, as.data.frame(predictor), minSize, 1:pop,
			vardis, 1, "Root", R, pval, maxdepth)
	class(tree) <- c("disstree", "list")
	tree$adjustment <- dissassoc(dissmatrix, disstreeleaf(tree), R=R)
	return(tree)
}

###########################
## Building node and finding predictor
###########################
DTNBuildNode <- function(dmat, pred, minSize, ind, vardis,
													depth, label, nbperm, pval, maxdepth) {
	node <- DTNInit(ind, vardis, depth, label)
	SCtot <- node$vardis*length(node$ind)
	SCres <- SCtot
	bestSpl <- NULL
	varnames <- colnames(pred)
	if (depth>=maxdepth) {
			return(node)
	}
	for (p in 1:ncol(pred)) {
		## cat("Checking", varnames[[p]], "...\n")
		spl <- DTNGroupFactorBinary(dmat, SCres, pred[, p], minSize, varnames[[p]], ind)
		if (!is.null(spl) && (is.null(bestSpl) || spl$SCres<bestSpl$SCres)) {
			bestSpl <- spl
			SCres <- spl$SCres
			## cat(varnames[[p]], " Ok", "\n")
		}
	}
	if (is.null(bestSpl)) {
		return(node)
	}
	if (nbperm>1) {
		spval <- DTNdissassoc(dmat, bestSpl$variable, ind, R=nbperm)
		## print(paste(label, bestSpl$varname, spval))
		if (spval>pval)return(node)
	}

	node$split <- bestSpl
	node$children <- list()
	node$children$left <- DTNBuildNode(dmat, as.data.frame(pred[bestSpl$variable, ]),
			minSize, ind[bestSpl$variable], vardis=bestSpl$lvar, depth=depth+1,
			label=paste(bestSpl$llabels, collapse="/"), nbperm, pval, maxdepth)
	node$children$right <- DTNBuildNode(dmat, as.data.frame(pred[!bestSpl$variable, ]),
			minSize, ind[!bestSpl$variable], vardis=bestSpl$rvar, depth=depth+1,
			label=paste(bestSpl$rlabels, collapse="/"), nbperm, pval, maxdepth)
	node$R2 <- 1-(SCres/SCtot)
	return(node)
}

###########################
## Find best binary partition
###########################
DTNGroupFactorBinary <- function(dissmatrix, currentSCres, pred, minSize, varname, ind) {
	totpop <- length(ind)
	grp <- factor(pred, ordered=(is.ordered(pred) || is.numeric(pred)))
	lgrp <- levels(grp)
	if (length(lgrp)<=1)return(NULL)
	grpint <- as.integer(grp)
	nbGrp <- length(lgrp)
	has.na <- FALSE
	llgrp <- lgrp
	## Here we add a group for missing values
	if (sum(is.na(grp))>0) {
		nbGrp <- length(lgrp)+1
		has.na <- TRUE
		llgrp[nbGrp] <- "<Missing>"
	}
	grpCond <- list()
	grpSize <- numeric(length=nbGrp)
	grpSize[] <- 0
	for (i in 1:length(lgrp)) {
		## on crée le groupe en question
		grpCond[[i]] <- (grpint==i)
		grpCond[[i]][is.na(grpCond[[i]])] <- FALSE
		grpSize[i] <- sum(grpCond[[i]])
	}
	## Treating missing values
	if (has.na) {
		grpCond[[nbGrp]] <- is.na(grp)
		grpSize[nbGrp] <- sum(grpCond[[nbGrp]])
	}
	inertiaMat <- matrix(0, nrow=nbGrp, ncol=nbGrp)
	for (i in 1:(nbGrp-1)) {
		grpindiv1 <- ind[grpCond[[i]]]
		for (j in (i+1):nbGrp) {
				grpindiv2 <- ind[grpCond[[j]]]
				r <- .Call("tmrinterinertia", dissmatrix, as.integer(grpindiv1),
						as.integer(grpindiv2), PACKAGE="TraMineR")
				## using only one half of the matrix
				inertiaMat[j, i] <- r
			}
		r <- .Call("tmrsubmatrixinertia", dissmatrix,
				as.integer(grpindiv1), PACKAGE="TraMineR")*sum(grpCond[[i]])
		inertiaMat[i, i] <- r
	}
	## FIXME This step is missing in the loop
	inertiaMat[nbGrp, nbGrp] <- .Call("tmrsubmatrixinertia", dissmatrix,
			as.integer(ind[grpCond[[nbGrp]]]), PACKAGE="TraMineR")*sum(grpCond[[nbGrp]])
	## Computing residuals
	SCres <- sum(diag(inertiaMat)/grpSize)
	if (SCres>currentSCres)return(NULL)
	## Fonction to comput inertia based on inertiaMat
	inertiaFunction <- function(inertiaMat, co, pop) {
		## Take care to add one
		## return(((sum(inertiaMat[co, co])+sum(diag(inertiaMat[co, co])))/2)/pop)
		## New way, inertiaMat is triangular -> we can just sum the matrix
		return(sum(inertiaMat[co, co])/pop)
	}
	bestSCres <- currentSCres
	bestRegroup <- NULL
	allgroups <- 1:(nbGrp)
	if (is.ordered(grp)) {
		maxGrp <- nbGrp-1
	} else {
		maxGrp <- ceiling(nbGrp/2)
	}
	for (p in 1:maxGrp) {
		if (is.ordered(grp)) {
			 combi <- list()
			 combi[[1]] <- 1:p
			 if (has.na) {
				 combi[[2]] <- c(1:p, nbGrp)
			 }
		} else {
			 combi <- combn(nbGrp, p, simplify=FALSE)
		}
		for (co in combi) {
			 popc <- sum(grpSize[co])
			 popothc <- totpop-popc
			 if (popc>minSize && popothc>minSize) {
					othc <- allgroups[!(allgroups %in% co)]
					ico <- inertiaFunction(inertiaMat, co, popc)
					iothc <- inertiaFunction(inertiaMat, othc, popothc)
					SCres <- ico+iothc

					if (SCres<bestSCres) {
						bestSCres <- SCres
						bestRegroup <- list(co=co, othc=othc, ico=ico,
								iothc=iothc, popc=popc, popothc=popothc)
					}
			 }
		}
	}
	if (is.null(bestRegroup)){
		return(NULL)
	}
	ret <- list(
			lpop=bestRegroup$popc,
			rpop=bestRegroup$popothc,
			lvar=bestRegroup$ico/bestRegroup$popc,
			rvar=bestRegroup$iothc/bestRegroup$popothc,
			SCres=bestSCres,
			varname=varname)
	ret$variable <- (grpint %in% bestRegroup$co)
	if (has.na) {
		ret$variable[is.na(grp)] <- (nbGrp %in% bestRegroup$co)
	}
	if (is.ordered(grp)) {
		if (has.na) {
			ret$llabels <- paste("<=", llgrp[max(bestRegroup$co[bestRegroup$co<nbGrp])], sep="")
			ret$rlabels=paste(">", llgrp[max(bestRegroup$co[bestRegroup$co<nbGrp])], sep="")
			if (nbGrp %in% bestRegroup$co) {
				ret$llabels <- paste(ret$llabels, ", ", llgrp[nbGrp], sep="")
			}
			else {
				ret$rlabels <- paste(ret$rlabels, ", ", llgrp[nbGrp], sep="")
			}
		}
		else {
			ret$llabels <- paste("<=", llgrp[max(bestRegroup$co)], sep="")
			ret$rlabels=paste(">", llgrp[max(bestRegroup$co)], sep="")
		}
	}
	else {
		ret$llabels=llgrp[bestRegroup$co]
		ret$rlabels=llgrp[bestRegroup$othc]
	}
	return(ret)
}

###########################
## Print method for disstree
###########################
print.disstree <- function(x, quote=FALSE, digits=3, ...) {
	cat("Dissimilarity tree\n")
	cat(paste("Global R2:", format(x$adjustment$stat$PseudoR2, digits =digits),"\n"))
	print(x$root, quote=quote, digits=digits, ...)

}
###########################
## Internal print method for disstree
###########################
print.DissTreeNode <- function(x, quote=FALSE, digits=3, ...) {
	gap <- character(x$depth)
	gap[] <- "      "
	string <- paste(paste(gap, collapse=" "), "|--", x$label, "[", length(x$ind), "] var:",
	format(x$vardis, digits =digits), collapse="")
	cat(paste(string,"\n"))
	if (!is.null(x$split)) {
		cat(paste(paste(gap, collapse=" "), " ", "|->", x$split$varname, " R2:",
			format(x$R2, digits=digits), "\n", collapse=""))
		for (i in x$children) {
			print(i, quote=quote, digits=digits, ...)
		}
	}
}



