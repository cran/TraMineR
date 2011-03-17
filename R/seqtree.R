seqtree <- function(formula, data=NULL, weighted=TRUE, minSize=0.05, maxdepth=5, R=1000, pval=0.01, weight.permutation="replicate", seqdist_arg=list(method="LCS", norm=TRUE), diss=NULL, squared=FALSE, first=NULL){
	formula.call <- formula
	seqdata <- eval(formula[[2]], data, parent.frame()) # to force evaluation
	if (!inherits(seqdata, "stslist")) {
		stop("Right hand of formula should be a stslist object (see seqdef)")
 	}
	if(is.null(diss)){
		seqdist_arg$seqdata <- seqdata
		diss <- do.call(seqdist, seqdist_arg)
	}
	if(weighted){
		weights <- attr(seqdata,"weights") 
	} 
	else {
		weights <- NULL
	}
	formula[[2]] <- diss
	tree <- disstree(formula=formula, data=data, weights=weights ,minSize=minSize, maxdepth=maxdepth, R=R, pval=pval, object=seqdata,weight.permutation=weight.permutation, squared=squared, first=first)
	tree$terms <- terms(formula.call)
	class(tree) <- c("seqtree", class(tree))
	return(tree)
}

print.seqtree <- function(x, digits = getOption("digits") - 2, medoid=TRUE, ...){
	stslistmedoid <- function(object, index) {
		x <- seqconc(object[index,], void=attr(object,"void"))
		x <- suppressMessages(seqformat(x, from='STS', to='SPS', compressed=TRUE))
		return(x)
	}
	print.disstree(x, digits = digits, medoid=medoid, medoid.fun=stslistmedoid,...)
}
