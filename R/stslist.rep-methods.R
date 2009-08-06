## ===========================
## Methods for stsstatd objects
## ===========================

print.stslist.rep <- function(x, ...) {
	criterion <- attr(x,"criterion")
	nbseq <- attr(x,"nbseq")
	rindex <- attr(x,"rindex")

	cat("\n [>] criterion:",criterion,"\n")
	cat(" [>]", nbseq,"sequence(s) in the original data set\n")
	cat(" [>]", nrow(x),"representative sequences\n")
	cat(" [>] overall quality:", rindex,"\n\n")
	NextMethod(x,...)
}

summary.stslist.rep <- function(object, ...) {
	criterion <- attr(object,"criterion")
	nbseq <- attr(object,"nbseq")
	rindex <- attr(object,"rindex")

	cat("\n [>] criterion:",criterion,"\n")
	cat(" [>]", nbseq,"sequence(s) in the original data set\n")
	cat(" [>]", nrow(object),"representative sequences\n")

	cat("\n [>] statistics for the representative set:\n\n")
	print(attr(object,"Quality"), digits=3, ...)
}
