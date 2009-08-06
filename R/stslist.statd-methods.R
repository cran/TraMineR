## ===========================
## Methods for stsstatd objects
## ===========================

print.stslist.statd <- function(x, digits=2, ...) {

	## computing max length of state label
	## to align valid states and entropy tables
	rl <- max(nchar(rownames(x$Frequencies)))
	## width <- max(nchar(colnames(x$Frequencies)), 4)	

	cat(rep(" ",rl),"[State frequencies]\n")
	print(x$Frequencies, digits=digits)

	VS <- t(as.matrix(x$ValidStates))
	rownames(VS) <- paste("N",rep(" ",rl-1),sep="")
	cat("\n", rep(" ",rl),"[Valid states]\n")
	print(VS, digits=digits)
	
	H <- t(as.matrix(x$Entropy))	
	rownames(H) <- paste("H",rep(" ",rl-1),sep="")
	cat("\n", rep(" ",rl),"[Entropy index]\n")
	print(H, digits=digits)
}

"[.stslist.statd" <- function(...) {
	stop(" [!] Operation not allowed", call.=FALSE)
} 
