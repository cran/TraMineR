## CHECK IF THE SEQUENCES ARE IN THE COMPRESSED FORMAT

seqfcheck <- function(seqdata) {
	seqdata <- as.matrix(seqdata)

	if (ncol(seqdata)==1) {
		if (length(grep("-",seqdata))>0) format <- "-"
		else if (length(grep(":",seqdata))>0) format <- ":"
		else format <- "?"
	}
	else {
		if (length(grep("-",seqdata))>0) format <- "-X"
		else format <- "X"
	}

	return(format)
}


