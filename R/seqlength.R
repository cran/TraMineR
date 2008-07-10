## ===================================
## Returns the length of the sequences 
## ===================================

seqlength <- function(seqdata) {
	if (!inherits(seqdata,"stslist")) {
		cat(" => data is not a sequence object, see function seqdef to create one\n")
		return()
	}
	
	missing <- attr(seqdata,"missing")

	if (!is.na(missing)) seqdata[seqdata==missing] <- NA

	if (seqdim(seqdata)[1]==1) {
		sl <- sum(!is.na(seqdata))
	} else {
		sl <- apply(seqdata,1, function(x) sum(!is.na(x)))
	}

	return(sl)
}
