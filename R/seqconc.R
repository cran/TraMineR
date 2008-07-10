## ========================================
## Concatenates vectors of states or events
## into character strings
## ========================================

sconc <- function(seq, sep) {
	seq <- seqtrunc(seq)
	cseq <- paste(seq, collapse=sep)
	}

seqconc <- function (data, var=NULL, sep="-", vname="Sequence") {

	seqdata <- seqxtract(data, var)

	if (seqdim(seqdata)[1]==1) cseq <- sconc(seqdata,sep)
	else cseq <- apply(seqdata, 1, sconc, sep)
	cseq <- as.matrix(cseq)

	## Rows and column names for the output
	rownames(cseq) <- paste("[",seq(1:length(cseq)),"]",sep="")
	colnames(cseq) <- vname

	return(cseq)
	}
