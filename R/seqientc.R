
seqientc <- function(seqdata, with.miss=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

	dss <- seqdss(seqdata, with.miss)
	ient <- seqient(seqdata, with.miss, norm=TRUE)
	seql <- seqlength(seqdata)

	seqndss <- seqlength(dss)

	ientc <- (seqndss/seql) * ient

	colnames(ientc) <- "Hw"

	return(ientc)
}
