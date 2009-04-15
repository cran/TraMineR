
summary.stslist <- function(object,...) {

	alphabet <- alphabet(object)
	nbstates <- length(alphabet)
	cpal <- cpal(object)
	labels <- attr(object,"labels")
	nr <- attr(object,"nr")
	void <- attr(object,"nr")

	nbseq <- seqdim(object)[1]
	seql <- seqlength(object)
	nuseq <- nrow(unique(object))

	cat(" [>] dimensionality of the sequence space:", (nbstates-1)*max(seql),"\n")
	cat(" [>]", nbseq, "sequences in the data set","\n")
	cat(" [>]", nuseq, "unique sequences in the data set","\n")
	cat(" [>] min/max sequence length:",min(seql),"/",max(seql),"\n")
	cat(" [>] alphabet:",paste(1:nbstates,alphabet,collapse=" ",sep="="),"\n")
	cat(" [>] colors:",paste(1:nbstates,cpal,collapse=" ",sep="="),"\n")
	cat(" [>] labels:",paste(1:nbstates,labels,collapse=" ",sep="="),"\n")
	cat(" [>] code for missing statuses:",nr,"\n")
	cat(" [>] code for void state:",void,"\n")
}

