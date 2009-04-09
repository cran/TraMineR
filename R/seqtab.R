## =========================
## Sequences frequency table
## =========================

seqtab <- function(seqdata, tlim=0, format="SPS") {

	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, see seqdef function to create one")

	## Eliminating empty sequences 
	seqdata <- seqdata[rowSums(seqdata!=attr(seqdata,"nr"))!=0,]

	if (seqfcheck(seqdata)=="-X") 
		warning("'-' character in states codes may cause invalid results")

	if (format=="SPS")	
		seqdata <- suppressMessages(seqformat(seqdata,from='STS', to='SPS', compressed=TRUE))
	else if (format=="STS")
		seqdata <- seqconc(seqdata)
	else 
		stop("Format must be one of: STS or SPS")

	fseq <- table(seqdata)

	if (tlim==0) tlim=length(fseq)
	Freq <- sort(fseq,decreasing=TRUE)
	Percent <- sort(fseq/sum(fseq),decreasing=TRUE)*100
	Percent <- round(Percent,2)
	
	return(data.frame(Freq,Percent)[1:tlim,])
}
	 
