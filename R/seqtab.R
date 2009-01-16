## =========================
## Sequences frequency table
## =========================

seqtab <- function(seqdata, tlim=0) {
	if (!inherits(seqdata,"stslist")) {
		stop("data is not a sequence object, see seqdef function to create one")
		}

	## Eliminating empty sequences 
	seqdata <- seqdata[rowSums(!is.na(seqdata))!=0,]

	if (seqfcheck(seqdata)=="-X") 
		warning("'-' character in states codes may cause invalid results")
	
	seqdata <- suppressMessages(seqformat(seqdata,from='STS', to='SPS', compressed=TRUE))
	fseq <- table(seqdata)

	if (tlim==0) tlim=length(fseq)
	Freq <- sort(fseq,decreasing=TRUE)
	Percent <- sort(fseq/sum(fseq),decreasing=TRUE)*100
	Percent <- round(Percent,2)
	
	return(data.frame(Freq,Percent)[1:tlim,])
	}
	 
