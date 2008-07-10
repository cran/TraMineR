## ========================================
## Defining a set of sequences as an object
## ========================================

seqdef <- function(data, var=NULL, informat='STS', stsep='-', alphabet=NULL, states=NULL, start=1, missing=NA, cnames=NULL, cpal=NULL, labels=NULL) {
	seqdata <- seqxtract(data, var, data.frame=TRUE)

	if (informat=='STS') {
		sf <- seqfcheck(seqdata)
		if (sf %in% c("-",":")) seqdata <- seqdecomp(seqdata,sep=sf)
		else if (sf=="-X") 
			message(" [!] found '-' character in states codes, not recommended.") 
	}
	else if (informat %in% c("SPS","SPS2","SPELL")) {
		seqdata <- seqdecomp(seqformat(seqdata,from=informat,to='STS',stsep=stsep))
		## if (is.null(cnames)) cnames <- colnames(seqdata)
	}

	seqdata <- as.data.frame(seqdata)

	if (!is.na(missing)) seqdata[seqdata==missing] <- NA
 
	## DEFINING THE CLASS AND SOME ATTRIBUTES
	class(seqdata) <- c("stslist","data.frame")

	if (is.null(start)) start <- 1
	attr(seqdata,"start") <- start
	attr(seqdata,"missing") <- missing

	## DEFINING THE ALPHABET
	if (is.null(alphabet)) plevels <- seqstatl(seqdata)
	else plevels <- alphabet

	if (is.null(states)) {
		message(" [>] distinct states appearing in the data: ",paste(seqstatl(seqdata),collapse="/"))
		attr(seqdata,"alphabet") <- plevels
		}
	else {
		## plevels <- states
		attr(seqdata,"alphabet") <- states
	}

	alphabet <- attr(seqdata,"alphabet")
	nbstates <- length(alphabet)

	if (nbstates==1) stop("alphabet has only one state!")

	message(" [>] alphabet: ",paste(1:nbstates,alphabet,collapse=" ",sep="="))

	## COLOR PALETTE AND STATES LABELS ATTRIBUTES
	if (is.null(cpal)) {
		if (nbstates==2) attr(seqdata,"cpal") <- brewer.pal(3,"Accent")[1:2]
		else if (nbstates<=8) attr(seqdata,"cpal") <- brewer.pal(nbstates,"Accent")
		else if (nbstates>8 & nbstates<=12) attr(seqdata,"cpal") <- brewer.pal(nbstates,"Set3")
	}
	else attr(seqdata,"cpal") <- cpal

	if (!is.null(labels)) attr(seqdata,"labels") <- labels
	else attr(seqdata,"labels") <- alphabet
	
	## ========================
	nbseq <- seqdim(seqdata)[1]
	seql <- seqlength(seqdata)

	message(" [>] ", nbseq, " sequences in the data set")
	message(" [>] min/max sequence length: ",min(seql),"/",max(seql))
	
	for (i in 1:seqdim(seqdata)[2]) {
		seqdata[,i] <- factor(seqdata[,i], 
			levels=plevels,
			labels=attr(seqdata,"alphabet"))
		}

	## COLUMNS NAMES
	if (!is.null(cnames)) colnames(seqdata) <- cnames
	else {
		cntmp <- colnames(seqdata)
		if (is.null(cntmp)) 
			colnames(seqdata) <- paste("T",start:(max(seql)+start-1),sep="")
		else if (all(is.na(cntmp)!=TRUE)) 
			colnames(seqdata) <- cntmp
		else colnames(seqdata) <- paste("T",start:(max(seql)+start-1),sep="")
	}

	## ROW NAMES
	rownames(seqdata) <- paste("[",1:seqdim(seqdata)[1],"]",sep="")
	
	return(seqdata)
}



