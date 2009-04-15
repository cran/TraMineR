## =============================================
## Representative sequence of a set of sequences
## =============================================

seqrep <- function(seqdata, method='modseq', dist.matrix, dist.method='LCS', norm=FALSE, 
	indel=1, sm, with.miss = FALSE) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one", call.=FALSE)

	slength <- ncol(seqdata)
	statelist <- alphabet(seqdata)	

	## State distribution
	freq <- seqstatd(seqdata)$Frequencies

	## ==============
	## Modal sequence
	## ==============
	if (method=="modseq") {
		
		## Finding the most frequent sequence	
		tab <- seqtab(seqdata, tlim=1, format="STS")

		## Converting to a sequence object
		seqp <- suppressMessages(seqdecomp(row.names(tab)))
		seqp <- suppressMessages(seqdef(seqp, alphabet=alphabet(seqdata), labels=attr(seqdata,"labels")))
	}  
	## ====================
	## Modal state sequence
	## ====================
	else if (method=="smode") {
		message(" [!] Method still under development")

		seqp <- matrix(nrow=1, ncol=slength)
		colnames(seqp) <- colnames(seqdata)

		## Constructing the transversal modal sequence
		for (i in 1:slength) {
			smax <- which(freq[,i]==max(freq[,i]))[1]
			seqp[,i] <- statelist[smax]
		}
		seqp <- suppressMessages(seqdef(seqp, alphabet=alphabet(seqdata), labels=attr(seqdata,"labels")))
	} 
	## ==================
	## Max score sequence
	## ==================
	else if (method=="mscore") {
		message(" [!] Method still under development")

		## Computing the scores
		mscore <- apply(seqdata,1, TraMineR.mscore, slength, statelist, freq)

		ctype <- which(mscore==max(mscore))
	} 
	## =================
	## Min dist sequence
	## =================
	else if (method=="dist") {

		if (missing(dist.matrix) || is.null(dist.matrix))
			dist.matrix <- suppressMessages(
				seqdist(seqdata, method=dist.method, norm=norm, 
					indel=indel, sm=sm, with.miss=with.miss)
				)

		## Searching for the minimum sum of distances	
		sumdist <- apply(dist.matrix, 1, sum)
		ctype <- which(sumdist==min(sumdist))
	}
	## =================
	## Max prob sequence
	## =================
	else if (method=="prob") {

		logprob <- seqlogp(seqdata)
		
		## Searching for the maximum prob	
		ctype <- which(logprob==min(logprob))
	} else 
		stop("Unknown method")
	
	## ============
	## Final object
	## ============
	
	## Occurence of the representative sequence
	if (method %in% c("dist", "mscore", "prob")) {
		nds <- nrow(unique(seqdata[ctype,]))

		message(" [>] found ", nds, " distinct representative sequence(s)")
		if (nds>1) 
			message(" [>] choosing the first sequence")
		
		seqp <- seqdata[ctype[1],]
	} else if (method %in% c("modseq","smode")) {
		## Index(es) of the representative sequence 	
		ctype <- seqfind(seqp, seqdata)
	}

	message(" [>] ", length(ctype), " representative sequence(s) present in the data")

	## Distance to representative sequence
	if (missing(dist.matrix) || is.null(dist.matrix) || (method=="smode" && length(ctype)==0)) {
		message(" [>] Computing distances to representative sequence using ",dist.method," metric")
		dist.repseq <- suppressMessages(
				seqdist(seqdata, method=dist.method, refseq=seqp, 
					norm=norm, indel=indel, sm=sm, with.miss=with.miss)
				)
	} else
		dist.repseq <- dist.matrix[, ctype[1]]

	## State frequencies of the representative sequence
	for (i in 1:slength) {
		fstate <- which(seqp[,i]==alphabet(seqdata))
		freq[-fstate,i] <- 0
	}

	## Final object
	if (method %in% c("dist","modseq","smode")) {
		repseq <- list(seqp, ctype, freq, as.vector(dist.repseq))
		names(repseq) <- c("Sequence", "Index", "Frequencies", "Distances")
	} else if (method=="mscore") {
		repseq <- list(seqp, ctype, mscore, freq, as.vector(dist.repseq))
		names(repseq) <- c("Sequence", "Index", "Scores", "Frequencies", "Distances")
	} else if (method=="prob") {
		repseq <- list(seqp, ctype, logprob, freq, as.vector(dist.repseq))
		names(repseq) <- c("Sequence", "Index", "LogProb", "Frequencies", "Distances")
	}	



	return(repseq)
}	
	
