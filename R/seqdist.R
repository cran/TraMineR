## ====================================================
## Computing distances between sequences
## Available metrics (method):
## OM = optimal matching
## LCP = Longest Common Prefix (Elzinga)
## LCS = Longest Common Subsequence (Elzinga)
## ====================================================

seqdist <- function(seqdata, method, refseq=NULL, 
	norm=FALSE, indel=1, sm=NA,
	with.miss=FALSE, full.matrix=TRUE) {

	## ======
	## CHECKS
	## ======
	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, use 'seqdef' function to create one")

	metlist <- c("OM","LCP", "LCS", "LCPinv")
	if (!method %in% metlist) 
		stop("Method must be one of: ", paste(metlist,collapse=" "))

	n <- seqdim(seqdata)[1]
	alphabet <- attr(seqdata,"alphabet")
	alphsize <- length(alphabet)
	message(" [>] ",n," sequences with ", alphsize,
		" distinct events/states (", paste(alphabet,collapse="/"),")")

	if (with.miss) {
		alphabet <- c(alphabet,attr(seqdata,"nr"))
		alphsize <- length(alphabet)
		message(" [>] including missing value as additional state" )
	}
	else 
		if (any(seqdata==attr(seqdata,"nr")))
			stop("found missing values in sequences, please set 'with.miss' option to nevertheless compute distances")

	## Checking if substitution cost matrix contains values for each state
	if (method=="OM") 
		if (nrow(sm)!=alphsize | ncol(sm)!=alphsize)
			stop("size of substitution cost matrix must be ",alphsize,"x", alphsize)

	## Reference sequence
	if (!is.null(refseq)) {
		if (refseq==0) {
			mfseq <- row.names(seqtab(seqdata,tlim=1))
			message(" [>] using most frequent sequence as reference: ", mfseq)

			idxmfseq <- suppressMessages(which(seqformat(seqdata, 
				from='STS', to='SPS', compressed=TRUE)==mfseq))
			message(" [>] most frequent sequence appears ", length(idxmfseq), " times")

			compseq <- seqdata[idxmfseq[1],]
		}
		if (refseq>0) {
			compseq <- seqdata[refseq,]
			message(" [>] using sequence ",refseq, seqformat(compseq,from="STS",to="SPS")," as reference")
		}
		lcompseq <- sum(!is.na(compseq))
		distmat <- FALSE
	} 
	else
		distmat <- TRUE

	## ==============
	## Preparing data
	## ==============
	seqdata <- seqnum(seqdata, with.miss=with.miss)

	## Selecting distinct sequences only and saving the indexes 
	dseq <- unique(seqdata)
	mcorr <- match(seqconc(seqdata),seqconc(dseq))

	nd <- seqdim(dseq)[1]
	message(" [>] ", nd," distinct sequences")

	l <- ncol(dseq)
	slength <- seqlength(dseq)

	dseq <- seqasnum(dseq, with.miss=with.miss)

	message(" [>] min/max sequence length: ",min(slength),"/",max(slength))

	debut <- Sys.time()

	message(" [>] computing distances using ",method, appendLF=FALSE)
	if (norm==TRUE) 
		message(" normalized metric ...",appendLF =FALSE) 
	else 
		message(" metric ...",appendLF =FALSE)
	
	## Function and arguments
	if (distmat==FALSE) {
		m <- vector(mode="numeric", length=nd)
		compseq <- seqasnum(seqnum(compseq))

		if (method=="OM") {
			for (i in 1:nd) 
				m[i] <- levenshtein(dseq[i,], slength[i], compseq, lcompseq, indel,sm,alphsize,norm)
			} 
		else	if (method=="LCP") {
			for (i in 1:nd) 
				m[i] <- LCPdist(dseq[i,],slength[i],compseq,lcompseq,norm)
			}
		else if (method=="LCS") {
			for (i in 1:nd) 
				m[i] <- LCSdist(dseq[i,],slength[i],compseq,lcompseq,norm)
			}

		## Constructing the final distance vector
		mcorr <- match(seqconc(seqdata),seqconc(dseq))
		distances <- m[mcorr]
		names(distances) <- NULL
	}
	else {
		magicSeq <- order(mcorr)
		magicIndex <- c(unique(rank(mcorr, ties.method="min")), nrow(seqdata)+1)-1

		if (method=="OM") {
			## One for OM, 2 for LCP
   			disttype <- as.integer(1)
		}
		else if (method=="LCS") {
			disttype <- as.integer(1)
			if(norm)norm<-as.integer(2)
			sm <- suppressMessages(seqsubm(seqdata, method="CONSTANT", cval=2, with.miss=with.miss, miss.cost=2))
			indel <- 1
			disttype <- as.integer(1)
#			print(norm)
		}		
		else if (method=="LCP") {
   		disttype <- as.integer(2) ##One for OM, 2 for LCP
			sm <- 0
			indel <- 0
		}
		else if (method=="LCPinv") {
   		disttype <- as.integer(3) ##One for OM, 2 for LCP
			sm <- 0
			indel <- 0
		}
 	
		distances <- .Call("cstringdistance",
			as.integer(dseq),
			as.integer(dim(dseq)),
			as.integer(slength), 
			as.double(indel),
			as.integer(alphsize),
			as.double(sm),
			as.integer(norm),
			as.integer(magicIndex),
			as.integer(magicSeq),
			disttype,
			PACKAGE="TraMineR")

		## Setting some attributes for the dist object
 		class(distances) <- "dist"
		attr(distances,"Size") <- length(magicSeq)
 		attr(distances,"method") <- method
 		attr(distances, "Labels") <- dimnames(seqdata)[[1]]
 		attr(distances, "Diag") <- FALSE
 		attr(distances, "Upper") <- FALSE
	}


	fin <- Sys.time()
	message(" (",round(difftime(fin,debut,units="mins"),2)," minutes)")

	if (full.matrix) return(as.matrix(distances))
	else return(distances)
}

