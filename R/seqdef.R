## ========================================
## Defining a set of sequences as an object
## ========================================

seqdef <- function(data, var, informat="STS", stsep="-", 
	alphabet, states, start=1, 
	left=NA, right="DEL", gaps=NA, missing=NA, void="%", nr="*",
	cnames, cpal, missing.color="darkgrey", labels, ...) {

	## ===================
	## Extracting the data
	## ===================
	seqdata <- seqxtract(data, var, data.frame=TRUE)

	## ============
	## INPUT FORMAT
	## ============
	if (informat=="STS") {
		sf <- seqfcheck(seqdata)
		if (sf %in% c("-",":")) seqdata <- seqdecomp(seqdata,sep=sf)
		else if (sf=="-X") 
			message(" [!] found '-' character in states codes, not recommended.") 
	}
	else if (informat %in% c("SPS","SPELL")) {
		seqdata <- seqformat(seqdata,from=informat,to='STS',stsep=stsep, ...)
		## if (is.null(cnames)) cnames <- colnames(seqdata)
	}

	cntmp <- colnames(seqdata)

	if (!is.na(missing)) seqdata[seqdata==missing] <- NA
	message(" [>] missing values in input file coded as: ", missing) 

	## =====================
	## DEFINING THE ALPHABET
	## =====================
	statl <- seqstatl(seqdata)
	if (missing(alphabet)) plevels <- statl
	else plevels <- alphabet

	## ===================
	## PREPARING THE DATA
	## ===================
	seqdata <- as.matrix(seqdata)
	if (any(is.na(seqdata)))
		seqdata <- seqprep(seqdata, left=left, gaps=gaps, right=right, void=void, nr=nr)

	## ======================================
	## DEFINING THE CLASS AND SOME ATTRIBUTES
	## ======================================
	seqdata <- as.data.frame(seqdata)
	class(seqdata) <- c("stslist","data.frame")

	if (is.null(start)) start <- 1
	attr(seqdata,"start") <- start
	attr(seqdata,"missing") <- missing
	attr(seqdata,"void") <- void
	attr(seqdata,"nr") <- nr

	## 
	if (missing(states)) {
		message(" [>] distinct states appearing in the data: ",paste(statl,collapse="/"))
		attr(seqdata,"alphabet") <- plevels
		}
	else {
		## plevels <- states
		attr(seqdata,"alphabet") <- states
	}

	alphabet <- attr(seqdata,"alphabet")
	nbstates <- length(alphabet)

	if (nbstates==1) stop("alphabet has only one state!")
	if (nbstates>12 && missing(cpal)) 
		stop("Can not attribute automatic color palete, number of states too high. Please specify one using 'cpal' option!")

	message(" [>] alphabet: ",paste(1:nbstates,alphabet,collapse=" ",sep="="))

	## ==========================================
	## COLOR PALETTE AND STATES LABELS ATTRIBUTES
	## ==========================================
	if (missing(cpal)) {
		if (nbstates==2) attr(seqdata,"cpal") <- brewer.pal(3,"Accent")[1:2]
		else if (nbstates<=8) attr(seqdata,"cpal") <- brewer.pal(nbstates,"Accent")
		else if (nbstates>8 & nbstates<=12) attr(seqdata,"cpal") <- brewer.pal(nbstates,"Set3")
	}
	else attr(seqdata,"cpal") <- cpal

	attr(seqdata,"missing.color") <- missing.color

	if (!missing(labels)) attr(seqdata,"labels") <- labels
	else attr(seqdata,"labels") <- alphabet
	labels <- attr(seqdata,"labels")

	message(" [>] labels: ",paste(1:nbstates,labels,collapse=" ",sep="="))
	
	## =====
	## Stats
	## =====
	nbseq <- nrow(seqdata)
	seql <- seqlength(seqdata)

	message(" [>] ", nbseq, " sequences in the data set")
	message(" [>] min/max sequence length: ",min(seql),"/",max(seql))

	## ==================================
	## Converting each column to a factor
	## ==================================
	for (i in 1:ncol(seqdata)) {
		seqdata[,i] <- factor(seqdata[,i], 
			levels=c(plevels,nr,void),
			labels=c(attr(seqdata,"alphabet"),nr,void))
		}
	
	## ======================
	## Rows and columns names
	## ======================
	if (!missing(cnames)) colnames(seqdata) <- cnames
	else {
		if (is.null(cntmp)) 
			colnames(seqdata) <- paste("T",start:(max(seql)+start-1),sep="")
		else if (all(is.na(cntmp)!=TRUE)) 
			colnames(seqdata) <- cntmp
		else colnames(seqdata) <- paste("T",start:(max(seql)+start-1),sep="")
	}

	rownames(seqdata) <- paste("[",1:nbseq,"]",sep="")
	
	## ======================
	## Returning the object
	## ======================
	return(seqdata)
}



