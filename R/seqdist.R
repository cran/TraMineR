## ====================================================
## Computing distances between sequences
## Available metrics (method):
## OM = optimal matching
## LCP = Longest Common Prefix (Elzinga)
## LCS = Longest Common Subsequence (Elzinga)
## ====================================================

seqdist <- function(seqdata, method, refseq=NULL, norm=FALSE, 
	indel=1, sm=NA,	with.missing=FALSE, full.matrix=TRUE) {
	debut <- Sys.time()
	## Checking correct arguments
	if (!inherits(seqdata,"stslist")) {
		stop(" [!] data is not a state sequence object, use 'seqdef' function to create one", call.=FALSE)
	}
	metlist <- c("OM","LCP", "LCS", "RLCP", "DHD", "HAM")
	if (missing(method)) {
		stop(" [!] You should specify a method to compute the distances. It must be one of: ", paste(metlist,collapse=" "))
	}
	if (!method %in% metlist) {
		stop(" [!] Method must be one of: ", paste(metlist,collapse=" "), call.=FALSE)
	}
	## Taking care of correct normalization settings
	if (is.logical(norm)) {
		if (norm) {
		## Normalize using Elzinga for LCP, LCS, RLCP
			if (method %in% c("LCP", "LCS", "RLCP")) {
				norm <- 2
			} else {## Normalize using Abbott for OM, HAM, and DHD
				norm <- 1
			}
      	} else {
      		norm <- 0
		}
	} else if (is.character(norm)) {
		## Using normalization name
		## Cast to integer for c code and normdist function
		## Match return the position, removing 1 to start at zero
		normIndice <- match(norm, c("none", "maxlength", "gmean", "maxdist", "YujianBo")) -1

		if (is.na(normIndice)) {  ##Not found
			stop(" [!] Unknown distance normalization method ", norm)
		}
		norm <- normIndice
	} else {
		stop(" [!] Unknown distance normalization method ", norm)
	}
	## Checking missing values
	if (!with.missing && any(seqdata==attr(seqdata,"nr"))) {
		stop("found missing values in sequences, please set 'with.missing=TRUE' to nevertheless compute distances")
	}
	if (method == "OM" && is.character(sm)) {
		if (sm == "TRATE" || sm =="CONSTANT") {
			sm <- seqsubm(seqdata, method=sm, cval=2, with.missing=with.missing, miss.cost=2)
		} else {
			stop(" [!] Unknown method ", sm, " to compute substitution costs")
		}
	}
	## Checking methods that are treated the same
	methodname <- method
	if (method == "LCS") {
		method <- "OM"
		sm <- suppressMessages(seqsubm(seqdata, method="CONSTANT", cval=2, with.missing=with.missing, miss.cost=2))
		indel <- 1
	} else if (method == "HAM") {
		method <- "DHD"
		sm <- seqsubm(seqdata, "CONSTANT", cval=1, with.missing=with.missing,
			miss.cost=1, time.varying=TRUE)
	}

	## =====================
	## Base information
	## =====================
	
	n <- nrow(seqdata)
	alphabet <- attr(seqdata,"alphabet")
	alphsize <- length(alphabet)
	message(" [>] ",n," sequences with ", alphsize, " distinct events/states")
	## Gaps in sequences
	if (with.missing) {
		alphabet <- c(alphabet,attr(seqdata,"nr"))
		alphsize <- length(alphabet)
		message(" [>] including missing value as additional state" )
	}

	## ===========================
	## Checking correct size of sm
	## ===========================
	
	## Checking if substitution cost matrix contains values for each state
	## and if the triangle inequality is respected
	if (method=="OM") {
		if (nrow(sm)!=alphsize | ncol(sm)!=alphsize) {
			stop(" [!] size of substitution cost matrix must be ", alphsize,"x", alphsize)
		}
		if (any(sm<0)) {
			stop(" [!] Negative substitution costs are not allowed")
		}
		if (any(diag(sm)!=0)) {
			stop(" [!] All element on the diagonal of sm (substitution cost) should be equal to zero")
		}
		if (indel <= 0) {
			stop(" [!] indel cost should be positive")
		}
		triangleineq <- checktriangleineq(sm, warn=FALSE, indices=TRUE)
		## triangleineq contain a vector of problematic indices.
		if (!is.logical(triangleineq)) {
			warning("The substitution cost matrix doesn't respect the triangle inequality.\n",
        			" At least, substitution cost between indices ",triangleineq[1]," and ",triangleineq[2],
        			" does not respect the triangle inequality. It costs less to first transform ",
        			triangleineq[1], " into ",triangleineq[3])
		}
		if (any(sm>2*indel)) {
			warning("Some substitution cost are greater that two times the indel cost.",
				" Such substitution cost will thus never be used.")
		}
	}
	## Checking if substitution cost matrix contains values for each state
	## and if the triangle inequality is respected
	if(method== "DHD"){
		## User entered substitution cost
		if(is.array(sm)){
			## checking correct dimension
			smdim <- dim(sm)
			if(!is.array(sm) || sum(smdim==c(alphsize, alphsize, ncol(seqdata)))!=3){
				stop(" [!] size of substitution cost matrix must be ", alphsize,"x", alphsize, "x", ncol(seqdata))
			}
		}
		else {
			sm <- seqsubm(seqdata, "TRATE", cval=1, with.missing=with.missing,
					miss.cost=1, time.varying=TRUE)
		}
	}


	## ==============
	## Preparing data
	## ==============
	seqdata <- seqnum(seqdata, with.missing=with.missing)

	## Selecting distinct sequences only and saving the indexes
	dseq <- unique(seqdata)
	mcorr <- match(seqconc(seqdata), seqconc(dseq))

	nd <- nrow(dseq)
	message(" [>] ", nd," distinct sequences")

	slength <- seqlength(dseq)

	dseq <- seqasnum(dseq, with.missing=with.missing)

	message(" [>] min/max sequence length: ",min(slength),"/",max(slength))

	## ===================================
	## Preparing data for Hamming distance
	## ===================================
	if (method=="DHD") {
		if(length(unique(slength))>1) {
			## Hamming is not defined for sequence of different length
			stop(methodname, " distance can only be computed between sequences of equal length")
		}
		## Here we use the indel parameter for the C function
		## But it contain the maximum possible cost of the hamming distance
		indel <- 0
		for (i in 1:max(slength)) {
			indel <- indel + max(sm[,,i])
		}
	}
	message(" [>] computing distances using ", methodname,
		ifelse(norm!=0," normalized", ""), " metric")
	
	## Function and arguments
	if (!missing(refseq) && !is.null(refseq)) {
		distances <- TraMineR.seqdist.refseq(seqdata, method, refseq,
			norm, indel, sm, alphsize, nd, dseq, slength, mcorr)	
	}
	else {
		distances <- TraMineR.seqdist.all(seqdata, method, 
			norm, indel, sm, alphsize, nd, dseq, slength, mcorr)
	}

	fin <- Sys.time()
	message(" [>] Total time: ", format(round(fin-debut, 3)))

	if (full.matrix && inherits(distances, "dist")) {
		return(dist2matrix(distances))
	}
	else {
		return(distances)
	}
}
