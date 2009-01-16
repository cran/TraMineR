## SEQUENCE TRUNCATION

TraMineR.trunc <- function(seq, left="DEL", right="DEL", gaps="DEL", 
	neutral="#", missing=NA, void="%") {

	sl <- length(seq)
	
	if (is.na(missing)) {
		na.pos <- which(is.na(seq))
		notna.pos <- which(!is.na(seq))
	}
	else {
		na.pos <- which(seq==missing)
		notna.pos <- which(seq!=missing)	
	}

	if (length(na.pos)==sl)
		return(seq)
	else {
		if (any(na.pos < min(notna.pos)))
			lc <- max(na.pos[na.pos < min(notna.pos)])
		else lc=0
		rc <- max(notna.pos)+1
		mm <- na.pos[na.pos > lc+1 & na.pos < rc-1]

		seq.trunc <- seq
	
		if (!is.na(left) & lc>0) {
			if (left=="DEL") seq.trunc[1:lc] <- void
			else if (left=="NEUTRAL") seq.trunc[1:lc] <- neutral
		}

		if (!is.na(right) & rc<=sl) {
			if (right=="DEL") seq.trunc[rc:sl] <- void
			else if (right=="NEUTRAL") seq.trunc[rc:sl] <- neutral
		}

		if (!is.na(gaps) & length(mm>0)) {
			if (gaps=="DEL") seq.trunc[mm] <- void
			else if (gaps=="NEUTRAL") seq.trunc[mm] <- neutral
		}
	
		ndel <- sum(seq.trunc==void, na.rm=TRUE)
	
		if (ndel>0) {
			seq.trunc <- seq.trunc[seq.trunc!=void]
			seq.trunc <- c(seq.trunc,rep(void,ndel))
			if (!is.null(names(seq)))
				names(seq.trunc) <- names(seq)
			else 
				names(seq.trunc) <- paste("[",1:sl,"]",sep="")
		}

		return(seq.trunc)
	}
} 
