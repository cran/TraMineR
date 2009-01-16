## Levenshtein distance between seq1 and seq2

levenshtein <- function(seq1, l1, seq2, l2, indel, sm, alphsize, norm) {
	
	dist <- .C("cLEVEN", as.integer(seq1), as.integer(seq2), 
		as.double(c(l1,l2,indel,alphsize)), as.double(sm), result = as.double(0), PACKAGE="TraMineR")$result
	if (norm) 
		return(dist/max(l1, l2))
	return(dist)
	}
