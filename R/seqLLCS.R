## ========================================
## Length of the longest common subsequence
## ========================================

seqLLCS <- function(seq1,seq2) {

	if (!inherits(seq1,"stslist") | !inherits(seq2,"stslist")) 
		stop("sequences must be sequence objects")

	l1 <- seqlength(seq1)
	l2 <- seqlength(seq2)

	result <- .C("cLCS", as.integer(seq1), as.integer(seq2), as.double(c(l1,l2)), result = as.integer(0), PACKAGE="TraMineR")$result
	
	return(result)
}
