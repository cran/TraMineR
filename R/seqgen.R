
seqgen <- function(n, length, alphabet, p=NULL) {

	nbstat <- length(alphabet)

	## FROM http://www.cs.chalmers.se/~dubhashi/Courses/BioAlg/Markov%20chains%20and%20DNA%20sequence%20analysis%20in%20R.htm

	m <- sample( alphabet, length, replace=TRUE, p )

	for (i in 1:(n-1)) m <- rbind(m, sample( alphabet, length, replace=TRUE, p ))

	colnames(m) <- paste("[",1:length,"]", sep="")
	rownames(m) <- 1:n

	return(m)
}
		
