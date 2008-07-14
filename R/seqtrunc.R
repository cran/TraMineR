## Compression de sequences (elimination des NA)

seqtrunc <- function(seq) 
	{
		idx <- !is.na(seq)
		ns <- sum(idx)
		return(seq[1:ns])
	} 
