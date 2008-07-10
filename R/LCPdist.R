## LCP distance between 2 sequences

LCPdist <- function(seq1,l1, seq2,l2, norm) {

	if (norm==TRUE) 
		dist <- 1-(seqLCP(seq1[1:l1],seq2[1:l2])/sqrt(l1*l2))
	else 
		dist <- l1+l2-2*seqLCP(seq1[1:l1],seq2[1:l2]) 

	return(dist)

	}
	
	

				
				
				
