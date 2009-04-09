## Compute the entropy of a distribution

entropy <- function(distrib)
	{
		distrib <- distrib[distrib!=0]
		p <- distrib/sum(distrib)
		e <- -sum(p*log(p))
		return(e)
	}
