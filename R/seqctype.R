## =================================
## Centrotype of a set of sequences
## =================================

seqctype <- function(seqdata, dist.matrix=NA, method="modeseq") {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one")
	
	slength <- dim(seqdata)[2]
	nbseq <- dim(seqdata)[1]

	if (method=="modeseq") {
		freq <- seqstatd(seqdata)$Frequencies

		statelist <- alphabet(seqdata)	

		seqp <- seqdata[1,]
		pfreq <- freq

		## Constructing the transversal modal sequence
		for (i in 1:slength) {
			smax <- which(freq[,i]==max(freq[,i]))[1]
			pfreq[-smax,i] <- 0
			seqp[,i] <- statelist[smax]
		}

		profile <- list(seqp,pfreq)
		names(profile) <- c("Profile","Frequencies")
	}

	else if (method=="mscore") {
		freq <- seqstatd(seqdata)$Frequencies

		statelist <- alphabet(seqdata)

		if (!length(dist.matrix)>1)
			dist.matrix <- seqdist(seqdata, method="LCS")

		mscore <- vector(mode="integer", length=nbseq)

		## Computing the scores
		for (i in 1:nbseq)
			for (j in 1:slength) {
				score <- freq[seqdata[i,j]==statelist, j]
				mscore[i] <- mscore[i]+score
			}

		ctype <- which(mscore==max(mscore))

		nds <- nrow(unique(seqdata[ctype,]))

		message(" [>] ", nds, " distinct sequence(s) having min. sum of distances")
		if (nds>1) 
			message(" [>] choosing the first sequence")
		
		ctype <- ctype[1]

		seqp <- seqdata[ctype,]

		pfreq <- freq

		for (i in 1:slength) {
			fstate <- which(seqp[,i]==alphabet(seqdata))
			pfreq[-fstate,i] <- 0
		}

		dist.matrix <- dist.matrix[, ctype]

		index <- seqfind(seqp, seqdata)

		profile <- list(seqp, index, mscore, pfreq, dist.matrix)
		names(profile) <- c("Profile", "Index", "Scores", "Frequencies", "Distances")
	}


	else if (method=="dist") {
		
		if (!length(dist.matrix)>1)
			dist.matrix <- seqdist(seqdata, method="LCS")

		sumdist <- apply(dist.matrix, 1, sum)
		ctype <- which(sumdist==min(sumdist))
		
		nds <- nrow(unique(seqdata[ctype,]))

		message(" [>] ", nds, " distinct sequence(s) having min. sum of distances")
		if (nds>1) 
			message(" [>] choosing the first sequence")
		
		ctype <- ctype[1]

		dist.matrix <- dist.matrix[, ctype]
		
		seqp <- seqdata[ctype,]

		index <- seqfind(seqp, seqdata)

		freq <- seqstatd(seqdata)$Frequencies

		pfreq <- freq

		for (i in 1:slength) {
			fstate <- which(seqp[,i]==alphabet(seqdata))
			pfreq[-fstate,i] <- 0
		}

		profile <- list(seqp,index,pfreq,dist.matrix)
		names(profile) <- c("Profile","Index","Frequencies","Distances")
	}	
	

	return(profile)
}	
	
