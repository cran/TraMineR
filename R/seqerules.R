
# j => i
seqerules <- function(res, minsup, sortv=NULL, decreasing=FALSE, wcount=FALSE) {
	
	
	results <- list()
	matevent <- seqeapplysub(res, rules=TRUE)
	matcount <- seqeapplysub(res, method="presence")
	k<-1
	for(i in 1:nrow(matevent)) {
		
		for(j in 1:ncol(matevent)) {
			if(matevent[i,j]!=0 && as.character(res$subseq[i,])!=as.character(res$subseq[j,])) {
				#	print(i)
				#	print(j)
				#	if (unlist(strsplit(as.character(res$subseq[j,]),"-")) != unlist(strsplit(as.character(res$subseq[i,]),"-"))[-1]) {
				prelist <- unlist(strsplit(as.character(res$subseq[j,]),"-"))
				conlist <- unlist(strsplit(as.character(res$subseq[i,]),"-"))
				if((conlist[length(conlist)] != prelist[length(prelist)]) && (prelist[1]==conlist[1])) {
					
					pre <- as.character(res$subseq[j,])
					conlist2 <- conlist[(length(prelist)+1):length(conlist)]
					# Pour éviter le cas où pre = A-B, con=A-Z-B-C et rule = A-B=>B-C
					if(prelist[length(prelist)] != conlist2[1]) {
						
						sizeofr <- length(prelist) + length(conlist2)
						results$size <- c(results$size, sizeofr)
						con <- paste(conlist2, collapse="-")
						
						precon <- res$subseq[i,]
						
						#conC <- 	unlist(strsplit(as.character(res$subseq[i,]), "-"))[-1]
						#conC <- res$subseq[i,]
						results$rules<-c(results$rules,paste(pre, "=>", con))
						
						#results$counta <- c(results$counta,res$data[j,2])
						#results$countb <- c(results$countb,res$data[as.character(res$subseq)==conC,2])
						#results$countab <- c(results$countab,res$data[i,2])
						ntotal <- length(res$seqe)
						results$a[k] <- pre
						
						c_a <- sum(matcount[,pre])
						c_b <- sum(matcount[,con])
						c_ab <- sum(matcount[,as.character(res$subseq[i,])])
						p_a <- c_a/ntotal
						p_b <- c_b/ntotal
						p_ab <- c_ab/ntotal
						
						results$counta[k] <- c_a
						results$b[k] <- con
						results$countb[k] <- c_b
						results$ab[k] <- as.character(res$subseq[i,])
						#results$countab <- c(results$countab, c_ab)
						results$countab[k] <- c_ab
						results$dbi[k] <- i
						results$dbj[k] <- j
						confidence <- p_ab/p_a
						lift <- confidence/p_b
#print("ok4.5")
						
						
						results$conf[k] <- confidence
						results$lift[k] <- lift
						
						# calcul du lift standardisé
						#lambda = max((p_a+p_b-1)/(p_a*p_b), (4*minsup)/((1+minsup)^2), minsup/(p_a*p_b), confidence/p_b)
						lambdahaut <- max((p_a+p_b-1), (1/ntotal))
						lambda = lambdahaut/(p_a*p_b)
						upsilon = 1/(max(p_a,p_b))
#						print(lambda)
#						print(upsilon)
						standardlift = (lift - lambda) / (upsilon - lambda)
						results$standardlift[k] <- standardlift
						
						
						icstat <- implicativestat(matcount[,pre], matcount[,con], type="indice")[2,2]
						results$icstat[k] <- icstat
						results$icstatp[k] <- pnorm(-icstat) 
						
						#print("ok5")
						k <- k+1
						
					}
				}
			} 
		}
	}
	#names(results$icstat)<-"Implicate" 
	#resdata <- data.frame("Rules" = results$rules, "Count prequel" = results$counta, "Count sequel" = results$countb, "Count rules" = results$countab, "Rule size" = results$size, "Debug i" = results$dbi, "Debug j" = results$dbj)
	if (wcount) {
		resdata <- data.frame("Rules" = results$rules, "Count prequel" = results$counta, "Count sequel" = results$countb, "Count rules" = results$countab[], "Rule size" = results$size, "Conf" = results$conf, "Lift" = results$lift, "Standardlift" = results$standardlift, "ImplicStat" = results$icstat, "p-value" = results$icstatp)
	}
	else { 
		resdata <- data.frame("Rules" = results$rules, "Rule size" = results$size, "Conf" = results$conf, "Lift" = results$lift, "Standardlift" = results$standardlift, "ImplicStat" = results$icstat, "p-value" = results$icstatp)
	}
	
	if(!is.null(sortv)) {
		or <- order(results[[sortv]],decreasing=decreasing)
		return(resdata[or,])
	}
	
	return(resdata)
	
}


