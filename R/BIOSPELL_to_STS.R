## ================================
## Convert from SPELL to STS format
## ================================

BIOSPELL_to_STS <- function(seqdata, id=1, begin=2, end=3, status=4, 
	process=TRUE, pdata=NULL, pvar=NULL, 
	limit=100, overwrite=TRUE, fillblanks=NULL, 
	tmin=NULL, tmax=NULL) {

	## if overwrite=TRUE, the latest spell overwrite the one before, if set to FALSE, the earlier is kept
	begincolumn <- seqdata[,begin]
	endcolumn <- seqdata[,end]
	frmoption <- NULL

	if (!is.null(pdata) & !is.null(pvar)) pdata <- pdata[, pvar]

#	if(process==TRUE && mean(begincolumn[!is.na(begincolumn) & begincolumn > 0]) > 1900 && is.null(pdata)) {
#		stop("Option mismatch : you want age sequences but you don't have ages in the begin column. You need to set birth years with pdata or specify process=FALSE")
#	}
#	if(process==FALSE && mean(begincolumn[!is.na(begincolumn) & begincolumn > 0]) < 100 && is.null(pdata)) {
#		stop("Option mismatch : you want year sequences, but you don't have years in the begin column. You need to set birth years with pdata or specify process=TRUE")
#	}
#	# option redondante : on a des années, on veut des années, et on donne des années de naissance en +
#	if(process==FALSE && mean(begincolumn[!is.na(begincolumn) & begincolumn > 0]) > 1900 && !is.null(pdata)) {
#		pdata <- NULL
#		warning("You want year sequences and you have years at the beginning and at the end of episodes. Ignoring the pdata option") 
#	}

	if (process==TRUE) {
		if (!is.null(pdata)) frmoption <- "year2age"
		else frmoption <- "age2age"
	}

	if (process==FALSE) {
		if (is.null(pdata)) frmoption <- "year2year"
		else frmoption <- "age2year"
	}

	## =========================
	## creation of the dataframe
	## =========================
	if(process==FALSE) {
		if(is.null(tmin) || is.null(tmax)) {
			if(frmoption=="year2year") {
				tmin <- min(begincolumn[!is.na(begincolumn) & begincolumn > 0])
				tmax <- max(endcolumn[!is.na(endcolumn) & endcolumn > 0])
			}
			if(frmoption=="age2year") {
				tmin <- min(pdata[,2][!is.na(pdata[,2]) & pdata[,2] > 0])
				max1 <- max(endcolumn[!is.na(endcolumn) & endcolumn > 0])
				max2 <- max(pdata[,2][!is.na(pdata[,2]) & pdata[,2] > 0])
				tmax <- max1 + max2
			}
			message(paste(" [>] time axis:", tmin, "->", tmax))
			if(is.na(tmin) || is.na(tmax)) {
				stop("Could not find the minimum or maximum year automatically, please use tmin/tmax options")
			}
		}	
		limit <- (tmax - tmin) + 1
		year <- tmin
		names.seqresult <- NULL
		for(i in 1:limit) {
			names.seqresult <- c(names.seqresult, (paste("y", year, sep="")))
			year <- year+1
		}
		#names.seqresult <- c(names.seqresult, "id")
	}
	else {
		#names.seqresult <- c((paste("a", seq(1:limit), sep="")),"id")
		names.seqresult <- c((paste("a", seq(1:limit), sep="")))
	}
	
	#seqresult <- matrix(nrow=1, ncol=limit+1)
	seqresult <- matrix(nrow=1, ncol=limit)
	seqresult <- as.data.frame(seqresult)
	if (is.factor(seqdata[,status])) { 
		for (k in 1:(limit)) { 
			seqresult[,k] <- factor(seqresult[,k], levels=levels(seqdata[,status]), labels=levels(seqdata[,status])) 
		} 
	}
	names(seqresult) <- names.seqresult
	## ================================
	## end of creation of the dataframe
	## ================================
	
	# on récupère la liste des individus
	lid <- unique(seqdata[,id])
	nbseq <- length(lid)

	# si un dataframe avec les années de naissances a été donné en argument, on récupère les ID et les années
	if(frmoption=="year2age" || frmoption=="age2year") {
		birthyrid1 <- pdata[,1]
		birthyr1 <- pdata[,2]
        }
	######################
	# individual loop
	for (i in 1:nbseq) {
		spell <- seqdata[seqdata[,id]==lid[i],]
		# number of spell for individual i
		idxmax <- nrow(spell)
		# we check if the first episode looks normal (starting age/year > 0) (problèmes avec un fichier du panel)
		if(frmoption=="age2age" || frmoption=="year2year" || frmoption=="age2year") {
			age1 <- spell[1,begin]
		}
		# if we need to convert years to ages, we need the birthyear
		if(frmoption=="year2age") {
			birthy <- birthyr1[birthyrid1==lid[i]]
			age1 <- spell[1,begin] - birthy
		}
		# if we convert from ages to years, we need the birthyear, but don't need to substract it to the time of beginning
		if(frmoption=="age2year") {
			birthy <- birthyr1[birthyrid1==lid[i]]
			#print("birthyr")
			#print(birthy)
			age1 <- spell[1,begin]
		}
		if (is.na(age1)) { age1 <- -1 }
	
		# we fill the line with NAs
		#seqresult[i,1:(limit+1)] <- c(rep(NA,(limit+1)))
		seqresult[i,1:(limit)] <- c(rep(NA,(limit)))
		#seqresult[i,limit+1] <- lid[i]
		
	      	if (age1 >= 0) {
			if (idxmax>0) {
				# by default, the most recent episode erases the one before
				if(overwrite==TRUE) {      
					spelllist <- 1:idxmax
				}
				# if we want the opposite, we just go from the last to the first episode
				else {
					spelllist <- idxmax:1
				}
				# for each spell
				for (j in spelllist) {
					#####################
					# definition of starting and ending point of the spell
					
					# spell are allready in age format, and we want age sequences
					if(frmoption=="age2age") {
						sstart <- spell[j,begin]
						sstop <- spell[j,end]
					}
					if(frmoption=="age2year") {
						sstart <- (birthy-tmin) +  (spell[j,begin]) +1
						sstop <- (birthy-tmin) +  (spell[j,end]) +1
					}
					# spell are in year format, and we want year sequences
					if(frmoption=="year2year") {
						sstart <- (spell[j,begin] - tmin)+1
						sstop <- (spell[j,end] - tmin)+1
						#print(sstart)
						#print(sstop)
					} 
					
					# spell are in year format, and we want age sequences
					if(frmoption=="year2age") {
						sstart <- spell[j,begin] - birthy
						sstop <- spell[j,end] - birthy
					}
					
					#######################
					# fillblanks option
					# if fillblanks is not null, the gaps between episodes is filled with its value
					if (!is.null(fillblanks)) {
						if (j>1) {
							# for every episode after the first one, we check if there is a gap between the one before and this one
							if(frmoption=="age2age" || frmoption=="age2year") {
								previousend <- spell[j-1,end]
							}
							if(frmoption=="year2age") {
								previousend <- spell[j-1,end] - birthy
							}
							if(frmoption=="year2year") {
								previousend <- (spell[j-1, end] - tmin)+1
							}
							if (sstart != previousend && sstart != (previousend+1)) {
								dur <- sstart - (previousend+1)
								if (dur>0 & (sstart-1 < limit) && sstart > 0 && spell[j-1,end] > 0) {
									seqresult[i,(previousend+1):(sstart-1)] <- rep(fillblanks, dur)
								}
							}

						}
					}
					#
					#########################
					
					#########################
					# conversion from episode to subsequence
					dur <- sstop - sstart
					# we check if all values look normal
					if (dur >= 0 && sstop > 0 && sstart >= 0) {
						state <- spell[j,status]
						if (!is.na(state)) {
							# if dur == 0, it means the individual stays in the state only one year
							if (dur == 0 && (sstop < limit) ) {
								seqresult[i,sstart] <- state
                                                        }
							else {
								if(sstop <= limit) {
									# if the sequence begins at age 0, we delete the first state
									if (sstart==0) { 
										sstart <- sstart+1 
										dur <- dur -1
										} 
									seqresult[i,sstart:sstop] <- rep(state, dur+1) 
								}
							}
					        }
					}
				}
			 }
		}
	}	
	return(seqresult)
}

