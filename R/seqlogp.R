seqlogp <- function(seqdata, transition="age") {
	## Liste des taux de transitions par age
	sl <- seqlength(seqdata)
	maxage <- max(sl)
	agedtr <- vector(mode="list", length=maxage)

	## On ajoute 1 pour que les codes correspondent aux index R (commence à 1)
	seqdatanum <- TraMineR:::seqasnum(seqdata)+1

	## Frequence du premier état
	firstfreq <- seqstatd(seqdata, digits=NULL)$Frequencies[, 1]

	if (transition=="global") {
		tr <- seqtrate(seqdata)
		agedtr[[1]] <- firstfreq
		for (i in 2:maxage) {
			agedtr[[i]] <- tr
		}
	}
	else if (transition=="age") {
		agedtr[[1]] <- firstfreq
		maxlevel <- max(seqdatanum, na.rm=TRUE)
		for (i in 2:maxage) {
			tr <- prop.table(table(seqdatanum[, i-1], seqdatanum[, i]), 1)
			## Attention table supprime les états qui n'apparaissent pas
			mytr <- matrix(0, nrow=maxlevel, ncol=maxlevel)
			mytr[as.integer(dimnames(tr)[[1]]), as.integer(dimnames(tr)[[2]])] <- tr
			agedtr[[i]] <- mytr
		}
	}
	else if (transition=="static") {
		## On crée quand même une matrice de transition (qui ne dépend pas de l'état précédant)
		## On peut ainsi utiliser le même algorithme
		agedtr[[1]] <- firstfreq
		freqs <- seqstatd(seqdata)$Frequencies
		for (i in 2:maxage) {
			maxlevel <- max(seqdatanum)
			mytr <- matrix(0, nrow=maxlevel, ncol=maxlevel)


			for (j in 1:length(freqs[, i])) {
					mytr[, j] <- freqs[j, i]
			}
			agedtr[[i]] <- mytr
		}
	}
	else {
		stop("Unknow method to compute transition rate")
	}
	loglik <- numeric(length=(nrow(seqdata)))
	loglik[] <- 0
	for (i in 1:nrow(seqdatanum)) {
		p <- agedtr[[1]][seqdatanum[i, 1]]
		loglik[[i]] <- -log(p)
		if (sl[i]>1) {
			for (j in 2:sl[i]) {
				tr <- agedtr[[j]]
				p <- tr[seqdatanum[i, j-1], seqdatanum[i, j]]
				loglik[[i]] <- loglik[[i]] -log(p)
			}
		}
	}
	return(loglik)
}
