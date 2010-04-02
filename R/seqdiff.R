###########################
## Locate the difference in sequence between groups
###########################
seqdiff <- function(seqdata, group, cmprange=c(0, 1),
     seqdist_arg=list(method="LCS", norm=TRUE)) {

	if (!inherits(seqdata, "stslist")) {
			stop("seqdata should be a stslist, see seqdef")
	}
	slenE <- ncol(seqdata)
	ret <- NULL
	startAt <- 1
	#Range where we compare
	totrange <- max(startAt, 1-cmprange[1]):min(slenE, slenE-cmprange[2])
	for (i in totrange) {
		gc()
		srange=c((i+cmprange[1]):(i+cmprange[2]))
		if (inherits(group, "stslist")) {
			cmpbase <- group[, i]
		} else {
			cmpbase <- group
		}
		#Getting complete case on range and group var
		subseq <- seqdata[, srange]
		seqok <- complete.cases(cmpbase, subseq)
		#Computing distance on range
		seqdist_arg$seqdata <- subseq[seqok, ]
		sdist <- suppressMessages(do.call(seqdist, args=seqdist_arg))
		tmp <- dissassoc(sdist, cmpbase[seqok], R=1)
		if (is.null(ret)) {
			ret <- list()
			ret$stat <- tmp$stat[c(1, 2, 4)]
			ret$variance <- tmp$groups$variance
		}
		else {
			ret$stat <- rbind(ret$stat, tmp$stat[c(1, 2, 4)])
			ret$variance <- rbind(ret$variance, tmp$groups$variance)
			colnames(ret$variance) <- rownames(tmp$groups)
		}
	}
	rownames(ret$stat) <- colnames(seqdata)[totrange]
	rownames(ret$variance) <- colnames(seqdata)[totrange]
	class(ret) <- "seqdiff"
	return(ret)
}

###########################
## Print method for seqdiff
###########################
print.seqdiff <- function(x, ...) {
	message("\nStatistics:")
	print(x$stat, ...)
	message("\nVariances:")
	print(x$variance, ...)
}

###########################
## Plot method for seqdiff
###########################
plot.seqdiff <- function(x, stat="PseudoR2", type="l", ylab=stat, xlab="",
	legendposition="top", ylim=NULL, ...) {

	if (stat=="Variance" || stat=="Residuals") {
		nbstates=ncol(x$variance)
		if (nbstates <= 8) cpal <- brewer.pal(nbstates, "Accent")
		 else if (nbstates > 8 & nbstates <= 12) cpal <- brewer.pal(nbstates, "Set3")

		if (stat=="Residuals") {
			toplot <- x$variance*(1-x$stat$PseudoR2)
		} else {
			toplot <- x$variance
		}
		if (is.null(ylim)) {
			ylim=c(min(toplot), max(toplot))
		}
		plot(1:nrow(x$variance), x$variance[, ncol(x$variance)], type=type, ylab=ylab, xlab=xlab, xaxt="n", col=cpal[nbstates], ylim=ylim, ...)

		for (i in 1:(ncol(x$variance)-1)) {
			lines(toplot[, i], type=type, col=cpal[i], ...)
		}
		legend(legendposition, fill = cpal, legend = colnames(x$variance))
		axis(1, at=1:nrow(x$variance), labels=rownames(x$stat) )
	}
	else if (stat %in% c("PseudoR2", "PseudoT", "PseudoF")) {
		plot(x$stat[, stat], type=type, ylab=ylab, xlab=xlab, xaxt="n", ...)
		axis(1, at=1:nrow(x$variance), labels=rownames(x$stat) )
	}
	else {
		stop("Unknow value for the 'stat' argument")
	}
}
