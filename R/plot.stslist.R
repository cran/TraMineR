## =============================
## PLot of STS sequence objects
## =============================

plot.stslist <- function(x, group=NULL, type="i", title=NULL, tlim=NULL, 
	pbarw=FALSE, sortv=NULL, 
	method="modseq", dist.matrix=NULL,
	cpal=NULL, missing.color=NA, 
	ylab, yaxis=NULL, axes="all", xtlab=NULL, cex.plot=1,
	withlegend="auto", ltext=NULL, cex.legend=1, 
	use.layout=(!is.null(group) | withlegend!=FALSE), legend.prop=NA, rows=NA, cols=NA, ...) {

	if (!inherits(x,"stslist"))
		stop("data is not a sequence object, use seqdef function to create one")

	## Extracting some sequence characteristics
	statl <- attr(x,"alphabet")
	nr <- attr(x,"nr")
	nbstat <- length(statl)

	## ================================
	## Setting color palette and labels
	## ================================
	if (is.null(ltext)) ltext <- attr(x,"labels")

	if (is.null(missing.color)) missing.color <- attr(x,"missing.color") 

	if (is.null(cpal)) cpal <- attr(x,"cpal")

	if (is.null(xtlab)) xtlab <- colnames(x)

	## Adding an entry for missing in the legend
	if (any(x==nr)) {
		cpal <- c(cpal,missing.color)
		ltext <- c(ltext,"missing")
		statl <- c(statl,nr)
		nbstat <- nbstat+1
		}

	## ==============================
	## Preparing if group is not null
	## ==============================
	if (!is.null(group)) {
		## Eliminate the unused levels
		group <- factor(group)
		nplot <- length(levels(group))
		gindex <- vector("list",nplot)
				
		for (s in 1:nplot)
			gindex[[s]] <- which(group==levels(group)[s])

		## Title of each plot
		if (!is.null(title)) 
			title <- paste(title,"-",levels(group))
		else 
			title <- levels(group)
	}
	else {
		nplot <- 1
		gindex <- vector("list",1)
		gindex[[1]] <- 1:nrow(x)
	}

	## ===================
	## Defining the layout
	## ===================
	if (use.layout | !is.null(group) ) {
		## Saving graphical parameters
		savepar <- par(no.readonly = TRUE)

		lout <- TraMineR.setlayout(nplot, rows, cols, withlegend, axes, legend.prop)
	  	layout(lout$laymat, heights=lout$heights, widths=lout$widths)
		axisp <- lout$axisp
		legpos <- lout$legpos
	}
	else {
		if(axes!=FALSE) axisp <- 1
		else axisp <- 0
		legpos <- NULL
	}

	## =======
	## Ploting
	## =======
	for (np in 1:nplot) {

		## Selecting sub sample for x
		## according to 'group'
		subdata <- x[gindex[[np]],]

		## Selecting sub sample for sort variable
		## according to 'group'
		if (!is.null(sortv))
			subsort <- sortv[gindex[[np]]]
		else 
			subsort <- NULL
	
		## State distribution plot
		if (type=="d")
			TraMineR.dplot(subdata, np, title[np], cpal, ylab, yaxis, axisp, xtlab, cex.plot, ...)
		## Sequence frequency plot
		else if (type=="f")
			TraMineR.fplot(subdata, np, title=title[np], tlim=tlim, cpal, pbarw, ylab, yaxis, axisp, xtlab, 					cex.plot, ...)
		## Sequence index plot
		else if (type=="i")
			TraMineR.iplot(subdata, np, title=title[np], tlim=tlim, sortv=subsort, 
				cpal, ylab, yaxis, axisp, xtlab, cex.plot, ...)
		## Mean times
		else if (type=="mt")
			TraMineR.mtplot(subdata, np, title=title[np], cpal, ylab, axisp, xtlab, cex.plot, ...)
		## Representative sequence
		else if (type=="r") {
			## Selecting distances according to group
			if (!is.null(dist.matrix)) {
				subdist <- dist.matrix[gindex[[np]],gindex[[np]]]
				dmax <- max(dist.matrix)
			}
			else {
				subdist <- NULL
				dmax=NULL
			}

			TraMineR.rplot(subdata, np, title=title[np], method=method, dist.matrix=subdist, 
				pbarw=pbarw, dmax=dmax, cpal, ylab, axisp, xtlab, cex.plot, ...)
		}
	}	

	## Plotting the legend
	if (!is.null(legpos)) 
		TraMineR.legend(legpos, ltext, cpal, cex=cex.legend)

	## Restoring graphical parameters
	if (use.layout | !is.null(group) ) 
		par(savepar)
}
