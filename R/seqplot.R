## ====================================================
## Generic function for plotting state sequence objects
## ====================================================

seqplot <- function(seqdata, group = NULL, type = "i", main = NULL, cpal = NULL,
  missing.color = NULL, ylab = NULL, yaxis = TRUE, axes = "all", xtlab = NULL,
  cex.axis = 1, with.legend = "auto", ltext = NULL, cex.legend = 1,
  use.layout = (!is.null(group) | with.legend != FALSE), legend.prop = NA,
  rows = NA, cols = NA, title, cex.plot, withlegend, ...) {

  TraMineR.check.depr.args(alist(main = title, cex.axis = cex.plot, with.legend = withlegend))

	if (!inherits(seqdata,"stslist"))
		stop(call.=FALSE, "seqplot: data is not a state sequence object, use seqdef function to create one")

	## Storing original optional arguments list
	oolist <- list(...)

  	if ("sortv" %in% names(oolist)) {sortv <- oolist[["sortv"]]}
  	leg.ncol <- if ("ncol" %in% names(oolist)) { oolist[["ncol"]] } else { NULL }
    oolist <- oolist[names(oolist) != "ncol"]
    if ("tlim" %in% names(oolist)) {
      oolist[["idxs"]] <- oolist[["tlim"]]
      msg.warn("'tlim' deprecated, use 'idxs' instead!")
      oolist <- oolist[names(oolist) != "tlim"]
    }
    barlab <- if ("bar.labels" %in% names(oolist)) { as.matrix(oolist[["bar.labels"]]) } else { NULL }

    diss <- NULL
  	if ("diss" %in% names(oolist)) {
      diss <- oolist[["diss"]]
    }
  	else if ("dist.matrix" %in% names(oolist)) {
      diss <- oolist[["dist.matrix"]]
      oolist[["diss"]] <- diss
    } # FIXME dist.matrix is deprecated


  if (type == "pc") { # modification of Reto Bürgin 16.08.2012
    oolist <- append(oolist, list(group = group, rows = rows, cols = cols))
    group <- NULL
  }

  if (type == "r") { # stuff moved here by GR 17.01.2018
    ## For type="r" each group should have at least 2 cases
    grp <- group
    if (is.null(grp)) grp <- rep(1,nrow(seqdata))
    if (any(xtabs( ~ group(grp)) < 2))
      stop("For type = 'r', each group must have 2 or more cases. At least one group has only 1.", call.=FALSE)

		if (is.null(diss)) {## (! "diss" %in% names(oolist)  && ! "dist.matrix" %in% names(oolist))){
      if (! "method" %in% names(oolist)){
			  stop("For type = 'r', you must provide a distance matrix or a method to compute it", call.=FALSE)
      } else {
        #msg("seqplot calls seqdist")
        oolist[["seqdata"]] <- seqdata
        dlist <- names(formals(seqdist))
        diss <- do.call(seqdist, args = oolist[names(oolist) %in% dlist])
        oolist[["diss"]] <- diss
        # removing seqdist args
        oolist <- oolist[!(names(oolist) %in% dlist[dlist != "weighted"])]
      }
    }
		#else { ## GR: should also be checked on the seqdist outcome
			if (inherits(diss, "dist")) {
      				diss <- dist2matrix(diss)
			}
    #}
    ## Setting unique Max theoretical distance for all groups
		if (!"dmax" %in% names(oolist)) {
			dmax <- max(diss)
			oolist <- c(oolist,list(dmax=dmax))
		}
  }



	## ==============================
	## Preparing if group is not null
	## ==============================

	if (!is.null(group)) {
          group <- group(group)

          ## Check length
          if (length(group)!=nrow(seqdata))
            stop(call.=FALSE, "group must contain one value for each row in the sequence object")

          nplot <- length(levels(group))
          gindex <- vector("list",nplot)

          if (type=="mt" & !is.null(barlab)){
            if (!(ncol(barlab) %in% c(1,nplot)) )
            stop(call.=FALSE, "When a matrix, bar.labels should have one column per group")
          }

          for (s in 1:nplot)
            gindex[[s]] <- which(group==levels(group)[s])

          ## Title of each plot
          if (!is.null(main))
            main <- paste(main,"-",levels(group))
          else
            main <- levels(group)
	} else { # single group
          nplot <- 1
          gindex <- vector("list",1)
          gindex[[1]] <- 1:nrow(seqdata)
	}

	## ===================
	## Defining the layout
	## ===================
	if (type=="Ht" | type =="pc") { with.legend=FALSE }

	## IF xaxis argument is provided
	## it interferes with axes argument
	if ("xaxis" %in% names(oolist)) {
		tmpxaxis <- oolist[["xaxis"]]
		if (tmpxaxis==TRUE) {axes="all"}
		else if (tmpxaxis==FALSE) {axes=FALSE}
		oolist <- oolist[!names(oolist) %in% "xaxis"]
	}

	if (use.layout | !is.null(group) ) {
		## Saving graphical parameters
		savepar <- par(no.readonly = TRUE)

		lout <- TraMineR.setlayout(nplot, rows, cols, with.legend, axes, legend.prop)
	  	layout(lout$laymat, heights=lout$heights, widths=lout$widths)

		## Should axis be plotted or not ?
		xaxis <- 1:nplot==lout$axisp

		legpos <- lout$legpos
	}
	else {
		if (axes!=FALSE) {xaxis <- TRUE}
		else {xaxis <- FALSE}
		legpos <- NULL
	}

	## ========
	## Plotting
	## ========
	for (np in 1:nplot) {
		## Storing ... arguments in a list
		olist <- oolist

		plist <- list(main=main[np], cpal=cpal, missing.color=missing.color,
			ylab=ylab, yaxis=yaxis, xaxis=xaxis[np],
			xtlab=xtlab, cex.axis=cex.axis)

		## Selecting sub sample for x
		## according to 'group'
		subdata <- seqdata[gindex[[np]],]

		## State distribution plot or Entropy index
		if (type=="d" || type=="Ht") {
			f <- seqstatd
			plist <- c(list(type=type), plist)

			plist <- plist[!names(plist) %in% "missing.color"]

			## Removing the 'cpal' argument which is not used
			## in Entropy index plots
			if (type=="Ht") {plist <- plist[!names(plist) %in% "cpal"]}
		}
		## Sequence frequency plot
		else if (type=="f") {
			with.missing <- TRUE
			f <- seqtab
		}
		## Sequence index plot
		else if (type=="i" || type=="I") {
			f <- function(seqdata) {return(seqdata)}
			with.missing <- TRUE

			## Selecting sub sample for sort variable
			## according to 'group'
			if ("sortv" %in% names(olist)) {
				if (!length(sortv)==1) {
					olist[["sortv"]] <- sortv[gindex[[np]]]
				}
			}

			if (type=="I") {
				if (!"idxs" %in% names(olist)) {olist <- c(olist, list(idxs=0))}
				if (!"space" %in% names(olist)) {olist <- c(olist, list(space=0))}
				if (!"border" %in% names(olist)) {olist <- c(olist, list(border=NA))}
			}
		}
		## Mean times
		else if (type=="mt") {
      f <- seqmeant
      if (!is.null(barlab)) {
        if (ncol(barlab)==1)
          olist[["bar.labels"]] <- as.vector(barlab)
        else
          olist[["bar.labels"]] <- as.vector(barlab[,np])
      }
    }
		## Modal states
		else if (type=="ms") {
			f <- seqmodst
		}
		## Representative sequence
		else if (type=="r") {
			f <- seqrep
			with.missing <- TRUE

			## Removing unused arguments
			plist <- plist[!names(plist) %in% "yaxis"]

			## Selecting distances according to group
			# fixed for deprecated dist.matrix

			olist[["diss"]] <- diss[gindex[[np]],gindex[[np]]]

    # FIXME dist.matrix is deprecated
    olist <- olist[names(olist) != "dist.matrix"]

    }
    else if (type == "pc") { # modification of Reto Bürgin 16.08.2012

                        plist[["main"]] <- main
                        olist <- c(olist, plist)
                        olist[["plot"]] <- FALSE
                        f <- seqpcplot
                        olist <- olist[names(olist) %in% names(formals(f))]
                        plist <- list()
                      }
		else
			stop("Unknown 'type' argument.")

		## Calling appropriate function and plotting
		flist <- names(formals(f))

		if ("with.missing" %in% names(olist)) {
			with.missing <- olist[["with.missing"]]
		} else if ("with.missing" %in% flist) {
			with.missing <- formals(f)$with.missing
		}

		## Xlim when plotting individual sequences
		if (type %in% c("i", "I", "f")) {
			if (!"xlim" %in% names(olist)) {
				olist <- c(olist, list(xlim=c(0, ncol(seqdata))))
			}
		}

		match.args <- names(olist) %in% flist
		fargs <- olist[match.args]
		fargs <- c(list(seqdata=subdata), fargs)
    #msg(paste("do.call(",f, fargs,")"))
		res <- do.call(f, args=fargs)

		olist <- olist[!match.args]
    ## suppress non plot arguments if necessary
    olist <- olist[!names(olist) %in% c("with.missing")]
    if (!(type %in% c("i","I"))) olist <- olist[!(names(olist) %in% c("sortv","weighted"))]
    if (type != "r") olist <- olist[!(names(olist) %in% c("dmax","stats"))]

		plist <- c(list(x=res), plist, olist)
		do.call(plot, args=plist)
	}

	## Plotting the legend
	if (!is.null(legpos)) {
		## Extracting some sequence characteristics
		nr <- attr(seqdata,"nr")

		if (is.null(ltext)) ltext <- attr(seqdata,"labels")

		if (is.null(missing.color)) missing.color <- attr(seqdata,"missing.color")

		if (is.null(cpal)) cpal <- attr(seqdata,"cpal")

		density <- if ("density" %in% names(oolist)) { oolist[["density"]] } else { NULL }
		angle <- if ("angle" %in% names(oolist)) { oolist[["angle"]] } else { NULL }

		## Adding an entry for missing in the legend
		if (with.missing & any(seqdata==nr)) {
			cpal <- c(cpal,missing.color)
			ltext <- c(ltext,"missing")
		## statl <- c(statl,nr)
		## nbstat <- nbstat+1
		}

		TraMineR.legend(legpos, ltext, cpal, cex=cex.legend, density=density, angle=angle, leg.ncol=leg.ncol)
	}

	## Restoring graphical parameters
	if (use.layout | !is.null(group)) {par(savepar)}
}
