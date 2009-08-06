###########################
## Internal plot method
###########################
DTNseqplot <- function(ind, seqdata, sortv, plottype, seqplot.title.cex,
												seqplot.use.title=TRUE, seqplot.label.loc="main",
												seqplot.node.loc="main", seqplot.split.loc="sub", ...) {
	if (seqplot.use.title) {
		top <- (as.integer(seqplot.label.loc=="main")+as.integer(seqplot.node.loc=="main")+as.integer(seqplot.split.loc=="main"))*seqplot.title.cex
		bottom <- (as.integer(seqplot.label.loc=="sub")+as.integer(seqplot.node.loc=="sub")+as.integer(seqplot.split.loc=="sub"))*seqplot.title.cex
		left <- (as.integer(seqplot.label.loc=="ylab")+as.integer(seqplot.node.loc=="ylab")+as.integer(seqplot.split.loc=="ylab"))*seqplot.title.cex
		par(mar=c(bottom, left, top, 0), font.sub=2, mgp=c(0, 0, 0))
	}
	if (!is.null(sortv)){
		seqplot(seqdata[ind, ], sortv=sortv[ind], ...)
	}
	else {
		seqplot(seqdata[ind, ], ...)
	}
}

###########################
## Shortcut to build sequences tree
###########################
seqtree2dot <- function(tree, filename, seqdata, imgLeafOnly=FALSE, sortv=NULL, ...) {
	disstree2dot(tree, filename, imagedata=NULL, seqdata=seqdata, title.cex=3,
			sortv=sortv, imagefunc=DTNseqplot, use.title=TRUE, label.loc="main",
			node.loc="main", split.loc="sub", seqplot.use.title=TRUE,
			seqplot.label.loc="main", seqplot.node.loc="main", seqplot.split.loc="sub",
			seqplot.title.cex=3, ...)
}


disstree2dotp <- function(tree, filename, imagedata=NULL, imgLeafOnly=FALSE,
													imagefunc=plot, ...){
	disstree2dot(tree, filename, imagedata=NULL, title.cex=3, imagefunc=imagefunc,
			use.title=TRUE, label.loc="main", node.loc="main", split.loc="sub",
			seqplot.use.title=TRUE, seqplot.label.loc="main", seqplot.node.loc="main",
			seqplot.split.loc="sub", seqplot.title.cex=3, ...)

}


###########################
## Generate dot file
###########################
disstree2dot <- function(tree, filename, digits=3, imagefunc=NULL, imagedata=NULL,
												imgLeafOnly=FALSE, devicefunc="jpeg", imageext="jpg",
												device.arg=list(), use.title=TRUE, label.loc="main",
												node.loc="main", split.loc="sub", title.cex=1, ...) {
	dotfile <- paste(filename, ".dot", sep="")
	node <- tree$root
	cat("digraph distree{\n", file=dotfile)
	if (!is.null(imagefunc) && is.null(imagedata)) {
		imagedata <- as.data.frame(node$ind)
	}
	DTN2DotInternal(preced=filename, node=node, pos="none", digits=digits,
			imagefunc=imagefunc, imagedata=imagedata, imgLeafOnly=imgLeafOnly,
			dotfile=dotfile, devicefunc=devicefunc, imageext=imageext,
			device.arg=device.arg, use.title=use.title, label.loc=label.loc,
			node.loc=node.loc, split.loc=split.loc, title.cex=title.cex, ...)

	DTN2DotInternalRelation(preced=filename, node=node, pos="none", dotfile=dotfile)
	cat("}\n", file=dotfile, append=TRUE)

}
DTN2DotInternal <- function(preced, node, pos, digits, imagefunc, imagedata,
														imgLeafOnly, dotfile, devicefunc, imageext, device.arg,
														use.title, label.loc=label.loc, node.loc, split.loc,
														title.cex, ...) {
	nodename <- paste(preced, "_", pos, sep="")
	stringcontentnode <- paste("Size: ", length(node$ind), " Var: ",
			format(node$vardis, digits =digits), sep="")
	if (!is.null(node$split)) {
		stringcontentsplit <- paste("Split: ", node$split$varname, " R2:",
				format(node$R2, digits=digits), sep="")
	} else {
		stringcontentsplit <- ""
	}
	if (!is.null(imagefunc) && (!imgLeafOnly || is.null(node$split))) {
		device.arg$file <- paste(nodename, imageext, sep=".")
		do.call(devicefunc, device.arg)
		imagefunc(imagedata[node$ind, ], ...)
		if (use.title) {
			title.arg <- list(main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
						line = NA, outer = FALSE, cex.main=title.cex,	cex.sub=title.cex,
						cex.lab=title.cex)
			title.arg[[node.loc]] <- stringcontentnode
			title.arg[[split.loc]] <- stringcontentsplit
			title.arg[[label.loc]] <- node$label
			if (label.loc==node.loc) {
				title.arg[[label.loc]] <- paste(node$label, stringcontentnode, sep="\n")
			}
			if (label.loc==split.loc) {
				title.arg[[label.loc]] <- paste(title.arg[[label.loc]], stringcontentsplit, sep="\n")
			}
			if (node.loc==split.loc) {
				if (label.loc!=split.loc) {
					title.arg[[node.loc]] <- paste(stringcontentnode, stringcontentsplit, sep="\n")
				}
			}
			do.call("title", title.arg)
		}
		dev.off()
		imgstr <- paste(" image=\"", nodename, ".jpg\", imagescale=true, ", sep="")
	}
	else {
		imgstr <- ""
	}
	if (!use.title) {
		str <- paste("\"", nodename, "\"[shape=box, ", imgstr, " label=", "\"",
				stringcontentnode, stringcontentsplit, "\"", sep="")
	}
	else {
		str <- paste("\"", nodename, "\"[shape=box, ", imgstr, " label=\" \"", sep="")
	}
	cat(paste(str, "];\n", sep=""), file=dotfile, append=TRUE)
	if (!is.null(node$split)) {
		DTN2DotInternal(preced=nodename, node=node$children$left, pos="left", digits=digits,
				imagefunc=imagefunc, imagedata=imagedata, imgLeafOnly=imgLeafOnly, dotfile=dotfile,
				devicefunc=devicefunc, imageext=imageext, device.arg=device.arg, use.title=use.title,
				label.loc=label.loc, node.loc=node.loc, split.loc=split.loc, title.cex=title.cex, ...)
		DTN2DotInternal(preced=nodename, node=node$children$right, pos="right", digits=digits,
				imagefunc=imagefunc, imagedata=imagedata, imgLeafOnly=imgLeafOnly, dotfile=dotfile,
				devicefunc=devicefunc, imageext=imageext, device.arg=device.arg, use.title=use.title,
				label.loc=label.loc, node.loc=node.loc, split.loc=split.loc, title.cex=title.cex, ...)
	}
}

DTN2DotInternalRelation <- function(preced, node, pos, dotfile){
	nodename <- paste(preced, "_", pos, sep="")
	if (node$depth!=1) {
		cat(paste("\"", preced, "\"->", "\"", nodename, "\";\n", sep=""), file=dotfile, append=TRUE)
	}
	if (!is.null(node$split)) {
		DTN2DotInternalRelation(preced=nodename, node=node$children$left, pos="left", dotfile=dotfile)
		DTN2DotInternalRelation(preced=nodename, node=node$children$right, pos="right", dotfile=dotfile)
	}
	

}