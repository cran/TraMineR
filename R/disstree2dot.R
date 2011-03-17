###########################
## Internal plot method
###########################
DTNplotfunc <- function(imagedata, plotfunc, plotfunc.title.cex,
				plotfunc.use.title=TRUE, plotfunc.label.loc="main",
				plotfunc.node.loc="main", plotfunc.split.loc="sub", ...) {
	if (plotfunc.use.title) {
		top <- (as.integer(plotfunc.label.loc=="main")+as.integer(plotfunc.node.loc=="main")+as.integer(plotfunc.split.loc=="main"))*plotfunc.title.cex
		bottom <- (as.integer(plotfunc.label.loc=="sub")+as.integer(plotfunc.node.loc=="sub")+as.integer(plotfunc.split.loc=="sub"))*plotfunc.title.cex
		left <- (as.integer(plotfunc.label.loc=="ylab")+as.integer(plotfunc.node.loc=="ylab")+as.integer(plotfunc.split.loc=="ylab"))*plotfunc.title.cex
		par(mar=c(bottom, left, top, 0), font.sub=2, mgp=c(0, 0, 0))
		#on.exit(par(oldpar))
	}
	plotfunc(imagedata, ...)
}

DTNseqplot <- function(ind, seqdata, sortv=NULL, dist.matrix=NULL, ...) {
	if (!is.null(sortv)){
		seqplot(seqdata[ind, ], sortv=sortv[ind], ...)
	} else if(!is.null(dist.matrix)){
		seqplot(seqdata[ind, ], dist.matrix=dist.matrix[ind,ind], ...)
	}
	else {
		seqplot(seqdata[ind, ], ...)
	}
}

seqtreedisplay <- function(tree, filename=NULL, seqdata=tree$info$object, imgLeafOnly=FALSE, sortv=NULL, dist.matrix=NULL, title.cex=3, withlegend="auto", legend.fontsize=title.cex, axes=FALSE, imageformat="png", withquality=TRUE, legendtext=NULL, showtree=TRUE, ...) {
	actualdir <- getwd()
	tmpdir <- tempdir()
	tmpseqtree <- basename(tempfile(pattern="tmpseqtree"))
	setwd(tmpdir)
	if(withquality){
		
		rowStat <- function(x, n){
			formd <- function (x){
				return(format(x, digits =getOption("digits")-2))
			}
			star <- function(val){
				return(as.character(cut(val, c(0, 0.001, 0.01, 0.05, 1), labels=c("***", "**", "*",""))))
			}
			return(paste("<TR><TD ALIGN=\"LEFT\">",n,":</TD><TD ALIGN=\"LEFT\">", formd(x$info$adjustment$stat[n,1]),
						 "</TD><TD ALIGN=\"LEFT\">", star(x$info$adjustment$stat[n,2]),"</TD></TR>", sep=""))
		}
		legendtext <- paste("<FONT POINT-SIZE=\"",round(title.cex*11,0),"\"><TABLE BORDER=\"0\" CELLPADDING=\"5\"><TR><TD COLSPAN=\"3\">Global quality</TD></TR>",
							rowStat(tree, "Pseudo F"), rowStat(tree, "Pseudo R2"), rowStat(tree, "Levene"),"</TABLE></FONT>", sep="")
	}
	if(imageformat!="jpg"){
		seqtree2dot(tree=tree, filename=tmpseqtree, seqdata=seqdata, imgLeafOnly=imgLeafOnly,
			sortv=sortv, dist.matrix=dist.matrix, title.cex=title.cex, withlegend=withlegend,
			legend.fontsize=legend.fontsize, axes=axes, devicefunc="png", imageext="png", legendtext=legendtext, ...)
	}
	else {
		seqtree2dot(tree=tree, filename=tmpseqtree, seqdata=seqdata, imgLeafOnly=imgLeafOnly,
			sortv=sortv, dist.matrix=dist.matrix, title.cex=title.cex, withlegend=withlegend,
			legend.fontsize=legend.fontsize, axes=axes, legendtext=legendtext, ...)
	}
	
	myshellrun <- function(cmd, ...) {
		if (.Platform$OS.type=="windows") {
			return(system(paste(Sys.getenv("COMSPEC"),"/c", cmd),...))
		}
		else {
			return(system(cmd,...))
		}
	}
	if(imageformat!="jpg"){
		dotval <- myshellrun(paste("dot -Tpng -o", tmpseqtree,".png ", tmpseqtree,".dot", sep=""))
	}else{
		dotval <- myshellrun(paste("dot -Tjpg -o", tmpseqtree,".jpg ", tmpseqtree,".dot", sep=""))
	}
	if (dotval==1) {
		stop("You should install GraphViz to use this function: see http://www.graphviz.org")
	}
	
	if (!(imageformat %in% c("jpg", "png"))) {
		imagickval <- myshellrun(paste("convert ", tmpseqtree,".png ",tmpseqtree,".", imageformat, sep=""))
		if (imagickval == 1) {
			stop("To use another format than jpeg or png, you should install ImageMagick: see http://www.imagemagick.org")
		}
	}
    if (showtree) {
    	if (.Platform$OS.type=="windows") {
    		myshellrun(paste("start ", tmpseqtree, ".", imageformat, sep=""), wait=FALSE)
    	}
    	else {
    		myshellrun(paste("display ", tmpseqtree, ".", imageformat, sep=""), wait=FALSE)
    	}
    }
	
	setwd(actualdir)
	if(!is.null(filename)){
		file.copy(file.path(tmpdir, paste(tmpseqtree, imageformat, sep=".")), filename)
	}
	return(invisible())
}

###########################
## Shortcut to build sequences tree
###########################
seqtree2dot <- function(tree, filename, seqdata=tree$info$object, imgLeafOnly=FALSE, sortv=NULL,dist.matrix=NULL, title.cex=3, withlegend="auto", legend.fontsize=title.cex, axes=FALSE, ...) {
	legendimage <- NULL
	arguments <- list(...)
	if(withlegend!=FALSE) {
		devicefunc <- ifelse(is.null(arguments[["devicefunc"]]), "jpeg", arguments[["devicefunc"]])
		imageext <- ifelse(is.null(arguments[["imageext"]]), "jpg", arguments[["imageext"]])
		if (is.null(arguments[["device.arg"]])) {
			device.arg <- list()
		}
		else {
			device.arg <- arguments[["device.arg"]]
		}
		legendimage <- paste(filename,"legend", imageext, sep=".")
		device.arg$file <- legendimage
		
		do.call(devicefunc, device.arg)
		plot.new()
		seqlegend(seqdata, fontsize=legend.fontsize, title="Legend", position="center",  bty="n")
		dev.off()
	}
	disstree2dotp(tree, filename, imagedata=NULL, imgLeafOnly=imgLeafOnly, seqdata=seqdata, title.cex=title.cex,
			sortv=sortv,dist.matrix=dist.matrix, imagefunc=DTNseqplot, withlegend=FALSE, axes=axes, legendimage=legendimage, ...)
}


disstree2dotp <- function(tree, filename, imagedata=NULL, imgLeafOnly=FALSE,
						imagefunc=plot, title.cex=3, ...){
	disstree2dot(tree, filename, imagedata=imagedata, imgLeafOnly=imgLeafOnly, title.cex=title.cex, imagefunc=DTNplotfunc, plotfunc=imagefunc,
			use.title=TRUE, label.loc="main", node.loc="main", split.loc="sub",
			plotfunc.use.title=TRUE, plotfunc.label.loc="main", plotfunc.node.loc="main",
			plotfunc.split.loc="sub", plotfunc.title.cex=3, ...)

}


###########################
## Generate dot file
###########################

###########################
## Generate dot file
###########################
disstree2dot <- function(tree, filename, digits=3, imagefunc=NULL, imagedata=NULL,
							imgLeafOnly=FALSE, devicefunc="jpeg", imageext="jpg",
							device.arg=list(), use.title=TRUE, label.loc="main",
							node.loc="main", split.loc="sub", title.cex=1, legendtext=NULL, legendimage=NULL, showdepth=FALSE, ...) {
	dotfile <- paste(filename, ".dot", sep="")
	node <- tree$root
	cat("digraph distree{\n", file=dotfile)
	if (!is.null(legendimage)) {
		str <- paste("\"node_legendimage\"[shape=box, image=\"", legendimage,"\", imagescale=true, label=\" \"", sep="")
		cat(paste(str, "];\n", sep=""), file=dotfile, append=TRUE)
	}
	
	if (!is.null(imagefunc) && is.null(imagedata)) {
		imagedata <- as.data.frame(node$info$ind)
	}
	formd <- function (x){
		return(format(x, digits =digits))
	}
	nodedepthranking <- list()
	DTN2DotInternal <- function(preced, node, pos, label) {
		nodename <- paste(preced, "_", pos, sep="")
		nodedepthranking[[nodename]] <<- node$info$depth
		stringcontentnode <- paste("n: ", formd(node$info$n), " s2: ",
				formd(node$info$vardis), sep="")
		if (!is.null(node$split)) {
			stringcontentsplit <- paste("Split: ", colnames(tree$data)[node$split$varindex], " R2:",
					formd(node$split$info$R2), sep="")
		} else {
			stringcontentsplit <- ""
		}
		if (!is.null(imagefunc) && (!imgLeafOnly || is.null(node$split))) {
			device.arg$file <- paste(nodename, imageext, sep=".")
			do.call(devicefunc, device.arg)
			if (!is.null(ncol(imagedata))) {
				imagefunc(imagedata[node$info$ind, ], ...)
			} else {
			  	imagefunc(imagedata[node$info$ind], ...)
			}
			if (use.title) {
				title.arg <- list(main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
							line = NA, outer = FALSE, cex.main=title.cex,	cex.sub=title.cex,
							cex.lab=title.cex)
				title.arg[[node.loc]] <- stringcontentnode
				title.arg[[split.loc]] <- stringcontentsplit
				title.arg[[label.loc]] <- label
				if (label.loc==node.loc) {
					title.arg[[label.loc]] <- paste(label, stringcontentnode, sep="\n")
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
			imgstr <- paste(" image=\"", nodename,".", imageext,"\", imagescale=true, ", sep="")
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
		split_print <- function(sp){
			str_split <- character(2)
			if (!is.null(sp$breaks)) {
				str_split[1] <- paste("<=", formd(sp$breaks))
				str_split[2] <- paste(">", formd(sp$breaks))
			}
			else {
				## lgrp <- levels(factor(tree$data[,sp$varindex]))
				str_split[1] <- paste("[", paste(sp$labels[sp$index==1], collapse=","),"]")
				str_split[2] <- paste("[", paste(sp$labels[sp$index==2], collapse=","),"]")
			}
			if(!is.null(sp$naGroup)){
				str_split[sp$naGroup] <- paste(str_split[sp$naGroup], "with NA")
			}
			return(str_split)
		}
		if (!is.null(node$split)) {
			str_split <- split_print(node$split)
			pos <- c("left", "right")
			for (i in 1:2) {
				DTN2DotInternal(preced=nodename, node=node$kids[[i]], pos=pos[i], label=str_split[i])
			}
		}
	}

	DTN2DotInternalRelation <- function(preced, node, pos){
		nodename <- paste(preced, "_", pos, sep="")
		if (node$info$depth!=1) {
			cat(paste("\"", preced, "\"->", "\"", nodename, "\";\n", sep=""), file=dotfile, append=TRUE)
		}
		if (!is.null(node$split)) {
			pos <- c("left", "right")
			for (i in 1:2) {
				DTN2DotInternalRelation(preced=nodename, node=node$kids[[i]], pos=pos[i])
			}
		}
		

	}
	
	DTN2DotInternal(preced=filename, node=node, pos="none", label="Root")
	if (!is.null(legendtext)) {
		str <- paste("\"node_legendtext\"[shape=box, label=<", legendtext, ">", sep="")
		cat(paste(str, "];\n", sep=""), file=dotfile, append=TRUE)
	}
	if(showdepth){
		fontsize <- paste("[shape=box, fontsize=", title.cex*20, "];\n", sep="")
		unikrank <- unlist(unique(nodedepthranking))
		cat(paste(unikrank, collapse=fontsize),fontsize, file=dotfile, append=TRUE)
		nodename <- names(nodedepthranking)
		for(rank in unikrank){
			cat("{ rank=same ;", rank, ";", paste(nodename[nodedepthranking==rank], collapse="; "), ";}\n", file=dotfile, append=TRUE)
		}
		cat(paste(sort(unikrank), collapse=" -> "), ";\n", file=dotfile, append=TRUE)
	}
	DTN2DotInternalRelation(preced=filename, node=node, pos="none")
	
	
	cat("}\n", file=dotfile, append=TRUE)

} 