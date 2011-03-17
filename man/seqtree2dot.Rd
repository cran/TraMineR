\name{seqtree2dot}
%\alias{disstree2dotp}
\alias{seqtree2dot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical representation of a dissimilarity tree}
\description{
This function offers shortcuts to generate a "dot" file and associated images files that can be used in GraphViz to get a graphical representation of the tree.
}
\usage{
seqtree2dot(tree, filename, seqdata = tree$info$object,
           imgLeafOnly = FALSE, sortv = NULL, dist.matrix = NULL,
           title.cex = 3, withlegend = "auto", 
           legend.fontsize = title.cex, axes = FALSE, ...)
  
%disstree2dotp(tree, filename, imagedata=NULL, imgLeafOnly=FALSE, imagefunc=plot, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{A tree object to be plotted as defined by \code{\link{disstree}}}
  \item{filename}{A filename, without extension, that will be used to generate image and dot files}
  \item{seqdata}{a sequence object as defined by the the \code{\link{seqdef}} function.}
  \item{imgLeafOnly}{If \code{TRUE}, only terminal node will be plotted}
  \item{sortv}{The name of an optional variable used to sort the data before plotting, see \code{\link{seqplot}}.}
  \item{dist.matrix}{The name of an optional dissimilarity matrix used to find representative sequences, \code{\link{seqrplot}}.}
  \item{title.cex}{Size of the font of the title of each plot. A value less than 1 decreases the font size, a value greater than 1 increases the font size. Defaults to 1.}
  \item{withlegend}{defines if and where the legend of the state colors is plotted. The default value \code{"auto"} sets the position of the legend automatically. Other possible value is \code{"right"}. Obsolete value \code{TRUE} is equivalent to \code{"auto"}.}
  \item{legend.fontsize}{Size of the font of the legend.}
  \item{axes}{if set to \code{"all"} (default value) x axes are drawn for each plot in the graphic. If set to \code{"bottom"} and \code{group} is used, axes are drawn only under the plots located at the bottom of the graphic area. If \code{FALSE}, no x axis is drawn.}
  \item{\dots}{other parameters that will be passed to \code{\link{seqplot}}}
}
\details{
	\code{\link{seqtreedisplay}} provides a much simpler interface to produce graphical reprentation of a \code{\link{seqtree}}.
	
	This function generates a "dot" file and one image file per node. For each node, it calls \code{\link{seqplot}} passing the selected lines of \code{seqdata} as argument. 
	You should at least specify the \code{type} of the plot to use (i.e. using \code{type="d"} for instance, see \code{\link{seqplot}} for more details).
	\code{seqtree2dot} is a shortcut for sequences objects using the plot function \code{\link{seqplot}}.
	
	
}
\value{
  Nothing but generates a file in the current working directory (see \code{\link{setwd}}).
}
\seealso{\code{\link{seqtreedisplay}} and \code{\link{seqtree}} for examples}

