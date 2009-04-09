\name{seqtree2dot}
%\alias{disstree2dotp}
\alias{seqtree2dot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical representation of a dissimilarity tree}
\description{
This function offers shortcuts to generate a "dot" file and associated images files that can be used in GraphViz to get a graphical representation of the tree.
}
\usage{
seqtree2dot(tree, filename, seqdata, imgLeafOnly=FALSE, sortv=NULL, ...)
  
%disstree2dotp(tree, filename, imagedata=NULL, imgLeafOnly=FALSE, imagefunc=plot, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{A tree object to be plotted as defined by \code{\link{disstree}}}
  \item{filename}{A filename, without extension, that will be used to generate image and dot files}
  \item{seqdata}{a sequence object as defined by the the \code{\link{seqdef}} function.}
  \item{sortv}{The name of an optional variable used to sort the data before plotting, see \code{\link{plot.stslist}}.}
  %\item{imagedata}{a \code{data.frame} that will be passed to imagefunc, see details}
  \item{imgLeafOnly}{If \code{TRUE}, only terminal node will be plotted}
  %\item{imagefunc}{A function to plot the individuals in a node, see details}
  \item{\dots}{other parameters that will be passed to \code{\link{plot.stslist}}}
}
\details{
	This function generates a "dot" file and one image file per node. For each node, it calls \code{\link{plot.stslist}} passing the selected lines of \code{seqdata} as argument.
	\code{seqtree2dot} is a shortcut for sequences objects using the plot function \code{\link{plot.stslist}}.
}
\value{
  Nothing but generates a file in the current working directory (see \code{\link{setwd}}).
}
\seealso{\code{\link{disstree}} for examples}

