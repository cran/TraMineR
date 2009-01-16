\name{disstree2dot}
\alias{disstree2dot}
\alias{seqtree2dot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical representation of a dissimilarity tree}
\description{
Generate a "dot" file and associated images files that can be used in GraphViz to get a graphical reprensetation of the tree.
}
\usage{
disstree2dot(tree, filename, digits = 3, 
  imagefunc = NULL, imagedata = NULL, imgLeafOnly = FALSE, 
  devicefunc = "jpeg", imageext = "jpg", device.arg = list(), 
  use.title = TRUE, label.loc = "main", node.loc = "main",
  split.loc = "sub", title.cex = 1, ...)

seqtree2dot(tree, filename, seqs, plottype="seqdplot", 
  imgLeafOnly=FALSE,sortv=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{The tree to be plotted}
  \item{filename}{A filename, without extension, that will be used to generate image and dot files}
  \item{seqs}{Sequences to be plotted}
  \item{plottype}{type of sequence plot, one of \code{\link{seqdplot}}, \code{\link{seqfplot}} and \code{\link{seqiplot}}}
  \item{sortv}{Used to sort sequence if \code{plottype} is equal to seqiplot}
  \item{digits}{Number of significant digits to plot}
  \item{imagefunc}{A function to plot the individuals in a node, see details}
  \item{imagedata}{a \code{data.frame} that will be passed to imagefunc, see details}
  \item{imgLeafOnly}{If \code{TRUE}, only terminal node will be plotted}
  \item{devicefunc}{A device function used, typically jpeg}
  \item{imageext}{extension for image files.}
  \item{device.arg}{Argument passed to \code{devicefunc}}
  \item{use.title}{If \code{TRUE}, node information will be printed using \code{title} command, see details}
  \item{label.loc}{Location of the node label}
  \item{node.loc}{Node content location}
  \item{split.loc}{Split information location}
  \item{title.cex}{\code{cex} applied to all title call (see \code{use.title}}
  \item{\dots}{other parameters that will be passed to \code{imagefunc}}
}
\details{
  This function call imagefunc using the selected lines of \code{imagedata} if present or a list of index if not.
}
\value{
  Nothing but generate a file in the current working directory (see \code{\link{setwd}}).
}
\seealso{\code{\link{disstree}} for example}

