\name{disstree2dot}
\alias{disstree2dot}
\alias{disstree2dotp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical representation of a dissimilarity tree}
\description{
Generate a "dot" file and associated images files that can be used in GraphViz to get a graphical representation of the tree.
}
\usage{
disstree2dot(tree, filename, digits = 3,
   imagefunc = NULL, imagedata = NULL, imgLeafOnly = FALSE,
   devicefunc = "jpeg", imageext = "jpg", device.arg = list(),
   use.title = TRUE, label.loc = "main", node.loc = "main",
   split.loc = "sub", title.cex = 1, legendtext=NULL,
   legendimage=NULL, showdepth = FALSE, ...)

disstree2dotp(tree, filename,
   imagedata=NULL, imgLeafOnly=FALSE,
   imagefunc=plot, title.cex = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{The tree to be plotted.}
  \item{filename}{A filename, without extension, that will be used to generate image and dot files.}
  \item{digits}{Number of significant digits to plot.}
  \item{imagefunc}{A function to plot the individuals in a node, see details.}
  \item{imagedata}{a \code{data.frame} that will be passed to \code{imagefunc}, see details.}
  \item{imgLeafOnly}{Logical: If \code{TRUE}, only terminal node will be plotted.}
  \item{devicefunc}{A device function, \code{"jpeg"} by default.}
  \item{imageext}{extension for image files.}
  \item{device.arg}{Argument passed to \code{devicefunc}.}
  \item{use.title}{Logical: If \code{TRUE}, node information will be printed using \code{\link{title}} command, see details.}
  \item{label.loc}{Location of the node label, see \code{\link{title}} for possible values.}
  \item{node.loc}{Node content location, see \code{\link{title}} for possible values.}
  \item{split.loc}{Split information location, see \code{\link{title}} for possible values.}
  \item{title.cex}{\code{cex} applied to all calls to \code{title} (see \code{use.title}).}
  \item{legendtext}{An optional text appearing in a distinct node.}
  \item{legendimage}{An optional image file appearing in a distinct node.}
  \item{showdepth}{Logical. If \code{TRUE}, information about depth of the tree is added to the plot.}
  \item{\dots}{other parameters that will be passed to \code{imagefunc}.}
}
\details{
  This function generates a "dot" file that can be used in GraphViz (\url{http://www.graphviz.org}).
It also generates one image per node through a call to
\code{imagefunc} passing the selected lines of \code{imagedata} if
present or otherwise a list of indexes (of individuals belonging to a
node).

if \code{use.title} is \code{TRUE}, \code{imagefunc} should take care
to leave enough space for the title.

\code{disstree2dotp} is a simplified interface of
\code{disstree2dot} which automatically leaves enough space for the title and subtitles. These functions are intended to be generic. 

See \code{\link{seqtreedisplay}} for a much simpler way to generate a graphical representation of a \code{\link{seqtree}}; that is a \code{disstree} built using a sequence object. 

}
\value{
  Nothing but generates a "dot" and several images files (on per node) in the current working directory (see \code{\link{setwd}}).
}

\seealso{\code{\link{seqtree}} and \code{\link{seqtreedisplay}}, \code{\link{disstree}} for an example.}
