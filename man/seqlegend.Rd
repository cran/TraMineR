\name{seqlegend}
\alias{seqlegend}
\title{Plots a legend for the states in a sequence object}
\description{
  Plots a legend for the states in a sequence object. Useful if several graphics are plotted together and only one legend is necessary. Unless specified by the user, the \emph{cpal} and \emph{labels} attributes of the sequence object are used for the colors and text appearing in the legend (see \code{\link{seqdef}}).
}
\usage{
seqlegend(seqdata, cpal=NULL, ltext=NULL, position="topleft", fontsize=1)
}
\arguments{
  \item{seqdata}{a sequence object created with the \code{\link{seqdef}} function.}
  \item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of distinct states. By default, the 'cpal' attribute of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}
  \item{ltext}{optional description of the states to appear in the legend. Must be a vector of character strings with number of elements equal to the number of distinct states. If unspecified, the 'label' attributes of the 'seqdata' sequence object is used (see \code{\link{seqdef}}). }
  \item{position}{the position of the legend in the graphic area. For accepted values, see \code{\link{legend}}. Default to 'topleft'.}
  \item{fontsize}{size of the font for the labels. A value less than 1 decreases the font size, a value greater than 1 increases the font size. Default to 1. }
}
%\seealso{}
\examples{
## Loading the 'actcal' example data set
## and defining a sequence object with 
## (activity statuses from jan. to dec. 2000)
## the data in columns 13 to 24
data(actcal)
actcal.seq <- seqdef(actcal,13:24,
	labels=c("> 37 hours", "19-36 hours", "1-18 hours", "no work"))

## Plotting the sequences frequency,
## the states distribution  
## and the legend
par(mfrow=c(2,2))
seqiplot(actcal.seq, tlim=0, withlegend=FALSE, withborder=FALSE)
seqfplot(actcal.seq, withlegend=FALSE)
seqdplot(actcal.seq, withlegend=FALSE)
seqlegend(actcal.seq)
}
\keyword{hplot}
