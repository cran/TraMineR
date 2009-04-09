\encoding{latin1}
\name{cpal}
\alias{cpal}
\alias{cpal<-}
\title{Get or set the color palette of a sequence object}
\description{
This function gets or sets the color palette of a sequence object, that is, the list of colors used to represent the states.
}
\details{In the plot functions provided for visualizing sequence objects, each color is associated to one state of the alphabet. The color palette is defined when creating a sequence object, either automatically using the \code{brewer.pal} function provided by the \code{RColorBrewer} package or by specifying a user defined color vector. The \code{cpal} function can be used to get or set the color palette of a previously defined sequence object.}
\usage{
cpal(seqdata)
cpal(seqdata) <- value
}
\arguments{
  \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
  \item{value}{a vector containing the colors, of the same length as the number of states in the alphabet. The colors can be passed as character strings representing color names such as returned by the \code{\link{colors}} function, as hexadecimal values or as RGB vectors using the \code{\link{rgb}} function. Each color is attributed to the corresponding state in the alphabet, the order being the one returned by the \code{\link{alphabet}}.}
}
\value{For 'cpal' a vector containing the colors.

For 'cpal<-' the updated sequence object.
}
\seealso{
 \code{\link{seqdef}}
}
\examples{
## Creating a sequence object with the columns 13 to 24 
## in the 'actcal' example data set
## The color palette is automatically set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## Retrieving the color palette
cpal(actcal.seq) 
seqiplot(actcal.seq)

## Setting a user defined color palette
cpal(actcal.seq) <- c("blue","red", "green", "yellow")
seqiplot(actcal.seq)
}
\keyword{attribute}