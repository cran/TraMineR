\name{seqdplot}
\alias{seqdplot}
\title{Graphic presenting the states frequencies}
\description{
  Graphic presenting the states frequencies for each time unit, as computed by the \code{\link{seqstatd}} function.
}
\usage{
seqdplot(seqdata, cpal=NULL, title=NULL, 
         withlegend=TRUE, withborder=TRUE, ltext=NULL, xtlab=NULL, bmar=1)
}
\arguments{
  \item{seqdata}{a sequence object created with the \code{\link{seqdef}} function.}
   \item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of distinct states. By default, the 'cpal' attribute of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}
  \item{title}{a character string giving the title for the graphic. Default to no title}
  \item{withlegend}{logical variable indicating wether a legend for the states colors is plotted. Default to TRUE}
  \item{withborder}{if TRUE, plots the borders of the bars. Default to TRUE.}
  \item{ltext}{optional description of the states to appear in the legend. Must be a vector of character strings with number of elements equal to the number of distinct states. If unspecified, the 'label' attributes of the 'seqdata' sequence object is used (see \code{\link{seqdef}}). }
  \item{xtlab}{optional labels for the x axis ticks labels. If unspecified, the column names of the 'seqdata' sequence object are used (see \code{\link{seqdef}}).}
 \item{bmar}{base margin for the top and bottom of the graphic, that is, the amount of white space surrouding it.}
}
\seealso{
 \code{\link{seqstatd}}
}
\examples{
## Loading the 'biofam' example data set
data(biofam)

## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30)
biofam.seq <- seqdef(biofam,10:25,informat='STS')

## Plotting the states distribution (family status from age 15 to 30)
seqdplot(biofam.seq)
}
\keyword{hplot}
