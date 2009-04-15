\name{seqmtplot}
\alias{seqmtplot}
\title{Graphic presenting the mean time spent in each state of the alphabet}

\description{Graphic presenting the mean time spent in each state of the alphabet}

\usage{
seqmtplot(seqdata, group=NULL, title=NULL, ...)
}

\arguments{
  \item{seqdata}{a sequence object as returned by the \code{\link{seqdef}} function.}

  \item{group}{Plots one plot for each level of the factor given as argument.}

  \item{title}{a character string giving the title for the graphic. Default to no title.}

  \item{...}{arguments to be passed to the \code{\link{plot.stslist}} function or other graphical parameters. Use \code{axes} for controling the display of axes, \code{withlegend} for controling the display of a legend for the plot, etc... }
}
\seealso{
 \code{\link{plot.stslist}, \link{seqstatd}}
}
\examples{
## Loading the 'mvad' example data set
data(mvad)

## Defining a sequence object with the data in columns 15 to 86
mvad.seq <- seqdef(mvad,15:86)

## Plotting the mean times, 
## for each level of the gcse5eq factor
seqmtplot(mvad.seq, group=mvad$gcse5eq)
}
\keyword{hplot}
