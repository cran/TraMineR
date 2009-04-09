\name{seqdplot}
\alias{seqdplot}
\title{Plot the sequence of state distributions}
\description{
Graphic presenting the states frequencies for each time point, as computed by the \code{\link{seqstatd}} function.
}
\usage{
seqdplot(seqdata, group=NULL, title=NULL, ...)
}
\arguments{
	\item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}

	\item{group}{plots one plot for each level of the factor given as argument.}

	\item{title}{a character string giving the title for the graphic. Default to no title.}

	\item{...}{arguments to be passed to the \code{\link{plot.stslist}} function or other graphical parameters. For example \code{axes} for controling the display of axes, \code{withlegend} for controling the display of a legend for the plot, etc ... }
}
\seealso{
 \code{\link{plot.stslist}, \link{seqstatd}}
}
\examples{
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)

## Plotting the states distribution (family status from age 15 to 30)
seqdplot(biofam.seq)

## Grouped by language
seqdplot(biofam.seq, group=biofam$plingu02, axes="bottom")
}
\keyword{hplot}
