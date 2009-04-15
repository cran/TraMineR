\name{seqfplot}
\alias{seqfplot}
\title{Graphic presenting the frequency of sequences}
\description{
  The sequences are presented as stacks of successive states, ordered by their relative frequency in the dataset.
}
\usage{
 seqfplot(seqdata, group=NULL, title=NULL, tlim=NULL, pbarw=FALSE, ...)
}
\arguments{
	\item{seqdata}{a sequence object created with the \code{\link{seqdef}} function.}

	\item{group}{Plots one plot for each level of the factor given as argument.}

	\item{title}{title for the graphic. Default to NULL.}

	\item{tlim}{if tlim>0, shows only the 'tlim' most frequent sequences. If \code{tlim=NULL} (default), the 10 most frequent sequences are displayed.}

	\item{pbarw}{if pbarw=TRUE, the width of the bars are proportional to the sequence frequency in the dataset. Default to FALSE.}

	\item{...}{arguments to be passed to the \code{\link{plot.stslist}} function or other graphical parameters. Use \code{axes} for controling the display of axes, \code{withlegend} for controling the display of a legend for the plot, etc... }
}
\seealso{
 \code{\link{plot.stslist}, \link{seqtab}}.
}
\examples{
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)

## Plot of the 10 most frequent sequences 
## with equal bar width
seqfplot(biofam.seq)

## The same plot with
## proportional bar width
## grouped by sex
seqfplot(biofam.seq, group=biofam$sex, pbarw=TRUE)


## The same plot without legend,
## axes and labels
seqfplot(biofam.seq, group=biofam$sex, pbarw=TRUE, 
	axes=FALSE, ylab=NA, withlegend=FALSE)
}
\keyword{hplot}
