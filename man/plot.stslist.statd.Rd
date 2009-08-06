\name{plot.stslist.statd}
\alias{plot.stslist.statd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for objects produced by the seqstatd function}
\description{
This is the plot method for output produced by the \code{\link{seqstatd}} function, i.e objects of class \emph{stslist.statd}.
}
\usage{
\method{plot}{stslist.statd}(x, type = "d", cpal = NULL, ylab = NULL, 
  yaxis = TRUE, xaxis = TRUE, xtlab = NULL, cex.plot = 1, space=0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{stslist.statd} as produced by the \code{\link{seqstatd}} function.}
  \item{type}{if \code{"d"} (default), a state distribution plot is produced. If \code{"Ht"} an entropy index plot is produced.}
  \item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of states in the alphabet. By default, the 'cpal' attribute of the \code{x} object is used.}
  \item{ylab}{an optional label for the y axis. If set to NA, no label is drawn.}
  \item{yaxis}{if TRUE or "cum", the y axis is plotted with a label showing the cumulated percentage frequency of the displayed sequences. If "pct", the percentage value for each sequence is displayed.}

  \item{xaxis}{if TRUE (default) the xaxis is plotted.}

  \item{xtlab}{optional labels for the x axis ticks. If unspecified, the \code{names} attribute of the input object is used.}

  \item{cex.plot}{expansion factor for setting the size of the font for the axis labels and names. The default value is 1. Values lesser than 1 will reduce the size of the font, values greater than 1 will increase the size.}

  \item{space}{the space between the stacked bars. Default to 0, i.e. no space.}

  \item{\dots}{further graphical parameters. For example \code{border=NA} to remove the bars borders, \code{space=0} to remove space between sequences. For more details about the graphical parameter arguments, see \code{barplot} and \code{par}.}
}
\details{This is the plot method for the output produced by the \code{\link{seqstatd}} function, i.e. objects of class \emph{stslist.statd}. If \code{type="d"}it produces a state distribution plot presenting the sequence of the states frequencies for each time point, as computed by the \code{\link{seqstatd}} function. If \code{type="Ht"}, the series of state distribution entropies is plotted.

This method is called by the generic \code{\link{seqplot}} function (if \code{type="d"} or \code{type="Ht"}) that produces more sophisticated plots, allowing grouping and automatic display of the states legend. The \code{seqdplot} and \code{seqHtplot} functions are shortcuts for calling \code{seqplot} with \code{type="d"} or \code{type="Ht"} respectively.
}

\examples{
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)

## State distribution    
biofam.statd <- seqstatd(biofam.seq)

## State distribution plot (default type="d" option)
plot(biofam.statd)

## Entropy index plot
plot(biofam.statd, type="Ht")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

