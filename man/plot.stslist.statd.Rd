\name{plot.stslist.statd}
\alias{plot.stslist.statd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for objects produced by the seqstatd function}

\description{
This is the plot method for output produced by the \code{\link{seqstatd}} function, i.e for objects of class \emph{stslist.statd}.
}
\usage{
\method{plot}{stslist.statd}(x, type = "d", cpal = NULL, ylab = NULL,
  yaxis = TRUE, xaxis = TRUE, xtlab = NULL, xtstep = NULL, tick.last = NULL,
  cex.axis = 1, space = 0, xlab = NULL, cex.plot, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{stslist.statd} as produced by the \code{\link{seqstatd}} function.}
  \item{type}{if \code{"d"} (default), a state distribution plot is produced. If \code{"Ht"} an entropy index plot is produced.}
  \item{cpal}{alternative color palette to be used for the states. If user specified, a vector of colors with number of elements equal to the number of states in the alphabet. By default, the 'cpal' attribute of the \code{x} object is used.}
  \item{ylab}{an optional label for the y axis. If set to NA, no label is drawn.}
  \item{yaxis}{if \code{TRUE} or \code{"cum"}, the y axis is plotted with a label showing the cumulated percentage frequency of the displayed sequences. If \code{"pct"}, the percentage value for each sequence is displayed.}

  \item{xaxis}{if \code{TRUE} (default) the x-axis is plotted.}

  \item{xtlab}{optional labels for the ticks of the x-axis. If unspecified, the \code{names} attribute of the input \code{x} object is used.}
  \item{xtstep}{optional interval at which the tick-marks and labels of the x-axis are displayed. For example, with xtstep=3 a tick-mark is drawn at position 1, 4, 7, etc... The display of the corresponding labels depends on the available space and is dealt with automatically. If unspecified, the \code{xtstep} attribute of the \code{x} object is used.}
  \item{tick.last}{Logical. Should a tick mark be enforced at the last position on the x-axis? If unspecified, the \code{tick.last} attribute of the \code{x} object is used.}

  \item{cex.axis}{Axis annotation magnification. See \code{\link{par}}.}

  \item{space}{the space between the stacked bars. Default is 0, i.e. no space.}

  \item{xlab}{Optional title for the x axis. See \code{\link{title}}.}

  \item{cex.plot}{Deprecated. Use \code{cex.axis} instead.}

  \item{\dots}{further graphical parameters such as \code{border=NA} to remove the borders of the bars. For more details about the graphical parameter arguments, see \code{barplot} and \code{par}.}
}
\details{This is the plot method for the output produced by the \code{\link{seqstatd}} function, i.e. for objects of class \emph{stslist.statd}. If \code{type="d"} it produces a state distribution plot presenting the sequence of the transversal state frequencies at each successive (time) position, as computed by the \code{\link{seqstatd}} function. With \code{type="Ht"}, the series of entropies of the transversal state distributions is plotted.

This method is called by the generic \code{\link{seqplot}} function (if \code{type="d"} or \code{type="Ht"}) that produces more sophisticated plots, allowing grouping and automatic display of the state color legend. The \code{seqdplot} and \code{seqHtplot} functions are shortcuts for calling \code{seqplot} with \code{type="d"} or \code{type="Ht"} respectively.
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
\keyword{Plot}

