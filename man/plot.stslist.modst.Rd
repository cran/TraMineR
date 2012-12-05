\name{plot.stslist.modst}
\alias{plot.stslist.modst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for modal state sequences}
\description{
Plot method for output produced by the seqmodst function, i.e objects of class stslist.modst.
}
\usage{
\method{plot}{stslist.modst}(x, cpal = NULL, ylab = NULL, yaxis = TRUE, xaxis = TRUE,
    xtlab = NULL, xtstep = NULL, cex.plot = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{stslist.modst} as produced by the \code{seqmodst} function.}
  \item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of states in the alphabet. By default, the 'cpal' attribute of the \code{x} object is used.}
  \item{ylab}{an optional label for the y axis. If set to NA, no label is drawn.}

  \item{yaxis}{if TRUE (default) the y axis is plotted.}

  \item{xaxis}{if TRUE (default) the x axis is plotted.}

  \item{xtlab}{optional labels for the x axis ticks. If unspecified, the \code{names} attribute of the \code{x} object is used.}

  \item{xtstep}{optional interval at which the tick-marks and labels of the x-axis are displayed. For example, with xtstep=3 a tick-mark is drawn at position 1, 4, 7, etc... The display of the corresponding labels depends on the available space and is dealt with automatically. If unspecified, the \code{xtstep} attribute of the \code{x} object is used.}

  \item{cex.plot}{expansion factor for setting the size of the font for the axis labels and names. The default value is 1. Values lesser than 1 will reduce the size of the font, values greater than 1 will increase the size.}

  \item{\dots}{further graphical parameters. For more details about the graphical parameter arguments, see \code{barplot} and \code{par}. }
}
\details{This is the plot method for the output produced by the \code{\link{seqmodst}} function, i.e. objects of class \emph{stslist.modst}. It produces a plot showing the sequence of modal states with bar width proportional to the state frequencies.

This method is called by the generic \code{\link{seqplot}} function (if \code{type="ms"}) that produces more sophisticated plots, allowing grouping and automatic display of the states legend. The \code{seqmsplot} function is a shortcut for calling \code{seqplot} with \code{type="ms"}.
}

\examples{
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)

## Modal state sequence
biofam.modst <- seqmodst(biofam.seq)
plot(biofam.modst)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Plot}
