\name{plot.stslist.freq}
\alias{plot.stslist.freq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for sequence frequency tables}

\description{
Plot method for output produced by the \code{seqtab} function, i.e objects of class \code{stslist.freq}.
}
\usage{
\method{plot}{stslist.freq}(x, cpal = NULL, missing.color = NULL, pbarw = TRUE,
  ylab = NULL, yaxis = TRUE, xaxis = TRUE,
  xtlab = NULL, cex.plot = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{stslist.freq} as produced by the \code{seqtab} function.}

  \item{cpal}{alternative color palette to be used for the states. If user specified, a vector of colors with number of elements equal to the number of states in the alphabet. By default, the 'cpal' attribute of the \code{x} object is used.}

  \item{missing.color}{alternative color for representing missing values inside the sequences. By default, this color is taken from the \code{missing.color} attribute of the \code{x} object.}

  \item{pbarw}{if \code{pbarw=TRUE} (default), the width of the bars are proportional to the sequence frequency in the dataset.}

  \item{ylab}{an optional label for the y axis. If set to \code{NA}, no label is drawn.}

  \item{yaxis}{if \code{TRUE} or \code{"cum"}, the y axis is plotted with a label showing the cumulated percentage frequency of the displayed sequences. If \code{"pct"}, the percentage value for each sequence is displayed.}

  \item{xaxis}{if \code{TRUE} (default) the x-axis is plotted.}

  \item{xtlab}{optional labels for the ticks of the x-axis. If unspecified, the \code{names} attribute of the \code{x} object is used.}

  \item{cex.plot}{expansion factor for setting the size of the font for the axis labels and names. The default value is 1. Values smaller than 1 will reduce the size of the font, values greater than 1 will increase the size.}

  \item{\dots}{further graphical parameters. For example \code{border=NA} to remove the bars borders, \code{space=0} to remove space between sequences. For more details about the graphical parameter arguments, see \code{barplot} and \code{par}. }
}
\details{This is the plot method for the output produced by the \code{\link{seqtab}} function, i.e. objects of class \emph{stslist.freq}. It produces a plot showing the sequences sorted bottom up according to their frequency in the data set.

This method is called by the generic \code{\link{seqplot}} function (if \code{type="f"}) that produces more sophisticated plots, allowing grouping and automatic display of the state color legend. The \code{seqfplot} function is a shortcut for calling \code{seqplot} with \code{type="f"}.
}

\examples{
## Loading the 'actcal' example data set
data(actcal)

## Defining a sequence object with data in columns 13 to 24
## (activity status from january to december 2000)
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal, 13:24, labels=actcal.lab)

## 10 most frequent sequences in the data
actcal.freq <- seqtab(actcal.seq)

## Plotting the object
plot(actcal.freq, main="Sequence frequencies - actcal data set")

## Plotting all the distinct sequences without borders
## and space between sequences
actcal.freq2 <- seqtab(actcal.seq, tlim=0)
plot(actcal.freq2, main="Sequence frequencies - actcal data set",
  border=NA, space=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
