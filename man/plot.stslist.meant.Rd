\name{plot.stslist.meant}
\alias{plot.stslist.meant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for objects produced by the seqmeant function}
\description{
This is the plot method for objects of class \emph{stslist.meant} produced by the \code{\link{seqmeant}} function.
}
\usage{
\method{plot}{stslist.meant}(x, cpal = NULL, ylab = NULL, yaxis = TRUE, xaxis = TRUE,
     cex.plot = 1, ylim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{stslist.meant} as produced by the \code{seqmeant} function.}
  \item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of states in the alphabet. By default, the '\code{cpal}' attribute of the '\code{seqdata}' sequence object is used (see \code{\link{seqdef}}).}
  \item{ylab}{an optional label for the y axis. If set to \code{NA}, no label is drawn.}
  \item{yaxis}{controls whether the y axis is plotted. Default is \code{TRUE}.}
  \item{xaxis}{if \code{TRUE} (default) the xaxis is plotted.}

  \item{cex.plot}{expansion factor for setting the size of the font for the axis labels and names. The default value is 1. Values lesser than 1 will reduce the size of the font, values greater than 1 will increase the size.}

  \item{ylim}{an optional vector setting the limits for the y axis. If \code{NULL} (default), limits are set to (0, max. sequence length).}

  \item{\dots}{further graphical parameters. For more details about the graphical parameter arguments, see \code{barplot} and \code{par}.}
}
\details{This is the plot method for the output produced by the \code{\link{seqmeant}} function, i.e., objects of class \emph{stslist.meant}. It produces a plot showing the mean times spent in each state of the alphabet.

This method is called by the generic \code{\link{seqplot}} function (if \code{type="mt"}) that produces more sophisticated plots, allowing grouping and automatic display of the states legend. The \code{\link{seqmtplot}} function is a shortcut for calling \code{seqplot} with \code{type="mt"}.
}

\examples{
## Loading the mvad data set and creating a sequence object
data(mvad)
mvad.labels <- c("employment", "further education", "higher education",
                    "joblessness", "school", "training")
mvad.scodes <- c("EM","FE","HE","JL","SC","TR")
mvad.seq <- seqdef(mvad, 15:86, states=mvad.scodes, labels=mvad.labels)

## Computing the mean times
mvad.meant <- seqmeant(mvad.seq)

## Plotting
plot(mvad.meant, main="Mean durations in each state of the alphabet")

## Changing the y axis limits
plot(mvad.meant, main="Mean durations in each state of the alphabet",
 ylim=c(0,40))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
