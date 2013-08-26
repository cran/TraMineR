\name{seqmeant}
\alias{seqmeant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mean durations in each state}
\description{
Compute the mean durations spent in each state of the alphabet for the set of sequences given as input.
}
\usage{
seqmeant(seqdata, weighted=TRUE, with.missing=FALSE, prop=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{a sequence object as defined by the \code{\link{seqdef}} function.}
  \item{weighted}{logical: if \code{TRUE}, the weights (weights attribute) attached to the sequence object are used for computing weighted mean durations.}
  \item{with.missing}{logical: if set to \code{TRUE}, cumulated durations are also computed for the missing status (gaps in the sequences). See \link{seqdef} on options for handling missing values when creating sequence objects.}
  \item{prop}{logical: if \code{TRUE}, proportions of time spent in each state are returned instead of absolute values. This option is especially useful when sequences contain missing states, since the sum of the state durations may not be the same for all sequences.}
}
%\details{}

\value{
An object of class \emph{stslist.meant}. There are \code{print} and \code{plot} methods for such objects.
}

\seealso{\code{\link{plot.stslist.meant}} for basic plots of \emph{stslist.meant} objects and \code{\link{seqmtplot}} (\code{\link{seqplot}} with \code{type="mt"}) argument for more sophisticated plots of the mean durations allowing grouping and legend.}

\author{Alexis Gabadinho}

\references{
Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical Software} \bold{40}(4), 1-37.
}

\examples{
## Defining a sequence object with columns 13 to 24
## in the actcal example data set
data(actcal)
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal,13:24,labels=actcal.lab)

## Computing the mean durations
seqmeant(actcal.seq)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Global characteristics}
