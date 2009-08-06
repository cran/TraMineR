\name{seqmodst}
\alias{seqmodst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sequence of modal states}
\description{
Sequence made of the modal state at each position.
}
\usage{
seqmodst(seqdata, dist = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{a state sequence object as defined by the \code{\link{seqdef}} function.}
  \item{dist}{experimental}
  \item{\dots}{experimental}
}
\details{
In case of multiple modal states at a given position, the first one is taken. Hence, the result may vary with the alphabet order.
}
\value{
an object of class \emph{stslist.modst}. This is actually a state sequence object (containing a single state sequence) with additional attributes, among which the \code{Frequencies} attribute containing the transversal frequency of each state in the sequence. There are print and plot methods for such objects. More sophisticated plots can be produced with the \code{seqplot} function.
}

\seealso{\code{\link{plot.stslist.modst}} for default plot method, \code{\link{seqplot}} for higher level plots.}

\examples{
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)

## Modal state sequence
seqmodst(biofam.seq)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
