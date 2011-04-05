\name{seqistatd}
\alias{seqistatd}
\title{States frequency for each individual sequence}
\description{
  Returns the state frequencies (total durations) for each sequence in the sequence object.
}
\usage{
 seqistatd(seqdata, with.missing=FALSE)
}
\arguments{
  \item{seqdata}{a sequence object (see \code{\link{seqdef}} function).}
  \item{with.missing}{if set to TRUE, total durations are also computed for the missing status (gaps in the sequences). See \link{seqdef} on options for handling missing values when creating sequence objects.} 
}

\references{
Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical Software} \bold{40}(4), 1-37.
}

%\seealso{}

\examples{
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
seqistatd(actcal.seq[1:10,])

## Example using "with.missing" argument
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)

seqistatd(ex1.seq)
seqistatd(ex1.seq, with.missing=TRUE)

}
\keyword{univar}
