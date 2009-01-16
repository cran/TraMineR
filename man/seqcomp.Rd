\name{seqcomp}
\alias{seqcomp}
\title{Compare two sequences}
\description{
Compare two sequences and return TRUE if they are equal and FALSE otherwise
}
\usage{
 seqcomp(x, y)
}
\arguments{
  \item{x}{a sequence object containing one sequence (typically a subscript of a main sequence object, see \code{\link{seqdef}})}
  \item{y}{a sequence object containing one sequence (typically a subscript of a main sequence object, see \code{\link{seqdef}})}
}

\value{TRUE if sequences are identical, FALSE otherwise}

\seealso{
 \code{}
}
\examples{
data(mvad)
mvad.shortlab <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad.seq <- seqdef(mvad, states=mvad.shortlab, 15:86)

## Comparing sequences 1 and 2 in mvad.seq
seqcomp(mvad.seq[1,],mvad.seq[2,])

## Comparing sequences 176 and 211 in mvad.seq
seqcomp(mvad.seq[176,],mvad.seq[211,])
}
\keyword{misc}
