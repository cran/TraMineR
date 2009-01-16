\name{seqfind}
\alias{seqfind}
\title{Find the occurrences of sequence(s) x in the set of sequences y}
\description{
Finds the occurrences of sequence(s) x in the set of sequences y. The function returns the indexes of sequence x in the y sequence object. 
}
\usage{
 seqfind(x, y)
}
\arguments{
  \item{x}{a sequence object containing one or more sequences.}
  \item{y}{a sequence object.}
}
\seealso{
 \code{}.
}
\value{index(es) of the occurence of sequence(s) x in the set of sequences y.}

\examples{
data(mvad)
mvad.shortlab <- c("EM", "FE", "HE", "JL", "SC", "TR")
mvad.seq <- seqdef(mvad, states=mvad.shortlab, 15:86)

## Finding occurrences of sequence 176 in mvad.seq
seqfind(mvad.seq[176,],mvad.seq)

## Finding occurrences of sequence 1 to 8 in mvad.seq
seqfind(mvad.seq[1:8,],mvad.seq)
}
\keyword{misc}
