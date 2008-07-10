\name{seqmpos}
\alias{seqmpos}
\title{Number of matching positions between two sequences.}
\description{
Returns the number of common elements, ie same states appearing at the same position in the two sequences.
}
\usage{
 seqmpos(iseq, jseq)
}
\arguments{
  \item{iseq}{a sequence from a sequence object or a vector of states.}
  \item{jseq}{a sequence from a sequence object or a vector of states.}
}
\seealso{
 \code{}.
}
\examples{
data(famform)
famform.seq <- seqdef(famform)

seqmpos(famform.seq[1,],famform.seq[2,])
seqmpos(famform.seq[2,],famform.seq[4,])
}
\keyword{misc}
