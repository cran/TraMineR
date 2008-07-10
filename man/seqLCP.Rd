\name{seqLCP}
\alias{seqLCP}
\title{Longest common prefix of two sequences.}
\description{
Returns the length of the longest common prefix of two sequences.
}
\usage{
seqLCP(iseq, jseq)
}
\arguments{
  \item{iseq}{a sequence from a sequence object or a vector of states.}
  \item{jseq}{a sequence from a sequence object or a vector of states.}
}
%\seealso{}
\examples{
data(famform)
famform.seq <- seqdef(famform)
## The LCP's length between sequences 1 and 2
## in the famform sequence object is 2
seqLCP(famform.seq[1,],famform.seq[2,])
}
\keyword{misc}
