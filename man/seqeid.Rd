\name{seqeid}
\alias{seqeid}
\title{Retrieve unique ids from an event sequence object.}
\description{
  Retrieve the unique ids from an event sequence object or from a list of event sequence object.
}
\usage{
seqeid(s)
}
\arguments{
\item{s}{An event sequence object (as created with \code{\link{seqecreate}}) or a list of event sequence objects}
}
%\seealso{}
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse)
seqeid(actcal.seqe)
}

\author{Matthias Studer (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Event sequences}
\keyword{Sequence-object attributes}
