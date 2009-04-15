\name{seqstatl}
\alias{seqstatl}
\title{List of distinct states or events (alphabet) in a sequence data set.}
\description{
  Returns a list containing distinct states or events found in a data frame or matrix containing sequence data, the alphabet.
}
\usage{
 seqstatl(data, var=NULL, format='STS')
}
\arguments{
  \item{data}{a data frame or matrix containing sequence data.}
  \item{var}{the list of columns containing the sequences. Default is NULL, i.e. all the columns. Whether the sequences are in the compressed (character strings) or extended format is automatically detected from the number of columns.}.
  \item{format}{the format of the sequence data set. One of 'STS', 'SPS', 'DSS'. Default is 'STS'. The \code{seqstatl} function uses the \code{\link{seqformat}} function to translate between formats when necessary.}
}
\seealso{
  \code{\link{seqformat}}
}
\examples{
data(actcal)
seqstatl(actcal,13:24)
}
\keyword{misc}
