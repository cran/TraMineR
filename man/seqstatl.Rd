\name{seqstatl}
\alias{seqstatl}
\title{List of distinct states or events (alphabet) for a sequence data set.}
\description{
  Returns a list containing distinct states or events found in a dataframe or matrix containing sequence data, the alphabet.
}
\usage{
 seqstatl(data, var=NULL, format='STS')
}
\arguments{
  \item{data}{a dataframe or matrix containing sequence data.}
  \item{var}{the list of columns containing the sequences. Defaut to NULL, ie all the columns. Whether the sequences are in the compressed (character strings) or extended format is automatically detected by counting the number of columns.}.
  \item{format}{the format of the sequence data set. One of 'STS', 'SPS', 'DSS'. Default to 'STS'. The seqstatl function will use the seqformat function to translate between formats if necessary.}
}
\seealso{
  \code{\link{seqformat}}
}
\examples{
data(actcal)
seqstatl(actcal,13:24)
}
\keyword{misc}
