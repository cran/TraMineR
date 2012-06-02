\name{seqdecomp}
\alias{seqdecomp}
\title{Convert a character string into a vector of states or events}
\description{
  For the moment, each character in the string will be considered to be one state or event = this function will not give accurate results if the character string representing the sequence contains events or states coded with more than one character.
}
\usage{
 seqdecomp(data, var=NULL, sep='-', miss="NA", vnames=NULL)
}
\arguments{
 \item{data}{a dataframe or matrix containing sequence data.}
  \item{var}{the list of columns containing the sequences. Default is \code{NULL}, ie all the columns. Whether the sequences are in the compressed (character strings) or extended format is automatically detected by counting the number of columns.}
  \item{sep}{the between states/events separator used in the input data set. Default to '\code{-}'.}
  \item{miss}{the symbol for missing values (if any) used in the input data set. Default is \code{NA}.}
  \item{vnames}{optional names for the column/variables of the output data set. Default is \code{NULL}.}
}
\seealso{
 \code{\link{seqconc}}.
}
\examples{
## Converts 'seq' into a vector of states of length 10
seq <- "A-A-A-A-B-B-B-C-C-C"
seqdecomp(seq)
}
\keyword{manip}
