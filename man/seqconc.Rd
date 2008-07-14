\name{seqconc}
\alias{seqconc}
\title{Concatenate vectors of states or events into a character string.}
\description{
Concatenate vectors of states or events into a character string. In the string, each state is separated by 'sep'. The missing values (NA's) are eliminated in the input sequences.
}
\usage{
 seqconc(data, var=NULL, sep="-",vname="Sequence")
}
\arguments{
  \item{data}{a dataframe or matrix containing sequence data.}
  \item{var}{the list of columns containing the sequences. Defaut to NULL, ie all the columns. Whether the sequences are in the compressed (character strings) or extended format is automatically detected by counting the number of columns.}
  \item{sep}{the character used as separator. By default, "-".}
  \item{vname}{an optional name for the variable containing the sequences. By default, "Sequence".}
}
\value{a vector of character strings, one for each row in the input data.}
\seealso{
 \code{\link{seqdecomp}}.
}
\examples{
data(actcal)
actcal.string <- seqconc(actcal,13:24)
}
\keyword{manip}
