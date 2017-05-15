\name{alphabet}
\alias{alphabet}
\alias{alphabet<-}
\title{Get or set the alphabet of a sequence object}
\description{
This function gets or sets the (short) labels associated to the states in the alphabet of a sequence object (the list of all possible states, some of which states may not appear in the data).
}
\details{A state sequence object --- created with the \code{\link{seqdef}} function --- stores sequences as a matrix where columns are factors. The levels of the factors are made of the alphabet as well as the codes for missing value and void elements. The alphabet function retrieves or sets the "alphabet" attribute of the sequence object.
The state names composing the alphabet are preferably short labels, since they are used for printing sequences. Longer labels for describing more precisely each state in legend are stored in the "labels" attribute of the sequence object.}
\usage{
alphabet(seqdata)
alphabet(seqdata) <- value
}
\arguments{
  \item{seqdata}{a state sequence object as defined with the \code{\link{seqdef}} function.}
 \item{value}{a character vector of the same length as the vector returned by the \code{alphabet} function, i.e. one label for each state in the alphabet.}
 }
\value{For '\code{alphabet}' a character vector containing the alphabet.

For '\code{alphabet <-}' the updated sequence object.
}
\seealso{
 \code{\link{seqdef}}
}
\examples{
## Creating a sequence object with the columns 13 to 24
## in the 'actcal' example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## Retrieving the alphabet
alphabet(actcal.seq)

## Setting the alphabet
alphabet(actcal.seq) <- c("FT", "PT", "LT", "NO")
}

\author{Alexis Gabadinho}

\keyword{Sequence-object attributes}
