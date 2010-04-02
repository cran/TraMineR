\encoding{latin1}
\name{seqfcheck}
\alias{seqfcheck}
\title{Check if sequences are in the compressed or extended format}
\description{Check wheter \emph{seqdata} contains sequences in the compressed format (as character strings, with states separated by a separator) or in the extended format (sequences are stored in a matrix, where each state is stored in a separate column of the matrix. For a more detailed description of the compressed and extended format, see \cite{Gabadinho, 2008}.
} 
\details{
Whether the sequence(s) are in compressed format is checked by counting the number of columns and searching for the \code{'-'} or \code{':'} separator. The function returns the separator if it has been found in the data. If the data contains more than one column, the data is supposed to be in the extended format, and 'X' is returned, unless some states codes contain the '-' character (eg. states coded with negative integer values), in which case \code{'-X'} is returned.

}
\usage{
seqfcheck(seqdata)
}
\arguments{
 \item{seqdata}{a vector, data frame or matrix containing sequence data.}
 }
\value{a character string coding the format of the sequence data, either \code{':'},\code{'-'},\code{'X'} or \code{'-X'}.} 
\seealso{
 \code{\link{seqconc}}, \code{\link{seqdecomp}}
}

\references{
Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2008). Mining Sequence Data in \code{R} with \code{TraMineR}: A user's guide. \emph{Department of Econometrics and Laboratory of Demography, University of Geneva}.
}

\examples{
## The sequences in the actcal data set
## are in the extended format
data(actcal)
head(actcal[,13:24])
seqfcheck(actcal[,13:24])

## The sequences in the famform data set
## are in the compressed format
data(famform)
famform
seqfcheck(famform)
}
\keyword{internal}
