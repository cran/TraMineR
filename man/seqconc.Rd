\encoding{latin1}
\name{seqconc}
\alias{seqconc}
\title{Concatenate vectors of states or events into a character string}
\description{
Concatenate vectors of states or events into a character string. In the string, each state is separated by 'sep'. The void elements in the input sequences are eliminated.
}
\usage{
seqconc(data, var=NULL, sep="-", vname="Sequence", void=NA)
}
\arguments{
  \item{data}{a dataframe or matrix containing sequence data.}

  \item{var}{the list of columns containing the sequences. Defaut to NULL, ie all the columns. Whether the sequences are in the compressed (character strings) or extended format is automatically detected by counting the number of columns.}

  \item{sep}{the character used as separator. By default, "-".}

  \item{vname}{an optional name for the variable containing the sequences. By default, "Sequence".}

	\item{void}{the code used for void elements appearing in the sequences (see \cite{Gabadinho et al. (2008)} for more details on missing values and void elements in sequences). Default to NA.}

}
\value{a vector of character strings, one for each row in the input data.}

\references{
Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2008). Mining Sequence Data in \code{R} with the \code{TraMineR} package: A user's guide. \emph{Department of Econometrics and Laboratory of Demography, University of Geneva}.
}

\seealso{
 \code{\link{seqdecomp}}.
}
\examples{
data(actcal)
actcal.string <- seqconc(actcal,13:24)
head(actcal.string)
}
\keyword{manip}
