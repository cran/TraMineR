\name{seqtransn}
\alias{seqtransn}
\title{Number of transitions in a sequence}
\description{
Computes the number of transitions (state changes) in each sequence of a sequence object.
}
\usage{
seqtransn(seqdata, with.missing=FALSE, norm=FALSE, pweight=FALSE)
}
\arguments{
	\item{seqdata}{a state sequence object as defined by the
\code{\link{seqdef}} function.}

	\item{with.missing}{logical: should non-void missing values be treated as a regular state? If \code{FALSE} (default) missing values are ignored.}

	\item{norm}{logical. If set as \code{TRUE}, the number of transitions is
divided by its theoretical maximum, the length of the sequence minus
1. When the length of the sequence is 1, the normalized value is set as 0.}

	\item{pweight}{logical. EXPERIMENTAL! If set as \code{TRUE}, the count of transition weights  
each transition by its probability of transition in the data.}
}

\details{A transition in a sequence is a state change between
time/position \eqn{t} and \eqn{t+1}. For example, the sequence
\code{"A-A-A-A-B-B-A-D-D-D"} contains 3 transitions. The maximum
number of transitions a sequence can contain is \eqn{\ell-1}{l-1}
where \eqn{\ell}{l} is the length of the sequence. The number of
transitions is obtained by subtracting 1 to the length of the the
Distinct Successive State (DSS) sequence. }

\value{a state sequence object containing the number of transitions
of each sequence in the object given as argument.}

\references{
Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical Software} \bold{40}(4), 1-37.
}

\seealso{
 \code{\link{seqdss}}.
}

\examples{
## Creating a sequence object from columns 13 to 24
## in the 'actcal' example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24)

## Computing the number of transitions
actcal.trans <- seqtransn(actcal.seq)

## Displaying the DSS for the first 10 sequences
actcal.trans[1:10]

## Example with with.missing argument
data(ex1)
ex1.seq <- seqdef(ex1, 1:13)

seqtransn(ex1.seq)
seqtransn(ex1.seq, with.missing=TRUE)
}

\author{Alexis Gabadinho (with Gilbert Ritschard for the help page)}

\keyword{Longitudinal characteristics}
