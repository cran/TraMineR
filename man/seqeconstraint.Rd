\name{seqeconstraint}
\alias{seqeconstraint}
\alias{print.seqeconstraint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Setting time constraints and the counting method}
\description{
  Function used to set time constraints and the counting method in  methods (\code{seqe..}) for event sequences such as \code{\link{seqefsub}} for searching frequent subsequences or \code{\link{seqeapplysub}} for checking occurrences of subsequences.
}
\usage{
seqeconstraint(maxGap = -1, windowSize = -1, ageMin = -1,
     ageMax = -1, ageMaxEnd = -1, countMethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxGap}{The maximum time gap between two events}
  \item{windowSize}{The maximum time span accepted for subsequences}
  \item{ageMin}{Minimal start time position allowed for subsequences. Ignored when equal to -1 (default).}
  \item{ageMax}{Maximal start time position allowed for subsequences. Ignored when equal to -1 (default).}
  \item{ageMaxEnd}{Maximal end time position allowed for subsequences. Ignored when equal to -1 (default).}
  \item{countMethod}{By default, subsequences are counted only
    one time by sequence (\code{'COBJ'} method). Alternative counting methods
    are \code{'CDIST_O'}, \code{'CWIN'},
    \code{'CMINWIN'} or \code{'CDIST'} respectively. See details.}
}

\details{
  \code{maxGap}, \code{windowSize}, \code{ageMin}, \code{ageMax} and
  \code{ageMaxEnd}. If so, two events should not be separated by more
  than \code{maxGap} and the whole subsequence should not exceed a
  \code{windowSize}  time span.
  The other parameters specify the start and end age of the subsequence,
  it should start between \code{ageMin} and \code{ageMax} and finish
  before \code{ageMaxEnd}.
  Parameters \code{ageMin}, \code{ageMax} and \code{ageMaxEnd} are
  interpreted as the number of positions (time units) from the beginning
  of the sequence.

  There are 5 options for the \code{countMethod} argument. (1) By default,
  the count is the number of sequences that contain the subsequence (\code{"COBJ"} method).
  Alternatives are (2) \code{"CDIST_O"} (counts all distinct occurrences in each sequence including possibly overlapping occurrences, i.e., occurrences sharing a same event occurrence), (3) \code{"CWIN"} (number of slidden windows of length \code{windowSize}  that contain an occurrence of the subsequence),
  (4) \code{"CMINWIN"} (number of minimal windows of occurrence) and (5) \code{"CDIST"} (distinct occurrences without event occurrences overlap). See
  references.
}

\value{
  A constraint object containing one item per constraint type.
}

\references{
  Joshi, Mahesh V., George Karypis, and Vipin Kumar (2001) A Universal Formulation of
  Sequential Patterns \emph{Proceedings of the KDD'2001 Workshop on
  Temporal Data Mining}, San Francisco.

  Ritschard, G., A. Gabadinho, N.S. Müller and M. Studer (2008), Mining event sequences: A social science perspective, \emph{International Journal of Data Mining, Modelling and Management, IJDMMM}, 1(1), 68-90.
  }

\seealso{\code{\link{seqefsub}}, \code{\link{seqeapplysub}}}

\author{Matthias Studer, Nicolas S. Müller and Reto Bürgin (alternative counting methods) (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Event sequences}
