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
    can be called by setting  this argument to \code{'CDIST_O'}, \code{'CWIN'},
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
  subsequences are counted by occurrence by object (\code{"COBJ"} method).
  Alternatives are (2) \code{"CDIST_O"} (distinct occurrences with possibility of
  event-timestamp overlap), (3) \code{"CWIN"} (one occurrence of span-window),
  (4) \code{"CMINWIN"} (number of minimal windows of occurrence) and (5) \code{"CDIST"} (distinct occurrences without event-timestamp overlap). See
  references.
}

\value{
  A constraint object containing one item per constraint type.
}

\references{
  Joshi, Mahesh V., George Karypis, and Vipin Kumar (2001) A Universal Formulation of
  Sequential Patterns \emph{Proceedings of the KDD'2001 Workshop on
  Temporal Data Mining}, San Francisco.
  }

\seealso{\code{\link{seqefsub}}, \code{\link{seqeapplysub}}}

