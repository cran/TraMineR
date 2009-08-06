\name{seqeconstraint}
\alias{seqeconstraint}
\alias{print.seqeconstraint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Setting time constraint}
\description{
  Function used to set time constraints in event sequence methods (seqe..) such as \code{\link{seqefsub}} for searching frequent subsequences or \code{\link{seqeapplysub}} for checking occurrences of subsequences.
}
\usage{
seqeconstraint(maxGap = -1, windowSize = -1, 
  ageMin = -1, ageMax = -1, ageMaxEnd = -1, countMethod = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxGap}{The maximum time gap between to events}
  \item{windowSize}{The maximum time span accepted for subsequences}
  \item{ageMin}{Minimal start time position allowed for subsequences. Ignored when equal to -1 (default).}
  \item{ageMax}{Maximal start time position allowed for subsequences. Ignored when equal to -1 (default).}
  \item{ageMaxEnd}{Maximal end time position allowed for subsequences. Ignored when equal to -1 (default).}
  \item{countMethod}{By default, subsequences are counted only one time by sequence. If set to 2, each occurence of the subsequence in a sequence is counted.}
}
\details{
 \code{maxGap}, \code{windowSize}, \code{ageMin}, \code{ageMax} and \code{ageMaxEnd}. If so, two events should not be separated by more than \code{maxGap} and the whole subsequence should not exceed a \code{windowSize}  time span.
  The other parameters specify the start and end age of the subsequence, it should start between \code{ageMin} and \code{ageMax} and finish before \code{ageMaxEnd}. 
  Parameters \code{ageMin}, \code{ageMax} and \code{ageMaxEnd} are interpreted as the number of positions (time units) from the beginning of the sequence.
}
\value{
  A constraint object containing one item per constraint type.
}
\seealso{\code{\link{seqefsub}}, \code{\link{seqeapplysub}}}

