\name{seqeconstraint}
\alias{seqeconstraint}
\alias{print.seqeconstraint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Setting time constraint}
\description{
  Function used to set time constraint in several event sequences methods such as \code{\link{seqefsub}}(searching frequent subsequences).
}
\usage{
seqeconstraint(maxGap = -1, windowSize = -1, 
  ageMin = -1, ageMax = -1, ageMaxEnd = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxGap}{The maximum time gap between to event}
  \item{windowSize}{The maximum window (subsequences time)}
  \item{ageMin}{Can be used to set a time period. If equal to -1 (default), it won't be considered.}
  \item{ageMax}{Can be used to set a time period. If equal to -1 (default), it won't be considered.}
  \item{ageMaxEnd}{Can be used to set a time period. If equal to -1 (default), it won't be considered.}
}
\details{
 \code{maxGap}, \code{windowSize}, \code{ageMin}, \code{ageMax} and \code{ageMaxEnd}. If so, two event should not be separated by more than \code{maxGap} and the whole subsequence should be included in a maximum time of \code{windowSize}.
  The other parameters specify the start and end age of the subsequence, it should start between \code{ageMin} and \code{ageMax} and finish before \code{ageMaxEnd}.
}
\value{
  A constraint object containing one item per constraint type.
}
\seealso{\code{\link{seqefsub}}}

