\name{seqefsub}
\alias{seqefsub}
%\alias{[.subseqelist}
%\alias{createsubseqelist}
%\alias{seqecreatesub}
\alias{is.subseqelist}
\alias{print.subseqelist}
\title{Searching for frequent subsequences}
\description{
  Returns the list of frequent subsequences satisfying the specified minimum support. Several time constraints can be set to restrict the search to specific time periods or subsequences durations.
}
\details{
  There are two usages of this function. The first is for searching subsequences satisfying a support condition.
  The support is counted per sequence and not per occurrence, i.e. when a sequence contains twice a same subsequence it is counted only once. The support can be set through \code{pMinSupport} as a percentage (between 0 and 1 and it will be rounded), or through minSupport as a number of sequences.
  Time constraints can also be imposed with the \code{constraint} argument, which must be the outcome of a call to the \code{\link{seqeconstraint}} function).

  The second possibility is for searching sequences that contain specified subsequences. This is done by passing the list of subsequences with the \code{strsubseq} argument. The subsequences must be formatted as the one used to display subsequences (see \code{\link{str.seqelist}}).
  Each group of events should be enclosed in parentheses () and separated with commas, and the succession of events should be denoted by a '-' that indicates a time gap.
  For instance "(FullTime)-(PartTime, Children)" stands for the subsequence "FullTime" followed by the group of the two simultaneously occurring events "PartTime" and "Children".

  Information about the sequences that contain the subsequences can then be obtained with the \code{\link{seqeapplysub}} function.

  Subsets of the returned \code{subseqelist} can be accessed with the [] operator (see example). There are print and plot methods for \code{subsequelist}.
}
\usage{
seqefsub(seq, strsubseq = NULL, minSupport = NULL,
    pMinSupport = NULL, constraint = seqeconstraint(), maxK = -1)
}

\arguments{
  \item{seq}{A list of event sequences}
  \item{strsubseq}{Can be used to look for specific subsequences. See details.}
  \item{minSupport}{The minimum support (in number of sequences)}
  \item{pMinSupport}{The minimum support (in percentage, will be rounded)}
  \item{constraint}{Time constraint object, i.e the result of a call to \code{\link{seqeconstraint}}}
  \item{maxK}{The maximum number of events allowed in a subsequence}
}
\value{
A \code{subseqelist} object which contain at least the following objects:
  \item{seqe}{The list of sequences in which the subsequences were searched (a \code{seqelist} event sequence object).}
  \item{subseq}{A list of subsequences (a \code{seqelist} event sequence object).}
  \item{data}{A data frame containing details (support, frequency, ...) about the subsequences}
  \item{constraint}{The constraint object used when searching the subsequences.}
  \item{type}{The type of search: 'frequent' or 'user'}
  %\item{count}{The number of sequences.}
}
\seealso{
See \code{\link{plot.subseqelist}} to plot the result.
 See \code{\link{seqecreate}} for creating event sequences. See \code{\link{seqeapplysub}} to count the number of occurrences of frequent subsequences in each sequence.
  See \code{\link{is.seqelist}} about \code{seqelist}.

}
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse)

##Searching for frequent subsequences, that is, appearing at least 20 times
fsubseq <- seqefsub(actcal.seqe, minSupport=20)
##The same using a percentage
fsubseq <- seqefsub(actcal.seqe, pMinSupport=0.01)
##Getting a string representation of subsequences
##Ten first subsequences
fsubseq[1:10]

##Using time constraints
##Looking for subsequence starting in summer (between june and september)
fsubseq <- seqefsub(actcal.seqe, minSupport=10,
  constraint=seqeconstraint(ageMin=6, ageMax=9))
fsubseq[1:10]

##Looking for subsequence contained in summer (between june and september)
fsubseq <- seqefsub(actcal.seqe, minSupport=10,
  constraint=seqeconstraint(ageMin=6, ageMax=9, ageMaxEnd=9))
fsubseq[1:10]

##Looking for subsequence enclosed in a 6 month period
## and with a maximum gap of 2 month
fsubseq <- seqefsub(actcal.seqe, minSupport=10,
  constraint=seqeconstraint(maxGap=2, windowSize=6))
fsubseq[1:10]
}
\keyword{nonparametric}
