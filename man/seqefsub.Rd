\name{seqefsub}
\alias{seqefsub}
\alias{[.subseqelist}
%\alias{createsubseqelist}
%\alias{seqecreatesub}
\alias{is.subseqelist}
\alias{print.subseqelist}
\title{Searching for frequent subsequences}
\description{
  Return a list of frequent subsequences satisfying a minimum support. Several time constraints can be set to restrict the search to specific time periods or subsequences durations.
}
\details{
  The support is counted per sequence and not per occurence. The support can be set through \code{pMinSupport} as a percentage (between 0 and 1 and it will be rounded), or throught minSupport as number of sequence.
  It is possible to specify time constraints using constraint argument (see \code{\link{seqeconstraint}} function).
  
  The \code{strsubseq} parameter can be used to loof for specific (user specified) subsequences. The format is the same as the one used to display subsequences (see \code{\link{str.seqelist}}).
  Each group of event should be enclosed in brace () and separeted with comma ,. A - can be used to specify a time gap.
  For instance "(FullTime)-(PartTime, Children)" will look for the subsequence "FullTime" event followed by the events "PartTime" and "Children" appearing at the same time.
}
\usage{
seqefsub(seq, strsubseq = NULL, minSupport = NULL, pMinSupport = NULL, 
  constraint = seqeconstraint(), maxK = -1)
}

\arguments{
  \item{seq}{A list of event sequences}
  \item{strsubseq}{Can be used to look for specific subsequences. See details.}
  \item{minSupport}{The minimum support (in number of sequence)}
  \item{pMinSupport}{The minimum support (in percentage, will be rounded)}
  \item{constraint}{Time constraint object, use \code{\link{seqeconstraint}}}
  \item{maxK}{The maximum number of event in a subsequence}
}
\value{
A \code{subseqelist} object which contain at least the following objects:
  \item{subseq}{A list of subsequence (event sequence object) as a \code{seqelist}.}
  \item{data}{A data.frame with subsequence specific information such as the support}
  \item{constraint}{The constraint object used to compute the subsequence.}
  \item{count}{The number of sequences.}

\code{subseqelist} can be selected using the [] operator (see example). There are print and plot methods for subsequelist
}
\seealso{
See \code{\link{plot.subseqelist}} to plot the result.
 See \code{\link{seqecreate}} for creating event sequences. See \code{\link{seqeapplysub}} to count the number of occurence of frequent subsequences in each sequence.
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
