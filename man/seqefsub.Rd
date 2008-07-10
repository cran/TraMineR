\name{seqefsub}
\alias{seqefsub}
\title{Searching for frequent subsequences}
\description{
  Return a list of frequent subsequences satisfying a minimum support. Several time constraints can be set to restrict the search to specific time periods or subsequences durations.
}
\details{
  The support is counted per sequence and not per occurence. The support can be set through \code{pMinSupport} as a percentage (between 0 and 1 and it will be rounded), or throught minSupport as number of sequence.
  It is possible to specify time constraints using \code{maxGap}, \code{windowSize}, \code{ageMin}, \code{ageMax} and \code{ageMaxEnd}. If so, two event should not be separated by more than \code{maxGap} and the whole subsequence should be included in a maximum time of \code{windowSize}.
  The other parameters specify the start and end age of the subsequence, it should start between \code{ageMin} and \code{ageMax} and finish before \code{ageMaxEnd}.


}
\usage{
seqefsub(seq,minSupport=NULL, pMinSupport=NULL,maxGap=-1, windowSize=-1, ageMin=-1, ageMax=-1,ageMaxEnd=-1, maxK=-1)
}

\arguments{
  \item{seq}{A list of event sequences}
   \item{minSupport}{The minimum support (in number of sequence)}
   \item{pMinSupport}{The minimum support (in percentage, will be rounded)}
  \item{maxGap}{The maximum time gap}
  \item{windowSize}{The maximum window (subsequences time}
      \item{ageMin}{Can be used to set a time period. If equal to -1 (default), it won't be considered.}
  \item{ageMax}{Can be used to set a time period. If equal to -1 (default), it won't be considered.}
  \item{ageMaxEnd}{Can be used to set a time period. If equal to -1 (default), it won't be considered.}
  \item{maxK}{The maximum number of event in a subsequence}
}
\value{
  \item{subseq}{A list of subsequence (event sequence object) as a \code{seqelist}.}
  \item{support}{A list with the support of each subsequence in number of occurrence.}
}
\seealso{
 See \code{\link{seqecreate}} for creating event sequences. See \code{\link{seqeapplysub}} to count the number of occurence of frequent subsequences in each sequence.
  See \code{\link{is.seqelist}} about \code{seqelist}.
}
\examples{
data(actcal.tse)
actcal.seqe<-seqecreate(actcal.tse$id,actcal.tse$time,actcal.tse$event)

#Searching for frequent subsequences, that is, appearing at least 20 times
fsubseq<-seqefsub(actcal.seqe,minSupport=20)
#The same using a percentage
fsubseq<-seqefsub(actcal.seqe,pMinSupport=0.01)
#Getting a string representation of subsequences
#Ten first subsequences
fsubseq$subseq[1:10]
#Support of this subsequences
fsubseq$support[1:10]


##Using time constraints
##Looking for subsequence starting in summer (between june and september)
fsubseq<-seqefsub(actcal.seqe,minSupport=10, ageMin=6, ageMax=9)
fsubseq$subseq[1:10]

##Looking for subsequence contained in summer (between june and september)
fsubseq<-seqefsub(actcal.seqe,minSupport=10, ageMin=6, ageMax=9, ageMaxEnd=9)
fsubseq$subseq[1:10]

##Looking for subsequence enclosed in a 6 month period and with a maximum gap of 2 month
fsubseq<-seqefsub(actcal.seqe,minSupport=10, maxGap=2, windowSize=6)
fsubseq$subseq[1:10]
}
\keyword{nonparametric}
