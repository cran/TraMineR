\name{seqeapplysub}
\alias{seqeapplysub}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Applying Subsequences to Event Sequences}
\description{
  Apply every subsequences \code{subseq} to each event sequences \code{seq} and compute the result of the given \code{method}.
}
\usage{
seqeapplysub(subseq, seq, method = "count", maxGap = -1, windowSize = -1, ageMin = -1, ageMax = -1, ageMaxEnd = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{subseq}{A list of subsequences}
 \item{seq}{A list of event sequences}
 \item{method}{Type of the result expected, should be one of "count", "presence" or "age"}
  \item{maxGap}{The maximum time gap between two group of event. If equal to -1 (default), it won't be considered.}
  \item{windowSize}{The maximum window (subsequence) time. If equal to -1 (default), it won't be considered.}
  \item{ageMin}{Can be used to set a time period. If equal to -1 (default), it won't be considered.}
  \item{ageMax}{Can be used to set a time period. If equal to -1 (default), it won't be considered.}
  \item{ageMaxEnd}{Can be used to set a time period. If equal to -1 (default), it won't be considered.}
}
\details{
  They are three methods implemented :
  \item{count}{Count the number of occurrence of a given subsequences in an event sequence.}
  \item{presence}{Returns one if the subsequence is present, 0 otherwise.}
  \item{age}{Returns the age of appearance of a subsequence in an event sequence. In case of multiple possiblities, the age of the first occurences is returned. When the subsequences is not in the sequence, -1 is returned.}

  It is possible to specify time constraints using \code{maxGap}, \code{windowSize}, \code{ageMin}, \code{ageMax} and \code{ageMaxEnd}. If so, two event should not be separated by more than \code{maxGap} and the whole subsequence should be included in a maximum time of \code{windowSize}.
  The other parameters specify the start and end age of the subsequence, it should start between \code{ageMin} and \code{ageMax} and finish before \code{ageMaxEnd}.

}
\value{
The return value is a matrix where every row correspond to a sequence (row names are set accordingly) and each column correspond to a subsequence (col names are set accordingly).
The matrix store the results of the count method.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~

\section{Warning}{Event sequence analysis module is still experimental}

\seealso{ See Also as \code{\link{seqecreate}} for more information on how to use event sequence analysis module}
\examples{
#loading data
data(actcal.tse)
#creating sequences
actcal.seqe<-seqecreate(actcal.tse$id,actcal.tse$time,actcal.tse$event)
##printing sequences
actcal.seqe[1:10]

#Looking for frequent subsequences
fsubseq<-seqefsub(actcal.seqe,pMinSupport=0.01)


#counting the number of occurence of each subsequence
msubcount<-seqeapplysub(fsubseq$subseq,actcal.seqe,method="count")
#First lines...
msubcount[1:10,1:10]
#Presence-absence of each subsequence
msubpres<-seqeapplysub(fsubseq$subseq,actcal.seqe,method="presence")
#First lines...
msubpres[1:10,1:10]

#age of first appearance of each subsequence
msubage<-seqeapplysub(fsubseq$subseq,actcal.seqe,method="age")
#First lines...
msubage[1:10,1:10]

#counting the number of occurence of each subsequence in summer
msubcount<-seqeapplysub(fsubseq$subseq,actcal.seqe,method="count",ageMin=6, ageMax=9, ageMaxEnd=9)
msubcount[1:10,1:10]
}
\keyword{nonparametric}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
