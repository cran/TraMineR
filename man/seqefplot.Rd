\name{seqefplot}
\alias{seqefplot}
\title{Plot frequencies of subsequences}
\description{
Plot frequencies of subsequences. If a group is specified, plot de difference of subsequence frequencies by the factor.
}
\usage{
seqefplot(subseq, seqe, group = NULL, mfrow=c(2,3),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subseq}{The subsequences to plot}
  \item{seqe}{The list of individual sequences}
  \item{group}{A factor or a vector describing a \code{group} }
  \item{mfrow}{When a group is used, a 2-length vector specifying the number of column and row (by default, only one row)}
%  \item{ylim}{A 2 length vector indicating minimum and maximum y values (by default to c(0,1))}
  \item{\dots}{arguments passed to \code{par} and \code{boxplot} commands}
}
%\details{}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
## loading data
data(actcal.tse)

## creating sequences
actcal.seqe <- seqecreate(actcal.tse$id,
	actcal.tse$time,
	actcal.tse$event)

## Looking for frequent subsequences
fsubseq <- seqefsub(actcal.seqe,pMinSupport=0.01)

## Frequence of first ten subsequences
seqefplot(fsubseq$subseq[1:10],actcal.seqe)

## Retrieve sex information

data(actcal)
## Plotting on 2 line 3 column
seqefplot(fsubseq$subseq[1:6],actcal.seqe, 
	group=actcal$sex, mfrow=c(2,3))
}
\keyword{hplot}

