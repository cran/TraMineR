\name{seqecmpgroup}
\alias{seqecmpgroup}
\title{Identifying discriminating subsequences}
\description{
  Identify and order the most discriminating subsequences according to a given statistical test.
}
\usage{
seqecmpgroup(subseq, seqe, group, method="bonferroni", p.valuelimit=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subseq}{The subsequences to order.}
  \item{seqe}{The list of event sequences}
  \item{group}{The group variable to discriminate}
  \item{method}{The test, bonferroni or chisq}
  \item{p.valuelimit}{Can be used to set a limit p-value (default to 0.05)}
}
\details{
The following test functions are implemented
  \code{bonferroni} returns a p-value with Bonferroni correction applied.
  \code{chisq} return the chisq value directly so decreasing should be set to TRUE.
}
\value{
  Returns a data.frame with
  \item{row name}{Name of parameter}
  \item{index}{Index of susequence}
  \item{stat}{Return value of the test}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(actcal.tse)
actcal.seqe<-seqecreate(actcal.tse$id,actcal.tse$time,actcal.tse$event)

#Searching for frequent subsequences, that is, appearing at least 20 times
fsubseq<-seqefsub(actcal.seqe,pMinSupport=0.01)

#searching for susbsequences discriminating the most men and women
data(actcal)
discr<-seqecmpgroup(fsubseq$subseq, actcal.seqe, group=actcal$sex, method="bonferroni")

#Plotting on 2 line 3 column the six most discriminating
seqefplot(fsubseq$subseq[discr$index[1:6]], actcal.seqe, group=actcal$sex, mfrow=c(2,3))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
