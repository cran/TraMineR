\name{seqecmpgroup}
\alias{seqecmpgroup}
\title{Identifying discriminating subsequences}
\description{
  Identify and order the most discriminating subsequences according to a given statistical test.
}
\usage{
seqecmpgroup(subseq, group, method="chisq", pvalue.limit=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subseq}{The subsequences to order. To get a list of subsequence see \code{\link{seqefsub}}}
  \item{group}{The group variable to discriminate}
  \item{method}{The test, bonferroni or chisq}
  \item{pvalue.limit}{Can be used to filter the results. Only subsequences with p-value lower than this parameter will be kept. If NULL all subsequences are returned (regardless there p-value).}
}
\details{
The following test functions are implemented
  \code{bonferroni} Chi squared test with Bonferroni correction applied.
  \code{chisq} Chi squared test.
}
\value{
  Returns a specific objet of type \code{subseqelistchisq} (subtype of \code{subseqelist}) that can be plotted and printed.
  \item{subseq}{List of subsequences}
  \item{labels}{labels of groups variables}
  \item{constraint}{time constraints usued to compute subsequences see \code{\link{seqeconstraint}}}
  \item{seqe}{The event sequences that were used to compute the test}
  \item{index}{Index of susequence}
  \item{data}{A \code{data.frame} with the support of each subsequence, index (original order), two column per group with frequencies and pearson residuals}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{See Also \code{\link{plot.subseqelistchisq}} to plot the results}
\examples{
data(actcal.tse)
actcal.seqe <- seqecreate(actcal.tse)

##Searching for frequent subsequences, that is, appearing at least 20 times
fsubseq <- seqefsub(actcal.seqe,pMinSupport=0.01)

##searching for susbsequences discriminating the most men and women
data(actcal)
discr <- seqecmpgroup(fsubseq, group=actcal$sex, method="bonferroni")
##Printing discriminating subsequences
print(discr)
##Plotting the six most discriminating subsequences
plot(discr[1:6])

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
