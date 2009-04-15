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
  \item{subseq}{A \code{subseqelist} object (list of subsequences) such as produced by \code{\link{seqefsub}}}
  \item{group}{Variable or factor defining the membership to the groups to discriminate}
  \item{method}{The required test, one of \code{bonferroni} or \code{chisq}}
  \item{pvalue.limit}{Can be used to filter the results. Only subsequences with a p-value lower than the value set for this parameter will be selected. If NULL all subsequences are returned (regardless their p-values).}
}
\details{
The following test functions are implemented
  \code{chisq} Pearson Independence Chi squared test.
  \code{bonferroni} Pearson Independence Chi squared test with Bonferroni correction.
}
\value{
  An objet of type \code{subseqelistchisq} (subtype of \code{subseqelist}) with the following elements 
  \item{subseq}{Sorted list of found discriminating subsequences}
  \item{seqe}{The event sequence object on which the tests were computed}
  \item{constraint}{time constraints used for searching the subsequences (see \code{\link{seqeconstraint}})}
  \item{labels}{levels (value labels) of the target group variable}
  \item{type}{Type of test used}
  \item{data}{A data frame with columns support, index (original order of the subsequence) and a pair of frequency and Pearson residual columns for each group}
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
fsubseq <- seqefsub(actcal.seqe, pMinSupport=0.01)

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
