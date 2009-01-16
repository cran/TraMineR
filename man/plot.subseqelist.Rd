\name{plot.subseqelist}
\alias{plot.subseqelist}
\title{Plot frequencies of subsequences}
\description{
Plot frequencies of subsequences.
}
\usage{
\method{plot}{subseqelist}(x, freq=NULL,cex=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The subsequences to plot}
  \item{freq}{The frequences to plot, support if NULL}
  \item{cex}{To control subsequences texts sizes. see \code{\link{par}}.}
 %  \item{ylim}{A 2 length vector indicating minimum and maximum y values (by default to c(0,1))}
  \item{\dots}{arguments passed to \code{\link{boxplot}}}
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
actcal.seqe <- seqecreate(actcal.tse)

## Looking for frequent subsequences
fsubseq <- seqefsub(actcal.seqe,pMinSupport=0.01)

## Frequence of first ten subsequences
plot(fsubseq[1:10], cex=2)
plot(fsubseq[1:10])

}
\keyword{hplot}

