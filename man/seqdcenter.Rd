\name{seqdcenter}
\alias{seqdcenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute distance to center of a group }
\description{
  Compute the distance to a group center according to a matrix of distance. The method use the decomposition of distance used by the ward criteria.
  In this case, distance to center equal the mean of the two by two distances.
}
\usage{
seqdcenter(distmatrix, group=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distmatrix}{A matrix of distance such as the one computed by \code{seqdist}}
  \item{group}{ If null, only one group is considered, otherwise group to compute center}
}
\value{
  A vector with the distance to center of group for each sequence
}
\examples{
## Define a sequence object
data(biofam)
bf.seq <- seqdef(biofam,10:25)
couts <- seqsubm(bf.seq, method="TRATE")

## compute distance
biofam.om <- seqdist(bf.seq,method="OM",indel=2,sm=couts)
#compute distance to center according to gender
distcenter <- seqdcenter(biofam.om, biofam$sex)

## ploting distribution of dist to center as a silhouette plot
## Ploting two plots
par(mfrow=c(1,2))
barplot(sort(distcenter[biofam$sex=="man"]), main="Man",border=NA,space=0,col="cyan")
barplot(sort(distcenter[biofam$sex=="woman"]), main="Woman",border=NA,space=0,col="cyan")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
