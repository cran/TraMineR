\name{disscenter}
\alias{disscenter}
\title{ Compute distance to center of a group}
\description{
 Compute the distance to a group center according to a matrix of distance. The method use the decomposition of distance used by the Ward criteria.

}
\usage{
disscenter(diss, group=NULL, medoids.index=FALSE, max.iter=20)
}
\arguments{
  \item{diss}{A matrix of dissimilarity such as the one computed by \code{seqdist}, or a "dist" object}
  \item{group}{ If null, only one group is considered, otherwise group to compute center}
%  \item{trim}{Percentage of individuals that are considered as outliers and included in computation of distance to center}
  \item{medoids.index}{If TRUE, return the index of the first encontered most central sequence. One index per group is returned.}
  \item{max.iter}{Maximum number of iteration to trim the outliers}
}
\value{
  A vector with the distance to center of group for each sequence, or a list of medoid index
}
\seealso{
   \code{\link{dissvar}} to compute pseudo variance using dissimilarities and for a basic introduction to concepts of pseudo variance analysis
   
   \code{\link{dissassoc}} to test association between dissimilarity and another variable
   
   \code{\link{disstree}} to analyse dissimilarities using induction trees
   
   \code{\link{dissreg}} to analyse dissimilarities in a way close to linear regression
}
\examples{
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])

## Building dissimilarities
mvad.lcs <- seqdist(mvad.seq, method="LCS")

## Compute distance to center according to group gcse5eq
dc <- disscenter(mvad.lcs,group=mvad$gcse5eq)

## Ploting distribution of dissimilarity  to center 
boxplot(dc~mvad$gcse5eq, col="cyan")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
