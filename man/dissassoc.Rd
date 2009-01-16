\name{dissassoc}
\alias{dissassoc}
\alias{hist.dissassoc}
\alias{print.dissassoc}
\encoding{latin1}
\title{Analysis of pseudo-variance based on dissimilarity measure}
\description{
Compute the pseudo variance (defined by a dissimilarity measure) explained by a categorical variabel.
}
\usage{
dissassoc(diss, group, R = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diss}{A dissimilarity matrix or a dist object (see \code{\link{dist}})}
  \item{group}{The group variable}
  \item{R}{Number of permutation to compute the pvalue. If equal to 1, no permutation test are performed.}
}
\details{
The association is based on a generalization of the principe of ANOVA to any kind of distance metric. The test return a pseudo R squared that can be interpred as a usual R squared. The statistical significance of the association is computed using permutation test.
This function also perform a test of pseudo-variance homogeneity (equality of variance) using a generalization of the T statistic.

There is a print method and hist method (to plot an histogramme of the significance value). 
}
\value{
Return an object of class \code{dissassoc} with the following componant:
  \item{groups}{A \code{data.frame} containing the number of case and the pseudo-variance of each group}
  \item{anova.table}{The pseudo ANOVA table}
  \item{stat}{The value of the statistics and their p-value}
  \item{perms}{The permutation object, see \code{\link{boot}}}
}
\references{
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009).
  Analyse de dissimilarités par arbre d'induction.
  \emph{Revue des Nouvelles Technologies de l'Information}, EGC'2009.
  
  Batagelj, V. (1988). Generalized ward and related clustering problems. In H. Bock (Ed.),
  \emph{Classification and related methods of data analysis}, pp. 67-74. North-Holland, Amsterdam.
  
  Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance.
\emph{Austral Ecology} 26, 32-46.
}
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }
\seealso{
  \code{\link{dissvar}} to compute pseudo variance using dissimilarities and for a basic introduction to concepts of pseudo variance analysis
  
  \code{\link{disstree}} to analyse dissimilarities using induction trees
  
  \code{\link{dissreg}} to analyse dissimilarities in a way close to linear regression
  
  \code{\link{disscenter}} to compute the distance of each object to its center of group using dissimilarities

}
\examples{
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])

## Building dissimilarities
mvad.lcs <- seqdist(mvad.seq, method="LCS")

## R=1 imply no permutation test
da <- dissassoc(mvad.lcs, group=mvad$gcse5eq, R=10)
print(da)
hist(da)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
