\name{dissassoc}
\alias{dissassoc}
\alias{hist.dissassoc}
\alias{print.dissassoc}
\encoding{latin1}
\title{Analysis of discrepancy based on dissimilarity measure}
\description{
Compute the discrepancy (defined by a dissimilarity measure) explained by a categorical variable.
}
\usage{
dissassoc(diss, group, R = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{diss}{A dissimilarity matrix or a dist object (see \code{\link{dist}})}
  \item{group}{The grouping variable}
  \item{R}{Number of permutations for computing the p-value. If equal to 1, no permutation test is performed.}
}
\details{
The association is based on a generalization of the ANOVA principle to any kind of distance metric. The test returns a pseudo R-squared that can be interpreted as a usual R-squared. The statistical significance of the association is computed by means of permutation tests.
This function also perform a test of discrepancy homogeneity (equality of variance) using a generalization of the T statistic.

There is a print method and hist method (to produce an histogram of the significance values).
}
\value{
Returns an object of class \code{dissassoc} with the following components:
  \item{groups}{A data frame containing the number of cases and the discrepancy of each group}
  \item{anova.table}{The pseudo ANOVA table}
  \item{stat}{The value of the statistics and their p-values}
  \item{perms}{The permutation object, see \code{\link{boot}}}
}
\references{
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009).
  Analyse de dissimilarités par arbre d'induction.
  \emph{Revue des Nouvelles Technologies de l'Information}, EGC'2009.

  Batagelj, V. (1988). Generalized Ward and related clustering problems. In H. Bock (Ed.),
  \emph{Classification and related methods of data analysis}, pp. 67-74. North-Holland, Amsterdam.

  Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance.
\emph{Austral Ecology} 26, 32-46.
}
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }
\seealso{
  \code{\link{dissvar}} to compute discrepancy using dissimilarities and for a basic introduction to concepts of discrepancy analysis

  \code{\link{disstree}} to analyse dissimilarities using induction trees

	\code{\link{dissmfac}} to perform multi-factor analysis of variance using dissimilarities

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
