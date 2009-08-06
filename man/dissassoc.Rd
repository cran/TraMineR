\encoding{latin1}
\name{dissassoc}
\alias{dissassoc}
\alias{hist.dissassoc}
\alias{print.dissassoc}
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
The association is based on a generalization of the ANOVA principle to any kind of distance metric. The test returns a 
pseudo R-squared that can be interpreted as a usual R-squared. The statistical significance of the association is computed 
by means of permutation tests.
This function also performs a test of discrepancy homogeneity (equality of variance) using a generalization of the T 
statistic.
\cr
There are \code{print} and \code{hist} methods (the latter producing an histogram of the significance values).
}
\value{
Returns an object of class \code{dissassoc} with the following components:
  \item{groups}{A data frame containing the number of cases and the discrepancy of each group}
  \item{anova.table}{The pseudo ANOVA table}
  \item{stat}{The value of the statistics and their p-values}
  \item{perms}{The permutation object, see \code{\link{boot}}}
}
\references{
  Studer, M., G. Ritschard, A. Gabadinho, and N. S. Müller (2009) 
  Discrepancy	analysis of complex objects using dissimilarities.
  In H. Briand, F. Guillet, G. Ritschard, and D. A. Zighed (Eds.),
  \emph{Advances in Knowledge Discovery and Management},
  Studies in Computational Intelligence. Berlin: Springer.
	
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009).
  Analyse de dissimilarités par arbre d'induction. In EGC 2009,
  \emph{Revue des Nouvelles Technologies de l'Information}, Vol. E-15, pp. 7--18.

  Batagelj, V. (1988) Generalized Ward and related clustering problems. In H. Bock (Ed.),
  \emph{Classification and related methods of data analysis}, Amsterdam: North-Holland, pp. 67--74.

  Anderson, M. J. (2001) A new method for non-parametric multivariate analysis of variance.
\emph{Austral Ecology} \bold{26}, 32--46.
}
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

\seealso{
\code{\link{dissvar}} to compute the pseudo variance from dissimilarities and for a basic introduction to concepts of 
pseudo variance analysis.
\cr
\code{\link{disstree}} for an induction tree analyse of objects characterized by a dissimilarity matrix.
\cr
\code{\link{disscenter}} to compute the distance of each object to its group center from pairwise dissimilarities.
 \cr
 \code{\link{dissmfac}} to perform multi-factor analysis of variance from pairwise dissimilarities.
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
