\name{dissvar}
\alias{dissvar}
\encoding{latin1}
\title{Dissimilarity based discrepancy}
\description{
  Compute the discrepancy from the pairwise dissimilarities between objects. The discrepancy is a measure of dispersion of the set of objects.
}
\usage{
dissvar(diss)
}
\arguments{
  \item{diss}{A dissimilarity matrix or a \code{dist} object (see \code{\link{dist}})}
}
\details{
The discrepancy is an extension of the concept of variance to other kind of objects for which we have a dissimilarity measure.
The discrepancy \eqn{s^2} is defined as:
\deqn{s^2=\frac{1}{2n^2}\sum_{i=1}^{n}\sum_{j=1}^{n}d_{ij}}{s^2=(1/(2n^2)) * sum sum d_ij}

\emph{Mathematical ground}:
In the Euclidean case, the sum of squares can be expressed as:
  \deqn{SS=\sum_{i=1}^{n}(y_i-\bar{y})^2=\frac{1}{2n}\sum_{i=1}^{n}\sum_{j=1}^{n}(y_i-y_j)^2}{SS= sum (y_i - y_mean)^2=(1/(2n)) * sum sum (y_i - y_j)^2}

  The concept of discrepancy generalizes the equation by allowing to replace the term \eqn{(y_i - y_j)^2} with any  measure of dissimilarity \eqn{d_{ij}}.
}

\value{
The pseudo variance.}

\references{
  Studer, M., G. Ritschard, A. Gabadinho, and N. S. Müller (2009) 
  Discrepancy	analysis of complex objects using dissimilarities.
  In H. Briand, F. Guillet, G. Ritschard, and D. A. Zighed (Eds.),
  \emph{Advances in Knowledge Discovery and Management},
  Studies in Computational Intelligence. Berlin: Springer.
	
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009)
  Analyse de dissimilarités par arbre d'induction. In EGC 2009,
  \emph{Revue des Nouvelles Technologies de l'Information}, Vol. E-15, pp. 7-18.

  Batagelj, V. (1988) Generalized ward and related clustering problems. In H. Bock (Ed.),
  \emph{Classification and related methods of data analysis}, Amsterdam: North-Holland, pp. 67-74.

  Anderson, M. J. (2001) A new method for non-parametric multivariate analysis of variance.
\emph{Austral Ecology} \bold{26}, 32-46.

}
\seealso{
   \code{\link{dissassoc}} to test association between objects represented by their dissimilarities and a covariate.
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

## Pseudo variance of the sequences
print(dissvar(mvad.lcs))
}
