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
The discrepancy \deqn{s^2} is defined as:
\deqn{s^2=\frac{1}{2n^2}\sum_{i=1}^{n}\sum_{j=1}^{n}d_{ij}}{s^2=(1/2n^2)* sum sum d_ij}

\emph{Mathematical ground}:
In the Euclidean case, the sum of square can be rewritten as:
  \deqn{SS=\sum_{i=1}^{n}(y_i-\bar{y})^2=\frac{1}{2n}\sum_{i=1}^{n}\sum_{j=1}^{n}(y_i-y_j)^2}{SS= sum (y_i-y_mean)^2=(1/2n)*sum sum (y_i-y_j)^2}

  The concept of discrepancy generalize the equation by replacing the term \deqn{(y_i-y_j)^2} by another measure of dissimilarities.

}
\value{
The pseudo variance.}
\references{

  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009).
  Analyse de dissimilarités par arbre d'induction.
  \emph{Revue des Nouvelles Technologies de l'Information}, EGC'2009.

  Batagelj, V. (1988). Generalized ward and related clustering problems. In H. Bock (Ed.),
  \emph{Classification and related methods of data analysis}, pp. 67-74. North-Holland, Amsterdam.

  Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance.
\emph{Austral Ecology} 26, 32-46.

}
\seealso{
   \code{\link{dissassoc}} to test association of objects characterized by their dissimilarities with some explanatory variable

   \code{\link{disstree}} to analyse dissimilarities with induction trees

   \code{\link{dissmfac}} to analyse dissimilarities in a way close to linear regression

   \code{\link{disscenter}} to compute the distance of each object to its group center from a dissimilarity matrix

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
