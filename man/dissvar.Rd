\name{dissvar}
\alias{dissvar}
\encoding{latin1}
\title{Dissimilarity based pseudo variance}
\description{
  Compute the pseudo variance based on the dissimilarity between object.  
}
\usage{
dissvar(diss)
}
\arguments{
  \item{diss}{A dissimilarity matrix or a dist object (see \code{\link{dist}})}
}
\details{
In the euclidian case, the sum of square can be rewritten as:
  \deqn{SS=\sum_{i=1}^{n}(y_i-\bar{y})^2=\frac{1}{2n}\sum_{i=1}^{n}\sum_{j=1}^{n}(y_i-y_j)^2}{SS= sum (y_i-y_mean)^2=(1/2n)*sum sum (y_i-y_j)^2}
  This function generalized this concept by allowing to use any kind of dissimilarity to compute a pseudo variance term.
}
\value{
Return the pseudo variance.}
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
   \code{\link{dissassoc}} to test association between dissimilarity and another variable
   
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

## Pseudo variance of the sequences
print(dissvar(mvad.lcs))
}
