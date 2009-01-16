\name{dissreg}
\alias{dissreg}
\alias{gower_matrix}
\alias{print.dissregression}
\encoding{latin1}
\title{Regression analysis of dissimilarity matrix}
\description{
  Perform a regression analysis of a dissimilarity matrix. 
}
\usage{
dissreg(formula, data, R = 1000, gower = FALSE, squared = TRUE,
 permutation = "dissmatrix")
}

\arguments{
  \item{formula}{A formula for the regression. The left hand side should be a dissimilarity matrix or a dist object.}
  \item{data}{data to search for variable in formula}
  \item{R}{Number of permutation to assess significance}
  \item{gower}{Is the dissimilarity matrix already a gower matrix ?}
  \item{squared}{should we square the dissimilarity matrix ?}
  \item{permutation}{if equal to \code{dissmatrix}, permutations are done on the dissimilarity matrix, else if equal to "model" permutation are done on variable matrix}
}
\details{
  This method is, in some way, a generalization of \code{dissassoc} in order to several explicative variables.
  
  The function is based on the program written for scipy (Python) by Ondrej Libiger and Matt Zapala. See references below for a full reference.
  
  This function compute the part of variance explained by a list of covariates using a decomposition of the variance explained.  
}
\value{
  Return a \code{dissregression} object with the following componant:
  \item{mreg}{The part of variance explained by each coefficient (comparing full model to model without the specified variable) and its significativity using permutation test}
  \item{call}{Function call}
  \item{perms}{Permutations values as a \code{boot} object}
  \item{perm_method}{Permutation method used to compute significance}
}
\references{
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009).
  Analyse de dissimilarités par arbre d'induction.
  \emph{Revue des Nouvelles Technologies de l'Information}, EGC'2009.
  
  Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance.
\emph{Austral Ecology} 26, 32-46.
  
  McArdle, B. H. et M. J. Anderson (2001). Fitting multivariate models to community data: A
comment on distance-based redundancy analysis. \emph{Ecology} 82(1), 290-297.

  Zapala, M. A. et N. J. Schork (2006). Multivariate regression analysis of distance matrices for
  testing associations between gene expression patterns and related variables. \emph{Proceedings of
the National Academy of Sciences of the United States of America} 103(51), 19430-19435.
}
\seealso{
\code{\link{dissvar}} to compute pseudo variance using dissimilarities and for a basic introduction to concepts of pseudo variance analysis

\code{\link{dissassoc}} to test association between dissimilarity and another variable

\code{\link{disstree}} to analyse dissimilarities using induction trees

\code{\link{disscenter}} to compute the distance of each object to its center of group using dissimilarities

}
\examples{
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])

## Building dissimilarities
mvad.lcs <- seqdist(mvad.seq, method="LCS")
print(dissreg(mvad.lcs ~ male + Grammar + funemp + 
	gcse5eq + fmpr + livboth,  data=mvad, R=10))
}
