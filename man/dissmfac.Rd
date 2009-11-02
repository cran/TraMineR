\name{dissmfac}
\alias{dissmfac}
\alias{dissreg}
\alias{gower_matrix}
\alias{print.dissregression}
\alias{print.dissmultifactor}
\encoding{latin1}
\title{Multi-factor ANOVA from a dissimilarity matrix}
\description{
  Perform a multi-factor analysis of variance from a dissimilarity matrix.
}
\usage{
dissmfac(formula, data, R = 1000, gower = FALSE, squared = TRUE,
 permutation = "dissmatrix")
}

\arguments{
  \item{formula}{A regression-like formula. The left hand side should be a dissimilarity matrix or a dist object.}
  \item{data}{data to search for variables in formula}
  \item{R}{Number of permutations to assess significance}
  \item{gower}{Logical: Is the dissimilarity matrix already a Gower matrix?}
  \item{squared}{Logical: should we square the dissimilarity matrix?}
  \item{permutation}{if equal to \code{dissmatrix}, permutations are done on the dissimilarity matrix, else if equal to "model" permutations are done on the variable matrix. Depending on the number of observation, "model" can be quicker.}
}
\details{
  This method is, in some way, a generalization of \code{\link{dissassoc}} that can account for several explanatory variables. This function compute the part of variance explained by a list of covariates using a decomposition of the discrepancy (variance) explained. This function is slower than \code{\link{dissassoc}} for one factor. More on that, the latter also perform a test of discrepancy homogeneity (equality of variance) using a generalization of the T statistic.

  The function is based on the program written for scipy (Python) by Ondrej Libiger and Matt Zapala. See Zapala and Schork (2006) for a full reference.
}

\value{
  A \code{dissmultifactor} object with the following components:
  \item{mfac}{The part of variance explained by each variable (comparing full model to model without the specified variable) and its significance using permutation test}
  \item{call}{Function call}
  \item{perms}{Permutation values as a \code{boot} object}
  \item{perm_method}{Permutation method used to compute significance}
}
\references{
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009)
  Discrepancy	analysis of complex objects using dissimilarities.
  In H. Briand, F. Guillet, G. Ritschard, and D. A. Zighed (Eds.),
  \emph{Advances in Knowledge Discovery and Management},
  Studies in Computational Intelligence. Berlin: Springer.
	
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009).
  Analyse de dissimilarités par arbre d'induction. In EGC 2009,
  \emph{Revue des Nouvelles Technologies de l'Information}, Vol. E-15, pp. 7-18.

  Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance.
\emph{Austral Ecology} 26, 32-46.

  McArdle, B. H. et M. J. Anderson (2001). Fitting multivariate models to community data: A
comment on distance-based redundancy analysis. \emph{Ecology} 82(1), 290-297.

  Zapala, M. A. et N. J. Schork (2006). Multivariate regression analysis of distance matrices for
  testing associations between gene expression patterns and related variables. \emph{Proceedings of
the National Academy of Sciences of the United States of America} 103(51), 19430-19435.
}
\seealso{
\code{\link{dissvar}} to compute the pseudo variance from dissimilarities and for a basic introduction to concepts of pseudo variance analysis.
   \cr
   \code{\link{dissassoc}} to test association between objects represented by their dissimilarities and a covariate.
\cr
\code{\link{disstree}} for an induction tree analyse of objects characterized by a dissimilarity matrix.
\cr
\code{\link{disscenter}} to compute the distance of each object to its group center from pairwise dissimilarities.
}
\examples{
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])

## Building dissimilarities
mvad.lcs <- seqdist(mvad.seq, method="LCS")
print(dissmfac(mvad.lcs ~ male + Grammar + funemp +
	gcse5eq + fmpr + livboth, data=mvad, R=10))
}
