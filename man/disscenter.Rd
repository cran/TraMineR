\name{disscenter}
\alias{disscenter}
\encoding{latin1}
\title{Compute distance to center of a group}
\description{
 Compute the dissimilarity between each objects and group centers according to a matrix of dissimilarity.
}
\usage{
disscenter(diss, group=NULL, medoids.index=NULL)
}
\arguments{
  \item{diss}{A matrix of dissimilarity such as the one computed by \code{\link{seqdist}}, or a \code{dist} object (see \code{\link{dist}}}
  \item{group}{If null, only one group is considered, otherwise group to compute center}
%  \item{trim}{Percentage of individuals that are considered as outliers and included in computation of distance to center}
  \item{medoids.index}{If NULL, return dissimilarity to center. If equal to "first", return the index of the first encountered most central sequence. One index per group is returned. If equal to "all", all medoids index are returned. If \code{group} is set, one list per group is returned.}
%  \item{max.iter}{Maximum number of iteration to trim the outliers}
}
\details{
  This method compute the dissimilarity between a given objects and its group center. The group center may not belong to the space formed by the objects (in the same way, the average do not belong to a space formed by discrete measure).
  This distance can also be understood as the contribution to the discrepancy (see \code{\link{dissvar}}).
  The dissimilarity between a given objects and its group center may be negative if the dissimilarity measure do not respect the triangle inequality.
  
  It can be shown that this dissimilarity is equal to \cite{Batagelj, 1988}:
    \deqn{d_{x\tilde{g}}=\frac{1}{n}\big(\sum_{i=1}^{n}d_{xi}-SS\big)}{d_(xg)=1/n *(sum d_xi - SS)}
  Where SS is the sum of square (see \code{\link{dissvar}}).
}
\value{
  A vector with the dissimilarity to center of group for each sequence, or a list of medoid index
}

\references{

  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009).
  Analyse de dissimilarités par arbre d'induction.
  \emph{Revue des Nouvelles Technologies de l'Information}, EGC'2009.
  
  Batagelj, V. (1988). Generalized ward and related clustering problems. In H. Bock (Ed.),
  \emph{Classification and related methods of data analysis}, pp. 67-74. North-Holland, Amsterdam.
  
}

\seealso{
   \code{\link{dissvar}} to compute discrepancy using dissimilarities and for a basic introduction to concepts of discrepancy analysis
   
   \code{\link{dissassoc}} to test association between dissimilarity and another variable
   
   \code{\link{disstree}} to analyse dissimilarities using induction trees
   
   \code{\link{dissmfac}} to perform multi-factor analysis of variance using dissimilarities
}
\examples{
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])

## Building dissimilarities
mvad.lcs <- seqdist(mvad.seq, method="LCS")

## Compute distance to center according to group gcse5eq
dc <- disscenter(mvad.lcs, group=mvad$gcse5eq)

## Ploting distribution of dissimilarity  to center 
boxplot(dc~mvad$gcse5eq, col="cyan")

## Retrieving index of the first medoids, one per group
dc <- disscenter(mvad.lcs, group=mvad$Grammar, medoids.index="first")
print(dc)

## Retrieving index of all medoids in each group
dc <- disscenter(mvad.lcs, group=mvad$Grammar, medoids.index="all")
print(dc)
}

