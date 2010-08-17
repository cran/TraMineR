\name{disstree}
\alias{disstree}
\alias{print.disstree}
\encoding{latin1}
\title{Dissimilarity Tree}
\description{
  Tree structured discrepancy analysis of non-measurable objects described by their pairwise dissimilarities.
}
\usage{
disstree(formula, data= NULL, minSize = 0.05, maxdepth = 5,
   R = 1000, pval = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula where the left hand side is a dissimilarity matrix and the right hand specifies the candidate partitioning variables to partition the population}
  \item{data}{a data frame where arguments in formula will be searched}
  \item{minSize}{minimum number of cases in a node, will be treated as a proportion if less than 1.}
  \item{maxdepth}{maximum depth of the tree}
  \item{R}{Number of permutations used to assess the significance of the split.}
  \item{pval}{Maximum p-value, in percent}
}
\details{
  The procedure iteratively splits the data. At each step, the procedure selects the variable and split that explains the greatest part of the discrepancy, 
  i.e. the split for which we get the highest pseudo R2. 
  The significance of the retained split is assessed through a permutation test.
}
\value{
An object of class \code{disstree} that contains the following components:
  \item{root}{A node object (see below), root of the tree}
  \item{adjustment}{A \code{\link{dissassoc}} object} providing global statistics for tree. 
  \item{formula}{The formula used to generate the tree}

Each node object contains itself the following components:
  \item{split}{Selected predictor, NULL for terminal nodes}
  \item{vardis}{Node discrepancy, see \code{\link{dissvar}}}
  \item{children}{Child nodes, NULL for terminal nodes}
  \item{ind}{Index of individuals in this node}
  \item{depth}{Depth of the node, starting from root node}
  \item{label}{Node label}
  \item{R2}{R squared of the split, NULL for terminal nodes}
}
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
  \emph{Classification and related methods of data analysis}, Amsterdam: Norht-Holland, pp. 67-74.
  
  Anderson, M. J. (2001) A new method for non-parametric multivariate analysis of variance.
\emph{Austral Ecology} \bold{26}, 32-46.
  
  Piccarreta, R. et F. C. Billari (2007) Clustering work and family trajectories by using a divisive
algorithm. \emph{Journal of the Royal Statistical Society A} \bold{170}(4), 1061--1078.
}
%\author{}
\seealso{
   \code{\link{seqtree2dot}} to generate graphic representation of \code{disstree} objects when analyzing state sequences.
   \cr
   \code{\link{disstree2dot}} is a more general interface to generate such representation.
   \cr
   \code{\link{dissvar}} to compute discrepancy using dissimilarities and for a basic introduction to discrepancy analysis.
   \cr
   \code{\link{dissassoc}} to test association between objects represented by their dissimilarities and a covariate.
   \cr
	 \code{\link{dissmfac}} to perform multi-factor analysis of variance from pairwise dissimilarities.
\cr
\code{\link{disscenter}} to compute the distance of each object to its group center from pairwise dissimilarities.

}
\examples{
data(mvad)

## Defining a state sequence object
mvad.seq <- seqdef(mvad[, 17:86])

## Computing dissimilarities
mvad.lcs <- seqdist(mvad.seq, method="LCS")
dt <- disstree(mvad.lcs~ male + Grammar + funemp + gcse5eq + fmpr + livboth, 
    data=mvad, R = 10)
print(dt)

## Using simplified interface to generate a file for GraphViz
seqtree2dot(dt, "mvadseqtree", seqdata=mvad.seq, type="d",
	border=NA, withlegend=FALSE, axes=FALSE, ylab="", yaxis=FALSE)

## Generating a file for GraphViz
disstree2dot(dt, "mvadtree", imagefunc=seqdplot, imagedata=mvad.seq, 
	## Additional parameters passed to seqdplot
	withlegend=FALSE, axes=FALSE, ylab="")
  
## Second method, using a specific function
myplotfunction <- function(individuals, seqs, mds,...) {
	par(font.sub=2, mar=c(3,0,6,0), mgp=c(0,0,0))

	## using mds to order sequence in seqiplot
	mds <- cmdscale(seqdist(seqs[individuals,], method="LCS"),k=1)
	seqiplot(seqs[individuals,], sortv=mds,...)
	}
 
## Generating a file for GraphViz
## If imagedata is not set, index of individuals are sent to imagefunc
disstree2dot(dt, "mvadtree", imagefunc=myplotfunction, title.cex=3,
	## additional parameters passed to myplotfunction
	seqs=mvad.seq, mds=mvad.mds,
	## additional parameters passed to seqiplot (through myplotfunction)
	withlegend=FALSE, axes=FALSE,tlim=0,space=0, ylab="", border=NA)

## To run GraphViz (dot) from R and generate an "svg" file
## shell("dot -Tsvg -O mvadtree.dot")
}
