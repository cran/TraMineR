\name{disstree}
\alias{disstree}
\alias{print.disstree}
\encoding{latin1}
\title{Dissimilarity Tree}
\description{
  Analyse non-measurable objects described through a set of dissimilarity by recursively partionning the population.
}
\usage{
disstree(formula, data= NULL, minSize = 0.05,maxdepth = 5,
   R = 1000, pval = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ A formula where de left hand side is a dissimilarity matrix, the right hand side should be a list of candidate variable to partion the population}
  \item{data}{a data.frame where arguments in formula can be identified}
  \item{minSize}{minimum number of observation in a node, in percentage if less than 1.}
  \item{maxdepth}{maximum depth of the tree}
  \item{R}{Number of permutation used to assess significativity of a partition.}
  \item{pval}{Maximum p-value, in percent}
}
\details{
  At each step, this procedure choose the variable that explains the biggest part of the pseudo variance to partition the population. It assess the significance of the choosen variable by performing a permutation test.
}
\value{
Return an object of class \code{disstree}, a list with the following component:
  \item{node}{A tree object (see below)}
  \item{adjustement}{global adjustement of the tree}
The node object is a list with the following component:
  \item{split}{Choosen predictor, NULL  for terminal nodes}
  \item{vardis}{Node pseudo variance, see \code{\link{dissvar}}}
  \item{children}{Child node, NULL for terminal nodes}
  \item{ind}{Index of individuals in this node}
  \item{depth}{Depth of the node, starting from root node}
  \item{label}{Label of this node}
  \item{R2}{R squared of the split, NULL for terminal nodes}
}
\references{
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009).
  Analyse de dissimilarités par arbre d'induction.
  \emph{Revue des Nouvelles Technologies de l'Information}, EGC'2009.
  
  Batagelj, V. (1988). Generalized ward and related clustering problems. In H. Bock (Ed.),
  \emph{Classification and related methods of data analysis}, pp. 67-74. North-Holland, Amsterdam.
  
  Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance.
\emph{Austral Ecology} 26, 32-46.
  
  Piccarreta, R. et F. C. Billari (2007). Clustering work and family trajectories by using a divisive
algorithm. \emph{Journal of the Royal Statistical Society A} 170(4), 1061-1078.
}
%\author{}
\seealso{
   \code{\link{dissvar}} to compute pseudo variance using dissimilarities and for a basic introduction to concepts of pseudo variance analysis
   
   \code{\link{dissassoc}} to test association between dissimilarity and another variable
   
   \code{\link{dissreg}} to analyse dissimilarities in a way close to linear regression
   
   \code{\link{disscenter}} to compute the distance of each object to its center of group using dissimilarities

}
\examples{
data(mvad)

## Defining a state sequence object
mvad.seq <- seqdef(mvad[, 17:86])

## Building dissimilarities
mvad.lcs <- seqdist(mvad.seq, method="LCS")
dt <- disstree(mvad.lcs~ male + Grammar + funemp + gcse5eq + fmpr + livboth, 
    data=mvad, R = 10000)
print(dt)

## Compute quality of the tree
print(dissassoc(mvad.lcs, disstreeleaf(dt), R=1))
  
## Using simplified interface to generate a file for GraphViz
seqtree2dot(dt, "mvadseqtree", seqs=mvad.seq, plottype="seqdplot", 
	border=NA, withlegend=FALSE)

## Generating a file for GraphViz
disstree2dot(dt, "mvadtree", imagefunc=seqdplot, imagedata=mvad.seq, 
	## Additionnal parameters passed to seqdplot
	withlegend=FALSE, axes=FALSE)
  
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
	## additionnal parameters passed to myplotfunction
	seqs=mvad.seq, mds=mvad.mds,
	## additionnal parameters passed to seqiplot (through myplotfunction)
	withlegend=FALSE, axes=FALSE,tlim=0,space=0, ylab="", border=NA)

## To run GraphViz (dot) from R
## shell("dot -Tsvg -O mvadtree.dot")
}
