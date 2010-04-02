\encoding{latin1}
\name{dissrep}
\alias{dissrep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extracting sets of representative objects using a dissimilarity matrix}
\description{
The function extracts a set of representative objects that exhibits the key features of the whole data set, the goal being to get easy sounded interpretation of the latter. The user can set either the desired coverage level (the proportion of objects having a representative in their neighborhood) or the desired number of representatives.
}
\usage{
dissrep(diss, criterion="density",
    score=NULL, decreasing=TRUE,
    trep=0.25, nrep=NULL, tsim=0.1, dmax=NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
 \item{diss}{A dissimilarity matrix or a \code{dist} object (see \code{\link{dist}})}

  \item{criterion}{the representativeness criterion for sorting the
candidate list. One of \code{"freq"} (frequency), \code{"density"}
(neighborhood density) or \code{"dist"} (centrality). An optional
vector containing the scores for sorting the candidate objects may
also be provided. See below and details.}

  \item{score}{an optional vector containing the representativeness
scores used for sorting the objects in the candidate list. The length
of the vector must be equal to the number of rows/columns in the
distance matrix, i.e the number of objects.}

  \item{decreasing}{if a score vector is provided, indicates whether
the objects in the candidate list must be sorted in ascending or
decreasing order of this score. The first object in the candidate
list is supposed to be the most representative.}

  \item{trep}{controls the size of the representative set by setting
the desired coverage level, i.e the proportion of objects having a
representative in their neighborhood. Neighborhood radius is defined
by \code{tsim}.}

  \item{nrep}{number of representatives. If \code{NULL} (default),
\code{trep} argument is used to control the size of the
representative set.}

 \item{tsim}{threshold for setting the redundancy and neighborhood
radius. Defined as a percentage of the maximum (theoretical)
distance. Defaults to 0.1 (10\%). Object $y$ is considered as
redundant to/in the neighborhood of object $x$ if the distance from
$y$ to $x$ is less than \code{tsim}*\code{dmax}. The neighborhood
diameter is thus twice this threshold. }

  \item{dmax}{maximum theoretical distance. Redundancy and neighborhood
diameters are defined as a proportion of this maximum theoretical
distance. If \code{NULL}, it is derived from the distance matrix.}
 }

\details{
 The representative set is obtained by an heuristic that
first builds a sorted list of candidates using a representativeness
score and then eliminates redundancy. The available criterions for
sorting the candidate list are: \emph{sequence frequency},
\emph{neighborhood density}, \emph{centrality}. Other user defined
sorting criterions can be provided using the \code{score} argument.

The \emph{frequency} criterion uses the frequencies as
representativeness score. The frequency of an object in the data is
computed as the number of other objects with whom the dissimilarity
is equal to 0. The more frequent an object the more representative it
is supposed to be. Hence, objects are sorted in decreasing frequency
order. Indeed, this criterion is the neighborhood (see below)
criterion with the neighborhood diameter set to 0.

The \emph{neighborhood density} criterion uses the
number---density---of objects in the neighborhood of each candidate.
This requires indeed to set the neighborhood diameter. We suggest to
set it as a given proportion of the maximal (theoretical) distance
between two objects. Candidates are sorted in decreasing density
order.

The \emph{centrality} criterion uses the sum of distances to all
other objects, i.e. the centrality as a representativeness criterion.
The smallest the sum, the most representative the candidate.

For more details, see \cite{Gabadinho et al., 2009}.
}

\value{An object of class \code{diss.rep}. This is a vector containing
the indexes of the representative objects with the following additional attributes:
  \item{Scores}{a vector with the representative score of each object given the chosen criterion.}
  \item{Distances}{a matrix with the distance of each object to its nearest representative.}
  \item{Statistics}{contains several quality measures for each representative
in the set: number of objects attributed to the representative,
number of object in the representatives neighborhood, mean distance
to the representative.}
  \item{Quality}{overall quality measure.}

Print and summary methods are available.
}

\references{
 Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009).
Summarizing Sets of Categorical Sequences, In International
Conference on Knowledge Discovery and Information Retrieval, Madeira,
6-8 October, INSTICC.}

\seealso{\code{\link{seqrep}}, \code{\link{plot.stslist.rep}}}

\examples{
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)

## Computing the distance matrix
costs <- seqsubm(biofam.seq, method="TRATE")
biofam.om <- seqdist(biofam.seq, method="OM", sm=costs)

## Representative set using the neighborhood density criterion
biofam.rep <- dissrep(biofam.om)
biofam.rep
summary(biofam.rep)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
