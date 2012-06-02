\encoding{latin1}
\name{seqrep}
\alias{seqrep}
\title{Extracting sets of representative sequences}

\description{ The function attempts to find an optimal (as small as
possible for a required coverage) set of representative
sequences that exhibits the key features of the whole sequence data
set, the goal being to easily get sounded interpretation of the set of sequences.
}

\usage{ seqrep(seqdata, criterion="density", score=NULL,
    decreasing=TRUE, trep=0.25, nrep=NULL,
    tsim=0.1, dmax=NULL, dist.matrix=NULL, weighted=TRUE, ...)
}

\arguments{
  \item{seqdata}{a state sequence object as defined by the \code{\link{seqdef}} function.}

  \item{criterion}{the representativeness criterion for sorting the candidate list. One of \code{"freq"} (sequence
      frequency), \code{"density"} (neighborhood density), \code{"mscore"} (mean state frequency), \code{"dist"}
      (centrality) and \code{"prob"} (sequence likelihood). See details.}

  \item{score}{an optional vector containing the representativeness scores used to sort the sequences in the candidate
      list. The length of the vector must be equal to the number of sequences in the sequence object.}

  \item{decreasing}{if a score vector is provided, indicates whether the objects in the candidate list must be sorted in
      ascending or descending order of this score. Default is TRUE, i.e. descending. The first object in the candidate list
      is then supposed to be the most representative. }

  \item{trep}{coverage threshold, i.e. minimum proportion of sequences that should have a representative in their
      neighborhood (neighborhood diameter is defined by \code{tsim}).}

 \item{nrep}{number of representative sequences. If \code{NULL} (default), the size of the representative set is
     controlled by \code{trep}.}

  \item{tsim}{threshold for setting the redundancy and neighborhood
radius. Defined as a percentage of the maximum (theoretical)
distance. Defaults to 0.1 (10\%). Sequence \eqn{y} is considered as
redundant to/in the neighborhood of sequence \eqn{x} if the distance
from \eqn{y} to \eqn{x} is less than \code{tsim}*\code{dmax}. The
neighborhood diameter is thus twice this threshold.}

  \item{dmax}{maximum theoretical distance. The neighborhood diameter is defined as a proportion of this
      maximum theoretical distance. If \code{NULL}, it is derived from the distance matrix.}

  \item{dist.matrix}{a matrix containing the pairwise distances between sequences in \code{seqdata}. If \code{NULL}, the
      matrix is computed by calling the \code{\link{seqdist}} function. In that case, optional arguments to be passed to
      the \code{seqdist} function (see ... hereafter) should also be provided.}

  \item{weighted}{logical: Should weights assigned to the state sequence object be accounted for? (See \code{\link{seqdef}}.)
  Set as \code{FALSE} to ignore the weights.}

  \item{...}{optional arguments to be passed to the \code{seqdist} function, mainly \code{dist.method} specifying the
      metric for computing the distance matrix, \code{norm} for normalizing the distances, \code{indel} and \code{sm} for
      indel and substitution costs when Optimal Matching metric is chosen. See \code{\link{seqdist}} manual page for
      details.}
}

\details{
 The representative set is obtained by an heuristic that
first builds a sorted list of candidates using a representativeness
score and then eliminates redundancy. The available criterions for
sorting the candidate list are: \emph{sequence frequency},
\emph{neighborhood density}, \emph{mean state frequency},
\emph{centrality} and \emph{sequence likelihood}.

The \emph{sequence frequency} criterion uses the sequence frequencies as representativeness score. The more frequent a
sequence the more representative it is supposed to be. Hence, sequences are sorted in decreasing frequency order.

The \emph{neighborhood density} criterion uses the
number---density---of sequences in the neighborhood of each candidate
sequence. This requires indeed to set the neighborhood diameter
\code{tsim}. We suggest to set it as a given proportion of the
maximal theoretical distance between two sequences. Sequences are
sorted in decreasing density order.

The \emph{mean state frequency} criterion is the mean value of the transversal frequencies of the successive states.
Let \eqn{s=s_{1}s_{2}\cdots s_{\ell}}{s=(s_1, s_2, ..., s_l)} be a sequence of length \eqn{\ell}{l} and \eqn{(f_{s_1},
f_{s_2}, \ldots, f_{s_\ell})}{f(s_1), f(s_2), ..., f(s_l)} the frequencies of the states at (time-)position \eqn{(t_1,
t_2,\ldots t_{\ell})}{t_1, t_2, ..., t_l}. The mean state frequency is the sum of the state frequencies divided by the
sequence length
\deqn{
MSF(s)=\frac{1}{\ell} \sum_{i=1}^{\ell} f_{s_{i}}
}{
	MSF(s)=1/l sum f(s_i)
}
%
The lower and upper boundaries of \eqn{MSF} are \eqn{0} and \eqn{1}. \eqn{MSF} is equal to \eqn{1} when all the sequences
in the set are the same, i.e. when there is a single distinct sequence. The most representative sequence is the one with
the highest score.

The \emph{centrality} criterion uses the sum of distances to all other sequences as a representativeness criterion. The
smallest the sum, the most representative the sequence.

The \emph{sequence likelihood} \eqn{P(s)} is defined as the product of the probability with which each of its observed
successive state is supposed to occur at its position.
Let  \eqn{s=s_{1}s_{2} \cdots s_{\ell}}{s_1, s_2, s_l}  be a sequence of length \eqn{\ell}{l}. Then
\deqn{
P(s)=P(s_{1},1) \cdot P(s_{2},2) \cdots P(s_{\ell},\ell)
}{
P(s)=P(s_1,1) * P(s_2,2) * ... * P(s_l,l)
}
with \eqn{P(s_{t},t)}{P(s_t,t)} the probability to observe state \eqn{s_t} at position \eqn{t}.
\cr
The question is how to determinate the state probabilities \eqn{P(s_{t},t)}{P(s_t,t)}. One commonly used method for
computing them is to postulate a Markov model, which can be of various order. The implemented criterion considers the
probabilities derived from the first order Markov model, that is each \eqn{P(s_{t},t)}{P(s_t,t)}, \eqn{t>1} is set to the
transition rate \eqn{p(s_t|s_{t-1})}{p(s_t)|s_t-1} estimated across sequences from the observations at positions \eqn{t}
and \eqn{t-1}. For \eqn{t=1}, we set \eqn{P(s_1,1)} to the observed frequency of the state \eqn{s_1} at position 1.

The likelihood \eqn{P(s)} being generally very small, we use
\eqn{-\log P(s)} as sorting criterion. The latter quantity is minimal
when \eqn{P(s)} is equal to 1, which leads to sort the sequences in
ascending order of their score.

For more details, see \cite{Gabadinho et al., 2009}.
}

\value{An object of class \code{stslist.rep}. This is actually a
state sequence object (containing a list of state sequences) with the
following additional attributes:
  \item{Scores}{a vector with the representative score of each sequence in the original set given the chosen criterion.}
  \item{Distances}{a matrix with the distance of each sequence to its nearest representative.}
  \item{Statistics}{contains several quality measures for each representative sequence in the set: number of sequences
      attributed to the representative, number of sequence in the representatives neighborhood, mean distance to the
      representative.}
  \item{Quality}{overall quality measure.}

Print,plot and summary methods are available. More elaborated plots
are produced by the \code{seqplot} function using the \code{type="r"}
argument, or the \code{seqrplot} alias. }

\references{Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009). Summarizing Sets of Categorical Sequences, In
International Conference on Knowledge Discovery and Information Retrieval, Madeira, 6-8 October, INSTICC.}

\seealso{\code{\link{seqplot}}, \code{\link{plot.stslist.rep}}}

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
biofam.rep <- seqrep(biofam.seq, dist.matrix=biofam.om, criterion="density")
biofam.rep
summary(biofam.rep)
plot(biofam.rep)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
