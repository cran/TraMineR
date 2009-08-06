\encoding{latin1}
\name{seqrep}
\alias{seqrep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extracting sets of representative sequences}
\description{
The function attempts to find an optimal (as small as possible while assuring a large coverage) set of representative sequences that exhibits the key features of the whole sequence data set, the goal being to get easy sounded interpretation of the latter.
}
\usage{
seqrep(seqdata, criterion = "freq", trep = 0.25, nrep = NULL, 
  tsim = 0.1, dmax = NULL, dist.matrix = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{a state sequence object as defined by the \code{\link{seqdef}} function.}
  \item{criterion}{the representativeness criterion for sorting the candidate list. One of \code{"freq"} (sequence frequency), \code{"density"} (neighborhood density), \code{"mscore"} (mean state frequency), \code{"dist"} (centrality) and \code{"prob"} (sequence likelihood). See details.}

  \item{trep}{parameter for controlling the size of the candidate list and thus the size of the representative set. The size of the candidate list is limited so that the cumulated frequency of the retained distinct candidates reaches a threshold proportion \code{trep} of the whole data set. Setting for instance \eqn{trep=25\%} ensures that at least \eqn{25\%} of the sequences will have a representative in their neighborhood and that the final representative set will contain at most \eqn{25\%} of the distinct sequences of the whole set. Thus \code{trep} defines also a minimum coverage level.}

 \item{nrep}{number of representative sequences. If \code{NULL} (default), trep is used to control the size of the representative set.}

  \item{tsim}{threshold for setting the redundancy and neighborhood limits. Defined as a percentage of the maximum theoretical distance. Defaults to 0.1 (10\%).}

  \item{dmax}{maximum theoretical distance. Redundancy and neighborhood diameters are defined as a proportion of this maximum theoretical distance. If \code{NULL}, it is derived from the distance matrix.}

  \item{dist.matrix}{a matrix containing the pairwise distances between sequences in \code{seqdata}. If \code{NULL}, the matrix is computed by calling the \code{\link{seqdist}} function. In that case, optional arguments to be passed to the \code{seqdist} function (see ... hereafter) should also be provided.}

  \item{...}{optional arguments to be passed to the \code{seqdist} function, mainly \code{dist.method} specifying the metric for computing the distance matrix, \code{norm} for normalizing the distances, \code{indel} and \code{sm} for indel and substitution costs when Optimal Matching metric is chosen. See \code{\link{seqdist}} manual page for details.}
}
\details{
The representative set is obtained by an heuristic that first builds a sorted list of candidates using a representativeness score and then eliminates redundancy. The available criterions for sorting the candidate list are: \emph{sequence frequency}, \emph{neighborhood density}, \emph{mean state frequency}, \emph{centrality} and \emph{sequence likelihood}.

The \emph{sequence frequency} criterion uses the sequence frequencies as representativeness score. The more frequent a sequence the more representative it is supposed to be. Hence, sequences are sorted in decreasing frequency order.

The \emph{neighborhood density} criterion uses the number --- the density --- of sequences in the neighborhood of each candidate sequence. This requires indeed to set the neighborhood diameter. We suggest to set it as a given proportion of the maximal theoretical distance between two sequences. Sequences are sorted in decreasing density order.

The \emph{mean state frequency} criterion is the mean value of the transversal frequencies of the successive states.
Let \eqn{s=s_{1}s_{2}\cdots s_{\ell}}{s=(s_1, s_2, ..., s_l)} be a sequence of length \eqn{\ell}{l} and \eqn{(f_{s_1}, f_{s_2}, \ldots, f_{s_\ell})}{f(s_1), f(s_2), ..., f(s_l)} the frequencies of the states at (time-)position \eqn{(t_1, t_2,\ldots t_{\ell})}{t_1, t_2, ..., t_l}. The mean state frequency is the sum of the state frequencies divided by the sequence length
\deqn{
MSF(s)=\frac{1}{\ell} \sum_{i=1}^{\ell} f_{s_{i}}
}{
	MSF(s)=1/l sum f(s_i)
}
%
The lower and upper boundaries of \eqn{MSF} are \eqn{0} and \eqn{1}. \eqn{MSF} is equal to \eqn{1} when all the sequences in the set are the same, i.e. when there is a single distinct sequence. The most representative sequence is the one with the highest score.

The \emph{centrality} criterion uses the sum of distances to all other sequences, i.e. the centrality as a representativeness criterion. The smallest the sum, the most representative the sequence.

The \emph{sequence likelihood} \eqn{P(s)} is defined as the product of the probability with which each of its observed successive state is supposed to occur at its position.
Let  \eqn{s=s_{1}s_{2} \cdots s_{\ell}}{s_1, s_2, s_l}  be a sequence of length \eqn{\ell}{l}. Then
\deqn{
P(s)=P(s_{1},1) \cdot P(s_{2},2) \cdots P(s_{\ell},\ell)
}{
P(s)=P(s_1,1) * P(s_2,2) * ... * P(s_l,l)
}
with \eqn{P(s_{t},t)}{P(s_t,t)} the probability to observe state \eqn{s_t} at position \eqn{t}.
\cr
The question is how to determinate the state probabilities \eqn{P(s_{t},t)}{P(s_t,t)}. One commonly used method for computing them is to postulate a Markov model, which can be of various order. The implemented criterion considers the probabilities derived from the first order Markov model, that is each \eqn{P(s_{t},t)}{P(s_t,t)}, \eqn{t>1} is set to the transition rate \eqn{p(s_t|s_{t-1})}{p(s_t)|s_t-1} estimated across sequences from the observations at positions \eqn{t} and \eqn{t-1}. For \eqn{t=1}, we set \eqn{P(s_1,1)} to the observed frequency of the state \eqn{s_1} at position 1.

The likelihood \eqn{P(s)} being generally very small, we use \eqn{-\log P(s)} as sorting criterion. The latter quantity is minimal when \eqn{P(s)} is equal to 1, which leads to sort the sequences in ascending order of their score.

For more details, see \cite{Gabadinho et al., 2009}.
}
\value{An object of class \code{stslist.rep}. This is actually a state sequence object (containing a list of state sequences) with the following additional attributes:
  \item{Scores}{a vector with the representative score of each sequence in the original set given the chosen criterion.}
  \item{Distances}{a matrix with the distance of each sequence to its nearest representative.}
  \item{Quality}{contains several quality measures for each representative sequence in the set: number of sequences attributed to the representative, number of sequence in the representatives neighborhood, mean distance to the representative.}
  \item{rindex}{overall quality measure.}

Print, plot and summary methods are available. More elaborated plots are produced by the \code{seqplot} function using the \code{type="r"} argument, or the \code{seqrplot} alias.
}

\references{Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009). Summarizing Sets of Categorical Sequences, In International Conference on Knowledge Discovery and Information Retrieval, Madeira, 6-8 October, INSTICC.}

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

## Representative set using the frequency criterion
biofam.rep <- seqrep(biofam.seq, dist.matrix=biofam.om)
biofam.rep
summary(biofam.rep)

## Representative set using the neighborhood density criterion
seqrep(biofam.seq, dist.matrix=biofam.om, criterion="density")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
