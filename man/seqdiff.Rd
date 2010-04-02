\encoding{latin1}
\name{seqdiff}
\alias{seqdiff}
\alias{print.seqdiff}
\alias{plot.seqdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Decompose the difference between groups of sequences}
\description{
  Decompose the difference between groups of sequences
}
\usage{
seqdiff(seqdata, group, cmprange = c(0, 1),
        seqdist_arg=list(method="LCS",norm=TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{a state sequence object created with the \code{\link{seqdef}} function.}
  \item{group}{The group variable}
  \item{cmprange}{The range used to compare subsequences}
  \item{seqdist_arg}{argument passed directly to seqdist as a list}
}
\details{
  Analyses at each timestamp the sequence discrepancy within a sliding
  time window (of range defined by \code{cmprange})
  that is explained by the \code{group} variable.
  The method computes a distance matrix, using \code{\link{seqdist}} at
  each timestamp and then derives the explained
  discrepancy with \code{\link{dissassoc}}.

  There are print and plot methods for the result returned.
}

\value{
  A \code{seqdiff} object, with the following items:

  \item{stat}{A \code{data.frame} with three statistics (PseudoF, PseudoR2 and PseudoT)
  for each timestamp of the sequence,
      see \code{\link{dissassoc}}}

  \item{variance}{A \code{data.frame} with, at each time stamp, the discrepancy within each group defined by the
      \code{group} variable and for the whole population.}
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

 }
\seealso{\code{\link{dissassoc}} to analyse the association with the whole sequence}
\examples{
## Defining a state sequence object
data(mvad)
mvad.seq <- seqdef(mvad[, 17:86])

## Building dissimilarities
mvad.diff <- seqdiff(mvad.seq, group=mvad$gcse5eq)
print(mvad.diff)
plot(mvad.diff)
plot(mvad.diff, stat="Variance")
}
