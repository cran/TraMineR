\name{seqdiff}
\alias{seqdiff}
\alias{print.seqdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Position-wise discrepancy analysis between groups of sequences}
\description{
  The function analyses how the differences between groups of sequences evolve along the positions. It runs a sequence of  discrepancy analyses on sliding windows.
}
\usage{
seqdiff(seqdata, group, cmprange = c(0, 1),
  seqdist.args = list(method = "LCS", norm = TRUE), with.missing = FALSE,
  weighted = TRUE, squared = FALSE, seqdist_arg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqdata}{a state sequence object created with the \code{\link{seqdef}} function.}
  \item{group}{The group variable.}
  \item{cmprange}{The time range of the sliding window on which subsequences are compared.}
  \item{seqdist.args}{List of arguments passed to \code{\link{seqdist}} for computing the distances.}
  \item{with.missing}{Logical. If \code{TRUE}, missing values are considered as an additional state. If \code{FALSE} subsequences with missing values are removed from the analysis.}
  \item{weighted}{Logical. If \code{TRUE}, \code{seqdiff} uses the weights specified in \code{seqdata}.}
  \item{squared}{Logical. If \code{TRUE} the dissimilarities are squared for computing the discrepancy.}
  \item{seqdist_arg}{Deprecated. Use \code{seqdist.args} instead.}
}
\details{
  The function analyses how the part of discrepancy explained by the \code{group} variable evolves along the position axis. It runs successively discrepancy analyses within a sliding time-window of range \code{cmprange}). At
  each position, the method uses \code{\link{seqdist}} to compute a distance matrix over the time-window and then derives the explained discrepancy on that window with \code{\link{dissassoc}}.

  There are print and plot methods for the returned value.}

\value{
  A \code{seqdiff} object, with the following items:

  \item{stat}{A \code{data.frame} with three statistics (PseudoF, PseudoR2
  and PseudoT) for each time stamp of the sequence,
      see \code{\link{dissassoc}}}

  \item{discrepancy}{A \code{data.frame} with, at each time stamp, the discrepancy within each group defined by the \code{group} variable and for the whole population.}
}

\references{
  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2010)
  Discrepancy analysis of complex objects using dissimilarities.
  In F. Guillet, G. Ritschard, D. A. Zighed and H. Briand (Eds.),
  \emph{Advances in Knowledge Discovery and Management},
  Studies in Computational Intelligence, Volume 292, pp. 3-19. Berlin: Springer.

  Studer, M., G. Ritschard, A. Gabadinho and N. S. Müller (2009)
  Analyse de dissimilarités par arbre d'induction. In EGC 2009,
  \emph{Revue des Nouvelles Technologies de l'Information}, Vol. E-15, pp. 7-18.

 }
\seealso{\code{\link{dissassoc}} to analyse the association of the \code{group} variable with the whole sequence}
\examples{
## Define a state sequence object
data(mvad)
## First 12 months of first 200 trajectories
mvad.seq <- seqdef(mvad[1:200, 17:28])

## Position-wise discrepancy analysis
mvad.diff <- seqdiff(mvad.seq, group=mvad$gcse5eq[1:200])
print(mvad.diff)
plot(mvad.diff, stat=c("Pseudo R2", "Levene"), xtstep=6)
plot(mvad.diff, stat="discrepancy")
}

\author{Matthias Studer (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity-based analysis}
%\keywords{States Sequences}
