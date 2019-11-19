\name{seqipos}
\alias{seqipos}
%
\author{Gilbert Ritschard}
%
\title{Proportion of positive spells or states per sequence}
%
\description{
  Returns the proportion of positive states within each sequence, i.e. of positive spells if computed on the sequences of distinct successive states (DSS).
}
\usage{
  seqipos(seqdata, dss=TRUE, pos.states=NULL, neg.states=NULL, with.missing=FALSE)
}
\arguments{
  \item{seqdata}{a state sequence object (\code{stslist}) as returned by \code{\link[TraMineR]{seqdef}}.}

  \item{dss}{logical. Should the proportion be computed inside the DSS sequences?}

  \item{pos.states}{vector of positive states.}

  \item{neg.states}{vector of negative states.}

	\item{with.missing}{logical: should non-void missing values be treated as a regular state? If \code{FALSE} (default) missing values are ignored.}
}
\details{
  Letting \eqn{npos} be the number of positive states and \eqn{nneg} the number of negative states in the sequence, the function returns the value of \eqn{npos/(npos + nneg)} for each sequence.

  When \code{dss=TRUE} (default), the \eqn{npos} and \eqn{nneg} are respectively the sum of positive and negative spells (in a same state).

  When \code{neg.states=NULL}, states not listed on the \code{pos.states} argument are considered as negative and conversely when \code{pos.states=NULL}. When \code{with.missing=TRUE}, the missing state \code{nr} will be considered as positive if included in the \code{pos.states} list or as negative if in \code{neg.states}.  When \code{with.missing=FALSE} (default) missing states are ignored. However, when missing values appear within a spell, such as in AA*ABB, the DSS sequence will have two consecutive occurrences of the spell. For the example we would have AAB, i.e. if A is positive and B negative, a proportion of 2/3 of positive spells.

  When both \code{pos.states} and \code{neg.states} are provided, states of the alphabet (including the \code{nr} code when \code{with.missing=TRUE}) that belong nor to \code{pos.states} nor to \code{neg.states} are ignored.
}

\value{
  A numeric vector with the proportion of positive states (or spells) within each sequence.
}

\seealso{
\code{\link{seqindic}}
}
\examples{
data(ex1)
sx <- seqdef(ex1[,1:13], right="DEL")
nr <- attr(sx,'nr') ## code for missing values

seqipos(sx, pos.states=c("A","B"))
seqipos(sx, pos.states=c("A","B"), with.missing=TRUE)
## ignoring state "D"
seqipos(sx, pos.states=c("A","B"), neg.states=c("C",nr), with.missing=TRUE)

seqipos(sx, pos.states=c("A","B"), dss=FALSE)
seqipos(sx, pos.states=c("A","B",nr), dss=FALSE, with.missing=TRUE)

}
\keyword{Longitudinal characteristics}
