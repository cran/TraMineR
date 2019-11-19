\encoding{utf8}

\name{seqprecarity}

\alias{seqprecarity}
\alias{seqprecorr}
\alias{print.seqprec}
\alias{print.seqprecorr}

\title{Precarity index}

\description{
The precarity index returned by the function \code{seqprecarity} is a modified complexity index that takes account of the ordering (possibly only partial order) of the states and the precarity degree of the first state in the sequence. The signed proportions (transition-based weighted proportions of negative transitions) used for setting the default correction factor can be extracted per se by means of the function \code{seqprecorr}. (See details.)
}

\usage{
seqprecarity(seqdata, correction=NULL, otto=.2, a=1, b=1.2,
    stprec=NULL, method="TRATEDSS",
    state.order=alphabet(seqdata), state.equiv=NULL,
    with.missing=FALSE, ...)

seqprecorr(seqdata, state.order=alphabet(seqdata), state.equiv=NULL,
    penalized="BOTH", method="TRATEDSS", weight.type="ADD", stprec=NULL,
    with.missing=FALSE, border.effect = 10, tr.type)
}

\arguments{
  \item{seqdata}{a state sequence object (class \code{stslist}) as returned by the \code{\link[TraMineR]{seqdef}} function.}

  \item{correction}{Vector of non-negative correction factor values. If \code{NULL}, the correction factors are set as 1 minus the value returned by \code{seqprecorr}.}

  \item{otto}{Scalar in the range [0,1]. Trade-off weight between the precarity degree of the initial state and the corrected complexity. Default is \code{otto=.2}.}

  \item{a}{Non-negative real value. Exponent weight of the complexity. Default is 1.}

  \item{b}{Non-negative real value. Exponent weight of the correction factor. Default is 1.2. }

  \item{stprec}{Vector of state precarity degrees. If \code{NULL} the values are derived from the state order using \code{\link{seqprecstart}}. If not \code{NULL}, the values will be normalized by their maximum.}

  \item{state.order}{Vector of short state labels defining the order of the states. First the less precarious (most positive) state and then the other states in increasing precariousness order. States of the alphabet that are not included here (and are not equivalent to one of the listed state) define the non-comparable states.}

  \item{state.equiv}{List of state equivalence classes. Each class in the list is given as the vector of the short labels of the states forming the class. }

  \item{method}{One of \code{'FREQ'}, \code{'TRATE'}, \code{'TRATEDSS'} (default), \code{'RANK'}, or \code{'ONE'}. Method for determining transition weights. Weights based on transition probabilities: \code{'FREQ'} overall frequency of the transitions, \code{'TRATE'} transition probabilities, and \code{'TRATEDSS'} transition probabilites in the DSS sequences. \code{'RANK'} differences between state precarity degrees. \code{'ONE'} no weight.}

  \item{weight.type}{One of \code{'ADD'} (default), \code{'INV'}, or \code{'LOGINV'}. When \code{method} is one of \code{'FREQ'}, \code{'TRATE'} or \code{'TRATEDSS'}, how weights are derived from the transition probabilities: \code{'ADD'} additive (1-p), \code{'INV'} inverse (1/p), and \code{'LOGINV'} log of inverse. Ignored when any other \code{method} is selected.}

  \item{penalized}{One of \code{'BOTH'} (default), \code{'NEG'}, \code{'POS'}, or \code{'NO'}. What should be penalized or rewarded? \code{'NEG'} only negative transitions, \code{POS} only positive transitions (with negative penalization), \code{'BOTH'} penalize negative transitions and reward positive ones, \code{NO} no penalization. Can also be logical with \code{TRUE} equivalent to \code{'BOTH'} and \code{FALSE} equivalent to \code{'NO'}.}

  \item{with.missing}{Logical. Should the missing state be considered as an element of the alphabet?
    }

  \item{border.effect}{Value (strictly greater than 1) used to adjust estimated transition probabilities to avoid border effect. Default is 10. See details.}

  \item{tr.type}{Deprecated. Use \code{method} instead.}

%  \item{norm}{Logical. Should the weights be normalized by the maximal number of transitions.
%      Not a true normalization because, depending on the transition weights, their sum in a sequence can exceed the number of transitions.}

  \item{...}{Arguments passed to \code{seqprecorr}}
}


\details{The precarity index of a sequence \eqn{x} is based on the complexity index  (\cite{Gabadinho et al., 2010}) \eqn{c(x)} (See the \code{\link[TraMineR]{seqici}} function) and is defined as

 \deqn{prec(x) = \lambda a(x_1) + (1-\lambda) r(x)^\beta c(x)^\alpha}{prec(x) = lambda * a(x_1) + (1 - lambda)*r(x)^b*c(x)^a}

where \eqn{a(x_1)} is the degree of precarity (argument \code{stprec}) of the starting state in sequence \eqn{x}, \eqn{r(x)} is the correction factor (argument \code{correction}) for the sequence. The \eqn{\lambda}{lambda} parameter (argument \code{otto}) determines the trade-off between the importance of the precarity of the starting state and of the corrected complexity index. Parameters \eqn{\alpha}{a} and \eqn{\beta}{b} (argument \code{a} and \code{b}) are exponent weights of respectively the complexity and the correction factor.

When \code{correction = NULL} (default), \eqn{r(x)} is determined as \eqn{r(x) = 1 + q(x)}, where the penalty \eqn{q(x)} is the difference between the weighted proportion of negative and positive transitions. This penalty \eqn{q(x)} can be obtained separately by means of the function \code{seqprecorr}.
The argument \code{penalized} allows to chose between three strategies for computing \eqn{q(x)}: only penalizing negative weights (in which case \eqn{q(x)} is the proportion of negative transitions), only rewarding (with negative penalties) positive transitions, and applying both positive and negative penalties. The transitions can be weighted and the type of transition weights used is selected with the \code{method} argument.
For weights based on transition probabilities, the way how theses probabilites are transformed into weights is controlled with \code{weight.type}. To avoid a border effect, when any computed transition probability \eqn{p} is close from 1 (\eqn{p > 1 - .1/d}), all \eqn{p}'s are adjusted as \eqn{p - p/d}, where \eqn{d} is the \code{border.effect} parameter.

When equivalent classes are provided, the transitions between states of a same class are ignored and the class mean precarity degree is assigned to each state of the class (see \code{\link{seqprecstart}}).

As for non-comparable states (those not listed on the \code{state.order} argument and not equivalent to a listed state), the transition to or from such states is ignored and they all receive the mean value (.5) as state precarity degree.

When there are missing states in the sequences, set \code{with.missing = TRUE} to treat the missing state as an additional state. In that case the missing state will be considered as non-comparable unless you include the \code{nr} attribute of \code{seqdata} in \code{state.order} or \code{state.equiv}. With \code{with.missing = FALSE}, transitions to and from the missing state will just be ignored. However, when \code{with.missing = FALSE}, there will be no \code{stprec} value for the missing state and, as a result, the precarity value will be \code{NA} for sequences starting with a missing state. To drop these starting missing states, use \code{left='DEL'} in \code{\link{seqdef}} when defining the state sequence object.
}



\value{

For \code{seqprecarity}, an object of class \code{seqprec} with the value of the precarity index for each sequence. The returned object has an attribute \code{stprec} that contains the state precarity degree used at the starting position. The associated print method (\code{print.seqprec}) prints the state precarity values without the additional attribute.


For \code{seqprecorr}, an object of class \code{seqprecorr} with the weighted proportions \eqn{q(x)} and as additional attributes: \code{tr} the used transition weights; \code{signs} the transitions signs; \code{state.noncomp} the non-comparable states; and \code{state.order} the used state order. The associated print method (\code{print.seqprecorr}) prints the \eqn{q(x)} values without the additional attributes.

}

\seealso{
 \code{\link{seqici}}, \code{\link{seqprecstart}}, \code{\link{seqindic}}.
}

\references{
Ritschard, G., Bussi, M., and O'Reilly, J. (2018), "An index of precarity for
measuring early employment insecurity", in G. Ritschard, and M. Studer, \emph{Sequence Analysis and Related Approaches: Innovative Methods and Applications}, Series Life Course Research and Social
Policies, Vol. 10, pp 279-295. Cham: Springer.

Gabadinho, A., Ritschard, G., Studer, M. and Müller, N.S. (2010), "Indice de complexité pour le tri et la comparaison de séquences catégorielles", In Extraction et gestion des connaissances (EGC 2010), \emph{Revue des nouvelles technologies de l'information RNTI}. Vol. E-19, pp. 61-66.

}

\examples{
## Defining a sequence object with columns 13 to 24
## in the 'actcal' example data set
data(actcal)
actcal <- actcal[1:200,] ## Here, only a subset
actcal.seq <- seqdef(actcal[,13:24])

## precarity using the original state order
prec <- seqprecarity(actcal.seq)
ici <- seqici(actcal.seq) ## complexity

seqn <- seqformat(actcal.seq, to="SPS", compress=TRUE)
tab <- data.frame(seqn,ici,prec)
names(tab) <- c("seq","ici","prec")
head(tab)

## Assuming A and B as equivalent regarding precarity
prec2 <- seqprecarity(actcal.seq, state.equiv=list(c('A','B')))
tab <- cbind(tab,prec2)
names(tab)[ncol(tab)] <- "prec2"

## and letting C be non-comparable
prec3 <- seqprecarity(actcal.seq, state.order=c("A","B","D"),
  state.equiv=list(c('A','B')))
tab <- cbind(tab,prec3)
names(tab)[ncol(tab)] <- "prec3"
head(tab)

## Extracting the q(x) used for the correction factor (1-q(x))
q <- seqprecorr(actcal.seq, state.order=c("A","B","D"),
  state.equiv=list(c('A','B')))
## correction factor
corr.f <- 1 + q
## number of non neutral correction factors
length(corr.f[corr.f != 1])

## Precarity with rank-based transition weights
prec.rank <- seqprecarity(actcal.seq, method='RANK')

## Precarity in presence of missing values:
##  missing state treated as an additional state
data(ex1)
## by default right missings are dropped from the sequences
s.ex1 <- seqdef(ex1[,1:13])
seqprecarity(s.ex1, with.missing=TRUE)


}

\author{Gilbert Ritschard}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Longitudinal characteristics}
