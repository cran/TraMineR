\encoding{utf8}

\name{seqprecarity}

%\alias{seqidegrad}
\alias{seqprecarity}
%%\alias{seqprecarity.old}
\alias{seqprecorr}
\alias{print.seqprec}
\alias{print.seqprecorr}

\title{Precarity index}

\description{
%The degrading index returned by \code{seqidegrad} is a normalized transformation of the difference between the sum of downward and upward state changes (transitions).
The precarity index is essentially the complexity index corrected by the difference between the proportion of downward and upward transitions, and the undesirableness of the first state in the sequence.
}

\usage{
seqprecarity(seqdata, correction=NULL,
    otto=.2, a=1, b=1.2, stprec=NULL, method = "TRATEDSS",
    state.order=alphabet(seqdata, with.missing), state.equiv=NULL, 
    with.missing=FALSE, ...)

seqprecorr(seqdata, state.order=alphabet(seqdata, with.missing), state.equiv=NULL,
    penalized="BOTH", method="TRATEDSS", weight.type="ADD",
    stprec=NULL, with.missing=FALSE, border.effect=10, tr.type)
}

\arguments{
  \item{seqdata}{a state sequence object (class \code{stslist}) as returned by the \code{\link[TraMineR]{seqdef}} function.}
  
  %\item{type}{Integer: 1 or 2. \code{type = 1} is the precarity index of \cite{Ritschard et al (2018)} and \code{type = 2} is the modified index defined in \cite{Ritschard (2020)}.}

  \item{correction}{Vector of non-negative correction factor values. If \code{NULL}, the correction factors are set as 1 minus the value returned by \code{seqprecorr}.}

  \item{otto}{Scalar in the range [0,1]. Trade-off weight between the precarity degree of the initial state and the corrected complexity. Default is \code{otto=.5} for \code{seqprecarity} (.2 for \code{seqprecarity.old}.}

  \item{a}{Non-negative real value. Exponent weight of the complexity. Default is 1.}

  \item{b}{Non-negative real value. Exponent weight of the correction factor. Default is 1.2. }

  \item{stprec}{Vector of state undesirableness degrees. If \code{NULL} the values are derived from the state order using \code{\link{seqprecstart}}. If not \code{NULL}, the values in \code{stprec} should conform the order of the alphabet. 
      Use negative values to indicate non-comparable states. See details.}

  \item{state.order}{Vector of short state labels defining the order of the states. First the less precarious (most positive) state and then the other states in increasing precariousness order. States of the alphabet that are not included here (and are not equivalent to one of the listed state) define the non-comparable states.}

  \item{state.equiv}{List of state equivalence classes. Each class in the list is given as the vector of the short labels of the states forming the class. }

  \item{method}{One of \code{'FREQ'}, \code{'FREQ+'}, \code{'TRATE'}, \code{'TRATE+'}, \code{'TRATEDSS'} (default for \code{seqprecarity.old} and \code{seqprecorr}), \code{'TRATEDSS+'}, \code{'RANK'} (default for \code{seqprecarity} and \code{seqidegrad}), \code{'RANK+'}, or \code{'ONE'}. Method for determining transition weights. Weights based on transition probabilities: \code{'FREQ'} overall frequency of the transitions, \code{'TRATE'} transition probabilities, and \code{'TRATEDSS'} transition probabilites in the DSS sequences. \code{'RANK'} differences between state precarity degrees. \code{'ONE'} no weight. With the \code{+} form the returned penalty is adjusted by the mean transition weight in the sequence.}

  \item{weight.type}{One of \code{'ADD'} (default), \code{'INV'}, or \code{'LOGINV'}. When \code{method} is one of \code{'FREQ'}, \code{'TRATE'} or \code{'TRATEDSS'}, how weights are derived from the transition probabilities: \code{'ADD'} additive (1-p), \code{'INV'} inverse (1/p), and \code{'LOGINV'} log of inverse. Ignored when any other \code{method} is selected.}

  \item{penalized}{One of \code{'BOTH'} (default), \code{'NEG'}, \code{'POS'}, or \code{'NO'}. What should be penalized or rewarded? \code{'NEG'} only negative transitions, \code{POS} only positive transitions (with negative penalization), \code{'BOTH'} penalize negative transitions and reward positive ones, \code{NO} no penalization. Can also be logical with \code{TRUE} equivalent to \code{'BOTH'} and \code{FALSE} equivalent to \code{'NO'}.}

  \item{with.missing}{Logical. Should the missing state be considered as an element of the alphabet?
    }

  \item{border.effect}{Value (strictly greater than 1) used to adjust estimated transition probabilities to avoid border effect. Default is 10. See details.}

%  \item{degr}{Logical. Should the degradation index be used as correction factor? Default is \code{FALSE} for backward compatibility. See details.}

%  \item{start.integr}{Logical. Should precarity level of start state be adjusted by the integration index of the first spell? See details.}

%  \item{spell.integr}{Logical. Should transition weight be adjusted by the integration index of the following spell? See details.}


  \item{tr.type}{Deprecated. Use \code{method} instead.}

%  \item{norm}{Logical. Should the weights be normalized by the maximal number of transitions.
%      Not a true normalization because, depending on the transition weights, their sum in a sequence can exceed the number of transitions.}

  \item{...}{Arguments passed to \code{seqprecorr}}
}


\details{
The \code{seqprecorr} function returns the penalizing factor \eqn{q(x)}, i.e. the difference between the proportions of downward and upward transitions (state changes).

The argument \code{penalized} allows to chose between three strategies for computing \eqn{q(x)}: only penalizing negative weights (in which case \eqn{q(x)} is the proportion of negative transitions), only rewarding (with negative penalties) positive transitions, and applying both positive and negative penalties. The transitions can be weighted and the type of transition weights used is selected with the \code{method} argument. For weights based on transition probabilities, the way how theses probabilites are transformed into weights is controlled with \code{weight.type}. To avoid a border effect, when any computed transition probability \eqn{p} is close from 1 (\eqn{p > 1 - .1/d}), all \eqn{p}'s are adjusted as \eqn{p - p/d}, where \eqn{d} is the \code{border.effect} parameter. With \code{method="RANK"}, the weights are set as the differences between the to and from state undesirableness.


The precarity index of a sequence \eqn{x} is based on the complexity index (\cite{Gabadinho et al., 2010}) \eqn{c(x)} (See the \code{\link[TraMineR]{seqici}} function), and the  undesirableness \eqn{a(x_1)} of the starting state. It is defined as

 \deqn{prec(x) = \lambda a(x_1) + (1-\lambda) r(x)^\beta c(x)^\alpha}{prec(x) = lambda * a(x_1) + (1 - lambda)*r(x)^b*c(x)^a}

where \eqn{r(x)} is the correction factor (argument \code{correction}) for the sequence. The \eqn{\lambda}{lambda} parameter (argument \code{otto}) determines the trade-off between the importance of the undesirableness of the starting state and of the corrected complexity index. Parameters \eqn{\alpha}{a} and \eqn{\beta}{b} (argument \code{a} and \code{b}) are exponent weights of respectively the complexity and the correction factor.

When \code{correction = NULL} (default), \eqn{r(x)} is determined as \eqn{r(x) = 1 + q(x)}, where the penalty \eqn{q(x)} is the degrading index. %provided by \code{seqidegrad}. 
When \code{correction = NULL} and \code{type=2}, the correction is set as \eqn{r(x) = (1 + q(x))/2}.


%When \code{spell.integr=TRUE}, the weight of each transition is multiplied by the integration index of the following spell in the computation of the proportions of downward and upward transitions. When \code{start.integr=TRUE} the precarity degree of the start state is multiplied by the integration index of the first spell.

When \code{stprec} is a vector, negative values indicate non-comparable sates that receive each the mean positive undesirableness value. After this transformation, the vector is normalized such that the minimum is 0 and the maximum 1.

When equivalent classes are provided, the class mean undesirableness degree is assigned to each state of the class (see \code{\link{seqprecstart}}). For the count of transitions a same state value is assigned to all equivalent states.

Non-comparable states (those not listed on the \code{state.order} argument and not equivalent to a listed state) all receive the mean undesirableness value. For the count of transitions, transitions from and to non-comparable states are ignored and replaced by a transition between the states that immediately precede and follow a spell in non-comparable states.

When there are missing states in the sequences, set \code{with.missing = TRUE} to treat the missing state as an additional state. In that case the missing state will be considered as non-comparable unless you include the \code{nr} attribute of \code{seqdata} in \code{state.order} or \code{state.equiv}. With \code{with.missing = FALSE}, transitions to and from the missing state will just be ignored and the undesirableness value of the first valid state will be used as starting undesirableness. 
%To drop these starting missing states, use \code{left='DEL'} in \code{\link{seqdef}} when defining the state sequence object.
}



\value{

For \code{seqprecarity}, an object of class \code{seqprec} with the value of the precarity index for each sequence. The returned object has an attribute \code{stprec} that contains the state precarity degree used at the starting position. The associated print method (\code{print.seqprec}) prints the state precarity values without the additional attribute.


For \code{seqprecorr} %and \code{seqidegrad}, 
an object of class \code{seqprecorr} with the weighted proportions \eqn{q(x)} 
%for \code{seqprecorr} or the degradation index \eqn{(1-q(x))/2} for \code{seqidegrad}
and as additional attributes: \code{tr} the used transition weights; \code{signs} the transitions signs; \code{state.noncomp} the non-comparable states; and \code{state.order} the used state order. The associated print method (\code{print.seqprecorr}) prints the outcome values without the additional attributes.

}

\seealso{
 \code{\link{seqici}}, \code{\link{seqprecstart}}, \code{\link{seqindic}}.
}

\references{
Ritschard, G. (2020), "Measuring the nature of individual sequences", manuscript.

Ritschard, G., Bussi, M., and O'Reilly, J. (2018), "An index of precarity for
measuring early employment insecurity", in G. Ritschard, and M. Studer, \emph{Sequence Analysis and Related Approaches: Innovative Methods and Applications}, Series Life Course Research and Social
Policies, Vol. 10, pp 279-295. Cham: Springer.

Gabadinho, A., Ritschard, G., Studer, M. and M\"{u}ller, N.S. (2010), "Indice de complexit\'{e} pour le tri et la comparaison de s\'{e}quences cat\'{e}gorielles", In Extraction et gestion des connaissances (EGC 2010), \emph{Revue des nouvelles technologies de l'information RNTI}. Vol. E-19, pp. 61-66.

}

\examples{
## Defining a sequence object with columns 13 to 24
## in the 'actcal' example data set
data(actcal)
actcal <- actcal[1:20,] ## Here, only a subset
actcal.seq <- seqdef(actcal[,13:24], alphabet=c('A','B','C','D'))

## precarity using the original state order
prec <- seqprecarity(actcal.seq)
ici <- seqici(actcal.seq) ## complexity

seqn <- seqformat(actcal.seq, to="SPS", compress=TRUE)
tab <- data.frame(seqn,ici,prec)
names(tab) <- c("seq","ici","prec")
tab

## Assuming A and B as equivalent regarding precarity
prec2 <- seqprecarity(actcal.seq, state.equiv=list(c('A','B')))
tab <- cbind(tab,prec2)
names(tab)[ncol(tab)] <- "prec2"

## and letting C be non-comparable
prec3 <- seqprecarity(actcal.seq, state.order=c("A","B","D"),
  state.equiv=list(c('A','B')))
tab <- cbind(tab,prec3)
names(tab)[ncol(tab)] <- "prec3"


## Extracting the q(x) used for the correction factor (1-q(x))
q <- seqprecorr(actcal.seq, state.order=c("A","B","D"),
  state.equiv=list(c('A','B')))
## correction factor
corr.f <- (1 + q)
## number of non neutral correction factors
length(corr.f[corr.f != 1])

## Precarity with transition weights based on transition probabilities
prec.trdss <- seqprecarity(actcal.seq, method='TRATEDSS')

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
