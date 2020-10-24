\name{seqindic}
\alias{seqindic}
%
\author{Gilbert Ritschard}
%
\title{Sequence indicators}
%
\description{
  Table of per sequence values of selected indicators.
}
\usage{
seqindic(seqdata, indic=c("visited","trans","entr","cplx","turb2n"), with.missing=FALSE,
    ipos.args=list(), prec.args=list(), w=1)
}
\arguments{
  \item{seqdata}{a state sequence object (class \code{stslist}) as returned by \code{\link[TraMineR]{seqdef}}.}

  \item{indic}{vector of character strings. List of selected indicators among \code{"lgth"} (sequence length), \code{"nonm"} (number of non-missing elements), \code{"dlgth"}, (number of successive distinct states), \code{"visited"} (number of visited states), \code{"recu"} (recurrence index), \code{"trans"} (number of state changes), \code{"transp"} (number of state changes as a proportion of maximum possible transitions), \code{"meand"} (mean spell duration), \code{"dustd"} (duration standard deviation),
      %\code{"meand2"} (mean spell duration accounting for non-visited states), \code{"dustd2"} (duration standard deviation accounting for non-visited states), 
      \code{"entr"} (longitudinal normalized entropy), \code{"volat"} (objective volatility), \code{"nsubs"} (number of subsequences of the DSS sequence), \code{"cplx"} (complexity index),
      \code{"turbn"} (normalized turbulence), \code{"turb"} (turbulence), 
      %\code{"turbn2"} (normalized turbulence accounting for non-visited states), \code{"turb2"} (turbulence accounting for non-visited states), 
      \code{"ppos"} (proportion of positive states), \code{"nvolat"} (normative volatility), \code{"vpos"} volatility of positive-negative state sequences, \code{"inpos"} (integration potential),  \code{"prec"} (precarity index). Can also be \code{"all"}, \code{"basic"}, \code{"diversity"}, \code{"complexity"},  or \code{"binary"} that each select a subset of the indicators. Setting \code{indic="all"} selects all indicators but \code{"ppos"}, \code{"vpos"}, \code{"nvolat"}, \code{"inpos"}, and \code{"prec"} that require an additional argument list.}

	\item{with.missing}{logical: should non-void missing values be treated as a regular state? If \code{FALSE} (default) missing values are ignored.}

  \item{ipos.args}{list: when any of \code{"ppos"}, \code{"nvolat"}, \code{"vpos"}, or \code{"inpos"} is requested, the arguments passed to \code{seqipos} together with \code{seqdata} and \code{with.missing}. See \code{\link{seqipos}}}

  %\item{integr.args}{list: when \code{"integr"} requested, the arguments passed to \code{seqintegration} together with \code{seqdata} and \code{with.missing}. See \code{\link{seqintegration}}.}

  \item{prec.args}{list: when \code{"prec"} requested, the arguments passed to \code{seqprecarity} together with \code{seqdata} and \code{with.missing}. See \code{\link{seqprecarity}}.}

  \item{w}{real in range [0,1]: when \code{volat} weight coefficient passed to \code{\link{seqivolatility}}. Default is 1.}
}
\details{
  The number of visited states is the number of different elements in the sequence, e.g. 2 for \code{aababba}. The recurrence index \code{'recu'} is the average number of visits of visited states, i.e. \code{Dlgth/Visited}, the number of spells on the number of visited states.

  The sequence length, number of transitions, longitudinal entropy, volatility, duration standard deviation, complexity, turbulence, proportion of positive states, integrative potential, and precarity are computed respectively with functions \code{\link{seqlength}}, \code{\link{seqtransn}}, \code{\link{seqient}}, \code{\link{seqivolatility}}, \code{\link{seqivardur}}, \code{\link{seqici}}, \code{\link{seqST}}, \code{\link{seqipos}}, \code{\link{seqintegration}}, and \code{\link{seqprecarity}}. See corresponding help pages for details.

  The proportion of positive states (\code{'ppos'}) and the normative volatility (\code{'nvolat'}) are the proportions of positive states in respectively the original sequences and the DSS. They ignore the value of \code{dss} in the \code{ipos.args} list.

  The \code{with.missing} argument applies to all indicators but the length. \code{'lgth'} returns the length obtained with \code{with.missing=TRUE}, and \code{'nonm'} the length obtained with \code{with.missing=FALSE}.
}
\value{
  A data frame with the selected indicators. Names are:
  \cr \code{Lght}: Length of the sequence
  \cr \code{NonM}: Number of non missing elements
  \cr \code{Dlgth}: Length of the sequence of distinct successive states (DSS)
  \cr \code{Visited}: Number of visited states
  \cr \code{Recu}: Average number of visits to visited states
  \cr \code{Trans}: Number of state changes
  \cr \code{Transp}: Number of state changes as a proportion of maximum number of transitions
  \cr \code{Meand}: Mean spell duration
  \cr \code{Dustd}: Duration standard deviation  %\cr \code{Meand2}: Mean spell duration accounting for non visited states   %\cr \code{Dustd2}: Duration standard deviation accounting for non visited states
  \cr \code{Entr}: Longitudinal entropy
  \cr \code{Volat}: Objective volatility
  \cr \code{Nsubs}: Number of subsequences of the DSS sequence
  \cr \code{Cplx}: Complexity
  \cr \code{Turbn}: Normalized turbulence
  \cr \code{Turb}: Turbulence   %\cr \code{Turbn2}: Normalized turbulence accounting for non visited states    %\cr \code{Turb2}: Turbulence accounting for non visited states    %\cr \code{Integr}: Integrative capacity
  \cr \code{Ppos}: Proportion of positive states
  \cr \code{Nvolat}: Normative volatility
  \cr \code{Vpos}: Volatility of positive-negative state sequences
  \cr \code{Inpos}: Integrative capacity   %\cr \code{Bad}: Badness    %\cr \code{Degrad}: Degrading
  \cr \code{Prec}: Precarity   %\cr \code{Prec2}: Revised Precarity
}
\seealso{
\code{\link{seqlength}}, \code{\link{seqtransn}}, \code{\link{seqient}}, \code{\link{seqivolatility}}, \code{\link{seqivardur}}, \code{\link{seqici}},  \code{\link{seqST}}, \code{\link{seqipos}}, \code{\link{seqintegration}} and \code{\link{seqprecarity}}
}
\examples{
data(ex1)
sx <- seqdef(ex1[,1:13], right="DEL")
print(sx, format='SPS')
seqindic(sx, indic=c("lgth","nonm","visited","turbn","cplx"))
seqindic(sx, indic=c("lgth","nonm","visited","turbn","cplx"), with.missing=TRUE)
seqindic(sx, indic=c("lgth","dlgth","ppos","inpos","prec"), with.missing=TRUE,
          ipos.args=list(pos.states=c("A","B")),
          prec.args=list(state.order=c("A","B","C"), state.equiv=list(c("C","D"))))
seqindic(sx, indic=c("volat","binary"), ipos.args=list(pos.states=c("A","B")))
seqindic(sx, indic=c("basic","inpos"), ipos.args=list(pos.states="D"))
}
\keyword{Longitudinal characteristics}
