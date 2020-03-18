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
seqindic(seqdata, indic=c("visited","trans","entr","cplx"), with.missing=FALSE,
    ipos.args=list(), integr.args=list(), prec.args=list(), w=1)
}
\arguments{
  \item{seqdata}{a state sequence object (class \code{stslist}) as returned by \code{\link[TraMineR]{seqdef}}.}

  \item{indic}{vector of character strings. List of selected indicators among \code{"lgth"} (sequence length), \code{"nonm"} (number of non-missing elements), \code{"dlgth"}, (number of successive distinct states), \code{"visited"} (number of visited states), \code{"trans"} (number of state changes), \code{"transp"} (number of state changes as a proportion of maximum possible transitions), \code{"entr"} (longitudinal normalized entropy), \code{"volat"} (objective volatility), \code{"cplx"} (complexity index), \code{"turbn"} (normalized turbulence), \code{"turb"} (turbulence), \code{"integr"} (integrative potential), \code{"ppos"} (proportion of positive states), \code{"vpos"} volatility of positive-negative state sequences, \code{"inpos"} (potential to integrate a positive state"),  \code{"prec"} (precarity index). Setting \code{indic="all"} selects all indicators but \code{"integr"}, \code{"ppos"}, \code{"vpos"}, \code{"inpos"}, and \code{"prec"} that require an additional argument list.}

	\item{with.missing}{logical: should non-void missing values be treated as a regular state? If \code{FALSE} (default) missing values are ignored.}

  \item{ipos.args}{list: when any of \code{"ppos"}, \code{"vpos"}, or \code{"inpos"} is requested, the arguments passed to \code{seqipos} together with \code{seqdata} and \code{with.missing}. See \code{\link{seqipos}}}

  \item{integr.args}{list: when \code{"integr"} requested, the arguments passed to \code{seqintegration} together with \code{seqdata} and \code{with.missing}. See \code{\link{seqintegration}}}

  \item{prec.args}{list: when \code{"prec"} requested, the arguments passed to \code{seqprecarity} together with \code{seqdata} and \code{with.missing}. See \code{\link{seqprecarity}}}
      
  \item{w}{real in range [0,1]: when \code{volat} weight coefficient passed to \code{\link{seqivolatility}}. Default is 1.}
}
\details{
  The number of visited states is the number of different elements in the sequence, e.g. 2 for \code{aababba}.

  The sequence length, number of transitions, longitudinal entropy, complexity, turbulence, proportion of positive states, and precarity are computed respectively with functions \code{\link{seqlength}}, \code{\link{seqtransn}}, \code{\link{seqient}}, \code{\link{seqici}}, \code{\link{seqST}}, \code{\link{seqipos}}, and \code{\link{seqprecarity}}. See corresponding help pages for details.

  The \code{with.missing} argument applies to all indicators but the length. \code{lgth} returns the length obtained with \code{with.missing=TRUE}, and \code{nonm} the length obtained with \code{with.missing=FALSE}.
}
\value{
  A data frame with the selected indicators. Names are:
  \cr \code{Lght}: Length of the sequence
  \cr \code{NonM}: Number of non missing elements
  \cr \code{Dlgth}: Length of the sequence of distinct successive states (DSS)
  \cr \code{Visited}: Number of visited states
  \cr \code{Trans}: Number of state changes
  \cr \code{Transp}: Number of state changes as a proportion of maximum number of transitions
  \cr \code{Entr}: Longitudinal entropy
  \cr \code{Volat}: Objective volatility
  \cr \code{Cplx}: Complexity
  \cr \code{Turbn}: Normalized turbulence
  \cr \code{Turb}: Turbulence
  \cr \code{Integr}: Integrative capacity
  \cr \code{Ppos}: Proportion of positive states
  \cr \code{Vpos}: Volatility of positive-negative state sequences
  \cr \code{Inpos}: Potential to integrate a positive state
  \cr \code{Prec}: Precarity
}
\seealso{
\code{\link{seqlength}}, \code{\link{seqtransn}}, \code{\link{seqient}}, \code{\link{seqici}},  \code{\link{seqST}}, \code{\link{seqipos}}, and \code{\link{seqprecarity}}
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
seqindic(sx, indic=c("volat","integr"), integr.args=list(state="D"))
}
\keyword{Longitudinal characteristics}
