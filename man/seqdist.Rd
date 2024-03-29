% Author for TraMineR 2: Pierre-Alexandre Fonta (2016-2017)

\name{seqdist}
\alias{seqdist}
\title{Distances (dissimilarities) between sequences}
\description{
  Computes pairwise dissimilarities between sequences or dissimilarity from
  a reference sequence. Several dissimilarity measures can be chosen, including
  optimal matching (OM) and many of its variants, distance based on the count
  of common attributes, and distances between sequence state distributions.
}
\usage{
seqdist(seqdata, method, refseq = NULL, norm = "none", indel = "auto", sm = NULL,
  with.missing = FALSE, full.matrix = TRUE, kweights = rep(1.0, ncol(seqdata)),
  tpow = 1.0, expcost = 0.5, context, link = "mean", h = 0.5, nu,
  transindel = "constant", otto, previous = FALSE, add.column = TRUE,
  breaks = NULL, step = 1, overlap = FALSE, weighted = TRUE,
  global.pdotj = NULL, prox = NULL, check.max.size=TRUE)
}
\arguments{
  \item{seqdata}{
    State sequence object of class \code{stslist}.
    The sequence data to use.
    Use \code{\link{seqdef}} to create such an object.
  }
  \item{method}{
    String.
    The dissimilarity measure to use.
    It can be \code{"OM"}, \code{"OMloc"}, \code{"OMslen"}, \code{"OMspell"},
    \code{"OMstran"}, \code{"HAM"}, \code{"DHD"}, \code{"CHI2"}, \code{"EUCLID"},
    \code{"LCS"}, \code{"LCP"}, \code{"RLCP"}, \code{"NMS"}, \code{"NMSMST"},
    \code{"SVRspell"}, or \code{"TWED"}. See the Details section.
  }
  \item{refseq}{
    \code{NULL}, Integer, State Sequence Object, or List.
    Default: \code{NULL}.
    The baseline sequence to compute the distances from.

    When an integer, the index of a sequence in \code{seqdata} or \code{0} for the most frequent sequence.

    When a state sequence object, it must contain a single sequence and have the same
    alphabet as \code{seqdata}.

    When a list, it must be a list of two sets of indexes of \code{seqdata} rows.
  }
  \item{norm}{
    String.
    Default: \code{"none"}.
    The normalization to use when \code{method} is one of \code{"OM"},
    \code{"OMloc"}, \code{"OMslen"}, \code{"OMspell"},
    \code{"OMstran"}, \code{"TWED"}, \code{"HAM"}, \code{"DHD"}, \code{"LCS"},
    \code{"LCP"}, \code{"RLCP"}, \code{"CHI2"}, \code{"EUCLID"}.
    It can be \code{"none"}, \code{"auto"}, or, except for
    \code{"CHI2"} and \code{"EUCLID"}, \code{"maxlength"},
    \code{"gmean"}, \code{"maxdist"}, or \code{"YujianBo"}. \code{"auto"} is
    equivalent to \code{"maxlength"} when \code{method} is one of \code{"OM"},
    \code{"HAM"}, or \code{"DHD"}, to \code{"gmean"} when \code{method} is one
    of \code{"LCS"}, \code{"LCP"}, or \code{"RLCP"}, to \code{YujianBo} when
    \code{method} is one of \code{"OMloc"}, \code{"OMslen"}, \code{"OMspell"},
    \code{"OMstran"}, \code{"TWED"}. See the Details section.
    % Others: CHI2, EUCLID, OMloc, OMslen, OMspell, OMstran, NMS, NMSMST, SVRspell, TWED.
  }
  \item{indel}{
    Double, Vector of Doubles, or String.
    Default: \code{"auto"}.
    Insertion/deletion cost(s). Applies when \code{method} is one of \code{"OM"}, \code{"OMslen"}, \code{"OMspell"},
    or \code{"OMstran"}.

    The single state-independent insertion/deletion cost when a double.
    %, or \code{"TWED"}.

    The state-dependent insertion/deletion costs when a vector of doubles.
    The vector should contain an indel cost by state in the order of the alphabet.

    When \code{"auto"}, the indel is set as \code{max(sm)/2} when \code{sm} is
    a matrix and is computed by means of \code{\link{seqcost}} when \code{sm} is
    a string specifying a cost method.

    % Others: OMloc, HAM, DHD, CHI2, EUCLID, LCS, LCP, RLCP, NMS, NMSMST, SVRspell.
  }
  \item{sm}{
    \code{NULL}, Matrix, Array, or String. Substitution costs.
    Default: \code{NULL}.

    The substitution-cost matrix when a matrix and \code{method} is one of
    \code{"OM"}, \code{"OMloc"}, \code{"OMslen"}, \code{"OMspell"},
    \code{"OMstran"}, \code{"HAM"}, or \code{"TWED"}.

    The series of the substitution-cost matrices when an array and
    \code{method = "DHD"}. They are grouped in a 3-dimensional array with the
    third index referring to the position in the sequence.

    One of the strings \code{"CONSTANT"}, \code{"INDELS"}, \code{"INDELSLOG"},
    or \code{"TRATE"}. Designates a \code{\link{seqcost}} method
    to build \code{sm}. \code{"CONSTANT"} is not relevant for \code{"DHD"}.

    %for \code{"OM"}, \code{"OMloc"}, \code{"OMslen"}, \code{"OMspell"},
    %\code{"OMstran"}, \code{"HAM"}, and \code{"TWED"}, \code{"CONSTANT"} for
    %\code{"OM"} and \code{"HAM"}, \code{"TRATE"} for \code{"OM"}, \code{"HAM"},
    %and \code{"DHD"}.

    \code{sm} is mandatory when \code{method} is one of \code{"OM"},
    \code{"OMloc"}, \code{"OMslen"}, \code{"OMspell"}, \code{"OMstran"},
    or \code{"TWED"}.

    \code{sm} is autogenerated when \code{method} is one of \code{"HAM"} or
    \code{"DHD"} and \code{sm = NULL}. See the Details section.

    Note: With \code{method = "NMS"} or \code{method = "SVRspell"}, use
    \code{prox} instead.

    % Others: CHI2, EUCLID, LCS, LCP, RLCP, NMS, NMSMST, SVRspell.
  }
  \item{with.missing}{
    Logical.
    Default: \code{FALSE}.
    Should the non-deleted missing value be added to the alphabet as an additional
    state? If \code{FALSE} and \code{seqdata} or \code{refseq} contains such
    gaps, an error is raised.
  }
  \item{full.matrix}{
    Logical.
    Default: \code{TRUE}.
    When \code{refseq = NULL}, if \code{TRUE}, the full distance matrix is
    returned, if \code{FALSE}, an object of class \code{\link{dist}} is returned,
    that is, a vector containing only values from the lower triangle of the
    distance matrix. Objects of class \code{dist} are smaller and can be passed
    directly as arguments to most clustering functions.
  }
  \item{kweights}{
    Double or vector of doubles.
    Default: vector of \code{1}s.
    The weights applied to subsequences when \code{method} is one of \code{"NMS"},
    \code{"NMSMST"}, or \code{"SVRspell"}. It contains at position \eqn{k} the
    weight applied to the subsequences of length \eqn{k}. It must be positive. 
    Its length should be equal to the number of columns of \code{seqdata}. If shorter,
    longer subsequences are ignored. If a scalar, it is transformed into
    \code{rep(kweights,ncol(sedata))}.
  }
  \item{tpow}{
    Double.
    Default: \code{1.0}.
    The exponential weight of spell length when \code{method} is one of
    \code{"OMspell"}, \code{"NMSMST"}, or \code{"SVRspell"}.
  }
  \item{expcost}{
    Double.
    Default: \code{0.5}.
    The cost of spell length transformation when \code{method = "OMloc"} or
    \code{method = "OMspell"}. It must be positive. The exact interpretation is
    distance-dependent.
  }
  \item{context}{
    Double.
    Default: \code{1-2*expcost}.
    The cost of local insertion when \code{method = "OMloc"}. It must be positive.
  }
  \item{link}{
    String.
    Default: \code{"mean"}.
    The function used to compute substitution costs when \code{method = "OMslen"}.
    One of \code{"mean"} (arithmetic average) or \code{"gmean"} (geometric mean
    as in the original proposition of Halpin 2010).
  }
  \item{h}{
    Double.
    Default: \code{0.5}.
    It must be greater than or equal to 0.

    The exponential weight of spell length when \code{method = "OMslen"}.

    The gap penalty when \code{method = "TWED"}. It corresponds to the lambda
    in \cite{Halpin (2014), p 88}. It is usually chosen in the range [0,1]
  }
  \item{nu}{
    Double.
    Stiffness when \code{method = "TWED"}. It must be strictly greater than 0
    and is usually less than 1.
    See \cite{Halpin (2014), p 88}.
  }
  \item{transindel}{
    String.
    Default: \code{"constant"}.
    Method for computing transition indel costs when \code{method = "OMstran"}.
    One of \code{"constant"} (single indel of 1.0), \code{"subcost"} (based on
    substitution costs), or \code{"prob"} (based on transition probabilities).
  }
  \item{otto}{
    Double.
    The origin-transition trade-off weight when \code{method = "OMstran"}. It
    must be in [0, 1].
  }
  \item{previous}{
    Logical.
    Default: \code{FALSE}.
    When \code{method = "OMstran"}, should we also account for the transition
    from the previous state?
  }
  \item{add.column}{
    Logical.
    Default: \code{TRUE}.
    When \code{method = "OMstran"}, should the last column (and also the first
    column when \code{previous = TRUE}) be duplicated? When sequences have different
    lengths, should the last (first) valid state be duplicated.
  }
  \item{breaks}{
    \code{NULL}, List of pairs Integers.
    Default: \code{NULL}.
    The list of the possibly overlapping intervals when \code{method = "CHI2"}
    or \code{method = "EUCLID"}.
  }
  \item{step}{
    Integer.
    Default: \code{1}.
    The length of the intervals when \code{method = "CHI2"} or
    \code{method = "EUCLID"} and \code{breaks = NULL}. It must be positive.
    It must also be even when \code{overlap = TRUE}.
  }
  \item{overlap}{
    Logical.
    Default: \code{FALSE}.
    When \code{method = "CHI2"} or \code{method = "EUCLID"} and
    \code{breaks = NULL}, should the intervals overlap?
  }
  \item{weighted}{
    Logical.
    Default: \code{TRUE}.
    When \code{method} is \code{"CHI2"} or when \code{sm} is a string (method),
    should the distributions of the states account for the sequence weights
    in \code{seqdata}? See \code{\link{seqdef}}.
  }
  \item{global.pdotj}{
    Numerical vector, \code{"obs"}, or \code{NULL}.
    Default: \code{NULL}.
    Only for \code{method = "CHI2"}.
    The vector of state proportions to be used as marginal distribution. When \code{NULL}, the state distribution on the corresponding interval is used. When \code{"obs"}, the overall state distribution in \code{seqdata} is used for all intervals. When a vector of proportions, it is used as marginal distribution for all intervals.
  }
  \item{prox}{
    \code{NULL} or Matrix.
    Default: \code{NULL}.
    The matrix of state proximities when \code{method = "NMS"} or
    \code{method = "SVRspell"}.
  }
  \item{check.max.size}{
    Logical. Should \code{seqdist} stop when maximum allowed number of unique sequences is exceeded? Caution, setting \code{FALSE} may produce unexpected results or even crash R.
  }
}
\details{
  The \code{seqdist} function returns a matrix of distances between sequences
  or a vector of distances from the reference sequence when \code{refseq} is set.
  The available metrics (see \code{method} option) include:
  \itemize{
    \item{
      \emph{Edit distances}: optimal matching (\code{"OM"}), localized OM
      (\code{"OMloc"}), spell-length-sensitive OM (\code{"OMslen"}), OM of spell
      sequences (\code{"OMspell"}), OM of transition sequences (\code{"OMstran"}),
      Hamming (\code{"HAM"}), dynamic Hamming (\code{"DHD"}), and the time warp edit
      distance (\code{"TWED"}).
    }
    \item{
      \emph{Metrics based on counts of common attributes}: distance based on
      the longest common subsequence (\code{"LCS"}), on the longest common prefix
      (\code{"LCP"}), on the longest common suffix (\code{"RLCP"}), on the number
      of matching subsequences (\code{"NMS"}), on the number of matching
      subsequences weighted by the minimum shared time (\code{"NMSMST"}) and,
      the subsequence vectorial representation distance (\code{"SVRspell"}).
    }
    \item{
      \emph{Distances between state distributions}: Euclidean (\code{"EUCLID"}),
      Chi-squared (\code{"CHI2"}).
    }
  }

  See \cite{Studer and Ritschard (2014, 2016)} for a description and the comparison
  of the above dissimilarity measures except \code{"TWED"} for which we refer to
  \cite{Marteau (2009)} and \cite{Halpin (2014)}.

  Each method can be controlled with the following parameters:

  \tabular{ll}{
    method \tab parameters \cr
    ------------------ \tab ---------------------------------\cr
    \verb{OM} \tab \verb{sm, indel, norm} \cr
    \verb{OMloc} \tab \verb{sm, expcost, context, norm} \cr
    \verb{OMslen} \tab \verb{sm, indel, link, h, norm} \cr
    \verb{OMspell} \tab \verb{sm, indel, norm, tpow, expcost, norm} \cr
    \verb{OMstran} \tab \verb{sm, indel, transindel, otto, previous, add.column, norm} \cr
    \verb{HAM, DHD} \tab \verb{sm, norm}\cr
    \verb{CHI2} \tab \verb{breaks, step, overlap, norm, weighted, global.pdotj, norm}\cr
    \verb{EUCLID} \tab \verb{breaks, step, overlap, norm}\cr
    \verb{LCS, LCP, RLCP} \tab \verb{norm}\cr
    \verb{NMS} \tab \verb{prox, kweights}\cr
    \verb{NMSMST} \tab \verb{kweights, tpow}\cr
    \verb{SVRspell} \tab \verb{prox, kweights, tpow}\cr
    \verb{TWED} \tab \verb{sm, (indel), h, nu, norm} \cr
    ------------------ \tab ---------------------------------
  }

  \code{"LCS"} is \code{"OM"} with a substitution cost of 2 (\code{sm = "CONSTANT",
  cval = 2}) and an \code{indel} of \code{1.0}. \code{"HAM"} is \code{"OM"} without
  indels. \code{"DHD"} is \code{"HAM"} with specific substitution costs at each
  position.

  \code{"HAM"} and \code{"DHD"} apply only to sequences of equal length.

  For \code{"TWED"}, the (single) indel serves only for empty sequences.
  The distance to an empty sequence is set as \eqn{n*}\code{indel}, where \eqn{n} is
  the length of the non empty sequence. By default (\code{indel="auto"}), indel is set
  as \code{2 * max(sm) + nu + h}.

  %%Currently, \code{"OM"} works only with sequences of equal lengths.

  When \code{sm = NULL}, the substitution-cost matrix is automatically created
  for \code{"HAM"} with a single substitution cost of 1 and for \code{"DHD"} with
  the costs derived from the transition rates at the successive positions, i.e. with
  \code{sm = "TRATE"}.

  Some distances can optionally be normalized by means of the \code{norm} argument.
  Let \eqn{d} be the distance, \eqn{m} the maximum possible of the distance
  given the lengths \eqn{p} and \eqn{q} of the two sequences, and \eqn{k} the
  length of the longer sequence. Normalization \code{"maxlength"} is \eqn{d/k}
  (Abbott's normalization), \code{"gmean"} is \eqn{1-(m-d)/(p*q)^.5} (Elzinga's
  normalization), \code{"maxdist"} is \eqn{d/m}, and "YujianBo" is \eqn{2*d/(m+d)}.
  For more details, see  \cite{Gabadinho et al. (2009, 2011)}.
  Actually, to avoid negative outcomes, the length \eqn{p}, \eqn{q}, and \eqn{k} are
  set as (max) indel times the corresponding length. For some distances, \eqn{m} is
  only a possibly non-reachable upper bound.

  When \code{norm="auto"}, \code{"gmean"} is applied to \code{"LCS"},
  \code{"LCP"} and \code{"RLCP"} distances, \code{"maxlength"} is applied to \code{"OM"}, \code{"HAM"}
  and \code{"DHD"}, and the normalization "YujianBo" of \cite{Yujian and Bo (2007)} that preserves the
  triangle inequality is used in the other cases except \code{"CHI2"} and \code{"EUCLID"}.
  For the latter two, the square of the
  distances are normalized by the number of intervals and the maximal distance
  on each interval. Note that for 'CHI2' the maximal distance on each interval
  depends on the state distribution on the interval.

  When sequences contain gaps and the \code{left = NA}, \code{gaps = NA}, or \code{right = NA}
  option was passed to
  \code{\link{seqdef}} (i.e. when there are non deleted missing values), the
  \code{with.missing} argument should be set as \code{TRUE}. If left as
  \code{FALSE} the function stops when it encounters a gap. This is to make the
  user aware that there are gaps in the sequences. For methods that need an
  \code{sm} value, \code{seqdist} expects a substitution-cost matrix with a row
  and a column entry for the missing state (symbol defined with the \code{nr}
  option of \code{\link{seqdef}}). Substitution-cost matrices returned by
  \code{\link{seqcost}} (and so \code{\link{seqsubm}}) include these additional
  entries when the function is called with \code{with.missing = TRUE}. More
  details on how to compute distances with sequences containing gaps can be
  found in \cite{Gabadinho et al. (2009)}.
}
\value{
  When \code{refseq} is \code{NULL} (default), the whole matrix of pairwise
  distances between sequences or, if \code{full.matrix = FALSE},
  the corresponding \code{dist} object of pairwise distances between sequences.

  When \code{refseq} is a \code{list} of two sets of indexes, the matrix
  of distances from the first set of sequences (rows) to the second set (columns).

  Otherwise, a vector with distances between the sequences in the
  state sequence object and the reference sequence specified with \code{refseq}.
}
\references{
  Studer, M. and G. Ritschard (2016), "What matters in differences between life
  trajectories: A comparative review of sequence dissimilarity measures",
  \emph{Journal of the Royal Statistical Society, Series A}. \bold{179}(2),
  481-511, \doi{10.1111/rssa.12125}

  Studer, M. and G. Ritschard (2014). "A Comparative Review of Sequence
  Dissimilarity Measures". \emph{LIVES Working Papers}, \bold{33}. NCCR LIVES,
  Switzerland, \doi{10.12682/lives.2296-1658.2014.33}

  Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011). Analyzing and
  Visualizing State Sequences in R with TraMineR. \emph{Journal of Statistical
  Software} \bold{40}(4), 1--37.

  Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009). Mining
  Sequence Data in \code{R} with the \code{TraMineR} package: A user's guide.
  Department of Econometrics and Laboratory of Demography, University of Geneva

  Halpin, B. (2014). Three Narratives of Sequence Analysis, in Blanchard, P.,
  Bühlmann, F. and Gauthier, J.-A. (Eds.) \emph{Advances in Sequence Analysis:
  Theory, Method, Applications}, Vol 2 of Series \emph{Life Course Research and
  Social Policies}, pages 75--103, Heidelberg: Springer. \doi{10.1007/978-3-319-04969-4_5}

  Marteau, P.-F. (2009). Time Warp Edit Distances with Stiffness Adjustment for
  Time Series Matching. \emph{IEEE Transactions on Pattern Analysis and Machine
  Intelligence}, \bold{31}(2), 306--318. \doi{10.1109/TPAMI.2008.76}

  Yujian, L. and Bo, L. (2007). A normalized Levenshtein distance metric.
  \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence},
  \bold{29}(6), 1091--1095. \doi{10.1109/TPAMI.2007.1078}

  See also all references in \cite{Studer and Ritschard (2014, 2016)}
}
\author{
  Matthias Studer, Gilbert Ritschard, Pierre-Alexandre Fonta, Alexis Gabadinho, Nicolas S. Müller.
  .
}
\seealso{
  \code{\link{seqcost}}, \code{\link{seqsubm}}, \code{\link{seqdef}}, and \code{\link{seqMD}} for
  multidomain (multichannel) distances using the cost additive trick.
}
\examples{
## =========================
## Examples without missings
## =========================

## Defining a sequence object with columns 10 to 25
## of a subset of the 'biofam' data set
data(biofam)
biofam.seq <- seqdef(biofam[501:600, 10:25])

## OM distances using the vector of indels and substitution
## costs derived from the estimated state frequencies
costs <- seqcost(biofam.seq, method = "INDELSLOG")
biofam.om <- seqdist(biofam.seq, method = "OM",
                     indel = costs$indel, sm = costs$sm)

## OM between sequences of transitions
biofam.omstran <- seqdist(biofam.seq, method = "OMstran",
                     indel = costs$indel, sm = costs$sm,
                     otto=.3, transindel="subcost")

## Normalized LCP distances
biofam.lcp.n <- seqdist(biofam.seq, method = "LCP",
                        norm = "auto")

## Normalized LCS distances to the most frequent sequence
biofam.dref1 <- seqdist(biofam.seq, method = "LCS",
                        refseq = 0, norm = "auto")

## LCS distances to an external sequence
ref <- seqdef(as.matrix("(0,5)-(3,5)-(4,6)"), informat = "SPS",
              alphabet = alphabet(biofam.seq))
biofam.dref2 <- seqdist(biofam.seq, method = "LCS",
                        refseq = ref)

## LCS distances between two subsets of sequences
set1 <- 1:10
set2 <- 31:36
biofam.dref2 <- seqdist(biofam.seq, method = "LCS",
                        refseq = list(set1,set2))


## Chi-squared distance over the full observed timeframe
biofam.chi.full <- seqdist(biofam.seq, method = "CHI2",
                           step = max(seqlength(biofam.seq)))

## Chi-squared distance over successive overlapping
## intervals of length 4
biofam.chi.ostep <- seqdist(biofam.seq, method = "CHI2",
                            step = 4, overlap = TRUE)


## ======================
## Examples with missings
## ======================
data(ex1)
## Ignore empty row 7
ex1.seq <- seqdef(ex1[1:6, 1:13])

## OM with indel and substitution costs based on
## log of inverse state frequencies
costs.ex1 <- seqcost(ex1.seq, method = "INDELSLOG",
                     with.missing = TRUE)
ex1.om <- seqdist(ex1.seq, method = "OM",
                  indel = costs.ex1$indel, sm = costs.ex1$sm,
                  with.missing = TRUE)

## Localized OM
ex1.omloc <- seqdist(ex1.seq, method = "OMloc",
                     sm = costs.ex1$sm, expcost=.1, context = .4,
                     with.missing = TRUE)

## OMspell with a scalar indel
indel <- max(costs.ex1$indel)
## OM of spells
ex1.omspell <- seqdist(ex1.seq, method = "OMspell",
                       indel = indel, sm = costs.ex1$sm,
                       with.missing = TRUE)

## Distance based on number of matching subsequences
ex1.nms <- seqdist(ex1.seq, method = "NMS",
                   with.missing = TRUE)

## Using the sequence vectorial representation metric
costs.fut <- seqcost(ex1.seq, method = "FUTURE", lag = 4,
                     proximities = TRUE, with.missing = TRUE)
ex1.svr <- seqdist(ex1.seq, method = "SVRspell",
                   prox = costs.fut$prox, with.missing = TRUE)
}
\keyword{Dissimilarity measures}
