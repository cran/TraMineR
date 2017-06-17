% Author for TraMineR 2: Pierre-Alexandre Fonta (2016-2017)

\name{seqformat}
\alias{seqformat}
\title{Conversion between sequence formats}
\description{
  Convert a sequence data set from one format to another.
}
\usage{
seqformat(data, var = NULL, from, to, compress = FALSE, nrep = NULL, tevent,
  stsep = NULL, covar = NULL, SPS.in = list(xfix = "()", sdsep = ","),
  SPS.out = list(xfix = "()", sdsep = ","), id = 1, begin = 2, end = 3,
  status = 4, process = TRUE, pdata = NULL, pvar = NULL, limit = 100,
  overwrite = TRUE, fillblanks = NULL, tmin = NULL, tmax = NULL, missing = "*",
  with.missing = TRUE, compressed, nr)
}
\arguments{
  \item{data}{
    Data Frame, Matrix, or State Sequence Object.
    The data to use.

    A data frame or a matrix with sequence data in one or more columns when
    \code{from = "STS"} or \code{from = "SPS"}. If sequence data are in a single
    column, they are assumed to be in the compressed form (see \code{stsep}).

    A data frame with sequence data in one or more columns when \code{from = "SPELL"}.
    If sequence data are not in four columns with the order individual ID, spell
    start time, spell end time, and spell state status, use \code{var} or
    \code{id} / \code{begin} / \code{end} / \code{status}.

    A state sequence object when \code{from = "STS"} or \code{from} is not specified.
  }
  \item{var}{
    \code{NULL}, List of Integers or Strings.
    Default: \code{NULL}.
    The indexes or the names of the columns with the sequence data in \code{data}.
    If \code{NULL}, all columns are considered.
  }
  \item{from}{
    String.
    The format of the input sequence data.
    It can be \code{"STS"}, \code{"SPS"}, or \code{"SPELL"}. It is not needed
    if \code{data} is a state sequence object.
  }
  \item{to}{
    String.
    The format of the output data.
    It can be \code{"STS"}, \code{"DSS"}, \code{"SPS"}, \code{"SRS"}, \code{"TSE"},
    or \code{"SPELL"}.
  }
  \item{compress}{
    % to STS, DSS, SPS
    Logical.
    Default: \code{FALSE}.
    When \code{to = "STS"}, \code{to = "DSS"}, or \code{to = "SPS"}, should the
    sequences (row vector of states) be concatenated into strings? See
    \code{\link{seqconc}}.
   }
  \item{nrep}{
    % to SRS
    Integer.
    The number of shifted replications when \code{to = "SRS"}.
  }
  \item{tevent}{
    % to TSE
    Matrix.
    The transition-definition matrix when \code{to = "TSE"}.
    It should be of size \eqn{d * d} where \eqn{d} is the number of distinct
    states appearing in the sequences. The cell \eqn{(i,j)} lists the events
    associated with a transition from state \eqn{i} to state \eqn{j}. It can be
    created with \code{\link{seqetm}}.
  }
  \item{stsep}{
    % from STS, SPS
    \code{NULL}, Character.
    Default: \code{NULL}.
    The separator between states in the compressed form (strings) when
    \code{from = "STS"} or \code{from = "SPS"}.
    If \code{NULL}, \code{\link{seqfcheck}} is called for detecting automatically
    a separator among "-" and ":". Other separators must be specified explicitly.
    See \code{\link{seqdecomp}}.
  }
  \item{covar}{
    % to SRS
    List of Integers or Strings.
    The indexes or the names of additional columns in \code{data} to include as
    covariates in the output when \code{to = "SRS"}.
    The covariates are replicated across the shifted replicated rows.
  }
  \item{SPS.in}{
    % from SPS
    List.
    Default: \code{list(xfix = "()", sdsep = ",")}.
    The specifications for the state-duration couples in the input data when
    \code{from = "SPS"}.
    The first specification, \code{xfix}, specifies the prefix/suffix character.
    Use a two-character string if the prefix and the suffix differ. Use
    \code{xfix = ""} when no prefix/suffix are present. The second
    specification, \code{sdsep}, specifies the state/duration separator.
  }
	\item{SPS.out}{
    % to SPS
    List.
    Default: \code{list(xfix = "()", sdsep = ",")}.
    The specifications for the state-duration couples in the output data when
    \code{to = "SPS"}.
    See \code{SPS.in} above.
	}
	\item{id}{
    % from SPELL
    % to TSE
    \code{NULL}, Integer, String, List of Integers or Strings.
    Default: \code{1}.

    When \code{from = "SPELL"}, the index or the name of the column containing
    the individual IDs in \code{data} (after \code{var} filtering).

    When \code{to = "TSE"}, the index or the name of the column containing
    the individual IDs in \code{data} (after \code{var} filtering) or the unique
    individual IDs. If \code{id} is not manually specified, \code{id} is set as
    \code{NULL} for backward compatibility with TraMineR 1.8-13 behaviour. If
    \code{id} is manually or automatically set as \code{NULL}, the original
    individual IDs are ignored and replaced by the index of the sequence in the
    input data.

    When \code{from = "SPELL"} and \code{to = "TSE"}, the index or the name of
    the column containing the individual IDs in \code{data} (after \code{var}
    filtering). The TSE output will use the original individual IDs.
  }
	\item{begin}{
    % from SPELL
    Integer or String.
    Default: \code{2}.
    The index or the name of the column containing the spell start times in
    \code{data} (after \code{var} filtering) when \code{from = "SPELL"}.
  }
	\item{end}{
	  % from SPELL
    Integer or String.
    Default: \code{3}.
    The index or the name of the column containing the spell end times in
    \code{data} (after \code{var} filtering) when \code{from = "SPELL"}.
  }
	\item{status}{
		% from SPELL
    Integer or String.
    Default: \code{4}.
    The index or the name of the column containing the spell status in
    \code{data} (after \code{var} filtering) when \code{from = "SPELL"}.
	}
	\item{process}{
	  % from SPELL
    Logical.
    Default: \code{TRUE}.
    When \code{from = "SPELL"}, if \code{TRUE}, create sequences on a process
    time axis, if \code{FALSE}, create sequences on a calendar time axis.
    %If times are ages in \code{data}, only \code{process = FALSE} is currently
    %implemented. %% gr: I do not understand this and have to check.
    
    This \code{process} argument as well as the associated \code{pdata} and \code{pvar}
    arguments are intended for \code{data} containing spell data with calendar begin 
    and end times. When those times are ages, use \code{process = FALSE} with 
    \code{pdata=NULL} to use those ages as process times. Option \code{process = TRUE} 
    does currently not work for age times. 
  }
  \item{pdata}{
    % from SPELL
    % to SPELL
    \code{NULL}, \code{"auto"}, or Data Frame.
    Default: \code{NULL}.

    If \code{NULL}, the start and end times of each spell are supposed to be, if
    \code{process = TRUE}, ages, if \code{process = FALSE}, years when
    \code{from = "SPELL"}.

    If \code{"auto"}, ages are computed using the start time of the first spell
    of each individual as her/his birthdate when \code{from = "SPELL"} and
    \code{process = TRUE}.

    A data frame containing the ID and the birth time of the individuals when
    \code{from = "SPELL"} or \code{to = "SPELL"}. Use \code{pvar} to specify
    the column names. The ID is used to match the birth time of individuals
    with the sequence data. The birth time is the start time from which the time
    axis will be computed. It is used to compute \code{tmin} and to guess
    \code{tmax}, if there are \code{NULL}, when \code{from = "SPELL"} and
    \code{process = FALSE}.
  }
  \item{pvar}{
    % from SPELL
    % to SPELL
    List of Integers or Strings.
    The indexes or names of the columns of the data frame \code{pdata} that
    contain the ID and the birth time of the individuals in that order.
  }
  \item{limit}{
    % from SPELL
    Integer.
    Default: \code{100}.
    The maximum age of age sequences when \code{from = "SPELL"} and \code{process = TRUE}.
    Age sequences will be considered to start at 1 and to end at \code{limit}.
  }
	\item{overwrite}{
    % from SPELL
    Logical.
    Default: \code{TRUE}.
    When \code{from = "SPELL"}, if \code{TRUE}, the most recent episode
    overwrites the older one when they overlap each other, if \code{FALSE},
    in case of overlap, the most recent episode starts after the end of the
    previous one.
  }
	\item{fillblanks}{
    % from SPELL
    Character.
    The value to fill gaps between episodes when \code{from = "SPELL"}.
  }
	\item{tmin}{
    % from SPELL
    \code{NULL}, Integer.
    Default: \code{NULL}.
    The start time of the axis when \code{from = "SPELL"} and \code{process = FALSE}.
    If \code{NULL}, the value is the minimum of the spell start times
    (see \code{begin}) or the minimum of the birth time of the individuals
    (see \code{pdata} when it is a data frame and \code{process = FALSE}).
  }
	\item{tmax}{
    % from SPELL & !process
    \code{NULL}, Integer.
    Default: \code{NULL}.
    The end time of the axis when \code{from = "SPELL"} and \code{process = FALSE}.
    If \code{NULL}, the value is the maximum of the spell end times (see \code{end})
    or the sum of the maximum of the spell end times and of the maximum of the
    birth time of the individuals (see \code{pdata} when it is a data frame and
    \code{process = FALSE}).
  }
  \item{missing}{
    String.
    Default: \code{"*"}.
    The code for missing states in \code{data}.
    It will be replaced by \code{NA} in the output data. The code is obtained
    from the attribute \code{nr} when \code{data} is a state sequence object
    (see \code{\link{seqdef}}).
  }
  \item{with.missing}{
    % to SPELL
    Logical.
    Default: \code{TRUE}.
    When \code{to = "SPELL"}, should the spells of missing states be included?
  }
  \item{compressed}{
    Deprecated. Use \code{compress} instead.
  }
  \item{nr}{
    Deprecated. Use \code{missing} instead.
  }
}
\details{
  The \code{seqformat} function is used to convert data from one format to
  another. The input data is first converted into the STS format and then
  converted to the output format. Depending on input and output formats, some
  information can be lost in the conversion process. The output is a matrix or
  a data frame, NOT a sequence object to be passed to TraMineR functions for
  plotting and mining sequences (use the \code{\link{seqdef}} function for that).
  See \cite{Gabadinho et al. (2009)} and \cite{Ritschard et al. (2009)} for more
  details on longitudinal data formats and converting between them.
  
  When data are in \code{"SPELL"} format, the begin and end times expected
  to be positions in the sequences. Therefore they should be strictly positive integers.
}
\value{
  A data frame for \code{SRS}, \code{TSE}, and \code{SPELL}, a matrix otherwise.
}
\references{
  Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009). Mining
  Sequence Data in \code{R} with the \code{TraMineR} package: A user's guide.
  Department of Econometrics and Laboratory of Demography, University of Geneva.

  Ritschard, G., A. Gabadinho, M. Studer and N. S. Müller. Converting between
  various sequence representations. in Ras, Z. & Dardzinska, A. (ed.)
  \emph{Advances in Data Management}, Springer, 2009, 223, 155-175.
}
\author{
  Alexis Gabadinho, Pierre-Alexandre Fonta, Nicolas S. Müller, Matthias Studer
  (and Gilbert Ritschard for the help page).
}
\seealso{\code{\link{seqdef}}}
\examples{
## ========================================
## Examples with raw STS sequences as input
## ========================================

## Loading a data frame with sequence data in the columns 13 to 24
data(actcal)

## Converting to SPS format
actcal.SPS.A <- seqformat(actcal, 13:24, from = "STS", to = "SPS")
head(actcal.SPS.A)

## Converting to compressed SPS format with no
## prefix/suffix and with "/" as state/duration separator
actcal.SPS.B <- seqformat(actcal, 13:24, from = "STS", to = "SPS",
  compress = TRUE, SPS.out = list(xfix = "", sdsep = "/"))
head(actcal.SPS.B)

## Converting to compressed DSS format
actcal.DSS <- seqformat(actcal, 13:24, from = "STS", to = "DSS",
  compress = TRUE)
head(actcal.DSS)


## ==============================================
## Examples with a state sequence object as input
## ==============================================

## Loading a data frame with sequence data in the columns 10 to 25
data(biofam)

## Limiting the number of considered cases to the first 20
biofam <- biofam[1:20, ]

## Creating a state sequence object
biofam.labs <- c("Parent", "Left", "Married", "Left/Married",
  "Child", "Left/Child", "Left/Married/Child", "Divorced")
biofam.short.labs <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")
biofam.seq <- seqdef(biofam, 10:25, alphabet = 0:7,
  states = biofam.short.labs, labels = biofam.labs)

## Converting to SPELL format
bf.spell <- seqformat(biofam.seq, from = "STS", to = "SPELL",
  pdata = biofam, pvar = c("idhous", "birthyr"))
head(bf.spell)


## ======================================
## Examples with SPELL sequences as input
## ======================================

## Loading two data frames: bfspell20 and bfpdata20
## bfspell20 contains the first 20 biofam sequences in SPELL format
## bfpdata20 contains the IDs and the years at which the
## considered individuals were aged 15
data(bfspell)

## Converting to SPS format
bf.sts <- seqformat(bfspell20, from = "SPELL", to = "STS",
  id = "id", begin = "begin", end = "end", status = "states",
  process = TRUE, pdata = bfpdata20, pvar = c("id", "when15"),
  limit = 16)
names(bf.sts) <- paste0("a", 15:30)
head(bf.sts)
}
\keyword{Data handling}
