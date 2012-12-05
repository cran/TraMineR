\encoding{latin1}
\name{seqformat}
\alias{seqformat}
\title{Conversion between sequence formats}

\description{
 Convert a sequence data set from one format to another.
}

 \details{
  The \code{seqformat} function is used to convert data
from one format to another. The input data is first converted into
the STS format and then converted to the output format. Depending on
input and output formats, some information can be lost in the
conversion process. The output is a matrix, NOT a sequence object to
be passed to TraMineR functions for plotting and mining sequences
(use the \code{\link{seqdef}} function for that). See \cite{Gabadinho
et al. (2009)} and \cite{Ritschard et al. (2009)} for more details on
longitudinal data formats and converting between them.
 }

\usage{
seqformat(data, var=NULL, id=NULL,
     from, to, compressed=FALSE,
     nrep=NULL, tevent, stsep=NULL, covar=NULL,
     SPS.in=list(xfix="()", sdsep=","),
     SPS.out=list(xfix="()", sdsep=","),
     begin=NULL, end=NULL, status=NULL,
     process=TRUE, pdata=NULL, pvar=NULL,
     limit=100, overwrite=TRUE,
     fillblanks=NULL, tmin=NULL, tmax=NULL)
 }

\arguments{
  \item{data}{a data frame or matrix containing sequence data.}

  \item{var}{List of columns with the sequence data. Default is \code{NULL}, i.e., all columns.
Sequences are assumed to be in compressed form (character strings) when there is a single column and in extended form otherwise.}

   \item{id}{Column containing the 'id' of the sequences. Mandatory with \code{from="SPELL"}
in order to identify the spells of a same sequence.}

   \item{from}{Format of the input data. One of \code{"STS"}, \code{"SPS"}, \code{"SPELL"}.
If \code{data} is a sequence object, format is automatically set to \code{"STS"}.}

   \item{to}{Format for output data. One of \code{"STS"}, \code{"SPS"}, \code{"SRS"}, \code{"DSS"}, \code{"TSE"}.}

   \item{compressed}{Logical. Should \code{"STS"}, \code{"SPS"} or \code{"DSS"} output be compressed into character strings? Ignored for other output formats.}

   \item{nrep}{Number of shifted replications for output in \code{"SRS"} format.}

   \item{tevent}{Transition definition matrix for converting to time-stamped-event (\code{"TSE"}) format. Should be a matrix of size \eqn{d * d} where \eqn{d} is the number of distinct
states appearing in the sequences. In this matrix, the
cell \eqn{(i,j)} lists the events associated with a transition from state \eqn{i} to state \eqn{j}.}

  \item{stsep}{Separator character between successive elements in compressed (character strings) input data.
If \code{NULL} (default value), the \code{\link{seqfcheck}} function
is called for detecting automatically a separator among "-" and ":".
Other separators must be specified explicitly.}

  \item{covar}{When \code{from="STS"} or \code{from="SPS"}, additional column names to be included as covariates in the output data frame.
When \code{to="SRS"} the covariates are replicated across the shifted replicated rows. Default is \code{NULL}. Ignored when \code{from="SPELL"}. }

	\item{SPS.in}{List with the \code{xfix=} and \code{sdsep=} specifications for the state-duration couples in input data in \code{SPS} form. The first specification, \code{xfix}, specifies the prefix/suffix character (use a two-character string if the prefix and suffix differ and set \code{xfix=""} when no prefix/suffix are present). The second one, \code{sdsep}, specifies the state/duration separator.}

	\item{SPS.out}{List with the \code{xfix} and \code{sdsep} specifications for output in \code{SPS} format. (see argument \code{SPS.in} above.)}


	\item{begin}{When converting from \code{SPELL}, the column with
the beginning position of the spell}

	\item{end}{When converting from \code{SPELL}, the column with the
end position of the spell}

	\item{status}{When converting from \code{SPELL}, the column with
the status}

	\item{process}{Logical: When converting from
\code{SPELL}, should sequences be created on a process time axis? Default is \code{TRUE}. Set as
\code{FALSE} for creating sequences on a calendar time axis.}

	\item{pdata}{When converting from \code{SPELL} and
\code{process=TRUE}, either \code{NULL}, \code{"auto"} or the name of
the data frame containing the individual 'birth' time, that is, the
initial time from which the process time will be computed. If set as
\code{NULL} (default), the starting and ending time of each spell are
supposed to be ages. If set as \code{"auto"}, ages are computed using
the starting time of the first spell of each individual as her/his
birth date. If external birth dates are provided, the \code{pdata} data must
contain two columns: an \code{id} to match the birth time with
\code{SPELL} data and a '\code{birth}' time.}

	\item{pvar}{When \code{pdata} is a data frame, a vector of two names or numbers, the first one specifying the column with the individual 'id', and the second one the 'birth' time.}

	\item{limit}{When converting from \code{SPELL}, size of the
resulting data frame when creating age sequences (by default ranges from
age 1 to age 100)}

	\item{overwrite}{When converting from \code{SPELL}, if
\code{overwrite} is set to \code{TRUE}, the most recent episode
overwrites the older one when they overlap each other. If set to
\code{FALSE}, the most recent episode starts in case of overlap after the end of the
previous one.}

	\item{fillblanks}{When converting from \code{SPELL}, if
\code{fillblanks} is not \code{NULL}, gaps between episodes are
filled with any character given as argument.}

	\item{tmin}{When converting from \code{SPELL}, if sequences are
to be defined on a calendar time axis, it defines the starting time
of the axis. If set to \code{NULL}, the minimum time is taken from
the 'begin' column in the data.}

	\item{tmax}{When converting from \code{SPELL}, if year sequences
are wanted, defines the ending year. If set to
\code{NULL}, it is guessed from the data (not so accurately!).} }

\value{A data frame}

\references{
 Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller
(2009). Mining Sequence Data in \code{R} with the \code{TraMineR}
package: A user's guide. Department of Econometrics and Laboratory of
Demography, University of Geneva.

Ritschard, G., A. Gabadinho, M. Studer and N. S. Müller. Converting
between various sequence representations. in Ras, Z. & Dardzinska, A.
(ed.) \emph{Advances in Data Management}, Springer, 2009, 223,
155-175 }

\seealso{
 \code{\link{seqdef}}
}

\examples{
## Converting sequences into SPS format
data(actcal)
actcal.SPS.A <- seqformat(actcal,13:24, from="STS", to="SPS")
head(actcal.SPS.A)

## SPS (compressed) format with no prefix/suffix "/" as state/duration separator
actcal.SPS.B <- seqformat(actcal,13:24,
	from="STS", to="SPS", compressed=TRUE,
	SPS.out=list(xfix="", sdsep="/"))
head(actcal.SPS.B)

## Converting sequences into DSS (compressed) format
actcal.DSS <- seqformat(actcal,13:24,
	from="STS", to="DSS", compressed=TRUE)
head(actcal.DSS)
}


\author{Alexis Gabadinho, Nicolas S. Müller and Matthias Studer (with Gilbert Ritschard for the help page)}

\keyword{Data handling}
