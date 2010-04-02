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

  \item{var}{the list of columns containing the sequences. Default is \code{NULL}, i.e. all the columns.
Whether the sequences are in the compressed (character strings) or
extended format is automatically detected by counting the number of
columns.}

   \item{id}{column containing the identification numbers for the sequences.
When using \code{SPELL} format as input, this identification number
is mandatory, in order to identify all spells belonging to each
individual in the data set.}

   \item{from}{format of the original data. Available formats are: \code{"STS"}, \code{"SPS"}, \code{"SPELL"}.
If \code{data} is a sequence object, format is automatically set to
\code{"STS"}.}

   \item{to}{format of the output data. Available formats are:
\code{"STS"}, \code{"SPS"}, \code{"SRS"}, \code{"DSS"},
\code{"TSE"}.}

   \item{compressed}{if \code{TRUE} and output format is one of
\code{"STS"}, \code{"SPS"} or \code{"DSS"}, the output sequences are
compressed into character strings}

   \item{nrep}{number of previous states replicated, for the \code{"SRS"} format}

   \item{tevent}{when converting to time-stamped-event (\code{"TSE"}) format,
a matrix of size \eqn{d * d} where \eqn{d} is the number of distinct
states appearing in the sequences must be given. In this matrix, the
cell \eqn{(i,j)} contains all events associated with a transition
from state \eqn{i} to state \eqn{j}.}

  \item{stsep}{the character used as separator in the original data if input format is a vector of character strings.
If \code{NULL} (default value), the \code{\link{seqfcheck}} function
is called for detecting automatically a separator among "-" and ":".
Other separators must be specified explicitly.}

  \item{covar}{the list of columns containing associated covariates to be included in the output data frame.
If to=\code{"SRS"} is chosen, the covariates are replicated across
each row. Default is \code{NULL}. }

	\item{SPS.in}{a list with the characters used as prefix/suffix
and state/duration separator for each state duration couple if input
data contains sequences in \code{SPS} format. Set the \code{xfix}
element of the list to \code{""} if there are no pre-suf-fixes.}

	\item{SPS.out}{a list with the characters used as prefix/suffix
and state/duration separator to be used for each state duration
couple if output is in \code{SPS} format. Set the \code{xfix} element
of the list to \code{""} if there are no pre-suf-fixes.}

	\item{begin}{when converting from \code{SPELL}, the column with
the beginning position of the spell}

	\item{end}{when converting from \code{SPELL}, the column with the
end position of the spell}

	\item{status}{when converting from \code{SPELL}, the column with
the status}

	\item{process}{If \code{TRUE} (default) when converting from
\code{SPELL}, sequences are created on a process time axis. If set to
\code{FALSE}, they are created on a calendar time axis.}

	\item{pdata}{when converting from \code{SPELL} and
\code{process=TRUE}, either \code{NULL}, \code{"auto"} or the name of
the data frame containing the individual 'birth' time, that is, the
entering time from which the process time will be computed. If set to
\code{NULL} (default), the starting and ending time of each spell are
supposed to be ages. If set to \code{"auto"}, ages are computed using
the starting time of the first spell of each individual as her/his
birth date. If external birth dates are provided, the data must
contain two columns: an \code{id} to match the birth time with
\code{SPELL} data and a '\code{birth}' time.}

	\item{pvar}{names or numbers of the columns containing the
individual identification number and the 'birth' time in
\code{pdata}.}

	\item{limit}{when converting from \code{SPELL}, size of the
resulting dataframe when creating age sequences (by default goes from
age 1 to age 100)}

	\item{overwrite}{when converting from \code{SPELL}, if
\code{overwrite} is set to \code{TRUE}, the most recent episode
overwrites the older one if they overlap each other. If set to
\code{FALSE}, the most recent episode starts from the end of the
previous one.}

	\item{fillblanks}{when converting from \code{SPELL}, if
\code{fillblanks} is not \code{NULL}, gaps between episodes are
filled with any character given as argument.}

	\item{tmin}{when converting from \code{SPELL}, if sequences are
to be defined on a calendar time axis, it defines the starting time
of the axis. If set to \code{NULL}, the minimum time is taken from
the 'begin' column in the data.}

	\item{tmax}{when converting from \code{SPELL}, if year sequences
are wanted, defines the ending year of the dataframe. If set to
\code{NULL}, it is guessed from the data (not so accurately!).} }

\value{a data frame}

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

\keyword{manip}
