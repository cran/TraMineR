\encoding{latin1}
\name{seqformat}
\alias{seqformat}
\title{Translation between sequence formats}
\description{
Translate a sequence data set from one format to another.
}
\details{
The 'seqformat' function is used to convert data from one format to another. The input data is first converted into the STS format and then converted to the output format. Depending on input and output formats, some information can be lost due to the steps in the conversion process. The output is a matrix, NOT a sequence object to be passed to TraMineR functions for plotting and mining sequences (use the \code{\link{seqdef}} function therefore). 
See \cite{Gabadinho et al. (2008)} for more details on longitudinal data formats and translation between them.
}

\usage{
seqformat(data, var=NULL, id=NULL, 
	from, to, compressed=FALSE,
	nrep=NULL, tevent, stsep, covar=NULL,
	SPS.in=list(xfix="()", sdsep=","),
	SPS.out=list(xfix="()", sdsep=","),
	begin=NULL, end=NULL, status=NULL, 
	process=TRUE, pdata=NULL, pvar=NULL, limit=100, overwrite=TRUE, 
	fillblanks=NULL, tmin=NULL, tmax=NULL)
}

\arguments{
  \item{data}{a data frame or matrix containing sequence data.}

  \item{var}{the list of columns containing the sequences. Defaut to \code{NULL}, ie all the columns. 
Whether the sequences are in the compressed (character strings) or extended format is automatically detected by counting the number of columns.}

   \item{id}{column containing the identification numbers for the sequences. When using \code{SPELL} format as input, this identification number is mandatory, in order to identify all spells belonging to each individual in the data set.}

   \item{from}{format of the original data. Avalaible formats are: STS, SPS, SPELL. If \code{data} is a sequence object, format is automatically set to \code{STS}.}

   \item{to}{format of the output data. Avalaible formats are: STS, SPS, SRS, TSE}
	
	\item{compressed}{if TRUE and output format is one of STS, SPS or DSS, the output sequences are compressed into character strings} 

   \item{nrep}{number of previous states replicated, for the 'SRS' format}

   \item{tevent}{when converting to time-stamped-event (TSE) format, a matrix of size 'ns' * 'ns' where 'ns' is the number of distinct states appearing in the sequences must be given. In this matrix, the cell a,b contains all events associated with a transition from state a to state b.}

  \item{stsep}{the character used as separator in the original data if input format is separated character strings. By default a separator is searched for with the \code{\link{seqfcheck}} function.}

  \item{covar}{the list of columns containing associated covariates to be included in the output data frame. If to='SRS' is choosed, the covariates are replicated accross each row. Default to NULL. }

	\item{SPS.in}{a list with the characters used as prefix/suffix and state/duration separator for each state duration couple if input data contains sequences in SPS format. Set the xfix element of the list to "" if there are no pre-suf-fix.} 

	\item{SPS.out}{a list with the characters used as prefix/suffix and state/duration separator to be used for each state duration couple if output is in SPS format. Set the xfix element of the list to "" if there are no pre-suf-fix.} 

	\item{begin}{when converting from SPELL, the column with the beginning of the spell}

	\item{end}{when converting from SPELL, the column with the end of the spell}

	\item{status}{when converting from SPELL, the column with the status}

	\item{process}{when converting from SPELL, create sequences on a process time axis. If set to false, create sequences on a calendar time axis.}

	\item{pdata}{when converting from SPELL, name of the data frame containing the individual 'birth' time, that is, the entering time from which the process time will be computed. The data must contain two columns: an id to match the birth time with SPELL data and a 'birth' time.}

	\item{pvar}{names or numbers of the columns containing the individual identification number and the 'birth' time in pdata.}

	\item{limit}{when converting from SPELL, size of the resulting dataframe when creating age sequences (by default goes from age 1 to age 100)}

	\item{overwrite}{when converting from SPELL, if overwrite is set to TRUE, the most recent episode overwrites the older one if they overlap each other. If set to false, the most recent episode starts from the end of the previous one.}

	\item{fillblanks}{when converting from SPELL, if fillblanks is not NULL, gaps between episodes are filled with any character given as argument.} 

	\item{tmin}{when converting from SPELL, if sequences are to be defined on a calendar time axis is, defines the starting time of the axis. If set to NULL, the minimum time is taken from the 'begin' column in the data.}

	\item{tmax}{when converting from SPELL, if year sequences are wanted, defines the ending year of the dataframe. If set to NULL, it is guessed from the data (not very accurately).}
}

\value{a data frame}

\references{
Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2008). Mining Sequence Data in \code{R} with \code{TraMineR}: A user's guide. \emph{Department of Econometrics and Laboratory of Demography, University of Geneva}.
}

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

