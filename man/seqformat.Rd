\name{seqformat}
\alias{seqformat}
\title{Translation between sequence formats.}
\description{
Translate a sequence data set from one format to another.
}
\usage{
seqformat(data, var=NULL, id=NULL, from=NULL, to=NULL, 
          nrep=NULL, tevent=NULL, stsep="-", covar=NULL)
}
\arguments{
  \item{data}{a data frame or matrix containing sequence data.}
  \item{var}{the list of columns containing the sequences. Defaut to NULL, ie all the columns. 
Whether the sequences are in the compressed (character strings) or extended format is automatically detected by counting the number of columns.}
   \item{id}{identification number for the sequences. When using SPELL formt as input, this identification number is mandatory, in order to identify all spells belonging to each individual in the data set.}
   \item{from}{format of the original data. Avalaible formats are: STS, SPS, SPELL}
   \item{to}{format of the output data. Avalaible formats are: STS, SPS, SPELL, SRS, TSE}
   \item{nrep}{number of previous states replicated, for the 'SRS' format}
   \item{tevent}{when converting to time-stamped-event (TSE) format, a matrix of size 'ns' * 'ns' where 'ns' is the number of distinct states appearing in the sequences must be given. In this matrix, the cell a,b contains all events associated with a transition from state a to state b.} 
  \item{stsep}{the character used as separator in the original data if input format is separated character strings. By default, "-".}
  \item{covar}{the list of columns containing associated covariates to be included in the output data frame. If to='SRS' is choosed, the covariates are replicated accross each row. Default to NULL. }
}
\details{the 'seqformat' function is used to convert data from one format to another. The input data is first converted into the STS format and then converted to the output format. Depending on input and output formats, some information can be lost due to the steps in the conversion process. The output is a matrix, NOT a sequence object to be passed to TraMineR functions for plotting and mining sequences (use the \code{\link{seqdef}} function therefore). } 
\value{a data frame}
\seealso{
 \code{\link{seqdef}}
}
\examples{
data(actcal)
actcal.SPS <- seqformat(actcal,13:24,from="STS",to="SPS")
}
\keyword{manip}
