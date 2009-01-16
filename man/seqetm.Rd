\name{seqetm}
\alias{seqetm}
\title{Creating event transition matrix}
\description{
This function automatically create an event transition matrix to transform state sequences data in TSE format (transitions).
}
\usage{
seqetm(seq, method = "transition", use.labels = TRUE, 
  sep = ">", bp = "", ep = "end")
}
\arguments{
  \item{seq}{The state sequence to create transitions from}
  \item{method}{The method to use. Actually, "transition","state" and "period" are available.}
  \item{use.labels}{If true, transitions names are build from sequence labels instead of sequence alphabet.}
  \item{sep}{Separator to use between period in "transition" method.}
  \item{bp}{Prefix to signal begining of a period}
  \item{ep}{Prefix to signal end of a period}
}
\details{
  Three methods are implemented and can be set through the method parameter:
  'transition' automaticaly generates a distinct event for each transition and the begining of the sequence;
  'state' generates one event per state, in order to analyse state sequence using event sequences methods;
  'period' automaticaly generates two events for each transition: one for the end of the previous period, and the other one for the beginning of the new period.
}
\value{
The event matrix.
}

\seealso{See Also \code{\link{seqformat}} for conversion to TSE format and \code{\link{seqecreate}} for handling transition data.}
\examples{
## Creating a sequence object with the columns 13 to 24
## in the 'actcal' example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24,
        labels=c("FullTime", "PartTime", "LowPartTime", "NoWork"))
## Creating a transition matrix, one event per transition
seqetm(actcal.seq,method = "transition")

## Creating a transition matrix, two events per transition
seqetm(actcal.seq,method = "state")

## Creating a transition matrix, two events per transition
seqetm(actcal.seq,method = "period")

## Creating a transition matrix, considering the begining
## as a special event, two events per transition
seqetm(actcal.seq,method = "period", bp="begin")
}

\keyword{misc}
