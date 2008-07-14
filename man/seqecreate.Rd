\name{seqecreate}
\alias{seqecreate}
\title{Create event sequence objects.}
\description{
  Create event sequence objects from the given input.
}
\details{
  The input format consist in three list of the same size idpers, timestamp and one event. Several events can appear at the same time.
  This correspond to the TSE format described in \code{\link{seqformat}}.
  The result can then be used in other methods, such as \code{\link{seqefsub}} or \code{\link{seqeapplysub}}.
%  A subsequence can be created using seqecreatesub, in this case, timestamp are only used to detect simultaneous events (same time).
}
\usage{
  seqecreate(id, timestamp, event, endEvent=NULL)
%  seqecreatesub(timestamp, event)
}
\arguments{
  \item{id}{The sequence id as an integer.}
  \item{timestamp}{Dates of events as a double.}.
  \item{event}{The event that appeared at the specified date.}
  \item{endEvent}{If specified this event will be considered as a flag for the end of observation time (total length of event sequences).}
}
\seealso{
 \code{\link{seqformat}} for TSE format information
 \code{\link{seqefsub}} for searching frequent subsequences.
 \code{\link{seqeapplysub}} for counting subsequences and more.
 \code{\link{seqelength}} about length (observation time) of event sequences.
}
\examples{
## Loading data
data(actcal.tse)
## Creating sequences
actcal.seqe <- seqecreate(actcal.tse$id, actcal.tse$time, 
	actcal.tse$event)
##printing sequences
actcal.seqe[1:10]
}
\keyword{manip}
