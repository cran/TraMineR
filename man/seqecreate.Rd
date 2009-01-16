\name{seqecreate}
\alias{seqecreate}
\title{Create event sequence objects.}
\description{
  Create event sequence objects from the given input.
}
\details{
  There are several way to create events sequences. The first one is to use the TSE format (see \code{\link{seqformat}}). It consist in three list of the same size idpers, timestamp and event. Several events can appear at the same time.
  You can set these parameters using the data argument with a data.frame. Then the function will look to an "id", "timestamp" and "event" column in the data frame. You can also set this parameters separately using the argument with the same names.
  
  The other way is to use a state sequences. Each transition will be considered as an event. You can specify the event appearing in each transition using tevent argument in the same way as in \code{\link{seqformat}}. You can build one automatically using \code{\link{seqetm}}. If \code{tevent} is set to null, then each transition will be considered as a specific event.
  You can pass a sequence using the data argument.
  
  The result can then be used in other methods, such as \code{\link{seqefsub}} or \code{\link{seqeapplysub}}.
%  A subsequence can be created using seqecreatesub, in this case, timestamp are only used to detect simultaneous events (same time).
}
\usage{
  seqecreate(data = NULL, id = NULL, timestamp = NULL, event = NULL, 
    endEvent = NULL, tevent =  "transition")
%  seqecreatesub(timestamp, event)
}
\arguments{
  \item{data}{A state sequence or a data.frame }
  \item{id}{The sequence id as an integer.}
  \item{timestamp}{Dates of events as a double.}.
  \item{event}{The event that appeared at the specified date.}
  \item{endEvent}{If specified this event will be considered as a flag for the end of observation time (total length of event sequences).}
  \item{tevent}{If \code{data} is a state sequence either a transition matrix or a method to generate it (see \code{\link{seqetm}}}
}
\seealso{
 \code{\link{seqformat}} for TSE format information
 \code{\link{seqefsub}} for searching frequent subsequences.
 \code{\link{seqecmpgroup}} to search for discriminant subsequences.
 \code{\link{seqeapplysub}} for counting subsequences and more.
 \code{\link{seqelength}} about length (observation time) of event sequences.
}
\examples{
##Starting with states sequences
##Loading data
data(biofam)
## Creating state sequences
biofam.seq <- seqdef(biofam,10:25,informat='STS')
## Creating event sequences from biofam
biofam.seqe <- seqecreate(biofam.seq)

## Loading data
data(actcal.tse)
## Creating sequences
actcal.seqe <- seqecreate(id=actcal.tse$id, timestamp=actcal.tse$time,
	event=actcal.tse$event)
##printing sequences
actcal.seqe[1:10]
## Using the data argument
actcal.seqe <- seqecreate(data=actcal.tse)

}
\keyword{manip}
