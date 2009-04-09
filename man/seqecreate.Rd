\name{seqecreate}
\alias{seqecreate}
\title{Create event sequence objects.}
\description{
  Create an event sequence object from the given input.
}
\details{
  There are several ways to create an event sequence object. The first one is by providing the events in TSE format (see \code{\link{seqformat}}), i.e. by providing three paired lists: id, timestamp and event, such that each triplet (id, timestamp, event) defines the event that occurs at time timestamp for case id. Several events at the same time for a same id are allowed. The lists can be provided with the arguments \code{id}, \code{timestamp} and \code{event}. An alternative is by providing a data frame as \code{data} argument in which case the function takes the required information from the "id", "timestamp" and "event" columns of that data frame. 
  
  The other way is to pass a state sequence object (as \code{data} argument) and to perfom an automatic state-to-event conversion. The simplest way to make a conversion is by means of a predefined method (see \code{\link{seqetm}}), such as "transition" (one distinct event per possible transition), "state" (one event when entering a new state) and "period" (a pair of events, one start-state event and one end-state event for each found transition).  For a more customized conversion, you can specify a transition matrix in the same way as in \code{\link{seqformat}}. Function \code{\link{seqetm}} can help you in creating your transition matrix. 
   
  The resulting event sequence object can then be used in other 'seqe' methods, such as \code{\link{seqefsub}} or \code{\link{seqeapplysub}}.
%  A subsequence can be created using seqecreatesub, in this case, timestamp are only used to detect simultaneous events (same time).
}
\usage{
  seqecreate(data = NULL, id = NULL, timestamp = NULL, event = NULL, 
    endEvent = NULL, tevent =  "transition", use.labels=TRUE)
%  seqecreatesub(timestamp, event)
}
\arguments{
  \item{data}{A state sequence object (see \code{\link{seqdef}}) or a data frame }
  \item{id}{Concerned sequence id's (integer), that is the 'id' column of the TSE format (ignored if data argument is provided).}
  \item{timestamp}{Time (double) at which events occur, that is the 'timestamp' column of the TSE format (ignored if data argument is provided).}
  \item{event}{Events that occurred at the specified time stamps, that is the 'event' column of the TSE format (ignored if data argument is provided).}
  \item{endEvent}{If specified this event will be considered as a flag for the end of observation time (total length of event sequences).}
  \item{tevent}{If \code{data} is a state sequence object either a transition matrix or a method to generate it (see \code{\link{seqetm}})}
  \item{use.labels}{If TRUE, transitions names are built from state labels rather than from the sequence alphabet.}
}
\seealso{
 \code{\link{seqformat}} for converting between sequence formats,
 \code{\link{seqefsub}} for searching frequent subsequences,
 \code{\link{seqecmpgroup}} to search for discriminant subsequences,
 \code{\link{seqeapplysub}} for counting subsequence occurrences and more,
 \code{\link{seqelength}} about length (observation time) of event sequences,
 \code{\link{seqdef}} to create a state sequence object.
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
