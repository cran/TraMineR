\name{seqdef}
\alias{seqdef}
\alias{plot.stslist}
\alias{print.stslist}
\title{Create a sequence object}
\description{
  Create a sequence object to be passed to other functions provided by the TraMineR package.
  There are specific method for plotting and printing sequences objects.
}
\usage{
seqdef(data, var=NULL, informat='STS', stsep="-", 
	  alphabet=NULL, states=NULL, 
       start=1, missing=NA, cnames=NULL,
	  cpal=NULL, labels=NULL)
}
\arguments{
  \item{data}{a data frame or matrix containing sequence data.}
  \item{var}{the list of columns containing the sequences. Defaut to NULL, ie all the columns. Whether the sequences are in the compressed (successive states in a character string) or extended format is automatically detected.}
   \item{informat}{format of the original data. Default is 'STS'. Avalaible formats are: STS, SPS, SPELL. See TraMineR manual for a description of the formats.}
   \item{stsep}{the character used as separator in the original data if input format is successive states in a character string. By default, "-".}
  \item{alphabet}{optional vector containing the alphabet (the list of all possible states). Use this option if some states in the alphabet don't appear in the data or if yopu want to reorder the states in the alphabet. The specified vector MUST contain AT LEAST all the states appearing in the data as they are labelled, and some optional states not appearing in the data. If left to NULL, the alphabet is set to the distinct states appearing in the data, in the same order as returned by the \code{\link{seqstatl}} function.}
  \item{states}{an optional vector containing the labels for the states. Must have a length equal to the number of states in the data, and the labels must be sorted according to the output of the \code{\link{seqstatl}} function.}
  \item{start}{starting time. For instance, if your sequences begin at age 15, you can specify 15. At this stage, used only for labelling column names.}
   \item{missing}{the code for the missing values appearing in the data. If specified, all cells containing this value will be replaced by NA's, the internal R code for missing values. If 'missing' is not specified, cells containing NA's are considered to be missing values.}
   \item{cnames}{optional names for the columns composing the sequence data. Those names will be used by default in the graphics as axis labels. If not specified, names are taken from the original columns names in the data.}
 \item{cpal}{an optional color palette for representing the states in the graphics. If not specified, colors palette is created with the RColorBrewer package, using the "accent" palette. Note that the maximum number of colors in this palette is 8. If the number of states in the data is greater, you have to specify your own palette. The list of available colors is displayed by the \code{\link{colors}} function. You can also use alternatively some other palettes from the RColorBrewer package.
 }
 \item{labels}{labels for the states, to appear in the graphics' legend.}
}
\details{Subscripts applied to sequence objects (eg. \code{seq[,1:5] or seq[1:10,]}) returns a sequence object with preserved (alphabet, missing) and adapted attributes (start, column names), unless only one column is selected, in which case a factor is returned.}
\value{a sequence object that can be passed as argument to other functions such as plotting functions (seqdplot, seqiplot or seqfplot), functions to compute distances (seqdist), etc...} 
\seealso{
 \code{}
}
\examples{
## Creating a sequence object with the columns 13 to 24 
## in the 'actcal' example data set
data(actcal)
actcal.seq <- seqdef(actcal,13:24,
	labels=c("> 37 hours", "19-36 hours", "1-18 hours", "no work"))

## Displaying the first 10 rows of the sequence object
actcal.seq[1:10,]

## Displaying the first 10 rows of the sequence object
## in SPS format
print(actcal.seq[1:10,], format="SPS")

## Frequency plot for the monthes June to September
seqfplot(actcal.seq[,6:9])

## Re-ordering the alphabet
actcal.seq <- seqdef(actcal,13:24,alphabet=c("B","A","D","C"))
alphabet(actcal.seq)

## Adding a state not appearing in the data to the
## alphabet
actcal.seq <- seqdef(actcal,13:24,alphabet=c("A","B","C","D","E"))
alphabet(actcal.seq)

## Adding a state not appearing in the data to the
## alphabet and changing the states labels
actcal.seq <- seqdef(actcal,13:24,
  alphabet=c("A","B","C","D","E"),
  states=c("FT","PT","LT","NO","TR"))
alphabet(actcal.seq)
actcal.seq[1:10,]
}
\keyword{manip}
