\name{seqiplot}
\alias{seqiplot}
\title{Visualization of individual sequences.}
\description{
  Plots a graphic showing individual sequences as stacked bars. 
}
\usage{
seqiplot(seqdata, tlim=1:10, sortv=NULL, statl=NULL, 
   title=NULL, cpal=NULL, 
   withlegend=TRUE, withborder=TRUE, space=NULL, 
   ltext=NULL, slab=FALSE, xtlab=NULL, bmar=1)
}
\arguments{
  \item{seqdata}{a sequence object created with the \code{\link{seqdef}} function.}
 \item{tlim}{interval representing the range to be plotted. Default to 1:10, ie the first 10 rows. Set to 0 if you want to plot all sequences in 'seqdata'.}
 \item{sortv}{name of an optional variable used to sort the data before plotting}
  \item{statl}{list of states for the plot. Usefull for several plots to have the same colors representing the same states. If not specified, the state list for the plot is taken from seqdata.}
  \item{withlegend}{if TRUE, a legend is added to the plot. Default to FALSE.}
   \item{title}{title for the graphic. Default to NULL.}
 \item{cpal}{alternative collor palette for the states. Must be a vector of color names with number of elements equal to the number of states in seqdata. By default, colors are obtained with the RColorBrewer package (accent palette).}
  \item{withborder}{if TRUE, plots the borders of the bars cells. Default to TRUE.}
  \item{space}{sets the space between each sequence in the plot. Setting to 0 yiels no space and is usefull if plotting a large amount of sequences. Default to NULL, ie space value is automaticaly computed.}
  \item{ltext}{optional description of the states to appear in the legend. Must be a vector of character strings with number of elements equal to the number of distinct states. If unspecified, the 'label' attributes of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}
  \item{slab}{}
  \item{xtlab}{optional labels for the x axis ticks labels. If unspecified, the column names of the 'seqdata' sequence object are used (see \code{\link{seqdef}}).}
 \item{bmar}{base margin for the top and bottom of the graphic, that is, the amount of white space surrouding it.}
}
%\seealso{}
\examples{
data(actcal)
 
## Plotting the first 10 sequences in the 'actcal' example data set
seqiplot(seqdef(actcal,13:24))

## Plotting the all the sequences in the 'biofam' example data set
## using 'withborder=FALSE' and 'space=0' options to have a nicier
## plot
data(biofam)
seqiplot(seqdef(biofam,10:25),tlim=0,withborder=FALSE,space=0)
}
\keyword{hplot}
