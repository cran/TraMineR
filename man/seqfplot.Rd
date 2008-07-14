\name{seqfplot}
\alias{seqfplot}
\title{Graphic presenting the frequency of sequences}
\description{
  The sequences are presented as stacks of successive states, ordered by their relative frequency in the dataset.
}
\usage{
 seqfplot(seqdata, tlim=10, title=NULL, cpal=NULL, pbarw=FALSE,
	withlegend=TRUE, ltext=NULL, xtlab=NULL, bmar=1)
}
\arguments{
	\item{seqdata}{a sequence object created with the \code{\link{seqdef}} function.}
	\item{tlim}{if tlim>0, shows only the 'tlim' most frequent sequences. Default to 10.}
	\item{title}{title for the graphic. Default to NULL.}
	\item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of distinct states. By default, the 'cpal' attribute of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}
	\item{pbarw}{if pbarw=TRUE, the width of the bars are proportional to the sequence frequency in the dataset. Default to FALSE.}
  \item{withlegend}{if TRUE, a legend is plotted. Default to TRUE.}
  \item{ltext}{optional description of the states to appear in the legend. Must be a vector of character strings with number of elements equal to the number of distinct states. If unspecified, the 'label' attributes of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}
  \item{xtlab}{optional labels for the x axis ticks labels. If unspecified, the column names of the 'seqdata' sequence object are used (see \code{\link{seqdef}}).}
 \item{bmar}{base margin for the top and bottom of the graphic, that is, the amount of white space surrouding it.}
}
\seealso{
 \code{\link{seqtab}}.
}
\examples{
## Plot of the 10 most frequent sequences 
## with equal bar width
data(biofam)

biofam.seq <- seqdef(biofam,10:25) 
seqfplot(biofam.seq)

## The same plot with
## proportional bar width
seqfplot(biofam.seq, pbarw=TRUE)
}
\keyword{hplot}
