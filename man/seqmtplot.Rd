\name{seqmtplot}
\alias{seqmtplot}
\title{Graphic presenting the mean time spent in each state of the alphabet}

\description{Graphic presenting the mean time spent in each state of the alphabet.}

\usage{
seqmtplot(seqdata, group=NULL, title=NULL,
	cpal=NULL,  
	ylab, axes="all", cex.plot=1,
	withlegend="auto", ltext=NULL, cex.legend=1,  
	use.layout=(!is.null(group) | withlegend!=FALSE), 
	legend.prop=NA, rows=NA, cols=NA, ...)
}

\arguments{
  \item{seqdata}{a sequence object as returned by the \code{\link{seqdef}} function.}

  \item{group}{Plots one plot for each level of the factor given as argument.}

  \item{title}{a character string giving the title for the graphic. Default to no title.}

   \item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of distinct states. By default, the 'cpal' attribute of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}

	\item{ylab}{An optional label for the y axis. If set to NA, no label is drawn.} 

	\item{axes}{if set to "all" (default value) x axes are drawn for each plot in the graphic. If set to 'bottom' and \code{group} is used, axes are drawn only under the plots located at the bottom of the graphic area. If FALSE, no x axis is drawn.}

	\item{cex.plot}{expansion factor for setting the size of the font for the axis labels and names. The default value is 1. Values lesser than 1 will reduce the size of the font, values greater than 1 will increase the size.}

  \item{withlegend}{set wether a legend for the states colors is plotted. The default value 'auto' sets the position of the legend automatically. Other options are 'right'. Obsolete option 'TRUE' is identical to 'auto'.}

  \item{ltext}{optional description of the states to appear in the legend. Must be a vector of character strings with number of elements equal to the number of distinct states. If unspecified, the 'label' attributes of the 'seqdata' sequence object is used (see \code{\link{seqdef}}). }

	\item{cex.legend}{expansion factor for setting the size of the font for the labels in the legend. The default value is 1. Values lesser than 1 will reduce the size of the font, values greater than 1 will increase the size.}

 \item{use.layout}{if TRUE, layout is used to arrange plots when using the group option or plotting a legend. If layout is used, the standard 'par(mfrow=....)' for arranging plots will not work anymore. If withlegend is FALSE and group is NULL, layout is automatically deactivated and 'par(mfrow=....)' will work.}

	\item{legend.prop}{sets the proportion of the graphic area used for plotting the legend. Default value is set according to the place (bottom or right of the graphic area) where the legend is plotted. Values from 0 to 1.} 

 \item{rows,cols}{optional arguments to arrange plots if group is used.} 

 \item{...}{arguments to be passed to the plot function or other graphical parameters.}
}
\seealso{
 \code{\link{seqstatd}}
}
\examples{
## Loading the 'mvad' example data set
data(mvad)

## Defining a sequence object with the data in columns 15 to 86
mvad.seq <- seqdef(mvad,15:86)

## Plotting the mean times, 
## for each level of the gcse5eq factor
seqmtplot(mvad.seq, group=mvad$gcse5eq)
}
\keyword{hplot}
