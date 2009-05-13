\name{plot.stslist}
\alias{plot.stslist}
\title{Plot method for sequence objects of class 'stslist'}
\description{
This is the plot function for sequence objects of class 'stslist' created by the seqdef function.
}
\usage{
\method{plot}{stslist}(x, group=NULL, type="i", title=NULL, 
	tlim=NULL, pbarw=FALSE, sortv=NULL, 
	cpal=NULL, missing.color=NULL, 
	ylab, yaxis=NULL, axes="all", xtlab=NULL, cex.plot=1,
	withlegend="auto", ltext=NULL, cex.legend=1, 
	use.layout=(!is.null(group) | withlegend!=FALSE), 
	legend.prop=NA, rows=NA, cols=NA, ...)
}
\arguments{
	\item{x}{a sequence object created with the \code{\link{seqdef}} function.}

	\item{group}{Plots one plot for each level of the factor given as argument.}
	
	\item{type}{the type of the plot. Available types are \code{"d"} for state distribution plots, \code{"f"} for sequence frequency plots, \code{"i"} for sequence index plot, \code{"mt"} for mean times plots.}

	\item{title}{title for the graphic. Default to NULL.}

	\item{tlim}{if \code{type="f"} and \code{tlim>0}, shows only the 'tlim' most frequent sequences (default to 10). If \code{type="i"}, \code{tlim} is the indexes of the sequences to plot (default value is 1:10), for instance 20:50 to plot sequences 20 to 50, c(2,8,12,25) to plot sequences 2,8,12 and 25 in \code{seqdata}. If set to 0, all sequences in \code{seqdata} are plotted.}

	\item{pbarw}{if \code{type="f"} and \code{pbarw=TRUE}, the width of the bars are proportional to the sequence frequency in the dataset. Default to FALSE.}

	\item{sortv}{if \code{type="i"}, name of an optional variable used to sort the sequences before plotting.}

	\item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of distinct states. By default, the 'cpal' attribute of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}

	\item{missing.color}{alternative color for representing missing values inside the sequences. By default, this color is taken from the "missing.color" attribute of the sequence object being plotted.} 

	\item{ylab}{An optional label for the y axis. If set to NA, no label is drawn.}

	\item{yaxis}{Controls wheter a y axis is plotted. If left to \code{'NULL'}, a the default value is set depending on the plot type. If \code{type="d"} the default is \code{'TRUE'}. If \code{type="i"} the default is \code{'FALSE'}, if set to \code{'TRUE'}, the sequence indexes are displayed. }

	\item{axes}{if set to "all" (default value) x axes are drawn for each plot in the graphic. If set to "bottom" and \code{group} is used, axes are drawn only under the plots located at the bottom of the graphic area. If FALSE, no x axis is drawn.}

  \item{xtlab}{optional labels for the x axis ticks labels. If unspecified, the column names of the 'seqdata' sequence object are used (see \code{\link{seqdef}}).}

	\item{cex.plot}{expansion factor for setting the size of the font for the axis labels and names. The default value is 1. Values lesser than 1 will reduce the size of the font, values greater than 1 will increase the size.}

  \item{withlegend}{set wether a legend for the states colors is plotted. The default value 'auto' sets the position of the legend automatically. Other options are 'right'. Obsolete option 'TRUE' is identical to 'auto'.}

  \item{ltext}{optional description of the states to appear in the legend. Must be a vector of character strings with number of elements equal to the number of distinct states. If unspecified, the 'label' attributes of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}

	\item{cex.legend}{expansion factor for setting the size of the font for the labels in the legend. The default value is 1. Values lesser than 1 will reduce the size of the font, values greater than 1 will increase the size.}

 \item{use.layout}{if TRUE, layout is used to arrange plots when using the group option or plotting a legend. If layout is used, the standard 'par(mfrow=....)' for arranging plots will not work anymore. If withlegend is FALSE and group is NULL, layout is automatically deactivated and 'par(mfrow=....)' will work.}

 \item{legend.prop}{sets the proportion of the graphic area used for plotting the legend when use.layout=TRUE and withlegend=TRUE. Default value is set according to the place (bottom or right of the graphic area) where the legend is plotted. Values from 0 to 1.} 

 \item{rows,cols}{optional arguments to arrange plots when use.layout=TRUE.} 

 \item{...}{arguments to be passed to the plot function or other graphical parameters.}
}

\details{
The state distribution plot (\code{type="d"}) presents the sequence of the states frequencies for each time point, as computed by the \code{\link{seqstatd}} function.

In sequence frequency plots (\code{type="f"}), the sequences are rendered as stacks of successive states, ordered by their relative frequency in the dataset.

In sequence index plots (\code{type="i"}), individual sequences are rendered with stacked bars depicting the statuses
over time. The interest of such plots, known as index-plots, has for instance been stressed by \cite{Scherer (2001)}, \cite{Brzinsky-Fay et al. (2006)} and \cite{Gauthier (2007)}. Notice that such index plots for thousands of sequences result in very heavy graphic files if they are stored in PDF or POSTSCRIPT format. To reduce the size, the figures can be saved in png format by using \code{\link{png}} instead of \code{\link{postscript}} or \code{\link{pdf}}.

The mean time plot (\code{type="mt"}) displays the mean time spent in each state of the alphabet.
}

\seealso{
 \code{\link{seqdplot}, \link{seqfplot}, \link{seqiplot}, \link{seqmtplot}}.
}
\examples{
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)

## Plot of the 10 most frequent sequences 
## with bar width proportional to the frequency
plot(biofam.seq, type="f", pbarw=TRUE)

## State distribution plot
## in the actcal data set
## grouped by sex
data(actcal)
actcal.lab <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal,13:24,labels=actcal.lab)

plot(actcal.seq, type="d", group=actcal$sex)
}
\keyword{hplot}
