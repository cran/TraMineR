\encoding{latin1}
\name{seqiplot}
\alias{seqiplot}
\title{Visualization of individual sequences}
\description{
  Plots a graphic showing individual sequences as stacked bars. 
}

\details{The seqiplot function renders individual sequences with stacked bars depicting the statuses
over time. The interest of such plots, known as index-plots, has for instance been stressed by \cite{Scherer (2001)}, \cite{Brzinsky-Fay et al. (2006)} and \cite{Gauthier (2007)}. 
Notice that such index plots for thousands of sequences result in very heavy graphic files if they are stored in PDF or POSTSCRIPT format. To reduce the size, the figures can be saved in png format by using \code{\link{png}} instead of \code{\link{postscript}} or \code{\link{pdf}}.
}

\usage{
seqiplot(seqdata, group=NULL, tlim=1:10, sortv=NULL, title=NULL,
	cpal=NULL, missing.color=NULL,
	withborder=TRUE, 
	ylab, axes="all", xtlab=NULL, cex.plot=1,
	withlegend="auto", ltext=NULL, cex.legend=1, 
	use.layout=(!is.null(group) | withlegend!=FALSE), 
	legend.prop=NA, rows=NA, cols=NA, ...)
}
\arguments{
	\item{seqdata}{a sequence object as defined by the the \code{\link{seqdef}} function.}

	\item{group}{Plots one plot for each level of the factor given as argument.}

	\item{tlim}{interval representing the range to be plotted. Default to 1:10, ie the first 10 rows. Set to 0 if you want to plot all sequences in 'seqdata'.}

	\item{sortv}{name of an optional variable used to sort the data before plotting.}

	\item{title}{title for the graphic. Default to NULL.}

	\item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of distinct states. By default, the 'cpal' attribute of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}

	\item{missing.color}{alternative color for representing missing values inside the sequences. By default, this color is taken from the 'missing.color' attribute of the sequence object being plotted.} 

	\item{withborder}{Now obsolete. Standard option 'border=NA' can be used (see ... optional arguments below) to omit borders around the bars.}

	\item{ylab}{An optional label for the y axis. If set to NA, no label is drawn.} 

	\item{axes}{if set to "all" (default value) x axes are drawn for each plot in the graphic. If set to "bottom" and \code{group} is used, axes are drawn only under the plots located at the bottom of the graphic area. If FALSE, no x axis is drawn.}

	\item{xtlab}{optional labels for the x axis ticks labels. If unspecified, the column names of the 'seqdata' sequence object are used (see \code{\link{seqdef}}).}

	\item{cex.plot}{expansion factor for setting the size of the font for the axis labels and names. The default value is 1. Values lesser than 1 will reduce the size of the font, values greater than 1 will increase the size.}
	
  \item{withlegend}{set wether a legend for the states colors is plotted. The default value 'auto' sets the position of the legend automatically. Other options are 'right'. Obsolete option 'TRUE' is identical to 'auto'.}

	\item{ltext}{optional description of the states to appear in the legend. Must be a vector of character strings with number of elements equal to the number of distinct states. If unspecified, the 'label' attributes of the 'seqdata' sequence object is used (see \code{\link{seqdef}}).}

	\item{cex.legend}{expansion factor for setting the size of the font for the labels in the legend. The default value is 1. Values lesser than 1 will reduce the size of the font, values greater than 1 will increase the size.}

	\item{use.layout}{if TRUE, layout is used to arrange plots when using the group option or plotting a legend. If layout is used, the standard 'par(mfrow=....)' for arranging plots will not work anymore. If withlegend is FALSE and group is NULL, layout is automatically deactivated and 'par(mfrow=....)' will work.}

	\item{legend.prop}{sets the proportion of the graphic area used for plotting the legend. Default value is set according to the place (bottom or right of the graphic area) where the legend is plotted. Values from 0 to 1.} 

	\item{rows,cols}{optional arguments to arrange plots when layout is used.}

	\item{...}{arguments to be passed to the plot function or other graphical parameters. Use 'border=NA' to omit borders around the bars. Use 'space=' to set the space between each sequence in the plot. Setting to 0 yiels no space and is usefull if plotting a large amount of sequences. Default to NULL, ie space value is automaticaly computed.}
}
%\seealso{}

\references{
Brzinsky-Fay, C., U. Kohler, and M. Luniak (2006). Sequence analysis with Stata. \emph{The Stata Journal} 6 (4), 435-460.

Gauthier, J.-A. (2007). \emph{Empirical categorizations of social trajectories: A sequential view on the life course}. thèse, Université de Lausanne, Faculté des sciences sociales et politique (SSP), Lausanne.

Scherer, S. (2001). Early career patterns: A comparison of Great Britain and West Germany. \emph{European Sociological Review} 17 (2), 119-144.
}

\examples{
## Plotting the first 10 sequences in the 'actcal' example data set
data(actcal)
seqiplot(seqdef(actcal,13:24))

## Plotting the all the sequences in the 'biofam' example data set
## using 'border=NA' and 'space=0' options to have a nicier plot
data(biofam)
seqiplot(seqdef(biofam,10:25),tlim=0, border=NA, space=0)

## Plotting the all the sequences in the 'mvad' example data set
## by levels of the 'gcse5eq' factor
data(mvad)
seqiplot(seqdef(mvad,15:86), group=mvad$gcse5eq, 
tlim=0, border=NA, space=0)

## The same plot without legend, axes and 
## ylab
data(mvad)
seqiplot(seqdef(mvad,15:86), group=mvad$gcse5eq, 
tlim=0, border=NA, space=0,
axes=FALSE, ylab=NA)
}
\keyword{hplot}
