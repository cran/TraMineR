\encoding{latin1}
\name{seqiplot}
\alias{seqiplot}
\title{Visualization of individual sequences}
\description{Plots a graphic showing individual sequences as stacked bars}

\details{The seqiplot function renders individual sequences with stacked bars depicting the statuses
over time. The interest of such plots, known as index-plots, has for instance been stressed by \cite{Scherer (2001)}, \cite{Brzinsky-Fay et al. (2006)} and \cite{Gauthier (2007)}. 
Notice that such index plots for thousands of sequences result in very heavy graphic files if they are stored in PDF or POSTSCRIPT format. To reduce the size, the figures can be saved in png format by using \code{\link{png}} instead of \code{\link{postscript}} or \code{\link{pdf}}.
}

\usage{
seqiplot(seqdata, group=NULL, title=NULL, tlim=NULL, sortv=NULL,  ...)
}
\arguments{
	\item{seqdata}{a sequence object as defined by the the \code{\link{seqdef}} function.}

	\item{group}{Plots one plot for each level of the factor given as argument.}

	\item{title}{title for the graphic. Default to NULL.}

	\item{tlim}{interval representing the range to be plotted. If \code{NULL} (default) tlim is set to \code{1:10}, ie the first 10 sequences are displayed. Set to 0 if you want to plot all sequences in 'seqdata'.}

	\item{sortv}{name of an optional variable used to sort the data before plotting.}

	\item{...}{arguments to be passed to the \code{\link{plot.stslist}} function or other graphical parameters. Use 'border=NA' to omit borders around the bars. Use 'space=' to set the space between each sequence in the plot. Setting to 0 yields no space and is usefull if plotting a large amount of sequences. Default to NULL, ie space value is automaticaly computed. Use \code{axes} for controling the display of axes, \code{withlegend} for controling the display of a legend for the plot, etc...}
}
\seealso{\code{\link{plot.stslist}}
}

\references{
Brzinsky-Fay, C., U. Kohler, and M. Luniak (2006). Sequence analysis with Stata. \emph{The Stata Journal} 6 (4), 435-460.

Gauthier, J.-A. (2007). \emph{Empirical categorizations of social trajectories: A sequential view on the life course}. thèse, Université de Lausanne, Faculté des sciences sociales et politique (SSP), Lausanne.

Scherer, S. (2001). Early career patterns: A comparison of Great Britain and West Germany. \emph{European Sociological Review} 17 (2), 119-144.
}

\examples{
## Plotting the first 10 sequences in the 'actcal' example data set
data(actcal)
seqiplot(seqdef(actcal,13:24))

## Plotting all sequences of the 'biofam' example data set
## sorted by birth year (birthy variable in the biofam data set)
## using 'border=NA' and 'space=0' options to have a nicer plot
data(biofam)
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
"Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)

seqiplot(biofam.seq, tlim=0, border=NA, space=0, sortv=biofam$birthy)

## Plotting the sequences in the 'mvad' example data set
## grouped by levels of the 'gcse5eq' factor
data(mvad)
seqiplot(seqdef(mvad,15:86), group=mvad$gcse5eq, title="gcse5eq",
tlim=0, border=NA, space=0)
}
\keyword{hplot}
