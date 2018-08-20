\name{plot.stslist}
\alias{plot.stslist}
\title{Plot method for state sequence objects}
\description{
This is the plot method for state sequence objects of class \code{stslist} created by the \code{seqdef} function. It produces a sequence index plot.
}
\usage{
\method{plot}{stslist}(x, idxs = NULL, weighted = TRUE, sortv = NULL,
  cpal = NULL, missing.color = NULL, ylab = NULL,
  yaxis = TRUE, xaxis = TRUE, ytlab = NULL, ylas = 0,
  xtlab = NULL, xtstep = NULL, tick.last = NULL, cex.axis = 1,
  tlim, cex.plot, ...)
}
\arguments{
	\item{x}{A state sequence object created with the \code{\link{seqdef}} function.}

	\item{idxs}{Indexes of the sequences to be plotted (default value is \code{1:10}), for instance \code{20:50} to plot sequences 20 to 50, \code{c(2,8,12,25)} to plot sequences 2,8,12 and 25 in \code{seqdata}. If set to \code{0}, all sequences in \code{seqdata} are plotted.}

     \item{weighted}{Logical: Should the bar representing each sequence be proportional to its weight? Ignored when no weights are assigned to sequences (see \code{\link{seqdef}}.)}

	\item{sortv}{A sorting variable or a sort method (one of \code{"from.start"} or \code{"from.end"}). See details.}

	\item{cpal}{alternative color palette to use for the states. If user specified, a vector of colors with number of elements equal to the number of states in the alphabet. By default, the \code{cpal} attribute of the \code{seqdata} sequence object is used (see \code{\link{seqdef}}).}

	\item{missing.color}{alternative color for representing missing values inside the sequences. By default, this color is taken from the "missing.color" attribute of the \code{x} sequence  object.}

	\item{ylab}{An optional label for the y axis. If set to \code{NA}, no label is drawn.}

	\item{yaxis}{Controls whether the y axis is plotted or not. When set to \code{TRUE}, sequence indexes are displayed. }

	\item{xaxis}{if \code{TRUE} (default), the x (time) axis is plotted.}

	\item{ytlab}{the labels of the plotted sequences to display on the y axis. Default is the indexes of the sequences as defined by the \code{idxs} argument. Can be set to \code{"id"} for displaying the row names (id) of the sequences instead of their indexes; row names can be assigned to the sequence object with the  \code{id} argument of the \code{\link{seqdef}} function or afterwards with \code{\link{rownames}}. Otherwise \code{ytlab} can be set to a vector of length equal to the number of  sequences to be plotted.}

  \item{ylas}{sets the orientation of the sequence labels appearing on the y axis. Accepted values are the same as for the \code{las} standard option
  \cr
0: always parallel to the axis (default),
    \cr
1: always horizontal,
    \cr
2: always perpendicular to the axis,
    \cr
3: always vertical.
}

  \item{xtlab}{optional labels for the x axis ticks labels. If unspecified, the column names of the \code{seqdata} sequence object are used (see \code{\link{seqdef}}).}

  \item{xtstep}{optional interval at which the tick-marks and labels of the x-axis are displayed. For example, with \code{xtstep=3} a tick-mark is drawn at position 1, 4, 7, etc... The display of the corresponding labels depends on the available space and is dealt with automatically. If unspecified, the \code{xtstep} attribute of the \code{x} object is used.}

  \item{tick.last}{Logical. Should a tick mark be enforced at the last position on the x-axis? If unspecified, the \code{tick.last} attribute of the \code{x} object is used.}

  \item{cex.axis}{Axis annotation magnification. See \code{\link{par}}.}

  \item{tlim}{Deprecated. Use \code{idxs} instead.}
  \item{cex.plot}{Deprecated. Use \code{cex.axis} instead.}

 \item{...}{arguments to be passed to the plot function or other graphical parameters.}
}

\details{
This is the default plot method for state sequence objects (produced by the \code{\link{seqdef}} function), i.e., for objects of class \emph{stslist}. It produces a sequence index plot, where individual sequences are rendered with stacked bars depicting the states over time.
\cr
\cr

This method is called by the generic \code{\link{seqplot}} function (if \code{type="i"}). The latter produces more sophisticated plots, allowing grouping and automatic display of the state color legend. The \code{seqiplot} function is a shortcut for calling \code{seqplot} with \code{type="i"}.

When a \code{sortv} variable is provided to \code{seqiplot} or \code{seqIplot}, its values define the order in which the sequences are plotted. With \code{sortv = "from.start"}, sequence are sorted by the elements of the alphabet at the successive positions starting from the beginning of the sequences. The {"from.end"} method proceeds similarly, but backward from the last position.

The interest of sequence index plots has for instance been stressed by \cite{Scherer (2001)} and \cite{Brzinsky-Fay et al. (2006)}. Notice that such index plots for thousands of sequences result in very heavy graphic files if they are stored in PDF or POSTSCRIPT format. To reduce the size, we suggest saving the figures in bitmap format by using for instance \code{\link{png}} instead of \code{\link{postscript}} or \code{\link{pdf}}.
}

\seealso{\code{\link{seqplot}}
}

\examples{
## Defining a sequence object with the data in columns 10 to 25
## (family status from age 15 to 30) in the biofam data set
data(biofam)
biofam <- biofam[500:600,] ## using a subsample only
biofam.lab <- c("Parent", "Left", "Married", "Left+Marr",
    "Child", "Left+Child", "Left+Marr+Child", "Divorced")
biofam.seq <- seqdef(biofam, 10:25, labels=biofam.lab)

## Plot of the 10 most frequent sequences
## with bar width proportional to the frequency
plot(biofam.seq)

## Plotting the whole data set
## with no borders
plot(biofam.seq, idxs=0, space=0, border=NA)

## =======
## Weights
## =======
data(ex1)
ex1.seq <- seqdef(ex1, 1:13, weights=ex1$weights)
plot(ex1.seq)
plot(ex1.seq, weighted=FALSE)
}
\keyword{Plot}
