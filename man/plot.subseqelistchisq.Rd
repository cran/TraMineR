\name{plot.subseqelistchisq}
\alias{plot.subseqelistchisq}
\title{Ploting discriminant subsequences}
\description{
Plot the result of \code{\link{seqecmpgroup}}
}
\usage{
\method{plot}{subseqelistchisq}(x, ylim = "uniform", rows = NA, cols = NA, 
	residlevels = c(2, 4), 
	cpal = brewer.pal(1 + 2 * length(residlevels), "RdBu"), 
	legendcol = NULL, legend.cex = 1,ptype="freq", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The subsequences to plot.}
  \item{ylim}{if \code{"uniform"} all axis have same limits.}
  \item{rows}{Number of graphic rows}
  \item{cols}{Number of graphic columns}
  \item{residlevels}{The pearson residuals levels to colorize}
  \item{cpal}{Colors palette used to color the results}
  \item{legendcol}{Should the legend be printed in column ? If NULL the best position will be choosen}
  \item{legend.cex}{Scale parameters to text legend}
  \item{ptype}{If equal to "resid", pearson residuals are plotted instead of frequencies}
  \item{\dots}{additionnal parameters passed to \code{\link{barplot}}}
}
%\details{}
\value{
  nothing
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{seqecmpgroup}}}
