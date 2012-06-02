\name{plot.subseqelistchisq}
\alias{plot.subseqelistchisq}
\title{Plot discriminant subsequences}
\description{
Plot the result of \code{\link{seqecmpgroup}}
}
\usage{
\method{plot}{subseqelistchisq}(x, ylim = "uniform", rows = NA, cols = NA,
    residlevels = c(0.05,0.01),
    cpal = brewer.pal(1 + 2 * length(residlevels), "RdBu"),
    legendcol = NULL, legend.cex = 1, ptype="freq", 
    legend.title = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The subsequences to plot (a \code{subseqelist} object).}
  \item{ylim}{if \code{"uniform"} all axes have same limits.}
  \item{rows}{Number of graphic rows}
  \item{cols}{Number of graphic columns}
  \item{residlevels}{Significance levels used to colorize the Pearson residual}
  \item{cpal}{Color palette used to color the results}
  \item{legendcol}{When \code{TRUE} the legend is printed vertically, when \code{FALSE} it is printed horizontally. If \code{NULL} (default) the best position will be chosen.}
  \item{legend.cex}{Scale parameters for text legend.}
  \item{ptype}{If set to \code{"resid"}, Pearson residuals are plotted instead of frequencies}
  \item{legend.title}{Legend title.}  
  \item{\dots}{Additional parameters passed to \code{\link{barplot}}}
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
