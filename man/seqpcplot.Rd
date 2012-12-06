%% -*- RD:ANSI -*-
\name{seqpcplot}
\alias{seqpcplot}
\alias{plot.seqe}
\alias{plot.seqelist}
\title{Parallel coordinate plot for sequence data}
\description{
  A decorated parallel coordinate plot to render the order of the
  successive elements in sequences. The sequences are displayed as
  jittered frequency-weighted parallel lines.
  The plot is also embedded as the \code{type="pc"} option of the
  \code{\link{seqplot}} function and serves as plot
  method for \code{seqe} and \code{seqelist} objects.
}


\usage{

seqpcplot(seqdata, group = NULL, weights = NULL,
          cex = 1, lwd = 1/4, cpal = NULL, grid.scale = 1/5,
          ltype = "unique", embedding = "most-frequent",
          lorder = NULL , lcourse = "upwards",
          filter = NULL, hide.col = "grey80",
          alphabet = NULL, order.align = "first",
          title = NULL, xlab = NULL, ylab = NULL,
          xaxis = TRUE, yaxis = TRUE, axes = "all",
          xtlab = NULL, cex.plot = 1,
          rows = NA, cols = NA, plot = TRUE,
          output = !plot, ...)

}
\arguments{
  \item{seqdata}{The sequence data. Either an event sequence
    object of class \code{seqelist} (see
    \code{\link{seqecreate}}) or a state sequence object of class
    \code{stslist} (see \code{\link{seqdef}}).}

  \item{group}{a vector (numeric or factor) of group memberships
    of length equal the number of
    sequences. When specified, one plot is generated for each
    different membership value.}

  \item{weights}{a numeric vector of weights of length equal the number
  of sequences. Overrides weights in the \code{seqdata} object.}

  \item{cex}{expansion factor for the squared symbols.}

  \item{lwd}{expansion factor for line widths. The expansion is
    relative to the size of the squared symbols.}

  \item{cpal}{color palette vector for line coloring.}

  \item{grid.scale}{Expansion factor for the translation zones.}

  \item{ltype}{the type of sequence that is drawn. Either \code{"unique"}
    to render unique patterns or \code{"non-embeddable"} to render
    non-embeddable sequences.}

  \item{embedding}{The method for embedding sequences embeddable in
    multiple non-embeddable sequences. Either \code{"most-frequent"}
    (default) or \code{"uniformly"}. Relevant only with \code{ltype =
      "non-embeddable"}.}

  \item{lorder}{line ordering. Either \code{"background"} or
  \code{"foreground"}.}

  \item{lcourse}{Method to connect simultaneous elements with the
    preceding and following ones. Either \code{"upwards"} (default) or
    \code{"downwards"}.}

  \item{filter}{list of line coloring options. See details.}

  \item{hide.col}{Color for sequences filtered-out by the
    \code{filter} specification.}

  \item{alphabet}{a vector of response levels in the order they should
    appear on the y-axis. This argument is solely relevant for
    \code{seqelist} objects.}

  \item{order.align}{Method for converting timestamps to individual
  order positions. Either \code{"first"} (default) or
  \code{"last"}. \code{"first"} numbers the positions from the beginning
  while \code{"last"} numbers them from the end.}

  \item{title}{title for the graphic.}
  \item{xlab}{label for the x axis}
  \item{ylab}{label for the y axis}
  \item{xaxis}{logical: Should x-axis be plotted?}
  \item{yaxis}{logical: Should y-axis be plotted?}

  \item{axes}{if set to \code{"all"} (default value) x-axes are drawn
    for each plot in the graphic. If set to \code{"bottom"} and
    \code{group} is used, axes are drawn only under the plots at
    the bottom of the graphic area. If \code{FALSE}, no x-axis is drawn.}

  \item{xtlab}{labels for the x-axis ticks.}

  \item{cex.plot}{expansion factor for the size of the font for
    the axis labels and names. The default value is 1. Values lesser
    than 1 will reduce the size of the font, values greater than 1 will
    increase the size.}

  \item{rows,cols}{integers to arrange the plot panel design.}

  \item{plot}{logical. If \code{FALSE} nothing is plotted and an object
    of class \code{seqpcplot} is returned by default.}

  \item{output}{logical. Should object
    of class \code{seqpcplot} be returned?}

  \item{...}{arguments to be passed to other methods, such as graphical
    parameters (see \code{\link{par}}).}
}
\details{
  For plots by groups specified with the \code{group} argument, plotted
  line widths and point sizes reflect relative frequencies within
  group.

  The \code{filter} argument serves to specify filters to gray less
  interesting patterns. The filtered-out patterns are displayed in the
  \code{hide.col} color. The \code{filter} argument expects a list with
  at least elements \code{type} and \code{value}. The following types
  are implemented:

  Type \code{sequence}: colors a specific pattern, for example assign\cr
  \code{filter = list(type = "sequence", value = "(Leaving
  Home,Union)-(Child)")}.

  Type \code{subsequence}: colors patterns which include a specific
  subsequence, for example\cr \code{filter = list(type =
  "subsequence", value = "(Child)-(Marriage)")} .

  Type \code{value}: gradually colors the patterns according to the
  numeric vector (of length equal to the number of sequences) provided as
  \code{"value"} element in the list. You can give something like
  \code{filter = list(type = "value", value = c(0.2, 1, ...))} or
  provide the distances to the medoid as \code{value} vector for
  example.

  Type \code{function}: colors the patterns depending on the values
  returned by a [0,1] valued function on the frequency x of the
  pattern. Two native functions can be used: \code{"minfreq"} and
  \code{"cumfreq"}. Use \code{filter = list(type = "function", value =
  "minfreq", level = 0.05)} to color patterns with a support of at least
  5\% (within group). Or, use\cr\code{filter = list(type = "function",
  value = "cumfreq", level = 0.5)} to highlight the 50\% most frequent
  patterns (within group). Other user-specified functions can be
  provided by giving something like\cr\code{filter = list(type =
  "function", value = function(x, arg1, arg2) {return(x / max(x) * arg1
  / arg2)}, arg1 = 1, arg2 = 1)}. This latter function adjusts
  gradually the color intensity of patterns according to the frequency
  of the pattern.

}


\seealso{ \code{\link{seqplot}}, \code{\link{seqdef}}, \code{\link{seqecreate}}}

\references{
  Bürgin, R., G. Ritschard and E. Rousseaux (2012), Exploration graphique de donnéees séquentielles. In Atelier Fouille Visuelle de Données : méthodologie et évaluation, EGC 2012, Bordeaux. Association EGC: http://www.egc.asso.fr/sdoc-237-egc2012_atelier_fvd.pdf.
}

\examples{

## ================
## plot biofam data
## ================

data(biofam)
lab <- c("Parent","Left","Married","Left+Marr","Child","Left+Child",
         "Left+Marr+Child","Divorced")

## plot state sequences in STS representation
## ==========================================

## creating the weighted state sequence object.
biofam.seq <- seqdef(data = biofam[,10:25], labels = lab,
                     weights = biofam$wp00tbgs)

## select the first 20 weighted sequences (sum of weights = 18)
biofam.seq <- biofam.seq[1:20, ]

par(mar=c(4,8,2,2))
seqpcplot(seqdata = biofam.seq, order.align = "time")

## .. or
seqplot(seqdata = biofam.seq, type = "pc", order.align = "time")

## Distinct successive states (DSS)
## ==========================================

seqplot(seqdata = biofam.seq, type = "pc", order.align = "first")

## .. or (equivalently)

biofam.DSS <- seqdss(seqdata = biofam.seq) # prepare format
seqpcplot(seqdata = biofam.DSS)

## plot TSE data converted from state sequences
## ============================================

## conversion
biofam.TSE <- seqformat(data = biofam.seq, from = "STS", to = "TSE",
                        tevent = seqetm(biofam.seq, method = "state"))
biofam.TSE$event <- factor(biofam.TSE$event, levels = lab) # define alphabet
biofam.TSE$time <- biofam.TSE$time + 15 # correct age

seqpcplot(seqdata = biofam.TSE, order.align = "time")


## plot event sequences
## ====================

biofam.seqe <- seqecreate(biofam.seq, tevent = "state") # prepare data

## plot the time in the x-axis
seqpcplot(seqdata = biofam.seqe, order.align = "time", alphabet = lab)

## ordering of events
seqpcplot(seqdata = biofam.seqe, order.align = "first", alphabet = lab)

## ... or
plot(biofam.seqe, order.align = "first", alphabet = lab)

## additional arguments
## ====================

## non-embeddable sequences
seqpcplot(seqdata = biofam.seqe, ltype = "non-embeddable",
          order.align = "first", alphabet = lab)

## align on last event
par(mar=c(4,8,2,2))
seqpcplot(seqdata = biofam.seqe, order.align = "last", alphabet = lab)

## use group variables
seqpcplot(seqdata = biofam.seqe, group = biofam$sex[1:20],
          order.align = "first", alphabet = lab)

## color patterns (Parent)-(Married) and (Parent)-(Left+Marr+Child)
par(mfrow = c(1, 1))
seqpcplot(seqdata = biofam.seqe,
          filter = list(type = "sequence",
                          value=c("(Parent)-(Married)",
                                  "(Parent)-(Left+Marr+Child)")),
          alphabet = lab, order.align = "first")

## color subsequence pattern (Parent)-(Left)
seqpcplot(seqdata = biofam.seqe,
          filter = list(type = "subsequence",
                          value = "(Parent)-(Left)"),
          alphabet = lab, order.align = "first")

## color sequences over 5% (within group) (function method)
seqpcplot(seqdata = biofam.seqe,
          filter = list(type = "function",
                          value = "minfreq",
                          level = 0.1),
          alphabet = lab, order.align = "first")

## highlight the 50% most frequent sequences
seqpcplot(seqdata = biofam.seqe,
          filter = list(type = "function",
                          value = "cumfreq",
                          level = 0.5),
          alphabet = lab, order.align = "first")

}

\author{Reto Bürgin (with Gilbert Ritschard for the help page)}

\keyword{Plot}
\keyword{Event sequences} 
