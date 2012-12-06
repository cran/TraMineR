\name{TraMineR-package}
\alias{TraMineR-package}
\alias{TraMineR}
\docType{package}
\title{
Trajectory miner: a toolbox for exploring and rendering sequence data
}
\description{
This package is a toolbox for sequence manipulation, description, rendering and more generally the mining of sequence data in the field of social sciences. Though it is primarily intended for analyzing state or event sequences that describe life courses such as family formation histories or professional careers its features also apply to many other kinds of categorical sequence data. It accepts many different sequence representations as input and provides tools for translating sequences from one format to another. It offers several statistical functions for describing and rendering sequences, for computing distances between sequences with different metrics among which optimal matching, the longest common prefix and the longest common subsequence, and simple functions for extracting the most frequent subsequences and identifying the most discriminating ones among them. A user's guide can be found on TraMineR's web page.

TraMineR provides tools for both state sequences and event
sequences. The first step when using the package is to define a
state sequence object (with \code{\link{seqdef}}) if you want to
explore state sequences, and an event sequence object (with
\code{\link{seqecreate}}) if you are interested in event
sequencing.

State sequences are defined from a series of variables giving the
states at the successive positions, while event sequences are
defined from (vertical) time stamped event data. The package,
however, can handle many other different data organizations and
provides tools to help converting state sequences into event
sequences and vice versa.
}
\details{
\tabular{ll}{
Package: \tab TraMineR\cr
Version: \tab 1.8-5\cr
Date: \tab 2012-12-5\cr
Depends: \tab R (>= 2.8.1), RColorBrewer, boot\cr
Suggests: \tab cluster, xtable\cr
License: \tab GPL (>= 2)\cr
URL: \tab http://mephisto.unige.ch/traminer\cr
Encoding: \tab latin1\cr
}
Index:
\preformatted{
TraMineR-package        Trajectory miner: a toolbox for exploring and
                        rendering sequence data
TraMineR.checkupdates   Check for TraMineR updates
actcal                  Example data set: Activity calendar from the
                        Swiss Household Panel
actcal.tse              Example data set: Activity calendar from the
                        Swiss Household Panel (time stamped event
                        format)
alphabet                Get or set the alphabet of a sequence object
biofam                  Example data set: Family life states from the
                        Swiss Household Panel biographical survey
cpal                    Get or set the color palette of a sequence
                        object
dissassoc               Analysis of discrepancy from dissimilarity
                        measures
disscenter              Compute distances to the center of a group
dissmfac                Multi-factor ANOVA from a dissimilarity matrix
dissrep                 Extracting sets of representative objects using
                        a dissimilarity matrix
disstree                Dissimilarity Tree
disstree2dot            Graphical representation of a dissimilarity
                        tree
disstreeleaf            Terminal node membership
dissvar                 Dissimilarity based discrepancy
ex1                     Example data set with missing values and
                        weights
ex2                     Example data sets with weighted and unweighted
                        sequence data
famform                 Example data set: sequences of family formation
mvad                    Example data set: Transition from school to
                        work
plot.seqdiff            Plotting a seqdiff object.
plot.stslist            Plot method for state sequence objects
plot.stslist.freq       Plot method for sequence frequency tables
plot.stslist.meant      Plot method for objects produced by the
                        seqmeant function
plot.stslist.modst      Plot method for modal state sequences
plot.stslist.rep        Plot method for representative sequence sets
plot.stslist.statd      Plot method for objects produced by the
                        seqstatd function
plot.subseqelist        Plot frequencies of subsequences
plot.subseqelistchisq   Plot discriminant subsequences
read.tda.mdist          Read a distance matrix produced by TDA.
seqLLCP                 Compute the length of the longest common prefix
                        of two sequences
seqLLCS                 Compute the length of the longest common
                        subsequence of two sequences
seqST                   Sequences turbulence
seqalign                Computation details about a pairwise alignment
seqcomp                 Compare two state sequences
seqconc                 Concatenate vectors of states or events into a
                        character string
seqdecomp               Convert a character string into a vector of
                        states or events
seqdef                  Create a state sequence object
seqdiff                 Position-wise discrepancy analysis between
                        groups of sequences
seqdim                  Dimension of a set of sequences
seqdist                 Distances (dissimilarities) between sequences
seqdistmc               Multichannel distances between sequences
seqdss                  Extract distinct states sequence from a
                        sequence object
seqdur                  Extract state durations from a sequence object.
seqeapplysub            Checking for the presence of given event
                        subsequences
seqecmpgroup            Identifying discriminating subsequences
seqeconstraint          Setting time constraints and the counting
                        method
seqecontain             Check if sequence contains events
seqecreate              Create event sequence objects.
seqefsub                Searching for frequent subsequences
seqeid                  Retrieve unique ids from an event sequence
                        object.
seqelength              Lengths of event sequences
seqetm                  Create a transition-definition matrix
seqeweight              Setting or retrieving weights of an event
                        sequence object.
seqfind                 Indexes of state sequence(s) x in state
                        sequence object y
seqformat               Conversion between sequence formats
seqfpos                 Search for the first occurrence of a given
                        element in a sequence
seqgen                  Random sequences generation
seqici                  Complexity index of individual sequences
seqient                 Within sequence entropies
seqistatd               State frequencies in each individual sequence
seqlegend               Plot a legend for the states in a sequence
                        object
seqlength               Sequence length
seqlogp                 Logarithm of the probabilities of state
                        sequences
seqmeant                Mean durations in each state
seqmodst                Sequence of modal states
seqmpos                 Number of matching positions between two
                        sequences.
seqnum                  Transform into a sequence object with numerical
                        alphabet.
seqpcplot               Parallel coordinate plot for sequence data
seqplot                 Plot state sequence objects
seqpm                   Find patterns in sequences
seqrecode               Recoding state sequence objects and factors
seqrep                  Extracting sets of representative sequences
seqsep                  Adds separators to sequences stored as
                        character string
seqstatd                Sequence of transversal state distributions and
                        their entropies
seqstatf                State frequencies in the whole sequence data
                        set
seqstatl                List of distinct states or events (alphabet) in
                        a sequence data set.
seqsubm                 Create a substitution-cost matrix
seqsubsn                Number of distinct subsequences in a sequence.
seqtab                  Frequency table of the sequences
seqtransn               Number of transitions in a sequence
seqtrate                Compute transition rates between states
seqtree                 Tree structured analysis of a state sequence
                        object.
seqtreedisplay          Graphical rendering of a sequence regression
                        tree
stlab                   Get or set the state labels of a sequence
                        object
}
Further information is available in the following vignettes:
\tabular{ll}{
\code{TraMineR-state-sequence} \tab Analyzing and Visualizing State Sequences (source)\cr
}
}
\author{
Alexis Gabadinho <alexis.gabadinho@unige.ch>, Matthias Studer <matthias.studer@unige.ch>, Nicolas S. Müller <nsmuller@gmail.com>, Reto Bürgin <reto.buergin@unige.ch> and Gilbert Ritschard <gilbert.ritschard@unige.ch>.
\cr
Maintainer: Gilbert Ritschard <gilbert.ritschard@unige.ch>
}

\references{
{Gabadinho, A., G. Ritschard, N. S. Müller and M. Studer (2011).
Analyzing and Visualizing State Sequences in R with TraMineR.
\emph{Journal of Statistical Software} \bold{40}(4), 1-37.
}

{Gabadinho, A., G. Ritschard, M. Studer and N. S. Müller (2009).
Mining Sequence Data in \code{R} with the \code{TraMineR} package:
A users guide. Department of Econometrics and Laboratory of
Demography, University of Geneva }
}

\examples{
## load the mvad data
library(TraMineR)
data(mvad)

## create a state sequence object from columns 17 to 86
mvad.seq <- seqdef(mvad[,17:86])

## distribution plot by sex (male)
seqdplot(mvad.seq, group=mvad$male, border=NA)

## compute the LCS pairwise distance matrix
## among the first 10 sequences
mvad.lcs <- seqdist(mvad.seq[1:10,], method="LCS")
}

\keyword{package}

