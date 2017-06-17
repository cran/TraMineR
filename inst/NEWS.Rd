\name{NEWS}
\title{NEWS file for the TraMineR package}


\section{Changes in version 2.0-6 [2017-06-17]}{
  \subsection{Changes in existing functions:}{
    \itemize{
      \item
      \code{seqformat()}:
      \itemize{
        \item
        Redesigned version of the function with a new extensible, robust and documented
        implementation.
        \item
        Rewritten documentation with argument types, default values, scopes and
        detailed explanations.
        \item
        Added and generalized conversion to \code{"SPELL"} from TraMineRextras;
        added a \code{with.missing} argument specifically for this conversion.
        \item
        Added the possibility to pass directly the unique individual IDs (row names)
        of the input sequences with \code{id} when converting to \code{"TSE"}.
        \item
        Clarified the meaning of the \code{id} argument and changed its position in the
        list of arguments to reflect this.
        \item
        Changed the default value of the \code{id}, \code{begin}, \code{end}, and
        \code{status} arguments (now \code{"SPELL"} oriented).
        \item
        Renamed \code{compressed} argument as \code{compress} to avoid
        confusion (here it applies to the output, not to input data).
        \item
        Renamed \code{nr} argument as \code{missing} to avoid confusion (here
        it specifies the code to consider as missing values in input data,
        which is not the TraMineR internal code for missing values in state
        sequence objects).
        \item
        Adapted function calls to match the renamed argument names
        (\code{compressed}, \code{nr}) and the new input type checking (ie.
        single strings are deprecated).
        \item
        See the updated documentation of \code{seqformat()} for details.
      }
      \item
      \code{SPS_to_STS()}: renamed \code{nr} argument as \code{missing} to match
      \code{seqformat()} argument names renaming.
    }
  }
  \subsection{Bug fixes:}{
    \itemize{
      \item
      \code{seqformat()}:
      \itemize{
        \item
        From \code{"SPELL"} to \code{"TSE"}: output IDs were incorrect.
        \item
        From \code{"SPELL"} to \code{"STS"}: states not appearing in the
        data were dropped from the output column levels.
      }
    }
  }
  \subsection{Misc:}{
    \itemize{
      \item
      Added new internal helper functions:
      \itemize{
        \item
        \code{msg.warn0()}: same as the existing \code{msg.warn()} but without
        white space insertion.
        \item
        \code{is.positive.integers()}: check if an object is a vector of
        positive integers.
        \item
        \code{is.a.character()}: check if an object is a (unique) character.
        \item
        \code{is.a.string()}: check if an object is a string.
        \item
        \code{is.strings()}: check if an object is a vector of strings.
        \item
        \code{is.index()}: check if an object is a positive integer or a string.
        \item
        \code{is.indexes()}: check if an object is a vector of positive integers
        or strings.
        \item
        \code{checkindex()}: check if an object is a valid data frame or matrix
        index; otherwise, an error is raised and an information message is
        displayed.
        \item
        \code{checkindexes()}: check if an object is a vector of valid data
        frame or matrix indexes; otherwise, an error is raised and an information
        message is displayed.
      }
      \item
      Grouped internal \code{is.xxx} helper functions into a single
      TraMineR-is_helpers.R file.
      \item
      Generalized \code{msg.warn()} / \code{msg.warn0()} and \code{msg()} /
      \code{msg0()} code.
      \item
      Internal C/C++ function \code{tmrWeightedInertiaDist()} made public as the R function \code{TraMineRInternalWeightedInertiaDist()}. (Used by package \code{Weightedcluster})
    }
  }
}


\section{Changes in version 2.0-5 [2017-05-13]}{

\subsection{Note:}{
  \itemize{
    \item
    This is a major update of the CRAN version of TraMineR.
    \item
    Check also changes in versions 1.9-14, 2.0-0, 2.0-1, 2.0-2, 2.0-3, and 2.0-4 that have not been released on the CRAN.
    }
}

\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqtrate()}: now accepts a sequence object (\code{seqdata} argument) containing only one sequence.
    \item
    \code{seqdist()}: now accepts a reference sequence object (\code{refseq} argument)
    containing missing values while the main sequence object (\code{seqdata} argument) doesn't.
    }
}

\subsection{Misc:}{
  \itemize{
    \item
    Replaced default deprecated values of \code{norm} in
    \code{seqdistmc()} and \code{seqtree()}.

    \item
    Fixed issues with examples in documentation page of
    \code{seqtree()}.
    }
}
}
\section{Changes in unreleased version 2.0-4 [2017-04-13]}{


\subsection{Changes in existing functions:}{
  \itemize{
    \item
    Renamed several argument names to increase consistency within TraMineR and
    between TraMineR and R. The aim is also to have a common naming
    convention within TraMineR.
    A new internal function (\code{checkargs()})is used to guarantee
    backward compatibility. If the old argument name is used instead of the new
    one, a warning message with an explanation is displayed and the execution
    continues. If the new and old argument names are used together, an error
    message is displayed and the execution stops.
    The following functions have at least one renamed argument:
    \code{dissrep()}, \code{disstree()}, \code{disstree2dot()},
    \code{disstree2dotp()}, \code{disstreedisplay()}, \code{is.eseq()}, \code{is.seqelist()}, \code{seqdiff()},
    \code{seqeconstraint()}, \code{seqecontain()}, \code{seqecreate()}, \code{seqefsub()}, \code{seqeid()},
    \code{seqelength()}, \code{seqelength<-}, \code{seqetm()}, \code{seqeweight()}, \code{seqeweight<-},
    \code{seqlegend()}, \code{seqpcplot()}, \code{seqplot()}, \code{seqrep()}, \code{seqtab()}, \code{seqtrate()},
    \code{seqtree()}, \code{seqtree2dot()}, \code{seqtreedisplay()}, \code{seqeisweighted()},
    \code{plot.seqalign()}, \code{plot.seqdiff()}, \code{plot.stslist()}, \code{plot.stslist.freq()},
    \code{plot.stslist.meant()}, \code{plot.stslist.modst()}, \code{plot.stslist.rep()},
    \code{plot.stslist.statd()}, \code{plot.subseqelistchisq()}. \cr
    See the help page of each of
    these functions for the mapping between old and new argument names.

    \item
    Renamed function: \code{is.seqe()} was renamed as \code{is.eseq()}.
    }
}

\subsection{Misc:}{
  \itemize{
    \item
    Fixed an issue with an example in the documentation page of \code{seqpcplot()}.
    \item
    Fixed equation typing errors in \code{seqrep()} documentation.
    }
}
}

%% =========
\section{Version 2.0-3 [2017-04-06] (Unreleased version)}{

\subsection{Misc:}{
  \itemize{
    \item
    In \code{src/tmrsequence.cpp} from TraMineR 1.8-13:
    \itemize{
      \item
      Fixed two memory errors detected by Valgrind.
      \item
      Fixed a PROTECT error.
      }
    \item
    Fixed issues with examples in documentation pages of
    \code{dissmfacw}, \code{disstree}, \code{seqtree}, and \code{plot.stslist.meant}.
    }
}}

\section{Version 2.0-2 [2017-03-28] (Unreleased version)}{

\subsection{Misc:}{
  \itemize{
    \item
    Changed character encoding from latin1 (ISO-8859-1) to UTF-8.
    \item
    Normalized line endings: LF.
    }
}
}

\section{Version 2.0-1 [2017-03-27] (Unreleased version)}{

\subsection{Misc:}{
  \itemize{
    \item
    Removed the following unused functions: \code{vidx()}, \code{seqmatsaltt()}, \code{seqmathenikoff()}.
    \item
    Removed the functions deprecated in TraMineR 1.x: \code{seqesetlength()}, \code{dissreg()}, \code{dissmfac()}.
    }
}
}

\section{Version 2.0-0 [2017-03-27] (Unreleased version)}{

\subsection{New function:}{
  \itemize{
    \item
    \code{seqcost()}: Evolution of \code{seqsubm()} that offers different ways (CONSTANT,
    TRATE, INDELS, INDELSLOG, FUTURE, FEATURES) to determine indel and
    substitution costs (see the documentation of \code{seqcost()} for details). Unlike
    \code{seqsubm()}, \code{seqcost()} returns both the indel and the substitution costs.
    }
}


\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqdist()}:
      \itemize{
        \item
        New major version with many new features and a new extensible, robust
        and documented implementation (R code).
        \item
        New methods: localized OM (\code{OMloc}), spell length sensitive OM (\code{OMslen}),
        OM of spells (\code{OMspell}), OM of sequences of transitions (\code{OMstran}),
        Time Warp Edit Distance (\code{TWED}), Number of Matching Subsequences (\code{NMS}),
        Number of Matching Subsequences weighted by the Minimum Shared Time (\code{NMSMST}),
        Subsequence Vectorial Representation (\code{SVRspell}), Euclidean distance (\code{EUCLID}),
        Chi-squared distance (\code{CHI2}).

        \item
        New arguments: \code{kweights}, \code{tpow}, \code{expcost}, \code{context}, \code{link}, \code{h}, \code{nu},
        \code{transindel}, \code{otto}, \code{previous}, \code{add.column}, \code{breaks}, \code{step}, \code{overlap},
        \code{weighted}, \code{prox}.

        \item
        \code{sm}: value \code{"CONSTANT"} has been removed for DHD as it doesn't make sense
        and the values \code{"INDELS"} and \code{"INDELSLOG"} have been added (see \code{seqcost()} documentation).

        \item
        \code{norm}: value \code{TRUE} is replaced by \code{"auto"} and \code{FALSE} by \code{"none"}.

        \item
        See the updated documentation of \code{seqdist()} for details.
        }
    \item
    \code{seqsubm()}: This is now an alias for \code{seqcost(...)$sm}.
    }
}
}

\section{Version 1.8-14 [2017-03-20] (Unreleased version)}{

\subsection{Misc:}{
  \itemize{
    \item
      Registration of native routines for R 3.4.x.
      See \url{https://stat.ethz.ch/pipermail/r-devel/2017-February/073755.html}.
      }
}
}

\section{Version 1.8-13}{

\subsection{Misc:}{
  \itemize{
    \item
    Changes in C-code for the \code{seqefsub} function: replaced call
    to function round by a call to \code{fround} to comply with
    forthcoming changes in R 3.4.0 (request of Brian D. Ripley).
    The change is transparent for the user.
    }
}
}

\section{Version 1.8-12}{

\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqST()}: new argument \code{norm} to ask for a normalized
    turbulence index.
    \item
    \code{seqformat()}: The transformation now stops with an error message
    when the columns referenced with the begin and end
    argument contain a non integer value.
    }
}

\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{plot.stslist()}: an unnecessary warning occurred when a
    vector of labels was passed as \code{ytlab} argument.
    \item
    \code{seqdef()}: now accepts to create a state sequence object
    with an alphabet that has only one element. Fixed an error
    that occurred when there was only one state
    }
}

\subsection{New data examples:}{
\itemize{
    \item
    \code{bfspell}: a small data set with 20 sequences in SPELL format.
    }
}

\subsection{Misc:}{
  \itemize{
    \item
    updated \code{seqformat} help page: now includes an example of a
    transformation from SPELL to STS.
    \item
    fixed bad use of \code{extern "C" {}} in TraMineR.h
    (done by B. Ripley, CRAN version 1.8-11.1)
    }
}
}

\section{Version 1.8-11}{

\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqmeant()}: New \code{serr} argument. When \code{serr=TRUE},
    \code{seqmeant} computes
    the variance and standard deviation of the total times spent in the
    different states, and the standard error of the mean total times.
    \item
    \code{seqmtplot()}: When \code{serr=TRUE}, error bars are displayed in the
    mean time plot.
    \item
    \code{seqdist()}:
      \itemize{
        \item
        New error message when \code{sm=NA} with method \code{"OM"}.
        \item
        New error message when \code{refseq} is a state sequence object with
        an alphabet assigned to it different from that of \code{seqdata}.
        }
    }
}

\subsection{Misc:}{
  \itemize{
    \item
    updated \code{disstreedisplay} help page (\code{tree} argument).
    \item
    updated \code{seqdist} help page (\code{refseq} argument and example).
    \item
    updated CITATION file (new ref and fixed a doi argument).
  }
}
}

\section{Version 1.8-10}{

\subsection{Misc:}{
  \itemize{
    \item
    Updated help pages: \code{seqeconstraint}, \code{seqtree}, \code{distree2dot}.
    \item
    Added required basic packages to the import statements to comply with R v3.3 requirements.
    \item
    One additional exported alias to a TraMineR internal function:
    \code{TraMineRInternalSeqgbar}.
    }
}
}

\section{Version 1.8-9}{

\subsection{New function:}{
  \itemize{
    \item
    \code{seqpcfilter()}: convenience function to define the coloring
    filter options to be passed as \code{filter} argument to \code{seqpcplot()}.
    }
}

\subsection{Changes in existing functions:}{
\itemize{
  \item
  \code{seqpm()}: New \code{sep} argument to allow searching for string patterns when
    states are not labelled with single characters.
  \item
  \code{seqpcplot()}:
  \itemize{
    \item
    New argument \code{seed} to control the jittering.
    \item
    The \code{filter} argument can now simply be a scalar, in which case
    the \code{minfreq} filter is applied with this numeric value as threshold.
    See also the new function \code{seqpcfilter()}.
    \item
    New argument \code{missing} to control  whether and how to display missing
    values.
    }
  }
}
\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqformat()}: fixed error occurring when converting from STS to TSE
    with a \code{tevent} matrix containing empty strings (i.e. \code{""})
    \item
    \code{dissmfacw()}: reported F values now obtained by dividing the within
    discrepancy in the denominator by (\eqn{n-m}), where \eqn{n} is the sample size
    and \eqn{m} the total number of predictors (contrasts for categorical
    factors). Up to here (\eqn{n-m-1}) was mistakenly used
    (Reported by Vicente Ponsoda.)
    }
}
\subsection{Misc:}{
  \itemize{
    \item
    Two additional exported aliases to TraMineR internal functions:
    \code{TraMineRInternalSeqeage} and \code{TraMineRInternalLegend}.
    }
}
}

\section{Version 1.8-8}{

\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqecreate()}: an error is thrown when events are not grouped
    by id in inputted TSE data. (Reported by Nicolas Jay).
    This requirement is now specified in \code{seqecreate} help page.
    }
}
\subsection{Misc:}{
  \itemize{
    \item
    exported alias functions allowing other packages to access TraMineR
    internal functions (see \code{?TraMineRInternalLayout}).
    }
}
}

\section{Version 1.8-7}{

\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqpcplot()}: Suppressed unnecessary \code{output} argument. The
    \code{seqpcplot} object is automatically retrieved when using the
    assignment operator, e.g., \code{p <- seqepcplot(...)}.
    }
}
\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqformat()}: When converting from STS to TSE, an error was raised if the
    tevent matrix had empty strings (i.e. \code{""}). Now, this is considered as no event.
    \item
    \code{seqpcplot()}: Fixed error that appeared at the use of \code{"_end"} events.
    \item
    \code{seqpcplot()}: An error occurred when plotting a state sequence object
    (of class \code{stslist}) with a numeric \code{cnames} attribute.
    }
}
\subsection{Misc:}{
  \itemize{
    \item
    Required packages \code{RColorBrewer} and \code{boot} now listed as \code{"imports"}
    in DESCRIPTION and using \code{import} in NAMESPACE.
    }
}
}

\section{Version 1.8-6}{

\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqformat()}: new \code{nr} argument to specify the missing state
    symbol in SPS input.
    \item
    \code{disstreeleaf()}: new logical \code{label} argument to specify
    whether the leaf membership should be labelled with
    the classification rules.
    }
}
\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqtreedisplay()} and \code{disstreedisplay()}: GraphViz installer no longer
    adds GraphViz to the PATH environment variable. Therefore the two
    functions have been adapted to search for GraphViz. In case GraphViz
    would not be found, you can specify the GraphViz installation
    directory with the new \code{gvpath} argument.
    \item
    \code{disstree()} and \code{seqtree()}: removed the warning when \code{R} equals \code{0} or \code{1}  (no permutations).
    \item
    \code{seqformat()}: fixed a problem with missing states when converting
    from SPS to STS (see Changes in function above.)
    \item
    \code{seqpcplot()}: fixed issues with \code{which} argument of \code{plot.seqpcplot}.
    \item
    \code{seqpcplot()}: fixed issues with arguments \code{xlab} and \code{title}.
    \item
    \code{seqpcplot()}: replaced a warning message by an error message at
    failures in finding plot positions for sequences. The error
    message advices to modify the (currently hidden) \code{maxit} argument.
    Additionally, automatically generated subtitles are now hidden
    when the argument \code{title} is used.
    \item
    \code{as.character.seqelist()}, \code{print.seqelist()}: fixed an issue with time
    display in event sequences which was in scientific notation for numbers
    with more than 2 digits. The function now uses the R format function and
    thus accounts for global formatting options such as \code{options(digits=)}
    and/or \code{options(scipen=)}.
    \item
    \code{seqtrate()}: fixed error with sequence objects having only two columns.
    }
}
\subsection{New functionalities:}{
  \itemize{
    \item
    \code{disstree2dot()} and \code{disstree2dotp()} gain a new argument called \code{title.outer}.
    If \code{title.outer=TRUE}, the title is printed in the outer margins.
    }
}
}

\section{Version 1.8-5}{

\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqefsub()}: reported support did not properly account for weights.
    \item
    \code{seqtreedisplay()}: corrected a bug when using representative sequences and
    a \code{dist} object was passed to the \code{dist.matrix} argument (reported by
    Emanuela Struffolino).
    \item
    \code{seqLLCS()} and \code{seqLLCP()}: added a check on the argument. Both sequences
    should belong to state sequence objects with a common alphabet.
    \item
    \code{seqpcplot()}: small change in default lower \code{ylim}.
    }
}
\subsection{Help pages:}{
  \itemize{
    \item
    help pages updated with author and keyword fields.
    }
}
}

\section{Version 1.8-4}{

\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqrep()}, \code{dissrep()}, \code{seqrplot()}: wrong (unweighted) \code{"na"} values were returned;
    also quality measures \code{"MD"} (mean distance to representative) and \code{"V"}
    (discrepancy) were not computed properly when more than one representative
    selected (since version 1.8-2).
    }
}
\subsection{User invisible changes:}{
  \itemize{
    \item
    Added an internal function to fix an issue with an internal C level
    function when called from other packages.
    }
}
}

\section{Version 1.8-3}{

\subsection{Information pages:}{
  \itemize{
    \item
    Updated online help pages.
    \item
    Updated list of references returned by \code{citation("TraMineR")}.
    }
}
\subsection{New functionalities:}{
  \itemize{
    \item
    \code{seqpcplot()}: parallel coordinate plot for sequence data.
    }
}
\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqdss()}: adding long state labels to returned sequence object.
    \item
    \code{seqplot()}: new option \code{type = "pc"}.
    \item
    \code{seqplot()}, \code{seqdplot()}, \code{seqiplot()}, \code{seqfplot()}, \code{seqmsplot()}, \code{seqrplot()}:
    if \code{density} and/or \code{angle} are used to produce shading lines
    instead of solid colors, the legend is plotted using the same
    parameters and thus corresponds to the colors/shades used in the plot.
    \item
    \code{seqplot()}: if \code{group} argument is a factor, the plots are
    now ordered the same way as the factor levels.
    \item
    \code{seqplot()}, \code{plot.stslist()}: if \code{sortv} is a factor, the
    sequences are now sorted according to the order of the factor levels.
    \item
    \code{seqmodst()} and \code{plot.stslist.modst()}: changed name of
    attribute \code{occurences} of object returned by \code{seqmodst} to
    \code{occurrences}. Made resulting changes in \code{plot.stslsit.modst}.
    }
}
\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqdef()}: When selecting subsets of sequence objects using \code{rownames}
    instead of row indexes, the corresponding weights were not selected.
    Fixed by setting (row)names of weight vector as the sequence
    \code{rownames}. (Alexis Gabadinho)
    }
}

}

\section{Version 1.8-2}{

\subsection{Vignettes:}{
  \itemize{
    \item
    Slightly modified JSS article vol. 40(4) added as a vignette
    on state sequence analysis.
    }
}
\subsection{Misc:}{
  \itemize{
    \item
    added \code{ex2} data sets to test and illustrate the handling of weights,
    type \code{help(ex2)} for details.
    }
}
\subsection{New functionalities:}{
  \itemize{
    \item
    New faster interface between C code and R.
    \item
    \code{seqalign()} and associated \code{print} and \code{plot} methods to see computation details
    about the alignment of two state sequences.
    }
}
\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqtree()}, \code{disstree()}: speed improvements.
    \item
    \code{seqtm()}: gives a warning when state names or state labels contain a comma.
    \item
    \code{seqdef()}: changed the display of alphabet, state labels and long labels when
    creating a state sequence object.
    \item
    \code{seqistatd()} and \code{seqmeant()}: added \code{prop} argument to calculate
    proportions of time spent in
    each state instead of absolute values.
    \item
    \code{seqplot()} and aliases: \code{group} now also accepts as argument a list of
    variables/vectors and produces a plot for each combination of the
    values of the variables in the list.
    \item
    \code{dissrep()}, \code{seqrep()}, \code{seqrplot()}: now accounts for weights when present.
    \item
    \code{seqtrate()}, \code{seqsubm()}: added two arguments:
      \itemize{
        \item
        \code{lag}: compute transition rates from \eqn{t} to \eqn{t+}\code{lag}, set to one by default.
        \item
        \code{with.missing}: If \code{TRUE}, compute transition rates to and from missing values.
        }
    \item
    \code{seqtreedisplay()}: now overwrites previous file if filename is not \code{NULL}.
    Tree quality measures displayed with R code.
    \item
    \code{checktriangleineq()}: internal function to check triangle inequality is now
    in C, which allows checking much bigger distance matrices.
    \item
    \code{seqIplot()}, \code{seqiplot()} and \code{seqplot()} with type \code{"I"} or \code{"i"}: the \code{sortv} argument now also accepts a sorting method,
    namely one of \code{"from.start"} or \code{"from.end"}.
    See the help page \code{?plot.stslist} for explanation.
    \item
    \code{seqeconstraint()} and other \code{seqe...} functions for event sequences: support of subsequences can now be determined by means of any of
    Joshi's 5 counting methods (see the ref manual page). The method should be
    specified with \code{seqeconstraint()}.
    \item
    \code{seqeapplaysub()}: when \code{method=NULL} is specified (now the default), the count method assigned to the event sequence object is used. With \code{method="count"} CDIST_O (number of distinct occurrences) is used as previously.
    \item
    \code{seqrep()}: attribute \code{Index} of the returned object is now a vector instead of an object of class \code{dissrep}.
    }
}
\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqtreedisplay()}: was changing current directory when an error occurred in
    the plotting function.
    \item
    \code{dissrep()}, \code{seqrep()}: error when \code{nrep} cannot be reached (reported by M. Studer)
    \item
    \code{seqefsub()}: when using \code{strsubseq} argument, \code{countMethod} of \code{seqeconstraint} was not taken into account (Reported by Reto Bürgin).
    \item
    \code{seqeconstraint()}: added consistency checks to avoid misuses (Reported by Reto Bürgin).
    \item
    \code{print.seqelist()} and \code{as.character.seqelist()}: generated \code{segfault} when converting long event sequences to character (Reported by Pierre Molinier).
    \item
    \code{seqsubm()}: Very small rounding errors (\code{1e-16}) were sometimes leading to non symmetric substitution cost matrix (Reported by Alexandre Pollien).
    }
}
}

\section{Version 1.8-1}{

\subsection{Misc:}{
  \itemize{
    \item
    updated references in the citation file and manual pages to point to the newly published article in the \emph{Journal of Statistical Software}
    \item
    other references update in the manual pages.
    }
}

}

\section{Version 1.8}{

\subsection{New functionalities:}{
  \itemize{
    \item
    \code{seqtreedisplay()}: drawing a sequence regression tree.
    \item
    \code{seqtree()}: creating a sequence regression tree from a dissimilarity matrix.
    \item
    \code{seqrecode()}: recoding state sequences objects (i.e., merging states).
    \item
    weights are now supported by all dissimilarity analysis functions.
    \item
    weights can be assigned to event sequence objects and are
    supported by all related functions.
    \item
    \code{seqdef()}: \code{xtstep} option added to set step between displayed tick-marks and labels on the x-axis of state sequence plots.
    \item
    \code{seqplot()}: \code{xtstep} option added to state sequence plots.
    }
}
\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqformat()}: fixed problems with \code{fillblanks} argument when converting from SPELL to STS.
    \item
    \code{seqdist()} and \code{seqdistmc()}: \code{method="HAM"} did not account for the provided substitution cost matrix; it used 1 for all substitution costs (Reported by Florian Hertel).
    \item
    \code{seqdist()}: fixed a (possible) memory leak.
    \item
    \code{seqdss()}, \code{seqdur()}: fixed bad handling of missing values in several cases:
      \itemize{
        \item
        sequences finishing with missing values,
        \item
        sequences made of only one distinct state and missing values.
        }
    \item
    \code{seqiplot()}, \code{seqIplot()}, \code{seqfplot()}: changed the automatic setting of the x-axis length, to ensure identical lengths of the x-axis when the maximal sequence length differs between groups (reported by M. Studer).
    \item
    \code{seqplot()} and aliases: fixed error with \code{xaxis} argument.
    \item
    \code{seqtransn()}: the returned normalized number of transitions for
    sequences of length 1 was \code{NaN} (value of \code{transn.norm=0/seqlength-1}). Now set to 0.
    \item
    \code{seqici()}: returns now correctly 0 instead of \code{NaN} for sequences of length 1.
    }
}
\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqient()}: new option \code{base} for choosing the base of the logarithm used
    to compute the entropy.
    \item
    \code{seqdist()}: enhanced check of substitution cost matrix.
      \itemize{
        \item
        The function can now be cleanly interrupted by the user.
        \item
        Timing information now uses the processor time (instead of the elapsed time).
        }
    \item
    \code{seqformat()} when converting from "SPELL":
      \itemize{
        \item
        new error message when a start time is lower than 1 and/or an end time is smaller than the start time.
        \item
        warning message when start time of episode 1 is missing (sequence creation is skipped)
        \item
        warning message when start/end time of an episode is missing (episode is skipped and filled with NA's)
        }
    \item
    \code{disstree()} and \code{dissassoc()} have been entirely redesigned; objects created with the old function are no longer supported.
    \item
    \code{seqdss()} and \code{seqdur()}: the number of columns of the returned object is now set to the maximum DSS length rather than to the length of the original state sequence object.
    }
}
\subsection{Misc:}{
  \itemize{
    \item
    \code{ex1} example data set: contains now an additional sequence \code{'s7'}
        with only missing values.
    }
}
}

\section{Version 1.6-2}{

\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqdef()}: now checks whether all states encountered in the input data are present in an optional user provided alphabet (\code{state} argument)
    \item
    \code{seqefsub()}: The support of a subsequence with a total support of 1 was set to 0 (hence, this only applies if the minimum support is 1.) This is now corrected (Reported by Anna Hera).
    \item
    Corrected compilation problems under SOLARIS (Reported by Prof Brian Ripley).
    \item
    \code{seqformat()}: when converting from SPELL to STS, the \code{fillblank} argument was not used.
    }
}
}

\section{Version 1.6-1}{

\subsection{New functionalities:}{
  \itemize{
  \item
  Added a startup message when loading the TraMineR library.
  }
}
\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqtab()} and \code{seqfplot()}: \code{tlim} argument now allows to return any selected frequent sequences,
    in the same way as the \code{tlim} argument used in \code{plot.stslist()} and \code{seqiplot()} (Requested by G. Ritschard).
    For example, \code{tlim=3:6} returns the third, fourth, fifth and sixth most frequent sequences in the set. Default for \code{tlim} is now \code{1:10} instead of \code{10}.
    \item
    \code{seqsubsn()}: added detection of missing state in the sequences and computation of number of subsequences by adding missing state to the alphabet.
    \item
    \code{seqST()}: added detection of missing state in the sequences and computation of turbulence by adding missing state to the alphabet (Requested by G. Ritschard).
    \item
    \code{seqrplot} now accepts \sQuote{half} matrices (\code{dist} objects) as produced by \code{seqdist}  with the \code{full.matrix=FALSE} option (Requested by L. Lesnard).
    \item
    \code{seqiplot()} and \code{plot.stslist()}: new \code{ytlab} option allowing to display sequence labels
    on the Y-axis in sequence index plots (if set to \code{"id"}, the sequence ids are displayed).
    An additional \code{ylas} option sets the orientation of the labels (Requested by Andrew ??
    and P. Jeuniaux).
    \item
    \code{seqsubm()}:
    \itemize{
      \item
      When \code{method="TRATE"}, the substitution costs are now based on the value of \code{cval}: \eqn{SC(i,j) = cval -P(i,j) - P(j,i)} where \eqn{P(i,j)} is the transition rate from state \eqn{i} to \eqn{j}.
      \item
      added a new \code{transition} argument to use only transition from \code{"previous"} or \code{"next"} state
      instead of the default \code{"both"} when \code{time.varying = TRUE}.
      \item
      Now, by default \code{cval} equals 2, unless \code{transition} is set to \code{both} and
      \code{time.varying} is \code{TRUE} in which case \code{cval} equals 4.
      }
    \item
    \code{plot.stslist.meant()}, \code{seqmtplot()}, \code{plot.stslist.freq()}, \code{seqfplot()}, \code{plot.stslist.modst()},
    \code{seqmsplot()}: added display of weighted \eqn{n} instead of \eqn{n} in axis label if \code{weighted=TRUE}.
    \item
    \code{seqdplot()}: disabled plot of a legend for missing state if \code{with.missing=FALSE} (Requested by
      M. Studer and G. Ritschard).
    \item
    \code{seqplot()}, \code{plot.stslist()}, \code{plot.stslist.statd()}, \code{plot.stslist.modst()}, \code{plot.stslist.freq()},
    \code{plot.stslist.meant()}: when \code{weighted=TRUE}, weighted \code{n} displayed in the axis label rounded to 2 digits (Requested by M. Studer).
    }
}
\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{plot.stslist.modst()}, \code{seqmsplot()}: fixed bad display of missing states if \code{with.missing=TRUE}.
    \item
    \code{seqmeant()}: added color for missing state to color palette when using \code{with.missing=TRUE}.
    \item
    \code{dissrep()}, \code{seqrep()}, \code{seqrplot()}: fixed bad coverage statistics when \code{tsim} set to other than default value of 0.10 (Reported by G. Ritschard).
    \item
    \code{seqdist()}: Bad handling of missing values in DHD distance (wrong default substitution costs with missing values, set to one instead of four).
    \item
    \code{seqdist()}: Fixed an error in \code{seqdist} when using \code{refseq} with missing values (Reported by an anonymous user).
    }
}
}

\section{Version 1.6}{

\subsection{New functionalities:}{
  \itemize{
    \item
    \code{seqIplot()} and \code{seqplot(..., type="I")}: sequence index plot displaying all sequences with no space (\code{space=0}) and no border (\code{border=NA}) by default.
    \item
    \code{stlab()}: retrieving or setting the long state labels of a sequence object.
    \item
    \code{seqici()}: computes the complexity index, a composite measure of sequence complexity.
    \item
    \code{seqtransn()}: computes the number of transitions in a sequence.
    }
}
\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqrplot()} and \code{plot.stslist.rep()}: missing states are now correctly
    plotted (reported by G. Ritschard).
    \item
    \code{dissrep()} (called by \code{seqrplot()} and \code{plot.stslist.rep()}): when using coverage threshold, does no longer select one too many representatives.
    \item
    \code{seqdplot()} / \code{plot.stslist.statd()}: fixed bad coloring of missing states (reported by M. Studer).
    \item
    \code{seqstatd()} and \code{print.stslist.statd()}: fixed error when printing \code{seqstatd} output and length of the longest short state label was >2. (reported by G. Ritschard).
    \item
    \code{seqdss()}: (with default \code{with.missing=FALSE} value) DSS with identical
    successive states resulted when there were missing values between two identical states.
    For example, before: DSS of \code{"A-A-*-A-*-A"} was \code{"A-A-A"} and now \code{"A"}.
    }
}
\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{summary.stslist()}: added information in the output of the summary method
      for state sequence objects.
    \item
    \code{seqdist()}, \code{seqdistmc()}, \code{seqdss()}, \code{seqdur()},  \code{seqient()}, \code{seqistatd()}, \code{seqmpos()},
    \code{seqnum()}, \code{seqsubm()}:
      \code{with.miss} argument replaced by \code{with.missing} (obsolete \code{with.miss} argument
      still works for backward compatibility)
    \item
    \code{seqmeant()}: added \code{with.missing} argument to account for missing states.
    \item
    \code{seqdef()}:
    \itemize{
      \item
      new warning when one or several sequences contain only missing values.
      \item
      when no weights are provided:
      \itemize{
        \item
        \code{weights} attribute of the returned sequence object is now set to \code{NULL} instead to a vector of 1's;
        \item
        the \code{"[>] sum of weights"} message is suppressed;
        }
      \item
      new \code{Version} attribute with the number of the TraMineR version
        used for creating the sequence object;
      \item
      message \code{"[>] missing values in input file"} changed to
        \code{"[>] found missing values ('...') in sequence data"}
        and displayed only if missing values are found in the input data.
      }
    \item
    \code{seqIplot()}, seqiplot(), seqplot(..., type="I"), seqplot(..., type="i") and plot.stslist():
      new \code{weighted} argument. If set as \code{TRUE} sequence bar widths are set
        proportional to weights.
    \item
    \code{seqiplot()} and \code{plot.stslist()}: when the sequence object
    contains less than
      10 sequences and \code{tlim=NULL}, only the actual sequences are plotted
      without additional \sQuote{void sequences}.
    \item
    \code{seqmodst()}: added \code{weighted} and \code{with.missing} arguments.
    \item
    \code{seqtab()} and \code{seqstatd()}: attribute \code{nbseq} of the returned object is now
      the sum of weights (if \code{weighted=TRUE} and the sequence object has weights)
      instead of the number of sequences.
    \item
    \code{seqfplot()} / \code{plot.stslist.freq()} and \code{seqdplot()} / \code{plot.stslist.statd()}:
      when weights are used, i.e. when the \code{weighted=TRUE} argument is passed
      and the sequence object has weights,
      \itemize{
        \item
        the \sQuote{n=...} in the y axis legend now gives the sum of weights instead of
        the number of sequences;
        \item
        \sQuote{weighted} is specified in the y axis legend;
        }
    \item
    \code{seqfplot()} and \code{plot.stslist.freq()}: more precise positioning of the 0 of the y axis;
    \item
    \code{seqrplot()} and \code{plot.stslist.rep()}:
      \itemize{
        \item
        new \code{stats} option. If set to \code{FALSE}, statistics are not plotted;
        \item
        label of the yaxis changed.
        }
    \item
    \code{dissrep()}: name of the main argument changed from \code{dist.matrix} to \code{diss}  as in other \code{diss...} functions.
    }
}
}

\section{Version 1.4-1}{

\subsection{Bug fixes:}{
  \itemize{
    \item
    \code{seqecreate()}: Problems with handling simultaneous events when
    creating event sequences with data not previously sorted on
    \code{id}, \code{timestamp} and \code{event}.
    \item
    \code{seqrep()} with \code{"density"} criterion: neighborhood diameter is now
    correctly set to \code{trep*dmax} instead of \code{trep}.
    }
}
\subsection{New functionalities:}{
  \itemize{
    \item
    \code{dissrep()}: extracts a set of representative objects using a
    dissimilarity matrix. This function is used by \code{seqrep}.
    }
}
\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqrep()}:
      \itemize{
      \item
      much faster extraction of the representative set;
      \item
      default criterion is now \code{"density"} instead of \code{"frequency"};
      \item
      \code{trep} now sets a coverage threshold for the representative set
       rather than a size threshold for the candidate list;
      \item
      in the output object, name of the attribute containing
       statistics for the representative set changed from \code{Quality}
       to \code{Statistics} and that of the attribute containing the
       overall quality measure changed from \code{rindex} to \code{Quality}.
       }
    \item
    \code{seqplot()}: now checks that the length of the vector given as
    \code{group} argument matches the number of sequences.
    \item
    \code{plot.stslist()}: now checks that the length of the vector given as
     \code{sortv} argument matches the number of sequences.
    \item
    \code{seqformat()}: dramatic speed improvement in conversion from SPELL data.
    }
}
}

\section{Version 1.4}{

\subsection{New functionalities:}{
  \itemize{
    \item
    \code{TraMineR.checkupdates()}: check for TraMineR updates.
    \item
    \code{seqdistmc()}: computes multichannel distances.
    \item
    \code{seqmeant()}: computes mean duration in each state.
    \item
    \code{seqmodst()}: returns the sequence of modal states.
    \item
    \code{seqmsplot()}: for plotting the sequence of modal states. This
    function is a shortcut for \code{seqplot} with \code{type="ms"}, see below.
    \item
    \code{seqrep()}: extracts a set of representative sequences.
    \item
    \code{seqrplot()}: for producing representative sequence plots. This function is
    a shortcut for \code{seqplot} with \code{type="r"}, see below.
    \item
    \code{seqHtplot()}: for producing Entropy Index plots. This function is
    a shortcut for \code{seqplot} with \code{type="Ht"}, see below.
    \item
    \code{seqlogp()}: Computing logarithm of sequence probabilities.
    \item
    \code{seqdef()}: new \code{weights=} option for providing a vector of weights.
    \item
    \code{seqstatd()}, \code{seqtrate()}, \code{seqlogp()} and \code{seqtab()}: new option \code{weighted=TRUE}
    for using the weights when computing the statistics.
    \item
    \code{seqtrate()} and \code{seqsubm()}: new \code{time.varying} argument for computing
    position dependant transition rates or costs.
    \item
    \code{seqdist()}: two additional methods are now available for computing
    distances, namely \code{"HAM"} (Hamming distance) and \code{"DHD"} (Dynamic
    Hamming Distance).
    \item
    Output produced by \code{seqstatd()}, \code{seqtab()}, \code{seqmeant()}, \code{seqmodst()}, and
    \code{seqrep()} can now be plotted with their \code{plot()} dedicated methods (see
    new classes and methods below).
    }
}
\subsection{Important changes:}{
  \itemize{
    \item
    \code{seqplot()}: is now the generic function for plotting state sequence objects with
    \code{type} argument. Available types are \code{"d"} for state distribution plots,
    \code{"f"} for sequence frequency plots, \code{"Ht"} for entropy index plots,
    \code{"i"} for sequence index plots, \code{"ms"} for modal state sequence plots,
    \code{"mt"} for meant time plots, \code{"r"} for representative plots.
    This function replaces the previous generic \code{plot.stslist()} function.
    \item
    \code{plot.stslist()}: now produces only a sequence index plot (see new classes and
    methods below).
    }
}
\subsection{Changes in existing functions:}{
  \itemize{
    \item
    \code{seqfplot()}:

      \itemize{
        \item
        new \verb{yaxis} option: with \code{yaxis="cum"} (default) cumulated percentages
        are displayed, while with \code{yaxis="pct"} individual sequence percentages
        are shown.
        \item
        \code{pbarw=TRUE} is now the default for the \code{pbarw} argument.
        }

    \item
    \code{seqtab()}: the \code{format} argument that specifies the format of the sequences
    displayed as rownames is now set by default to the short SPS format,
    e.g. \code{TR/9-EM/63}.
    \item
    \code{seqiplot()}: sequence indexes are now displayed by default on the y axis.
    This can be disabled with \code{yaxis=FALSE}.
    }
}
\subsection{Fixed minor bugs in seqformat():}{
   \itemize{
      \item
      changes concern mainly the \code{from="SPELL"} and
      \code{from="SPS"} options.
   }
}
\subsection{New classes and methods:}{
  \itemize{
    \item
    new class \code{stslist.statd} for objects produced by the \code{seqstatd()} function
    and methods for printing and plotting such objects.
    \item
    new class \code{stslist.freq} for objects produced by \code{seqtab()} function and
    methods for printing and plotting such objects.
    \item
    new class \code{stslist.meant} for objects produced by \code{seqmeant()} function
    and methods for printing and plotting such objects.
    \item
    new class \code{stslist.modst} for objects produced by \code{seqmodst()} function
    and methods for printing and plotting such objects.
    }
}

}

\section{Version 1.2-1}{

\subsection{Changes in function arguments:}{
  \itemize{
    \item
    \code{seqdef()}: new \code{id} argument for setting the rownames of the sequence object.
    \item
    \code{disscenter()}: new \code{medoid.index} argument to get the indexes of all
    medoids (rather than only the first one).
    }
}
\subsection{Minor bugs fixed:}{
  \itemize{
    \item
    Plotting missing states with \code{seqiplot()} and \code{seqfplot()} functions.
    \item
    Sum of transition rates with sequences of different lengths not equal to 1.
    }
}
}

\section{Version 1.2}{

\subsection{Changes regarding plotting functions:}{
  \itemize{
    \item
    New generic function \code{plot.stlist()} with option \code{type=} for
    plotting state sequence objects of class \code{ststlist} created by
    the \code{seqdef()} function. Old functions \code{seqdplot()}, \code{seqfplot()},
    \code{seqiplot()}, \code{seqmtplot()} work as in the previous version but by
    calling \code{plot.stslist} with the appropriate \code{type=} option
    (types are \code{'d'}, \code{'f'}, \code{'i'} and \code{'mt'}). However, the order of the
    functions' arguments may have changed and this may cause problems
    if the names of the arguments were not explicitly specified in your
    scripts (which is inadvisable anyway).
    }
}
\subsection{Changes in the dissimilarity and discrepancy analysis diss module:}{
  \itemize{
    \item
    \code{dissreg()} is renamed to \code{dissmfac()} for multi-factor dissimilarity analysis.
    \item
    \code{disstree()}: great speed and memory improvements.
    }
}
\subsection{Speed and memory improvements:}{
  \itemize{
    \item
    \code{seqformat()}: conversion to TSE format is now much faster.
    \item
    \code{seqdur()}
    \item
    \code{seqST()}
    }
}
\subsection{Changes in the computation of distances between sequences:}{
  \itemize{
    \item
    \code{seqdist()} now checks if substitution costs respect the triangle inequality.
    When this is the case it ensures that the resulting dissimilarity matrix
    also respects the triangle inequality.
    \item
    New options for selecting the distance normalization method.
    \item
    New reversed LCP, i.e. longest common suffix method. (\code{method="RLCP"})
    }
}
\subsection{Various changes:}{
  \itemize{
    \item
    CITATION file added.
    \item
    \code{seqST()} caused an error when run with more than 12 states: bug fixed.
    \item
    Fixed other minor bugs.
    }
}

}

\section{Version 1.1}{

\subsection{New example data}{
  \itemize{
    \item
    \code{mvad} example data set added
    }
}
\subsection{Name changes of functions}{
  \itemize{
    \item
    \code{seqLLCP()} instead of old \code{seqLCP()}
    \item
    \code{seqLLCS()} instead of old \code{seqLCS()}
    }
}
\subsection{New module}{
New \code{diss} module for analysing a dissimilarity matrix (such as the one returned by \code{seqdist})
  \itemize{
    \item
    \code{dissassoc()}: Computes association with a factor;
    \item
    \code{dissreg()}: Regression analysis of a dissimilarity matrix;
    \item
    \code{disstree()}: Tree analysis of a dissimilarity matrix;
    \item
    \code{dissvar()}: Computes a pseudo-variance from a dissimilarity matrix.
    }
}

\subsection{Changes in graphic functions}{
 Changes in the graphics functions (seqiplot, seqfplot, seqdplot, ...):
  \itemize{
    \item
    New \code{group} option allowing to draw several plots for the levels
    of a factor within a single command;
    \item
    TraMineR now uses \code{layout} for controlling the position of the
    plots and the legend in the graphic area. This is not compatible
    with \code{par(mfrow=...)}. To use the standard \code{par(mfrow=...)} method,
    one must set the \code{use.layout=FALSE} option in the plot functions;
    \item
    \code{seqmtplot()}: new function that plots the mean time spent in each state;
    \item
    Option \code{withborder=FALSE} is now obsolete and replaced by the standard \code{border=NA} option.
    }
}

\subsection{Changes in event sequence (seqe) module:}{
  \itemize{
    \item
    New plot for exhibiting discriminant subsequences (\code{seqecmpgroup}).
    \item
    Overall syntax have been reviewed and is now much simpler.
    \item
    \code{seqefsub} now allows to search for user specified subsequences (\code{seqefsub}).
    \item
    Event subsequence lists now have specific plot and print method (\code{seqefsub}, \code{seqecmpgroup})
    \item
    \code{seqecreate} now accepts state sequences and performs automatic conversion (\code{seqecreate})
    \item
    Time constraints are now implemented separately and are stored with the results (\code{seqeconstraint})
    }
}
\subsection{Changes in the seqformat() function:}{
  \itemize{
    \item
    New options for importing SPELL formatted data.
    \item
    The STS internal and the output in STS, SPS or DSS formats are now by
    default in extended format (a matrix with one state per column) instead of
    compressed format (a character string).
    Use the \code{compressed=TRUE} option to get an output in the compressed format
    (sequences as character strings).
    \item
    \code{"SPS1"} and \code{"SPS2"} formats are now replaced by the generic \code{"SPS"} format with
    options \code{SPS.in} and \code{SPS.out} for defining the separator and surrounding
    characters used for specifying the state/duration couples.
    }
}
\subsection{Other new functions:}{
  \itemize{
    \item
    \code{seqgen()}: generates a random sequence.
    }
}
}
