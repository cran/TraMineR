\name{disstreeleaf}
\alias{disstreeleaf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Terminal node membership}
\description{
  Return a factor with the terminal node membership of each case.
}
\usage{
disstreeleaf(tree, label=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{The tree, a \code{disstree} or \code{DissTreeNode} object.}
  \item{label}{Logical. If \code{TRUE}, the returned leaf memberships are labelled with the corresponding classifications rules.}
}


\seealso{\code{\link{disstree}} for examples}
\author{Matthias Studer (with Gilbert Ritschard for the help page)}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dissimilarity-based analysis}
