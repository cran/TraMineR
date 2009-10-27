\encoding{latin1}
\name{biofam}
\docType{data}
\alias{biofam}
\title{Example data set: Family life states from the Swiss Household Panel biographical survey}
\description{
The \emph{biofam} data set was constructed by Müller et al. (2007)
from the data of the retrospective biographical survey carried out by the Swiss Household Panel (SHP) in 2002. 
The data set contains sequences of family life states from age 15 to 30 (sequence length is 16) and a series of covariates. The sequences are a sample of 2000 sequences of those created from the SHP biographical survey. It includes thus only individuals who were at least 30 years old at the time of the survey.  The \emph{biofam} data set describes thus family life courses of 2000 individuals born between 1909 and 1972.

The states numbered from 0 to 7 are defined from the combination of five basic states, namely Living with parents (Parent), Left home (Left), Married (Marr), Having Children (Child), Divorced:\cr
\cr
0 = "Parent" \cr
1 = "Left" \cr
2 = "Married" \cr
3 = "Left+Marr" \cr
4 = "Child" \cr
5 = "Left+Child" \cr
6 = "Left+Marr+Child" \cr
7 = "Divorced" \cr

The covariates are: \cr
\cr
sex \cr
birthyr  (birth year)\cr
nat\_1\_02 (first nationality) \cr
plingu02 (language of questionnaire) \cr
p02r01   (religion)\cr
p02r04   (religious participation) \cr
cspfaj   (father's social status) \cr
cspmoj   (mother's social status) \cr
\cr
Two additional weights variables are inserted for illustrative purpose ONLY (since biofam is a subsample of the original data, these weights are not adapted): \cr
\cr
wp00tbgp (weights inflating to the swiss population)\cr
wp00tbgs (weights keeping sample size) \cr
}
\usage{data(biofam)}
\format{A data frame with 2000 rows, 16 state variables, 1 id variable and 7 covariates and 2 weights variables.}
\source{Swiss Household Panel \url{www.swisspanel.ch}}
\references{
  Müller, N. S., M. Studer, G. Ritschard (2007).
  Classification de parcours de vie à l'aide de l'optimal matching.
  In \emph{XIVe Rencontre de la Société francophone de classification
  (SFC 2007), Paris, 5 - 7 septembre 2007}, pp.\  157--160.
}
\keyword{datasets}
