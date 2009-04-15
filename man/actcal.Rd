\name{actcal}
\docType{data}
\alias{actcal}
\title{Example data set: Activity calendar from the Swiss Household Panel}
\description{
This data set contains individual monthly activity statuses from
January to December 2000. It is a subsample of data collected by the Swiss Household Panel (SHP).

The state column (variable) names are `jan00', `feb00', etc...

There are four possible states:\cr
\cr
A = Full-time paid job (> 37 hours)\cr
B = Long part-time paid job (19-36 hours) \cr
C = Short part-time paid job (1-18 hours) \cr
D = Unemployed (no work) \cr

The data set contains also the following covariates:\cr
\cr
age00 (age in 2000) \cr
educat00 (education level) \cr
civsta00 (civil status) \cr
nbadul00 (number of adults in household) \cr
nbkid00 (number of children) \cr
aoldki00 (age of oldest kid) \cr
ayouki00 (age of youngest kid) \cr
region00 (residence region) \cr
com2.00 (residence commune type) \cr 
sex (sex of respondent)\cr
birthy (birth year) \cr
}
\usage{data(actcal)}
\format{A data frame with 2000 rows, 12 state variables, 1 id variable and 11 covariates.}
\source{Swiss Household Panel}
\references{\url{www.swisspanel.ch}}
\keyword{datasets}
