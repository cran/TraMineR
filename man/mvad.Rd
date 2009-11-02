\encoding{latin1}
\name{mvad}
\docType{data}
\alias{mvad}
\title{Example data set: Transition from school to work}
\description{
The data comes from a study by McVicar and Anyadike-Danes on transition from school to work. 
The data consist of static background characteristics and a time series sequence of 72 monthly 
labour market activities for each of 712 individuals in a cohort survey. 
The individuals were followed up from July 1993 to June 1999.

States are:\cr
\cr
employment  (EM) \cr
FE = further education (FE)\cr
HE = higher education (HE) \cr
joblessness  (JL) \cr
school (SC) \cr
training (TR)\cr
\cr
The data set contains also ids and sample weights as well as the following binary covariates:\cr
\cr
male\cr
catholic\cr
Belfast, N.Eastern, Southern, S.Eastern, Western (location of school, one of five Education and Library Board areas in Northern Ireland)\cr
Grammar (type of secondary education, 1=grammar school)\cr
funemp (father's employment status at time of survey, 1=father unemployed)\cr
gcse5eq (qualifications gained by the end of compulsory education, 1=5+ GCSEs at grades A-C, or equivalent)\cr  
fmpr (SOC code of father's current or most recent job, 1=SOC1 (professional, managerial or related))\cr
livboth (living arrangements at time of first sweep of survey (June 1995), 1=living with both parents)\cr 
}



\usage{data(mvad)}
\format{A data frame containing 712 rows, 72 state variables, 1 id variable and 13 covariates.}
\source{McVicar and Anyadike-Danes (2002)}
\references{
McVicar, Duncan and Anyadike-Danes, Michael (2002). Predicting Successful and Unsuccessful Transitions from School to Work by Using Sequence Methods, 
\emph{Journal of the Royal Statistical Society. Series A (Statistics in Society)}, 165, 2, pp. 317--334.
}
\keyword{datasets}
