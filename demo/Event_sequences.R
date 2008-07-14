library(TraMineR)
data(actcal)
transition <- seqetm(actcal[,13:24], method="transition")
transition[1,1:4] <- c("FullTime"         , "Decrease,PartTime",
     "Decrease,LowPartTime", "Stop")
transition[2,1:4] <- c("Increase,FullTime", "PartTime"         ,
     "Decrease,LowPartTime", "Stop")
transition[3,1:4] <- c("Increase,FullTime", "Increase,PartTime",
    "LowPartTime"         , "Stop")
transition[4,1:4] <- c("Start,FullTime"   , "Start,PartTime"   ,
    "Start,LowPartTime"   , "NoActivity")
transition
actcal.tse <- seqformat(actcal,var=13:24,from='STS',to='TSE', tevent=transition)


actcal.seqe <- seqecreate(id=actcal.tse$id,
time=actcal.tse$time, event=actcal.tse$event)

sl <- numeric()
sl[1:2000] <- 12
# All sequences are of length 12
seqesetlength(actcal.seqe,sl)
actcal.seqe[1:10]

fsubseq <- seqefsub(actcal.seqe,minSupport=100)
msubcount<-seqeapplysub(fsubseq$subseq, actcal.seqe, method="count")
#First lines...
msubcount[1:9,6:7]
## Using time constraints
## Searching subsequences starting in summer (between June and September)
fsubseq <- seqefsub(actcal.seqe, minSupport=10, ageMin=6, ageMax=9)
fsubseq$subseq[1:10]
## Searching subsequences occurring in summer (between June and September)
fsubseq <- seqefsub(actcal.seqe, minSupport=10, ageMin=6, ageMax=9,
                    ageMaxEnd=9)
fsubseq$subseq[1:10]
## Searching subsequences enclosed in a 6 months period
## and with a maximum gap of 2 months
fsubseq <- seqefsub(actcal.seqe, minSupport=10, maxGap=2, windowSize=6)
fsubseq$subseq[1:10]

## Looking for frequent subsequences
fsubseq <- seqefsub(actcal.seqe, pMinSupport=0.01)
## Frequences of 10 first subsequences
seqefplot(fsubseq$subseq[1:10], actcal.seqe, col="cyan")
seqefplot(fsubseq$subseq[1:6],actcal.seqe, group=actcal$sex, mfrow=c(2,3),col="cyan")
## looking for subsequence with FullTime
seqecontain(fsubseq$subseq, c("FullTime"))
## Looking for subsequences that are present in at least 1% (20) of the sequences
fsubseq <- seqefsub(actcal.seqe, pMinSupport=0.01)

## Looking for the discriminating subsequences for sex
discr <- seqecmpgroup(fsubseq$subseq, actcal.seqe, group=actcal$sex,
                    method="bonferroni")
## Plotting the six most discriminating subsequences in 2 x 4 format
seqefplot(fsubseq$subseq[discr$index[1:8]], actcal.seqe,
          group=actcal$sex, mfrow=c(2,4), col="cyan")
