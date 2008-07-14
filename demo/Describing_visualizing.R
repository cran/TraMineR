## ======================================
## Examples for the chapter
## 'Describing and visualizing sequences'
## in TraMineR User's Guide
## ======================================
library(TraMineR)

## Set 'graphdir' according to your system 
## and uncomment the pdf() and dev.off() commands 
## if you want to save the graphics as pdf files
graphdir <- "/****/*****/****/*****/TraMineR-Manual/Graphiques/"

## ---------------------------------------
## List of states present in sequence data
## ---------------------------------------
data(actcal)
stated <- c("> 37 hours", "19-36 hours", "1-18 hours", "no work")
actcal.seq <- seqdef(actcal,13:24,labels=stated)
alphabet(actcal.seq)

sp.ex1 <- rbind("(000,12)-(0W0,9)-(0WU,5)-(1WU,2)", 
	"(000,12)-(0W0,14)-(1WU,2)")
sp.ex1
seqstatl(sp.ex1, format='SPS')

## ------------------
## State distribution
## ------------------
data(biofam)
bfstates <- c("N/N/N/N", "Y/N/N/N", "N/Y/*/N", "Y/Y/N/N",
			"N/N/Y/N", "Y/N/Y/N", "Y/Y/Y/N", "*/*/*/Y")

biofam.seq <- seqdef(biofam,10:25,labels=bfstates)

## pdf(file=paste(graphdir,"seqdplot-biofam.pdf",sep=""))
seqdplot(biofam.seq)
## dev.off()

sd <- seqstatd(biofam.seq)
agelab <- paste("a",seq(15,30),sep="")

## pdf(file=paste(graphdir,"fg_biofam-entropy.pdf",sep=""), width=8, height=6, pointsize=14)
plot(sd$Entropy,
	main="Entropy of biofam state distribution by age",
	xlab="Age", ylab="Entropy", type="h", lwd=3.5, col="blue",
	xes=F, cex=1.3, frame.plot=T)
axis(1,labels=agelab, at=1:16)
axis(2, at=c(0.0,.2,.4,.6,.8))
## dev.off()

## --------------------
## Sequence frequencies
## --------------------

## pdf(file=paste(graphdir,"actcal-seqfplot.pdf",sep="")
seqfplot(actcal.seq)
## dev.off()

seqfplot(actcal.seq, pbarw=TRUE)  

data(biofam)
biofam.seq <- seqdef(biofam,10:25, labels=bfstates)

seqfplot(biofam.seq, pbarw=TRUE)

seqtab(actcal.seq, tlim=10)

seqtab(actcal.seq[,7:9], tlim=10)

## --------------------
## Sequence frequencies
## --------------------
seqtrate(actcal.seq)
tr <- seqtrate(actcal.seq)
tr

rowSums(seqtrate(actcal.seq))

seqiplot(actcal.seq)

seqiplot(actcal.seq, sortv=actcal$sex, tlim=0,
         withborder=FALSE, space=0)

## -----------------------------
## State frequencies by sequence
## -----------------------------
seqistatd(actcal.seq[1:10,])

statd <- seqistatd(actcal.seq)
apply(statd,2,mean)

barplot(apply(statd,2,mean), main="Mean time spent in each state",
        col="lightblue")

## -----------------------------
## State frequencies by sequence
## -----------------------------
print(actcal.seq[1:10,],"SPS")

seqdss(actcal.seq[1:10,])

seqdur(actcal.seq[1:10,])

## ---------------
## Sequence length
## ---------------
data(famform)
famform.seq <- seqdef(famform)
famform.seq
seqlength(famform.seq)

## --------------------------------
## Searching for given subsequences
## --------------------------------
seqpm(actcal.seq,"DAAD")

daad <- seqpm(actcal.seq,"DAAD")
actcal.seq[daad$MIndex,]

## -----------------------
## Within sequence entropy
## -----------------------
seqient(actcal.seq[1:10,])

s1 <- c("A","A","A","B","B","B","C","C","C","D","D","D")
s2 <- c("A","D","A","B","C","B","C","B","C","D","A","D")
s3 <- c("A","B","A","B","A","B","C","D","C","D","C","D")
ex1 <- rbind(s1,s2,s3)
ex1

ex1 <- seqdef(ex1)

seqistatd(ex1)

seqient(ex1)
seqient(ex1,norm=FALSE)

actcal.ient <- seqient(actcal.seq)
summary(actcal.ient)

hist(actcal.ient,col="cyan",
	main="Entropy for the sequences in the actcal data set",
	xlab="Entropy")

hist(seqient(actcal.seq),col="cyan",
	main="Entropy for the sequences in the actcal data set",
	xlab="Entropy")

max(actcal.ient)
which(actcal.ient==max(actcal.ient))
actcal.seq[1836,]

actcal[actcal.ient==max(actcal.ient),]

## Entropy for the biofam data set
biofam.ient <- seqient(biofam.seq)

hist(biofam.ient,col="cyan",
	xlab="Entropy",
	main="Entropy for the sequences in the biofam data set")

biofam <- data.frame(biofam, seqient(biofam.seq))
names(biofam)
summary(biofam$Entropy)

q1 <- quantile(biofam$Entropy,0.01)
q49 <- quantile(biofam$Entropy,0.49)
q51 <- quantile(biofam$Entropy,0.51)
q99 <- quantile(biofam$Entropy,0.99)

ient.min <- biofam.seq[biofam$Entropy<=q1,]
ient.med <- biofam.seq[biofam$Entropy>=q49 & biofam$Entropy<=q51,]
ient.max <- biofam.seq[biofam$Entropy>=q99,]

par(mfrow=c(2,2))
seqiplot(ient.min,
         title="10 seq. with low entropy",
         withlegend=FALSE)
seqiplot(ient.med,
         title="10 seq. with medium entropy",
         withlegend=FALSE)
seqiplot(ient.max,
         title="10 seq. with high entropy",
         withlegend=FALSE)
seqlegend(biofam.seq)

table(biofam$birthyr)

biofam <- data.frame(biofam,
	ageg=cut(biofam$birthy,c(1909,1918,1928,1938,1948,1958),
	label=c("1909-18","1919-28","1929-38","1939-48","1949-58"),include.lowest=TRUE))
table(biofam$ageg)

boxplot(Entropy ~ ageg,
        data=biofam,
        xlab="Birth cohort",
        ylab="Sequences entropy",
        col="cyan")

## -------------------
## Sequence turbulence
## -------------------
data(actcal)
actcal.seq <- seqdef(actcal,13:24)
actcal.seq[2,]
seqdss(actcal.seq[2,])

seqsubsn(actcal.seq[2,],DSS=FALSE)
seqsubsn(actcal.seq[2,],DSS=TRUE)

sp.ex1
sp.ex1 <- seqdef(sp.ex1,informat="SPS")
seqST(sp.ex1)

## biofam data set
biofam <- data.frame(biofam, seqST(biofam.seq))

summary(biofam$Turbulence)

hist(biofam$Turbulence,
     col="cyan",
     xlab="Turbulence",
     main="Turbulences for the sequences in the biofam data set")

max.turb <- max(biofam$Turbulence)
subset(biofam, Turbulence==max.turb)

biofam[biofam$Turbulence==max.turb,]

max.seq <- which(biofam$Turbulence==max.turb)
print(biofam.seq[max.seq,],format='SPS')

## Correlation with entropy
cor(biofam$Turbulence,biofam$Entropy)
cor(biofam$Turbulence,biofam$Entropy, method='spearman')

plot(biofam$Turbulence,biofam$Entropy,
     main="Turbulence vs. Entropy",
     xlab="Turbulence",
     ylab="Entropy")

## Low, medium and high turbulence
q1  <- quantile(biofam$Turbulence,0.01)
q49 <- quantile(biofam$Turbulence,0.49)
q51 <- quantile(biofam$Turbulence,0.51)
q99 <- quantile(biofam$Turbulence,0.99)

turb.min <- biofam.seq[biofam$Turbulence<=q1,]
turb.med <- biofam.seq[biofam$Turbulence>=q49 & biofam$Turbulence<=q51,]
turb.max <- biofam.seq[biofam$Turbulence>=q99,]

par(mfrow=c(2,2))
seqiplot(turb.min,
         title="10 seq. with low turbulence",
         withlegend=FALSE)
seqiplot(turb.med,
         title="10 seq. with medium turbulence",
         withlegend=FALSE)
seqiplot(turb.max,
         title="10 seq. with high turbulence",
         withlegend=FALSE)
seqlegend(biofam.seq) 



