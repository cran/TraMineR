## ============================
## Examples for the chapter
## 'Measuring similarities and
## distances between sequences'
## in TraMineR User's Guide
## ============================

## ----
## MPOS
## ----
data(famform)
famform.seq <- seqdef(famform)
famform.seq

seqmpos(famform.seq[1,],famform.seq[2,])
seqmpos(famform.seq[2,],famform.seq[4,])

## ---
## LCP
## ---
data(famform)
famform.seq <- seqdef(famform)
famform.seq

seqLCP(famform.seq[1,],famform.seq[2,])
seqLCP(famform.seq[3,],famform.seq[4,])
seqLCP(famform.seq[3,],famform.seq[5,])

seqdist(famform.seq,method="LCP")

seqdist(famform.seq,method="LCP",norm=TRUE)

1-seqdist(famform.seq,method="LCP",norm=TRUE)

## ---
## LCS
## ---
data(biofam)
biofam.seq <- seqdef(biofam,10:25)
biofam.lcs <- seqdist(biofam.seq,method="LCS")

## --
## OM
## --
couts <- seqsubm(biofam.seq,method="TRATE")
couts

biofam.om <- seqdist(biofam.seq, method="OM", indel=3, sm=couts)

object.size(biofam.om)/1024^2

round(biofam.om[1:10,1:10],1)

## ----------
## LCS <=> OM
## ----------
ccouts <- seqsubm(biofam.seq,method="CONSTANT",cval=2)
ccouts

biofam.om2 <- seqdist(biofam.seq,method="OM",indel=1,sm=ccouts)
biofam.om2[1:10,1:10]

all.equal(biofam.om2,biofam.lcs)


