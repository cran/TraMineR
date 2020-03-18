seqivolatility <- function(seqdata, type=1, w=.5, with.missing=FALSE){

	if (!inherits(seqdata,"stslist"))
		stop(" [!] data is NOT a sequence object, see seqdef function to create one")
  if (!type %in% c(1,2))
		stop(" [!] type should be 1  or 2!")
  if (w>1 | w<0)
		stop(" [!] w should be in the range [0, 1]!")


  alph <- alphabet(seqdata)
  ##void <- attr(seqdata,"void")
	nr <- attr(seqdata,"nr")
  if (with.missing)
    alph <- c(alph,nr)
  alph.size <- length(alph)


  transp <- suppressMessages(
      seqtransn(seqdata, with.missing=with.missing, norm=TRUE))
  sdist <- suppressMessages(
      seqistatd(seqdata, with.missing=with.missing))
  nvisit <- rowSums(sdist>0)
  if (type == 1) {
    ret <- ifelse(nvisit - 1 <= 0, 0, (nvisit - 1)/(alph.size -1))
    ret <- w * ret + (1-w) * transp
  } else if (type == 2) {
    ret <- w * nvisit/alph.size + (1-w) * transp
  }

  return(as.vector(ret))
}
