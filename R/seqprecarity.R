## aliases for simplified precarity

seqprecarity <- function(seqdata, correction=NULL,
    otto=.2, a=1, b=1.2, stprec=NULL,
    method = "TRATEDSS",
    state.order=alphabet(seqdata, with.missing), state.equiv=NULL, with.missing=FALSE,
    ...){

  prec <- seqprecarity.private(seqdata, type=1, correction=correction,
    otto=otto, a=a, b=b, stprec=stprec,
    method = method,
    state.order=state.order, state.equiv=state.equiv, with.missing=with.missing,
    ##degr = FALSE, start.integr=FALSE, spell.integr=FALSE,
    ##norm.trpen=FALSE,
    ...)

  return(prec)
}


seqprecarity.private <- function(seqdata, type=1, correction=NULL,
    otto=.2, a=switch(type,1,.5), b=switch(type,1.2,1-a), stprec=NULL,
    method = switch(type,"TRATEDSS","RANK"),
    state.order=alphabet(seqdata, with.missing), state.equiv=NULL, with.missing=FALSE,
    ##degr = FALSE, start.integr=FALSE, spell.integr=FALSE,
    ##norm.trpen=FALSE,
    ...){

  #if (is.null(otto)) otto <- switch(type,.2,.5)
  #if (is.null(a)) a <- switch(type,1,.5)
  #if (is.null(method)) method <- switch(type,"TRATEDSS","RANK")

  if (!(type %in% c(1,2)))
		stop(call.=FALSE, "seqprecarity: type should be 1 or 2")
  if (type==2){
    start.integr=TRUE
    spell.integr=TRUE
    tr.sum=TRUE
    #degr=TRUE
    if (a < 0 || a > 1) stop(call.=FALSE, "with type=2, a must be in [0,1]")
  } else {
    start.integr=FALSE
    spell.integr=FALSE
    tr.sum = FALSE
    #degr=FALSE
  }

	if (!inherits(seqdata,"stslist"))
		stop(call.=FALSE, "seqprecarity: data is not a state sequence object, use seqdef function to create one")
  if (!is.null(stprec) && length(stprec) != length(alphabet(seqdata, with.missing=with.missing)))
    stop(call.=FALSE, "seqprecarity: length(stprec) should equal length of alphabet")

  ##if (!with.missing && any(seqdata[,1] == attr(seqdata,"nr")))
  ##  message(" [>] At least one sequence starts with a missing value!
  ##         Set with.missing=TRUE or use left='DEL' in seqdef.")

  ##if(is.null(stprec)){
  ##  stprec <- seqprecstart(seqdata, state.order=state.order, state.equiv=state.equiv, with.missing=with.missing)
  ##} else {## normalize by maximum value and assign class mean value to members of equiv class
  stprec <- seqprecstart(seqdata, state.order=state.order, state.equiv=state.equiv,
                        stprec=stprec, with.missing=with.missing)
  ##}

  if (is.null(correction)){
    correction <- 1 + seqdegrad.private(seqdata, method=method, state.order=state.order,
                  state.equiv=state.equiv, stprec=stprec, with.missing=with.missing,
                  tr.sum=tr.sum,
                  ...)
    #if (degr) {
    #  correction <- correction #/2
      #b=1
      #a=0
    #}
  }

##  index of complexity
  if (a != 0){
    ici <- suppressMessages(seqici(seqdata, with.missing=with.missing))
    #if (type==2) ici <- 1+ici
  }
  else
    ici <- 1
  alphabet <- alphabet(seqdata, with.missing=with.missing)
  sdss <- seqdss(seqdata,with.missing=with.missing)
  lalph <- sapply(sdss[,1],'match',alphabet)
  #nr1 <- which(seqdata[,1]==attr(seqdata,'nr'))
  #if (with.missing && length(nr1)>0 && ncol(seqdata)>1) {
  #  lalph2 <- sapply(seqdss(seqdata,with.missing=with.missing)[,2],'match',alphabet)
  #  lalph[nr1] <- lalph2[nr1]
  #}
  integr1 <- rep(1,length(lalph))

  if (start.integr){
    Dur <- seqdur(seqdata, with.missing=with.missing)
    make.sps <- function(dur){
      sps <- paste0(1:length(dur),'/',dur)
      return(sps)
    }
    sps <- t(apply(Dur,1,make.sps))
    sps[is.na(Dur)] <- NA
    seqtmp <- suppressMessages(seqdef(sps, informat='SPS', SPS.in=list(xfix='',sdsep='/')))
    integr1 <- seqintegration(seqtmp, state=1, pow=0)
    ## for sequence starting with missing we consider the second spell
    #if (length(nr1>0)){
    #  integr2 <- seqintegration(seqtmp, state=2, pow=0)
    #  integr1[nr1] <- integr2[nr1]
    #}

  }


#  if (degr)
#    prec <- integr1 * stprec[lalph] + correction
#  else
  if (type==1){
    prec <- otto*(stprec[lalph]*integr1) + (1-otto) * ici^a * correction^b
  }
  else {
    minstprec <- function(dssrow,stprec,alphabet){
      if (any(dssrow %in% alphabet))
        best <- min(stprec[alphabet %in% dssrow])
      else
        best <- NA
      return(best)
    }

    prec <- stprec[lalph]*integr1 + (ici + (correction-1))
    #prec <- cbind(prec, rep(0,nrow(prec)))
    prec <- cbind(prec,apply(as.matrix(sdss),1,minstprec,stprec=stprec,alphabet=alphabet))
    prec[,1] <- apply(prec,1,max)
    prec[,2] <- rep(1,nrow(prec))
    prec[,1] <- apply(prec,1,min)
    prec <- prec[,1,drop=FALSE]
  }

	class(prec) <- c("seqprec","matrix")

  ##attr(prec,'correction') <- correction
  attr(prec,'stprec') <- stprec
  colnames(prec) <- switch(type,"prec","prec2")

  return(prec)
}

print.seqprec <- function(x, ...){
  names <- dimnames(x)
  attributes(x) <- NULL
  x <- as.matrix(x)
  dimnames(x) <- names
  print(x, ...)
}
