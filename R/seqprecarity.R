seqprecarity <- function(seqdata, correction=NULL,
    otto=.2, a=1, b=1.2, stprec=NULL, method = "TRATEDSS",
    state.order=alphabet(seqdata), state.equiv=NULL, with.missing=FALSE, ...){

	if (!inherits(seqdata,"stslist"))
		stop(call.=FALSE, "seqprecarity: data is not a state sequence object, use seqdef function to create one")
  if (!is.null(stprec) && length(stprec) != length(alphabet(seqdata)))
    stop(call.=FALSE, "seqprecarity: length(stprec) should equal length(alphabet)")

  if (!with.missing && any(seqdata[,1] == attr(seqdata,"nr")))
    message(" [>] At least one sequence starts with a missing value! \n       Set with.missing=TRUE or use left='DEL' in seqdef.")


  if(is.null(stprec)){
    stprec <- seqprecstart(seqdata, state.order=state.order, state.equiv=state.equiv, with.missing=with.missing)
  } else {## normalize by maximum value and assign class mean value to members of equiv class
    stprec <- seqprecstart(seqdata, state.order=state.order, state.equiv=state.equiv,
                        stprec=stprec, with.missing=with.missing)
  }

  if (is.null(correction)){
    correction <- 1 + seqprecorr(seqdata, method=method, state.order=state.order,
                  state.equiv=state.equiv, stprec=stprec, with.missing=with.missing, ...)
  }
##  index of complexity
  ici <- suppressMessages(seqici(seqdata, with.missing=with.missing))
  alphabet <- alphabet(seqdata)
  if (with.missing) alphabet <- c(alphabet, attr(seqdata,"nr"))
  lalph <- sapply(seqdata[,1],'match',alphabet)


  prec <- otto*stprec[lalph] + (1-otto) * ici^a * correction^b

	class(prec) <- c("seqprec","matrix")

  ##attr(prec,'correction') <- correction
  attr(prec,'stprec') <- stprec
  colnames(prec) <- "prec"

  return(prec)
}

print.seqprec <- function(x, ...){
  names <- dimnames(x)
  attributes(x) <- NULL
  x <- as.matrix(x)
  dimnames(x) <- names
  print(x, ...)
}
