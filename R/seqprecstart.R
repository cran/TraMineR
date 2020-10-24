seqprecstart <- function(seqdata, state.order=alphabet(seqdata, with.missing), state.equiv=NULL, stprec=NULL, with.missing=FALSE) {
  ## cost of starting state

  if (is.null(state.order) && is.null(stprec))
    stop("state.order and stprec cannot both be null!")

  if (!is.null(stprec)) { ## state.order set from stprec
    if (!is.null(state.order)) {
       msg.warn("state.order overridden by stprec order!")
    }
    ord <- order(stprec)
    ord <- ord[(1+sum(stprec<0)):length(ord)] ## only positive values
    state.order <- alphabet(seqdata,with.missing)[ord]
    stprec[stprec<0] <- mean(stprec[ord])
  }

  if (!is.null(state.order))
    state.order <- states.check(seqdata, state.order, state.equiv, with.missing=with.missing)

  step <- 1/(length(state.order)-1)

  alphabet <- alphabet(seqdata, with.missing=with.missing)

  state.noncomp <- NULL
  if (length(unique(state.order)) < length(alphabet)){
    inoncomp <- which(is.na(match(alphabet,unique(state.order))))
    state.noncomp <- alphabet[inoncomp]
    state.order.plus <- c(state.order, state.noncomp)
  }
  else {
    state.order.plus <- state.order
  }


  ord <- match(state.order.plus,alphabet)
  ordo <- match(alphabet,state.order.plus)

  if(is.null(stprec)){
    stprec <- seq(from=0, to=1, by=step)
    ## assign mean cost to non ranked states
    stprec <- c(stprec, rep(mean(stprec),length(state.noncomp)))
    stprec <- stprec[ordo]
  }
  else ## user provided stprec should conform order of the alphabet
    stprec <- (stprec-min(stprec))/(max(stprec)-min(stprec))

  ## assign the class mean cost to all states of a same equivalent class
  if(!is.null(state.equiv)){
    lc <- length(state.equiv)
    for (i in 1:lc) {
      iequiv <- match(state.equiv[[i]],alphabet)
      stprec[iequiv] <- mean(stprec[iequiv])
    }
  }

  return(stprec)
}
