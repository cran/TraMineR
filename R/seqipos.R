## proportion of positive spells (or states)

seqipos <- function(seqdata, dss=TRUE, pos.states=NULL, neg.states=NULL, with.missing=FALSE){

	if (!inherits(seqdata,"stslist"))
		stop(" [!] data is NOT a sequence object, see seqdef function to create one")
  if (is.null(pos.states) & is.null(neg.states))
		stop(" [!] at least one of pos.states and neg.states should be non null!")

  alph <- alphabet(seqdata)
  void <- attr(seqdata,"void")
	nr <- attr(seqdata,"nr")
  if (with.missing)
    alph <- c(alph,nr)

  if (!is.null(pos.states) & !all(pos.states %in% alph)){
    stop(" [!] invalid values in pos.states: ",paste(pos.states[!pos.states %in% alph], collapse=",
    "))
  }
  if (!is.null(pos.states) & !all(neg.states %in% alph)){
    stop(" [!] invalid values in neg.states: ",paste(neg.states[!neg.states %in% alph], collapse=",
    "))
  }

  if (is.null(pos.states)) pos.states <- alph[!alph %in% neg.states]
  if (is.null(neg.states)) neg.states <- alph[!alph %in% pos.states]

  if (length(pos.states)!=length(unique(pos.states)))
    stop(" [!] Multiple occurrences of same state in pos.states")
  if (length(neg.states)!=length(unique(neg.states)))
    stop(" [!] Multiple occurrences of same state in neg.states")

  recodes <- list("p"=pos.states, "m"=neg.states)

  if (dss)
    s <- seqdss(seqdata, with.missing = with.missing)
  else
    s <- seqdata

  sbinary <- seqrecode(s, recodes = recodes, otherwise=attr(s,'void'))

  npos <- rowSums(sbinary=="p")
  nneg <- rowSums(sbinary=="m")
  ratio <- npos/(nneg + npos)
  return(ratio)
}
