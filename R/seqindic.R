## Table of indicators

seqindic <- function(seqdata, indic=c("visited","trans","entr","cplx"), with.missing=FALSE,
              ipos.args=list(), integr.args=list(), prec.args=list(), w=1) {

	if (!inherits(seqdata,"stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

  indic.list <- c("lgth","nonm","dlgth","visited","trans",
    "transp","entr","volat","cplx","turb","turbn",
    "all","ppos","vpos","inpos","prec","integr","visit")

  indic.ipos <- c("ppos","vpos","inpos")
  if (any(indic.ipos %in% indic)){
    if(!is.null(ipos.args[["seqdata"]])) warning( "[!] seqdata argument in ipos.args will be overwritten!" )
    ipos.args[["seqdata"]] <- seqdata
    if(!is.null(ipos.args[["with.missing"]])) warning( "[!] with.missing argument in ipos.args will be overwritten!" )
    ipos.args[["with.missing"]] <- with.missing
    if(!is.null(ipos.args[["index"]])) warning( "[!] index argument in ipos.args will be overwritten!" )
  }

  if (!all(indic %in% indic.list)){
    stop("invalid values in indic: ",paste(indic[!indic %in% indic.list], collapse=", "))
  }

  ##if ("ipos" %in% indic) indic[indic=="ipos"] <- "ppos"
  if ("visit" %in% indic) indic[indic=="visit"] <- "visited"

  if (any(indic=="all")) {
    indic.all <- indic.list[-((length(indic.list)-6):length(indic.list))]
    if("ppos" %in% indic) indic.all <- c(indic.all,"ppos")
    if("vpos" %in% indic) indic.all <- c(indic.all,"vpos")
    if("inpos" %in% indic) indic.all <- c(indic.all,"inpos")
    if("prec" %in% indic) indic.all <- c(indic.all,"prec")
    if("integr" %in% indic) indic.all <- c(indic.all,"integr")
    indic <- indic.all
  }

  tab <- as.data.frame(rownames(seqdata))
  lab <- "id"

  if("lgth" %in% indic){
  ## Sequence length
    lgth <- suppressMessages(
      seqlength(seqdata, with.missing=TRUE))
    tab <- cbind(tab,lgth)
    lab <- c(lab,"Lgth")
  }
  if("nonm" %in% indic){
  ## Number of non-missing elements
    nonm <- suppressMessages(
      seqlength(seqdata, with.missing=FALSE))
    tab <- cbind(tab,nonm)
    lab <- c(lab,"NonM")
  }
  if("dlgth" %in% indic){
  ## Length of dss
	  dlgth <- suppressMessages(
      seqlength(seqdss(seqdata, with.missing=with.missing), with.missing=with.missing))
    tab <- cbind(tab,dlgth)
    lab <- c(lab,"Dlgth")
  }
  if("visited" %in% indic){
  ## Number of visited states
    sdist <- suppressMessages(
      seqistatd(seqdata, with.missing=with.missing))
    nvisit <- rowSums(sdist>0)
    tab <- cbind(tab,nvisit)
    lab <- c(lab,"Visited")
  }
  if("trans" %in% indic){
	## Number of state changes (transitions)
	  trans <- suppressMessages(
      seqtransn(seqdata, with.missing=with.missing, norm=FALSE))
    tab <- cbind(tab,trans)
    lab <- c(lab,"Trans")
  }
  if("transp" %in% indic){
	## Proportion of state changes
	  trans <- suppressMessages(
      seqtransn(seqdata, with.missing=with.missing, norm=TRUE))
    tab <- cbind(tab,trans)
    lab <- c(lab,"Transp")
  }
  if("entr" %in% indic){
	## Longitudinal Entropy
	  ient <- suppressMessages(
		  seqient(seqdata, with.missing=with.missing, norm=TRUE))
    tab <- cbind(tab,ient)
    lab <- c(lab,"Entr")
  }
  if("volat" %in% indic){
	## Longitudinal Entropy
	  volat <- suppressMessages(
		  seqivolatility(seqdata, with.missing=with.missing, w=w))
    tab <- cbind(tab,volat)
    lab <- c(lab,"Volat")
  }
  if("cplx" %in% indic){
	## Complexity
	  ici <- suppressMessages(
		  seqici(seqdata, with.missing=with.missing, silent=TRUE))
    tab <- cbind(tab,ici)
    lab <- c(lab,"Cplx")
  }
  if("turbn" %in% indic){
	## Normalized Turbulence
	  turb <- suppressMessages(seqST(seqdata, norm=TRUE, with.missing=with.missing, silent=TRUE))
    tab <- cbind(tab,turb)
    lab <- c(lab,"Turbn")
  }
  if("turb" %in% indic){
	## Turbulence
	  turb <- suppressMessages(seqST(seqdata, norm=FALSE, with.missing=with.missing, silent=TRUE))
    tab <- cbind(tab,turb)
    lab <- c(lab,"Turb")
  }

  if("integr" %in% indic){
  ## index of precarity
    if(!is.null(integr.args[["state"]])){
      if(!is.null(integr.args[["seqdata"]]))
        warning( "[!] seqdata argument in integr.args will be overwritten!" )
      integr.args[["seqdata"]] <- seqdata
      if(!is.null(integr.args[["with.missing"]]))
        warning( "[!] with.missing argument in integr.args will be overwritten!" )
      integr.args[["with.missing"]] <- with.missing
      integr <- do.call(seqintegration, args=integr.args)
      tab <- cbind(tab,integr)
      lab <- c(lab,"Integr")
    } else
        warning( "[!] Integr not computed because no state specified in integr.args!" )

  }


  if("ppos" %in% indic){
  ## Proportion of positive states
    ipos.args[["index"]] <- "share"
    ipos <- do.call(seqipos, args=ipos.args)
    tab <- cbind(tab,ipos)
    lab <- c(lab,"Ppos")
  }

  if("vpos" %in% indic){
  ## volatility of pos-neg sequences
    ipos.args[["index"]] <- "volatility"
    ipos <- do.call(seqipos, args=ipos.args)
    tab <- cbind(tab,ipos)
    lab <- c(lab,"Vpos")
  }
  if("inpos" %in% indic){
  ## Potential to integrate pos states
    ipos.args[["index"]] <- "integration"
    ipos <- do.call(seqipos, args=ipos.args)
    tab <- cbind(tab,ipos)
    lab <- c(lab,"Inpos")
  }

  if("prec" %in% indic){
  ## index of precarity
    if(!is.null(prec.args[["seqdata"]]))
      warning( "[!] seqdata argument given in prec.args is overwritten!" )
    prec.args[["seqdata"]] <- seqdata
    if(!is.null(prec.args[["with.missing"]]))
      warning( "[!] with.missing argument given in prec.args is overwritten!" )
    prec.args[["with.missing"]] <- with.missing
    prec <- do.call(seqprecarity, args=prec.args)
    tab <- cbind(tab,prec)
    lab <- c(lab,"Prec")
  }

  names(tab) <- lab
  tab <- tab[,-1, drop=FALSE]
	return(tab)
}
