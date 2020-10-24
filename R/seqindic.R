## Table of indicators

seqindic <- function(seqdata, indic=c("visited","trans","entr","cplx","turb2n"), with.missing=FALSE,
              ipos.args=list(), prec.args=list(), w=1) {

	if (!inherits(seqdata,"stslist"))
		msg.stop("data is NOT a sequence object, see seqdef function to create one")

#  indic.list <- c("lgth","nonm","dlgth","visited","recu","meand","dustd","meand2","dustd2",
#    "trans","transp","entr","volat","cplx","turb","turbn","turb2","turb2n",
#    "all","vpos","ppos","nvolat","inpos","prec","integr","visit","basic","diversity","complexity","binary","Develop")

  if ("Develop" %in% indic){
    basic.list <- c("lgth","nonm","dlgth","visited","recu","trans","transp","meand")
    diversity.list <- c("meand","dustd","meand2","dustd2", "entr","volat")
    complexity.list <- c("nsubs","cplx","turb","turbn","turb2","turb2n")
    binary.list <- c("ppos","nvolat","vpos","inpos")
    ranked.list <- c("bad","degrad","prec","prec2")
    group.list <- c("all","basic","diversity","complexity","binary","ranked")
    indic.list <- c(basic.list,diversity.list,complexity.list,binary.list,ranked.list)
    indic <- indic[indic != 'Develop']
  } else {
    basic.list <- c("lgth","nonm","dlgth","visited","recu","trans","transp","meand")
    diversity.list <- c("meand","dustd","entr","volat")
    complexity.list <- c("nsubs","cplx","turb","turbn")
    binary.list <- c("ppos","nvolat","vpos","inpos")
    ranked.list <- c("prec")
    group.list <- c("all","basic","diversity","complexity","binary","ranked")
    indic.list <- c(basic.list,diversity.list,complexity.list,binary.list,ranked.list)
  }

  if ("visit" %in% indic) indic[indic=="visit"] <- "visited"
  if ("integr" %in% indic) indic[indic=="integr"] <- "inpos"

  if (!all(indic %in% c(indic.list, group.list))){
    msg.stop("invalid values in indic: ", paste(indic[!indic %in% c(indic.list, group.list)], collapse=", "))
  }

  if (any(indic=="all")) {
    indic <- indic[indic!="all"]
    indic <- unique(c(indic, basic.list, diversity.list, complexity.list))
  }
  if (any(indic=="basic")){
    indic <- indic[indic!="basic"]
    indic <- unique(c(indic, basic.list))
  }
  if (any(indic=="diversity")){
    indic <- indic[indic!="diversity"]
    indic <- unique(c(indic, diversity.list))
  }
  if (any(indic=="complexity")){
    indic <- indic[indic!="complexity"]
    indic <- unique(c(indic, complexity.list))
  }
  if (any(indic=="binary")){
    indic <- indic[indic!="binary"]
    indic <- unique(c(indic, binary.list))
  }
  if (any(indic=="ranked")){
    indic <- indic[indic!="ranked"]
    indic <- unique(c(indic, ranked.list))
  }


  if (any(binary.list %in% indic)){
    if(length(ipos.args)==0)
      msg.stop("At least one of the selected indicators requires an non empty ipos.args")
    if(!is.null(ipos.args[["seqdata"]])) msg.warn("seqdata argument in ipos.args will be overwritten!" )
    ipos.args[["seqdata"]] <- seqdata
    if(!is.null(ipos.args[["with.missing"]])) msg.warn("with.missing argument in ipos.args will be overwritten!" )
    ipos.args[["with.missing"]] <- with.missing
    if(!is.null(ipos.args[["index"]])) msg.warn("index argument in ipos.args will be overwritten!" )
  }

  if (any(ranked.list %in% indic)){
    if(length(prec.args)==0)
      msg.warn("No prec.args list. Alphabet order used as state order.")
    if(!is.null(prec.args[["seqdata"]])) msg.warn("seqdata argument in prec.args will be overwritten!" )
    prec.args[["seqdata"]] <- seqdata
    if(!is.null(prec.args[["with.missing"]])) msg.warn("with.missing argument in prec.args will be overwritten!" )
    prec.args[["with.missing"]] <- with.missing
    if(!is.null(ipos.args[["type"]])) msg.warn("type argument in prec.args will be overwritten!" )
  }


  #if ("integr" %in% indic && length(integr.args)==0)
  #  msg.stop("'integr' requires a non empty integr.args!")

  #if ("prec" %in% indic && length(prec.args)==0)
  #  msg.warn("'prec' requested with empty prec.args!")


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
  if ("dlgth" %in% indic || "visited" %in% indic || "recu" %in% indic){
	  dlgth <- suppressMessages(
      seqlength(seqdss(seqdata, with.missing=with.missing), with.missing=with.missing))
    if("dlgth" %in% indic){
    ## Length of dss
      tab <- cbind(tab,dlgth)
      lab <- c(lab,"Dlgth")
    }
    sdist <- suppressMessages(
      seqistatd(seqdata, with.missing=with.missing))
    nvisit <- rowSums(sdist>0)
    if("visited" %in% indic){
    ## Number of visited states
      tab <- cbind(tab,nvisit)
      lab <- c(lab,"Visited")
    }
    recu <- dlgth/nvisit
    if("recu" %in% indic){
    ## Number of visited states
      tab <- cbind(tab,recu)
      lab <- c(lab,"Recu")
    }
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
  if(any(c("meand","dustd") %in% indic)){
    vardur <- suppressMessages(
		  seqivardur(seqdata, with.missing=with.missing, type=1))
    if("meand" %in% indic){
  	## Longitudinal Entropy
  	  meand <- attr(vardur,'meand')
      tab <- cbind(tab,meand)
      lab <- c(lab,"MeanD")
    }
    if("dustd" %in% indic){
  	## Longitudinal Entropy
  	  dustd <- sqrt(vardur)
      tab <- cbind(tab,dustd)
      lab <- c(lab,"Dustd")
    }
  }
  if(any(c("meand2","dustd2") %in% indic)){
    vardur <- suppressMessages(
		  seqivardur(seqdata, with.missing=with.missing, type=2))
    if("meand2" %in% indic){
  	## Longitudinal Entropy
  	  meand <- attr(vardur,'meand')
      tab <- cbind(tab,meand)
      lab <- c(lab,"MeanD2")
    }
    if("dustd2" %in% indic){
  	## Longitudinal Entropy
  	  dustd <- sqrt(vardur)
      tab <- cbind(tab,dustd)
      lab <- c(lab,"Dustd2")
    }
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
  if("nsubs" %in% indic){
  ## number of subsequences of the DSS
    nsubs <- seqsubsn(seqdata, with.missing=with.missing)
    tab <- cbind(tab,nsubs)
    lab <- c(lab,"Nsubs")
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
	  turb <- suppressMessages(seqST(seqdata, norm=TRUE, with.missing=with.missing, silent=TRUE, type=1))
    tab <- cbind(tab,turb)
    lab <- c(lab,"Turbn")
  }
  if("turb" %in% indic){
	## Turbulence
	  turb <- suppressMessages(seqST(seqdata, norm=FALSE, with.missing=with.missing, silent=TRUE, type=1))
    tab <- cbind(tab,turb)
    lab <- c(lab,"Turb")
  }
  if("turb2n" %in% indic){
	## Normalized Turbulence
	  turb <- suppressMessages(seqST(seqdata, norm=TRUE, with.missing=with.missing, silent=TRUE, type=2))
    tab <- cbind(tab,turb)
    lab <- c(lab,"Turb2n")
  }
  if("turb2" %in% indic){
	## Turbulence
	  turb <- suppressMessages(seqST(seqdata, norm=FALSE, with.missing=with.missing, silent=TRUE, type=2))
    tab <- cbind(tab,turb)
    lab <- c(lab,"Turb2")
  }

###   if("integr" %in% indic){
###   ## index of integration
###     if(!is.null(integr.args[["state"]])){
###       if(!is.null(integr.args[["seqdata"]]))
###         warning( "[!] seqdata argument in integr.args will be overwritten!" )
###       integr.args[["seqdata"]] <- seqdata
###       if(!is.null(integr.args[["with.missing"]]))
###         warning( "[!] with.missing argument in integr.args will be overwritten!" )
###       integr.args[["with.missing"]] <- with.missing
###       integr <- do.call(seqintegration, args=integr.args)
###       tab <- cbind(tab,integr)
###       lab <- c(lab,"Integr")
###     } else
###         warning( "[!] Integr not computed because no state specified in integr.args!" )
###
###   }

  ipos.dss <- NULL
  if(!is.null(ipos.args[["dss"]])) ipos.dss <- ipos.args[["dss"]]
  if("ppos" %in% indic){
  ## Proportion of positive states
    ipos.args[["index"]] <- "share"
    ipos.args[["dss"]] <- FALSE
    ipos <- do.call(seqipos, args=ipos.args)
    ipos.args[["dss"]] <- ipos.dss
    tab <- cbind(tab,ipos)
    lab <- c(lab,"Ppos")
  }
  if("nvolat" %in% indic){
  ## normative volatility
    ipos.args[["index"]] <- "share"
    ipos.args[["dss"]] <- TRUE
    ipos <- do.call(seqipos, args=ipos.args)
    ipos.args[["dss"]] <- ipos.dss
    tab <- cbind(tab,ipos)
    lab <- c(lab,"Nvolat")
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

  if("bad" %in% indic){
  ## index of precarity
    dlist <- unique(c(names(formals(seqibad)),names(formals(seqprecstart))))
    bad <- do.call(seqibad, args=prec.args[names(prec.args) %in% dlist])
    tab <- cbind(tab,bad)
    lab <- c(lab,"Bad")
  }

  if("degrad" %in% indic){
  ## index of precarity
    dlist <- names(formals(seqidegrad))
    degrad <- do.call(seqidegrad, args=prec.args[names(prec.args) %in% dlist])
    tab <- cbind(tab,degrad)
    lab <- c(lab,"Degrad")
  }

  if("prec" %in% indic){
  ## index of precarity
    dlist <- unique(c(names(formals(seqprecarity.private)),names(formals(seqdegrad.private))))
    prec.args[["type"]] <- 1
    prec <- do.call(seqprecarity.private, args=prec.args[names(prec.args) %in% dlist])
    tab <- cbind(tab,prec)
    lab <- c(lab,"Prec")
  }

  if("prec2" %in% indic){
  ## index of precarity
    dlist <- unique(c(names(formals(seqprecarity.private)),names(formals(seqdegrad.private))))
    prec.args[["type"]] <- 2
    prec <- do.call(seqprecarity.private, args=prec.args[names(prec.args) %in% dlist])
    tab <- cbind(tab,prec)
    lab <- c(lab,"Prec2")
  }


  names(tab) <- lab
  tab <- tab[,-1, drop=FALSE]
	return(tab)
}
