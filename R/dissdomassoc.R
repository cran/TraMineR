## From Joint Sequence Analysis, (R. Piccarreta, SMR, 2017)
##  measures of association between dimensions
##  portions of code inspired from the assoc.domains function of seqhandbook
##  Main diff with assoc.domains is that here the functions supports weights
##  It also is faster.

dissdomassoc <- function(domdiss, jointdiss = NULL, what = c("pearson","R2"),
      dnames=names(domdiss), weights=NULL, w.rank=FALSE) {

## domdiss: list of dissimilarities matrices or objects (one per channel)
## jointdiss: dissimilarity matrix or object between sequences of combined states
## what: requested measure of association

  if (!is.list(domdiss) || length(domdiss) < 2)
    stop("domdiss should be a list of at least two distance matrices or objects")

  if (!all(is.numeric(unlist(domdiss, use.names=FALSE))))
    stop("non numeric values in distance matrices")

  #whatlist <- c("pearson","spearman","kendall","R2","cronbach","cron.subsets","all")
  ## kendall takes too much time for large number of diss values
  whatlist <- c("pearson","spearman","R2","cronbach","cron.subsets","all")
  if (!all(what %in% whatlist))
    stop("bad what values, allowed are ", paste(whatlist, collapse=","))
  if ("all" %in% what) {
    what <- unique(c(what,whatlist[c(1,2,3,5)]))
    what <- what[what!="all"]
  }
  if ("R2" %in% what & !any(c("pearson","spearman") %in% what) ) {
    what <- c("pearson",what)
  }

  ## setting weights
  if (inherits(domdiss[[1]],'dist'))
    ncases <- attr(domdiss[[1]],'Size')
  else
    ncases <- nrow(domdiss[[1]])
  if (is.null(weights)) {
    weights <- rep(1, ncases)
    weighted=FALSE
    }
  else {
    if (length(weights) != ncases)
      stop("length of weights not equal to number of cases!")
    if (all(weights==1))
      weighted <- FALSE
    else weighted <- TRUE
  }
  ww <- as.numeric(as.dist(weights %*% t(weights)))
  #sww <- sqrt(ww)

  ndom <- length(domdiss)
  ndomv <- ndom - 1
  ## transforming into vector of distances
  distlist <- lapply(domdiss, function(x) as.numeric(as.dist(x)))
  dissmat <- matrix(unlist(distlist, use.names=FALSE), ncol=ndom, byrow=FALSE)
  if (is.null(dnames)) dnames <- paste0('dom.',1:ndom)
  colnames(dissmat) <- dnames

  if (!is.null(jointdiss)){
    if (!all(is.numeric(jointdiss)))
      stop("when not NULL, jointdiss must be a distance matrix or object")
    jointdiss <- as.numeric(as.dist(jointdiss))
    dissmat <- cbind(dissmat,jointdiss)
    ndomv <- c(rep(ndom-1,ndom),ndom)
  }
  dissmat <- apply(dissmat,2,scale)

  if ("spearman" %in% what) { ## we replace columns with weighted ranked
    rankmat <- apply(dissmat,2,rank)
    if (weighted & w.rank) {
      cat("\nPlease wait, computing weighted ranks may take a while ...")
      #dissmat.spear <- apply(dissmat,2,w.rank,w=ww)
      ## above proper solution is much slower than weighted.rank from cNORM
      dissmat.spear <- apply(dissmat,2,weighted.rank,weights=ww)
      ## weighted.rank returns NA for min and max ranks
      ## we replace these NAs with the non-weighted ranks
      dissmat.spear[is.na(dissmat.spear)] <- rankmat[is.na(dissmat.spear)]
    } else {
      dissmat.spear <- rankmat
    }
    #rm(rankmat)
    dissmat.spear <- apply(dissmat.spear,2,scale)
    ##dissmat.spearw <- sww * dissmat.spear
  }
  ##dissmatw <- sww * dissmat

  ## defining function
  rsquare.corr <- function(correlation, jointdiss, ndom, ndomv) {
    corr.tmp <- correlation
    if (!is.null(jointdiss)) {
      corr.tmp[1:ndom,ndom+1] <- 0
    }
    return((rowSums(corr.tmp^2)-1)/ndomv)
  }

  res <- list()

  if ("pearson" %in% what){ ## wtd.cor from weights package
    ##correlation <- cor(dissmatw, method='pearson')
    wcorr <- wtd.cor(dissmat, weight=ww)
    correlation <- wcorr[["correlation"]]
    res[["Pearson"]] <- correlation
    res[["p.pearson"]] <- wcorr[["p.value"]]
    if ("R2" %in% what) res[["Pearson.Rsquare"]] <- rsquare.corr(correlation, jointdiss, ndom, ndomv)
  }
  if ("spearman" %in% what){
    ##correlation <- cor(dissmat.spearw)
    wcorr <- wtd.cor(dissmat.spear, weight=ww)
    correlation <- wcorr[["correlation"]]
    res[["Spearman"]] <- correlation
    res[["p.spearman"]] <- wcorr[["p.value"]]
    if ("R2" %in% what) res[["Spearman.Rsquare"]] <- rsquare.corr(correlation, jointdiss, ndom, ndomv)
  }
  #if ("kendall" %in% what){  ## kendall takes too much time
  #  correlation <- cor(dissmat, method='kendall')
  #  res[["Kendall"]] <- correlation
  #  if ("R2" %in% what) res[["Kendall.Rsquare"]] <- rsquare.corr(correlation, jointdiss, ndom, ndomv)
  #}
  if (any(c("cronbach","cron.subsets") %in% what)){
    ## Cronbach alpha for all domains
    sigmatot <- var(rowSums(dissmat[,1:ndom]))
    chron <- (ndom/(ndom-1))*(1-ndom/sigmatot)
    names(chron) <- paste0('(',paste0(dnames,collapse=','),')')
    res[["Cronbach"]] <- chron
  }

  ## Cronbach alpha for each combination of domains
  if ("cron.subsets" %in% what) {
    cron.subset <- numeric()
    if (ndom>2) {
      for(p in (ndom-1):2) {
        set.start <- length(cron.subset) + 1
        sets <- combn(1:ndom, p, simplify=FALSE)
        set.end <- length(cron.subset) + length(sets)
        for(i in 1:length(sets)) {
          sigmatot <- var(rowSums(dissmat[,sets[[i]],drop=FALSE]))
          cronbach <- (p/(p-1))*(1-p/sigmatot)
          cron.subset <- c(cron.subset,cronbach)
        }
        names(cron.subset)[set.start:set.end] <- lapply(sets, function(x) paste0('(',paste0(dnames[x],collapse=','),')'))
      }
      res[["Cronbach.subsets"]] <- cron.subset
    }
    else
      message("Two or less domains, no subset possible for Cronbach alpha")
  }
  res[["dnames"]] <- dnames

  class(res) <- c(class(res), "ddomassoc")

  return(res)
}


summary.ddomassoc <- function(object, ...){
  dnames <- object[["dnames"]]
  ndom <- length(dnames)
  cnam <- NULL
  rnam <- NULL
  is.pearson <- !is.null(object[["Pearson"]])
  if (is.pearson) cnam <- c(cnam, "Pearson","p.Pearson")
  is.spearman <- !is.null(object[["Spearman"]])
  if (is.spearman) cnam <- c(cnam, "Spearman","p.Spearman")
  ncol <- length(cnam)
  tab <- matrix(NA,nrow=ndom*(ndom-1)/2, ncol=ncol)
  colnames(tab) <- cnam

  tabnames <- NULL
  k <- 0 ## res row counter
  dn.split <- strsplit(dnames," x ")
  for (d1 in 1:(ndom-1)) {
    for (d2 in (d1+1):ndom ) {
       ## cat("\n d1 = ",d1, " d2 = ", d2)
      if (!any(dn.split[[d1]] %in% dn.split[[d2]])) {
        k <- k+1
        if(is.pearson) {
          tab[k,"Pearson"]    <- object[["Pearson"]][dnames[d1],dnames[d2]]
          tab[k,"p.Pearson"] <- object[["p.pearson"]][dnames[d1],dnames[d2]]
        }
        if(is.spearman) {
          tab[k,"Spearman"] <- object[["Spearman"]][dnames[d1],dnames[d2]]
          tab[k,"p.Spearman"] <- object[["p.spearman"]][dnames[d1],dnames[d2]]
        }
        tabname <- paste(dnames[d1],dnames[d2],sep="_vs_")
        tabnames <- c(tabnames,tabname)
      }
    }
  }
  tab <- tab[1:k,]
  rownames(tab) <- tabnames
  return(tab)
}


##################

### w.rank <- function(x, w) {
###   ox <- order(x)
###   ## original order
###   oox <- rank(x,ties.method='first')
###   ooxt <- rank(x) ## with ties
###   ## sorted weights
###   sw <- w[ox]
###   ## cumsum normalized to range from 1 to length(x)
###   csw <- cumsum(sw)
###   rcsw <- (length(w)-1)*(csw-1)/(csw[length(csw)]-1) + 1
###
###   ## replace weighted rank of ties by mean weighted rank of the ties
###   for (i in unique(ooxt)){
###     rcsw[oox[ooxt==i]]<-mean(rcsw[oox[ooxt==i]])
###   }
###
###   ## setting original order
###   return(rcsw[oox])
### }
