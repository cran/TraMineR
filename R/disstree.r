#################
##DisTreeNode
#################
# Donnée accessible
# predictor : liste des prédicteurs
# dissmatrix : matrices des dissimilarités internes
#

#Donnée interne
# split : predicteur choisi(NULL pour noeux terminaux)
# vardis : variabilité interne
# children : noeud enfant (NULL pour noeux terminaux)
# ind: liste des index des individus du noeuds.
# depth: profondeur du noeud
# label: label du noeud (valeurs du predicteur
# R2: R2 du split,(NULL pour noeux terminaux)
#
DTNInit<-function(ind, vardis, depth, label){
  node<-list()
  node$label<-label
  node$children<-NULL
  node$split<-NULL
  node$ind<-ind
  node$depth<-depth
  node$vardis<-vardis
  class(node)<-c("DissTreeNode","list")
  return(node)
}
disstreeleaf<-function(tree){
  if(!inherits(tree,"disstree"))stop("tree should be a disstree object")
  categorie<-numeric(length(tree$root$ind))
  categorie[]<--1
  counter<-0
  catobject<-list(categorie=categorie,counter=counter)
  return(factor(DTNdisstreeleaf(tree$root,catobject)$categorie))
}
DTNdisstreeleaf<-function(node,co){
  if(is.null(node$children)){
    co$categorie[node$ind]<-(co$counter+1)
    catobject<-list(categorie=co$categorie,counter=co$counter+1)
    return(catobject)
  }else{
    co1<-DTNdisstreeleaf(node$children$left,co)
    return(DTNdisstreeleaf(node$children$right,co1))
  }
}
disstree<-function(formula,data=NULL,minSize=0.05,maxdepth=5, R=1000,pval=0.01){
  formula.call<-formula
  dissmatrix <- eval(formula[[2]], data, parent.frame()) # to force evaluation
  if(inherits(dissmatrix,"dist"))dissmatrix<-as.matrix(dissmatrix)
  formula[[2]] <- NULL
  #Model matrix from forumla
  predictor <- model.frame(formula, data, drop.unused.levels = TRUE,na.action=NULL)
##  print(summary(predictor))
  pop<-nrow(dissmatrix)
  if(minSize<1)minSize<- round(pop*minSize)
  if(pop!=nrow(predictor))stop("dissimilarity matrix and data should be of the same size")
  vardis<-sum(dissmatrix)/(2*pop*pop)
  tree <- list(formula=formula.call)
  tree$root<-DTNBuildNode(dissmatrix,as.data.frame(predictor),minSize,1:pop,vardis,1,"Root",R,pval, maxdepth)
  class(tree)<-c("disstree", "list")
  tree$adjustement <- dissassoc(dissmatrix,disstreeleaf(tree),R=R) 
  return(tree)
}
DTNBuildNode<-function(dmat,pred,minSize,ind, vardis, depth, label,nbperm,pval, maxdepth){
  node<-DTNInit(ind, vardis, depth, label)
  SCtot<-node$vardis*length(node$ind)
  SCres<-SCtot
#    print(class(pred))
#  print(SCtot)
  bestSpl<-NULL
  varnames<-colnames(pred)
  if(depth>=maxdepth){
      return(node)
  }
  for(p in 1:ncol(pred)){
#    cat("Checking",varnames[[p]],"...\n")
    spl<-DTNGroupFactorBinary(dmat,SCres,pred[,p],minSize,varnames[[p]])
    if(!is.null(spl)&&(is.null(bestSpl)||spl$SCres<bestSpl$SCres)){
      bestSpl<-spl
      SCres<-spl$SCres
#      cat(varnames[[p]]," Ok","\n")
    }else{
#      cat(varnames[[p]]," NULL","\n")
    }
    
  }
  if(is.null(bestSpl)){
    return(node)
  }
  if(nbperm>1){
    sig<-dissassoc(dmat,bestSpl$variable, R=nbperm)
    spval<-sig$stat$PseudoF_Pval
#    print(spval)
    if(spval>pval)return(node)
  }

  node$split<-bestSpl
  node$children<-list()
  node$children$left<-DTNBuildNode(dmat[bestSpl$variable,bestSpl$variable],as.data.frame(pred[bestSpl$variable,]),minSize,ind[bestSpl$variable], vardis=bestSpl$lvar, depth=depth+1, label=paste(bestSpl$llabels, collapse="/"),nbperm,pval,maxdepth)
  node$children$right<-DTNBuildNode(dmat[!bestSpl$variable,!bestSpl$variable],as.data.frame(pred[!bestSpl$variable,]),minSize,ind[!bestSpl$variable], vardis=bestSpl$rvar, depth=depth+1, label=paste(bestSpl$rlabels, collapse="/"),nbperm,pval,maxdepth)
  node$R2<-1-(SCres/SCtot)
#  print("split done")
  return(node)
}
DTNGroupFactorBinary<-function(dissmatrix, currentSCres,pred, minSize, varname){
  totpop<-nrow(dissmatrix)
  ind<-1:totpop
  grp<-factor(pred, ordered=(is.ordered(pred)||is.numeric(pred)))
  lgrp<-levels(grp)
  if(length(lgrp)==1)return(NULL)
  nbGrp<-length(lgrp)
  has.na<-FALSE
  llgrp<-lgrp
  #Here we add a group for missing values
  if(sum(is.na(grp))>0){
    nbGrp<-length(lgrp)+1
    has.na<-TRUE
    llgrp[nbGrp]<-"<Missing>"
  }

  grpCond<-list()
  grpSize<-numeric(length=nbGrp)
  grpSize[]<-0
  for(i in 1:length(lgrp)){
    #on crée le groupe en question
    grpCond[[i]]<-(grp==lgrp[i])
    grpCond[[i]][is.na(grpCond[[i]])]<-FALSE
    grpSize[i]<-sum(grpCond[[i]])
#    print(grpCond[[i]])
#    print(paste(i,grpSize[i]))
  }
#   print("Sizes")
#  print(grpSize)
  #Treating missing values
  if(has.na){
    grpCond[[nbGrp]]<-is.na(grp)
    grpSize[nbGrp]<-sum(grpCond[[nbGrp]])
  }
# print("Sizes2")
#  return(NULL)
#  print(SCres)
  inertiaMat<-matrix(0,nrow=nbGrp, ncol=nbGrp)
#  print(nbGrp)
  for(i in 1:(nbGrp-1)){
    grpindiv1<-ind[grpCond[[i]]]
#    print((i+1):nbGrp)
    for(j in (i+1):nbGrp){
        grpindiv2<-ind[grpCond[[j]]]
        r<-.Call("tmrinterinertia", dissmatrix, as.integer(grpindiv1),as.integer(grpindiv2), PACKAGE="TraMineR")
#        print(paste(i,j))
#        inertiaMat[i,j]<-r
        inertiaMat[j,i]<-r
      }
    r<-.Call("tmrsubmatrixinertia", dissmatrix, as.integer(grpindiv1), PACKAGE="TraMineR")*sum(grpCond[[i]])
#    r2<-.Call("tmrinterinertia", dissmatrix, as.integer(grpindiv1),as.integer(grpindiv1), PACKAGE="TraMineR")
#    print(r-r2)
    inertiaMat[i,i]<-r
  }
#  print(inertiaMat)
  #FIXME This step is missing in the loop
  inertiaMat[nbGrp,nbGrp]<-.Call("tmrsubmatrixinertia", dissmatrix, as.integer(ind[grpCond[[nbGrp]]]), PACKAGE="TraMineR")*sum(grpCond[[nbGrp]])
  #Computing residuals
#  print(inertiaMat)
  SCres<-sum(diag(inertiaMat)/grpSize)
# cat("Residuals:",SCres,"\n")
#    return(NULL)
#  print(diag(inertiaMat)/grpSize)
  if(SCres>currentSCres)return(NULL)
  #Fonction to comput inertia based on inertiaMat
  inertiaFunction<-function(inertiaMat,co,pop){
    #Take care to add one
    #return(((sum(inertiaMat[co,co])+sum(diag(inertiaMat[co,co])))/2)/pop)
    #New way, inertiaMat is triangular -> we can just sum the matrix
    return(sum(inertiaMat[co,co])/pop)
  }
#  print(sum(inertiaMat)/2+sum(inertia))
#  print(inertiaMat)
#  print(inertia)
  bestSCres<-currentSCres
  bestRegroup<-NULL
  allgroups<- 1:(nbGrp)
#  print(inertiaFunction(inertiaMat,inertia,allgroups,totpop))
  if(is.ordered(grp)){
    maxGrp<-nbGrp-1
  }else{
    maxGrp<-ceiling(nbGrp/2)
  }
  for(p in 1:maxGrp){
    if(is.ordered(grp)){
       combi<-list()
       combi[[1]]<-1:p
       if(has.na){
         combi[[2]]<-c(1:p,nbGrp)
       }
    }else{
       combi<-combn(nbGrp,p,simplify=FALSE)
    }
    for(co in combi){
       popc<-sum(grpSize[co])
    #  print(popc)
       popothc<-totpop-popc
       if(popc>minSize&&popothc>minSize){
          othc<-allgroups[!(allgroups %in% co)]
          ico<-inertiaFunction(inertiaMat,co,popc)
          iothc<-inertiaFunction(inertiaMat,othc,popothc)
          SCres<-ico+iothc
    #        print(paste(paste(lgrp[co], collapse=" "),"/",paste(lgrp[othc], collapse=" "),"=",SCres))
            #print(list(co=co,othc=othc,ico=ico,iothc=iothc,popc<-popc,popothc=popothc))
          if(SCres<bestSCres){
            bestSCres<-SCres
            bestRegroup<-list(co=co,othc=othc,ico=ico,iothc=iothc,popc=popc,popothc=popothc)
          }
       }
    }
  }
  if(is.null(bestRegroup))return(NULL)
   ret<-list(
    lpop=bestRegroup$popc,
    rpop=bestRegroup$popothc,
    lvar=bestRegroup$ico/bestRegroup$popc,
    rvar=bestRegroup$iothc/bestRegroup$popothc,
    SCres=bestSCres,
    varname=varname
   )
  ret$variable<-(grp %in% lgrp[bestRegroup$co])
  if(has.na){
    ret$variable[is.na(grp)]<-(nbGrp %in% bestRegroup$co)
  }
  if(is.ordered(grp)){
    if(has.na){
      ret$llabels<-paste("<=",llgrp[max(bestRegroup$co[bestRegroup$co<nbGrp])],sep="")
      ret$rlabels=paste(">",llgrp[max(bestRegroup$co[bestRegroup$co<nbGrp])],sep="")
      if(nbGrp %in% bestRegroup$co){
        ret$llabels<-paste(ret$llabels,",",llgrp[nbGrp],sep="")
      }else{
        ret$rlabels<-paste(ret$rlabels,",",llgrp[nbGrp],sep="")
      }
    }else{
      ret$llabels<-paste("<=",llgrp[max(bestRegroup$co)],sep="")
      ret$rlabels=paste(">",llgrp[max(bestRegroup$co)],sep="")
    }
  }else{
    ret$llabels=llgrp[bestRegroup$co]
    ret$rlabels=llgrp[bestRegroup$othc]
  }
  return(ret)
}
print.disstree<-function(x,quote=FALSE,digits=3,...){
  print("Dissimilarity tree", quote=quote,...)
  print(paste("Global R2:",format(x$adjustement$stat$PseudoR2,digits =digits)), quote=quote,...)
  print(x$root, quote=quote, digits=digits,...)
}
print.DissTreeNode<-function(x,quote=FALSE,digits=3,...){
  gap<-character(x$depth)
  gap[]<-"  "
  string<-paste(paste(gap,collapse=" "),"|--",x$label,"[",length(x$ind),"] var:",format(x$vardis,digits =digits),collapse="")
  print(string,quote=quote,...)
  if(!is.null(x$split)){
    print(paste(paste(gap,collapse=" ")," ","|->",x$split$varname," R2:",format(x$R2,digits=digits),collapse=""),quote=quote,...)
    for(i in x$children){
      print(i)
    }
  }
}
DTNseqplot<-function(ind, seqs, sortv,plotfunc,seqplot.title.cex,  seqplot.use.title=TRUE,seqplot.label.loc="main",seqplot.node.loc="main", seqplot.split.loc="sub",...){
  if(seqplot.use.title){
    top<- (as.integer(seqplot.label.loc=="main")+as.integer(seqplot.node.loc=="main")+as.integer(seqplot.split.loc=="main"))*seqplot.title.cex
    bottom <- (as.integer(seqplot.label.loc=="sub")+as.integer(seqplot.node.loc=="sub")+as.integer(seqplot.split.loc=="sub"))*seqplot.title.cex
    left<- (as.integer(seqplot.label.loc=="ylab")+as.integer(seqplot.node.loc=="ylab")+as.integer(seqplot.split.loc=="ylab"))*seqplot.title.cex
    par(mar=c(bottom,left,top,0), font.sub=2, mgp=c(0,0,0))
  }
  if(!is.null(sortv))
    plotfunc(seqs[ind,], sortv=sortv[ind],...)
  else
    plotfunc(seqs[ind,],...)
}
seqtree2dot<-function(tree, filename, seqs, plottype="seqdplot", imgLeafOnly=FALSE,sortv=NULL,...){
  if(plottype=="seqdplot"){
    plotfunc<-seqdplot
  }else if(plottype=="seqiplot"){
    plotfunc<-seqiplot
  }else if(plottype=="seqmtplot"){
    plotfunc<-seqmtplot
  }else if(plottype=="seqfplot"){
    plotfunc<-seqfplot
  }else{
    stop("Unknow plot type")
  }
  disstree2dot(tree, filename, imagedata=NULL,seqs=seqs, title.cex=3,sortv=sortv, imagefunc=DTNseqplot, plotfunc=plotfunc,
  use.title=TRUE,label.loc="main",node.loc="main", split.loc="sub",
  seqplot.use.title=TRUE,seqplot.label.loc="main",seqplot.node.loc="main", seqplot.split.loc="sub",seqplot.title.cex=3,...)
}

disstree2dot<-function(tree,filename,digits=3,imagefunc=NULL, imagedata=NULL,imgLeafOnly=FALSE,devicefunc="jpeg", imageext="jpg", device.arg=list(), use.title=TRUE,label.loc="main",node.loc="main", split.loc="sub",title.cex=1,...){
  dotfile<-paste(filename,".dot",sep="")
  node<-tree$root
  cat("digraph distree{\n",file=dotfile)
  if(!is.null(imagefunc)&&is.null(imagedata)){
    imagedata<-as.data.frame(node$ind)
  }
  DTN2DotInternal(def=TRUE,preced=filename,node=node,pos="none",digits=digits,
    imagefunc=imagefunc,imagedata=imagedata,imgLeafOnly=imgLeafOnly,dotfile=dotfile,devicefunc=devicefunc,imageext=imageext, device.arg=device.arg,use.title=use.title,label.loc=label.loc,node.loc=node.loc, split.loc=split.loc,title.cex,...)
  DTN2DotInternal(def=FALSE,preced=filename,node=node,pos="none",digits=digits,
    imagefunc=imagefunc,imagedata=imagedata,imgLeafOnly=imgLeafOnly,dotfile=dotfile,devicefunc=devicefunc,imageext=imageext, device.arg=device.arg,use.title=use.title,label.loc=label.loc,node.loc=node.loc, split.loc=split.loc,title.cex,...)
  cat("}\n",file=dotfile,append=TRUE)

}
DTN2DotInternal<-function(def,preced,node,pos,digits,imagefunc,imagedata,imgLeafOnly,dotfile,devicefunc, imageext, device.arg,use.title,label.loc=label.loc,node.loc,split.loc,title.cex,...){
  nodename<-paste(preced,"_",pos,sep="")
  if(def){
    stringcontentnode<- paste("Size: ",length(node$ind)," Var: ",format(node$vardis,digits =digits),sep="")
    if(!is.null(node$split)){
      stringcontentsplit<-paste("Split: ",node$split$varname," R2:",format(node$R2,digits=digits),sep="")
    }else { 
      stringcontentsplit<-""
    }
    if(!is.null(imagefunc)&&(!imgLeafOnly||is.null(node$split))){
      device.arg$file<-paste(nodename,imageext,sep=".")
      do.call(devicefunc,device.arg)
      imagefunc(imagedata[node$ind,],...)
      if(use.title){
        title.arg<- list(main = NULL, sub = NULL, xlab = NULL, ylab = NULL, line = NA, outer = FALSE, cex.main=title.cex,  cex.sub=title.cex, cex.lab=title.cex)
        title.arg[[node.loc]]<-stringcontentnode
        title.arg[[split.loc]] <- stringcontentsplit
        title.arg[[label.loc]]<-node$label
        if(label.loc==node.loc){
          title.arg[[label.loc]]<-paste(node$label,stringcontentnode, sep="\n")
        }
        if(label.loc==split.loc){
          title.arg[[label.loc]]<-paste(title.arg[[label.loc]],stringcontentsplit, sep="\n")
        }
        if(node.loc==split.loc){
          if(label.loc!=split.loc){
            title.arg[[node.loc]]<-paste(stringcontentnode,stringcontentsplit, sep="\n")
          }
        }
        do.call("title", title.arg)
      }
      dev.off()
      imgstr<-paste(" image=\"",nodename,".jpg\", imagescale=true,",sep="")
    }else{
      imgstr<-""
    }
    if(!use.title){
      str<-paste("\"",nodename,"\"[shape=box, ",imgstr," label=","\"",stringcontentnode, stringcontentsplit,"\"",sep="")
    }else {
      str<-paste("\"",nodename,"\"[shape=box, ",imgstr," label=\" \"", sep="")
    }
    cat(paste(str,"];\n",sep=""),file=dotfile,append=TRUE)

  }else{
    if(node$depth!=1)cat(paste("\"",preced,"\"->","\"",nodename,"\";\n",sep=""),file=dotfile,append=TRUE)#[label=","\"",node$label,"\"];\n",sep=""),...)
  }
  if(!is.null(node$split)){
    DTN2DotInternal(def=def,preced=nodename,node=node$children$left,pos="left",digits=digits,
      imagefunc=imagefunc,imagedata=imagedata,imgLeafOnly=imgLeafOnly,dotfile=dotfile,
      devicefunc=devicefunc,imageext=imageext, device.arg=device.arg,use.title=use.title,label.loc=label.loc,node.loc=node.loc, split.loc=split.loc,title.cex=title.cex,...)
    DTN2DotInternal(def=def,preced=nodename,node=node$children$right,pos="right",digits=digits,
      imagefunc=imagefunc,imagedata=imagedata,imgLeafOnly=imgLeafOnly,dotfile=dotfile,
      devicefunc=devicefunc,imageext=imageext, device.arg=device.arg,use.title=use.title,label.loc=label.loc,node.loc=node.loc, split.loc=split.loc,title.cex=title.cex,...)
  }
}


