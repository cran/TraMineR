#######################################################
###Compare groups and returns most discriminating func
#######################################################

#seqecmpgroup.survival <- function(time, event, var, stest) {
#	library(survival)
#	test <- survdiff(Surv(time, event) ~ var, rho=stest)
#	pval <- 1-pchisq(test$chisq,length(test$n)-1)

#	pval
#}

#seqecmpgroup.chisq<-function(group,seqp){
#  chi<- chisq.test(group,seqp)
#  return(chi$statistic)
#}

seqecmpgroup<-function(subseq, seqe, group, method="bonferroni", p.valuelimit=0.05){
  seqecmpgroup.chisq<-function(group,seqp,ntest, p.valuelimit=0.05){
    chi<- chisq.test(group,seqp)
    if(chi$p.value>p.valuelimit)return(NA)
    return(chi$statistic)
  }

  seqecmpgroup.chisq.bonferroni<-function(group,seqp,ntest, p.valuelimit=0.05){
    chi<- chisq.test(group,seqp)
    ret<-(1-(1-chi$p.value)^ntest)
    if(ret>p.valuelimit)return(NA)
    return(ret)
  }
  if(method=="bonferroni"){
    testfunc<-seqecmpgroup.chisq.bonferroni
    decreasing<-FALSE
    seqmatrix<-seqeapplysub(subseq, seqe, method="presence")
  }else if(method=="chisq"){
    testfunc<-seqecmpgroup.chisq
    decreasing<-TRUE
    seqmatrix<-seqeapplysub(subseq, seqe, method="presence")
  }
  
	res<-data.frame()
	ntest<-ncol(seqmatrix)
	index<-integer(ntest)
	stat<-numeric(ntest)
	for (i in 1:ncol(seqmatrix)) {
    index[i]<-i
    stat[i]<-testfunc(group,seqmatrix[,i],ntest=ntest,p.valuelimit)
	}
  cres<-order(as.double(stat),decreasing =decreasing)#[!is.na(as.double(res$stat))]
  index<-index[cres]
  stat<-stat[cres]
  res<-data.frame(index,stat)
  rownames(res)<-colnames(seqmatrix)[cres]
	return(res)
}

