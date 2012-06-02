seqedist<-function(seqe, idcost, vparam, interval=TRUE, norm=TRUE){
    norm<-as.integer(norm)
    interval<-as.integer(interval)
    return(.Call(TMR_tmrseqedist, seqe, idcost, vparam, norm, interval));
}