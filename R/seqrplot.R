## ============================================
## PLOT A REPRESENTATIVE SEQUENCE
## ============================================

seqrplot <- function(seqdata, group=NULL, title=NULL, method='dist', 
	dist.matrix=NULL, dist.method='LCS', norm=FALSE, indel=1, sm, with.miss = FALSE,
	pbarw=TRUE, entropy=FALSE, mline=FALSE, fline=FALSE, 
	dist.rep=TRUE, dist.center=TRUE,
	 ...) {

	plot(seqdata, group=group, type="r", title=title, 
	method=method, dist.matrix=dist.matrix, 
	## dist.method=dist.method, norm=norm, indel=indel, sm=sm, with.miss=with.miss,
	pbarw=pbarw, entropy=entropy, mline=mline, fline=fline, 
	dist.rep=dist.rep, dist.center=dist.center,
	 ...)

}
