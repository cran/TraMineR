## ================================
## PLot of the sequences frequency
## ================================

seqfplot <- function(seqdata, group=NULL, title=NULL, tlim=NULL, pbarw=FALSE, ...) {
	plot(seqdata, group=group, type="f", title=title, tlim=tlim, pbarw=pbarw, ...) 
}
