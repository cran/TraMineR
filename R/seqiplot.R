## =============================
## Plotting individual sequences
## =============================

seqiplot <- function(seqdata, group=NULL, title=NULL, tlim=NULL, sortv=NULL, ...) {
	plot(seqdata, group=group, type="i", title=title, tlim=tlim, sortv=sortv, ...)
}
