## ========================
## State distribution plot
## ========================

seqdplot <- function(seqdata, group=NULL, title=NULL, ...) {
	plot(seqdata, group=group, type="d", title=title, ...)
}
