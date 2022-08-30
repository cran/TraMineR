## =============================
## Plotting medoids of relative frequeny group of sequences
## =============================

seqrfplot <- function(seqdata, group = NULL, main = NULL, ...) {
	seqplot(seqdata, group=group, type="rf", main=main, ...)
}
