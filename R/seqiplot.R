## =============================
## Plotting individual sequences
## =============================

seqiplot <- function(seqdata, group = NULL, main = NULL, ...) {
	seqplot(seqdata, group=group, type="i", main=main, ...)
}
seqIplot <- function(seqdata, group = NULL, main = NULL, ...) {
	seqplot(seqdata, group=group, type="I", main=main, ...)
}
