###########################
## Compute discrepancy
###########################

dissvar <- function(diss) {
	if (inherits(diss, "dist")) {
		return(sum(diss)/(attr(diss, "Size")^2))
	} else if (is.matrix(diss)) {
		return(sum(diss)/(2*(nrow(diss)^2)))
	} else {
		stop("diss argument should be a dist object or a dissimilarity matrix")
	}
}