## ====================
## Plotting the legend
## ====================

TraMineR.legend <- function(pos, text, colors, cex=1) {

	nbstat <- length(text)

	## Saving graphical parameters
	savepar <- par(no.readonly = TRUE)

	## Computing some parameters for the legend's plotting
	leg.ncol <- if (pos=="bottom") ceiling(nbstat/2) else 1
		## leg.ncol <- if (round(nbstat/3,0)>1) ceiling(nbstat/2) else 2
		leg.inset <- -0.2 + ((2-leg.ncol)*0.025)

	par(mar = c(1, 1, 0.5, 1) + 0.1, xpd=FALSE)

	plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
		## legend(position, fill = cpal, legend = ltext, cex = fontsize)

	legend(pos,
		# inset=c(0,leg.inset),
		legend=text,
		fill=colors,
		ncol=leg.ncol,
		bty="o",
		cex=cex)

	## Restoring graphical parameters
	par(savepar)
}
