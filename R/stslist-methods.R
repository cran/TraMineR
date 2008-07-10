## ===========================
## Methods for stslist objects
## ===========================

print.stslist <- function(x,format='STS',extended=FALSE,...) {
	if (format=='STS') {
		if (extended==FALSE) {
			x <- seqconc(x)
			print(x,quote=FALSE)
		} else NextMethod("print")
	}
	if (format=='SPS') {
		x <- seqformat(x,from='STS',to='SPS') 
		print(x,quote=FALSE)
	}
}

plot.stslist <- function(x,...) {
	seqiplot(x)
}

"[.stslist" <- function(x,i,j,drop=FALSE) {
	## Specialized only for column subscript
	## If one column we keep the original data.frame method
	## Otherwise we copy attributes and update "start" value
	if (!missing(j) && length(j)>1) {
		## Storing the attributes
		a <- attr(x,"alphabet")
	     s <- attr(x,"start")
	     m <- attr(x,"missing")
		cpal <- attr(x,"cpal")
		lab <- attr(x,"labels")

		## Applying method
	     x <- NextMethod("[")
	
		## Redefining attributes
		attr(x,"alphabet") <- a
		attr(x,"missing") <- m
	     attr(x,"start") <- s-1+j[1]
		attr(x,"cpal") <- cpal
		attr(x,"labels") <- lab
		return(x)
    }
    NextMethod("[")
 }


## "[.stslist" <- function(x,...) {
## 	NextMethod("[")
## }

Math.stslist <- function(...){
 stop("Invalid operation on sequences")
}


