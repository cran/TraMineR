seqeconstraint<-function(maxGap=-1, windowSize=-1, ageMin=-1, ageMax=-1,ageMaxEnd=-1, countMethod=1){
  ret<-list()
  ret$maxGap <- maxGap
  ret$windowSize <- windowSize
  ret$ageMin <- ageMin
  ret$ageMax <- ageMax
  ret$ageMaxEnd <-ageMaxEnd
  ret$countMethod <- countMethod
  class(ret)<-"seqeconstraint"
  return(ret)
}

print.seqeconstraint<-function(x,...){
	z<-data.frame(Constraint=names(x),Value=as.numeric(x))
	z <- z[z$"Value"!=-1,]
	if (z[z$"Constraint"=="countMethod","Value"] == 1) {
		z[z$"Constraint"=="countMethod","Value"] <- "One by sequence"
	}
	if (z[z$"Constraint"=="countMethod","Value"] == 2) {
		z[z$"Constraint"=="countMethod","Value"] <-  "Several by sequence"
	}
	if(nrow(z) > 0) { print(z,row.names=FALSE,...) }
}