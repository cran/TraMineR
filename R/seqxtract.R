## EXTRACT SEQUENCES FROM A DATA SET 

seqxtract <- function(data, var=NULL, data.frame=FALSE) {
	
	## Extracting the sequences from the data set
	if (is.null(var)) seqdata <- data
	else seqdata <- subset(data,,var)

	if (data.frame==FALSE) seqdata <- as.matrix(seqdata)

	return(seqdata)

	}
