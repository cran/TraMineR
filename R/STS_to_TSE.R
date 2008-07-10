## ================================
## Convert from STS to TSE format
## ================================

STS_to_TSE <- function(seqdata,id=NULL,tevent) {
	nseq <- seqdim(seqdata)[1]
	slength <- seqdim(seqdata)[2]
	statl <- seqstatl(seqdata)

	cat(" => Converting",nseq, "sequences to TSE format, please wait\n")
	## cat("     Transition-to-events conversion matrix:\n")
	## print(tevent)

	if (is.null(id)) id <- 1:nseq

	trans <- data.frame(id=NULL, time=NULL, event=NULL)

	for (i in 1:nseq) {
    ## First status=> entrance event (diagonal of tevent)
   	s1 <- seqdata[i,1]
   	if(!is.na(tevent[s1,s1])){ #if NA, we don't generate an event
     	levent <- tevent[s1,s1]
    	levent <- strsplit(levent,",")[[1]]
    	for (k in 1:length(levent))
    		trans <- rbind(trans, data.frame(id=id[i],time=0,event=levent[k]))
  	} #end if
    ## Rest of the sequence
		for (j in 1:(slength-1)) {
			s1 <- seqdata[i,j] ## Status at time t
			s2 <- seqdata[i,j+1] ## Status at time t+1
			
			if (!is.na(s1) & !is.na(s2) & s1!="" & s2!="") {
				s1 <- as.character(s1) ## Status at time t
				s2 <- as.character(s2) ## Status at time t+1

				if (s2!=s1) {
					levent <- tevent[s1,s2]
					levent <- strsplit(levent,",")[[1]]
					for (k in 1:length(levent))			
						trans <- rbind(trans, data.frame(id=id[i],time=j,event=levent[k]))
				}
			}
		}
	}
	return(trans)
}
