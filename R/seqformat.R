## ============================================
## TRANSLATION BETWEEN SEQUENCE REPRESENTATIONS
## ============================================

seqformat <- function(data, var=NULL, id=NULL, from=NULL, to=NULL, nrep=NULL, tevent=NULL, stsep="-", covar=NULL) {

	## Checking the format
	listform <- c("STS","SPS","SPS2","SPELL","DSS","SRS","TSE")
	if (!from %in% listform)	
		stop(" [!] input format must be one of: ", listform)

	if (!to %in% listform) {
			stop(" [!] output format must be one of: ", listform, "Goodbye!")
			return()
		}

	if (!is.null(covar)) covariates <- subset(data,,covar)
	if (!is.null(id)) ident <- as.matrix(subset(data,,id))
	else ident <- NULL

	if (from=="STS") {
		## Extracting the sequences from the data set
		seqdata <- seqxtract(data, var)
		if (seqfcheck(seqdata) %in% c("X","-X")) seqdata <- seqconc(seqdata)

		trans <- seqdata
	}

	if (from=="SPS") {
		## Extracting the sequences from the data set
		seqdata <- seqxtract(data, var)
		if (seqfcheck(seqdata) %in% c("X","-X")) seqdata <- seqconc(seqdata)

		trans <- SPS_to_STS(seqdata,1,stsep)
	}

	if (from=="SPS2") {
		## Extracting the sequences from the data set
		seqdata <- seqxtract(data, var)
		if (seqfcheck(seqdata)=="X") seqdata <- seqconc(seqdata)

		trans <- SPS_to_STS(seqdata,2,stsep)
	}

	if (from=="SPELL") {
		## Extracting the sequences from the data set
		seqdata <- seqxtract(data, var, data.frame=TRUE)
		
		trans <- SPELL_to_STS(seqdata)	
		## ident <- unique(ident)
	}

	## ===============
	## INTERNAL FORMAT
	## ===============
	rm(seqdata)
	nbin <- seqdim(trans)[1]
	if (from != "STS") message(" [>] ", from," data converted into ",nbin," STS sequences")

	## =============
	## OUTPUT FORMAT
	## =============
	if (to=="SPS") {
		out <- STS_to_SPS(trans,1)
		nbout <- seqdim(out)[1]
		}

	if (to=="SPS2") {
		out <- STS_to_SPS(trans,2)
		nbout <- seqdim(out)[1]
	}

	## To Distinct-State-Sequence format
	if (to=="DSS") {
		out <- STS_to_DSS(trans)
		nbout <- seqdim(out)[1]
	}

	## STS
	if (to=="STS") {
		out <- trans
		nbout <- seqdim(out)[1]
		}

	if (to=="SRS") {
		out <- STS_to_SRS(trans,nrep)
		if (!is.null(covar)) out <- merge(out,data.frame(id=seq(1:nbin),covariates))
		nbout <- nrow(out)
	}

	if (to=="TSE") {
		trans <- seqdecomp(trans)
		out <- STS_to_TSE(trans, ident, tevent)
		nbout <- nrow(out)
		}

	if (to!="STS") message(" [>] STS sequences converted to ",nbout," ", to," seq./rows")
	return(out)

}
