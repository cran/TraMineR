## ============================================
## TRANSLATION BETWEEN SEQUENCE REPRESENTATIONS
## ============================================

seqformat <- function(data, var=NULL, id=NULL, 
	from=NULL, to=NULL, compressed=FALSE,
	nrep=NULL, tevent=NULL, stsep, covar=NULL,
	SPS.in=list(xfix="()", sdsep=","),
	SPS.out=list(xfix="()", sdsep=","),
	begin=NULL, end=NULL, status=NULL, 
	process=TRUE, pdata=NULL, pvar=NULL, limit=100, overwrite=TRUE, 
	fillblanks=NULL, tmin=NULL, tmax=NULL) {

	## Checking the format
	listform <- c("STS","SPS","SPELL","DSS","SRS","TSE")
	if (!from %in% listform | !to %in% listform)	
		stop("input and output formats must be one of: ", listform)

	if (!is.null(covar)) covariates <- subset(data,,covar)

	if (!is.null(id)) ident <- as.matrix(subset(data,,id))
	else ident <- NULL

	## ========
	## From STS
	## ========
	if (from=="STS") {
		## Extracting the sequences from the data set
		seqdata <- seqxtract(data, var)

		if (missing(stsep)) {
			sf <- seqfcheck(seqdata)
			if (sf %in% c("-",":")) seqdata <- seqdecomp(seqdata,sep=sf)
		}
		else seqdata <- seqdecomp(seqdata,sep=sf)

		trans <- seqdata
	}

	## ========
	## From SPS
	## ========
	if (from=="SPS") {
		## Extracting the sequences from the data set
		seqdata <- seqxtract(data, var)

		sf <- seqfcheck(seqdata)
		if (sf %in% c("-",":")) seqdata <- seqdecomp(seqdata,sep=sf)

		trans <- SPS_to_STS(seqdata, spsformat=SPS.in)
	}

	## ==========
	## From SPELL
	## ==========
	if (from=="SPELL") {
		if (!is.null(var)) data <- subset(data,,var)

		if (!is.null(id)) ident <- as.matrix(subset(data,,id))
		else ident <- as.matrix(subset(data,,1))

		if (!is.null(begin)) sp2 <- as.matrix(subset(data,,begin))
		else sp2 <- as.matrix(subset(data,,2))

		if (!is.null(end)) sp3 <- as.matrix(subset(data,,end))
		else sp3 <- as.matrix(subset(data,,3))

		if (!is.null(status)) sp4 <- as.matrix(subset(data,,status))
		else sp4 <- as.matrix(subset(data,,4))

		seqdata <- data.frame(ident, sp2, sp3, sp4)

		## Extracting the sequences from the data set
		## seqdata <- seqxtract(data, var, data.frame=TRUE)
		
		trans <- BIOSPELL_to_STS(seqdata=seqdata, 
			process=process, pdata=pdata, pvar=pvar, 
			limit=limit, overwrite=overwrite, fillblanks=fillblanks, 
			tmin=tmin, tmax=tmax)	
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
		out <- STS_to_SPS(seqdata=trans, spsformat=SPS.out)
		nbout <- seqdim(out)[1]

		if (compressed) out <- seqconc(out)
		}

	## To Distinct-State-Sequence format
	if (to=="DSS") {
		out <- STS_to_DSS(trans)
		nbout <- seqdim(out)[1]

		if (compressed) out <- seqconc(out)
	}

	## STS
	if (to=="STS") {
		out <- trans
		nbout <- seqdim(out)[1]

		if (compressed) out <- seqconc(out)
		}

	if (to=="SRS") {
		out <- STS_to_SRS(trans,nrep)
		if (!is.null(covar)) out <- merge(out,data.frame(id=seq(1:nbin),covariates))
		nbout <- nrow(out)
	}

	if (to=="TSE") {
		out <- STS_to_TSE(trans, ident, tevent)
		nbout <- nrow(out)
		}

	if (to!="STS") message(" [>] STS sequences converted to ",nbout," ", to," seq./rows")
	return(out)

}
