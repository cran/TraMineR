###Based on the program written for scipy (Python) by
###Ondrej Libiger and Matt Zapala
###
###Based on some of the methods presented in:
###McArdle, B.H. and Anderson, M.J. (2001).
###Fitting multivariate models to community data:
###a comment on distance-based redundancy analysis. Ecology, 82, 290-297.



gower_matrix <- function(diss, squared=TRUE) {
	n <- nrow(diss)
	# Creating centered (Gower's) matrix from the distance matrix:
	aux_centered_matrix <- diag(n) - tcrossprod(matrix(1, nrow=n))/n
	if (squared) {
		return(aux_centered_matrix %*% ((diss^2)/(-2)) %*% aux_centered_matrix) #to Optimise
	}
	else {
		return(aux_centered_matrix %*% (diss/(-2)) %*% aux_centered_matrix) #Optimise
	}
}

dissreg <- function(formula, data, R=1000, gower=FALSE,
										squared=TRUE, permutation="dissmatrix") {

	warning("dissreg function is deprecated. It has been renamed dissmfac.")
	return(dissmfac(formula, data, R, gower, squared=TRUE, permutation))
}

dissmfac <- function(formula, data, R=1000, gower=FALSE,
										squared=TRUE, permutation="dissmatrix") {
	
		#To ensure dissmatrix can be a dist object as well
		dissmatrix <- as.matrix(eval(formula[[2]], data, parent.frame()))
		formula[[2]] <- NULL
		#Model matrix from forumla
		predictor_frame <- model.frame(formula, data, drop.unused.levels = TRUE)
		seqok <- complete.cases(predictor_frame)
		#To make sure unused levels after removing NA values are actually removed
		predictor_frame <- model.frame(formula, data[seqok, ], drop.unused.levels = TRUE)
		predictor_terms <- attr(predictor_frame, "terms")
		var_complete_name <- c(attr(predictor_terms, "term.labels"), "Total")
		predictor <- model.matrix(formula, predictor_frame)
		if (!gower) {
				g_matrix <- gower_matrix(dissmatrix[seqok, seqok], squared)
		} else {
			g_matrix <- dissmatrix[seqok, seqok]
		}
		n <- nrow(g_matrix)

		var_list <- attr(predictor, "assign")
		
		#Building Hat matrix list (one per model)
		#decomposition
		var_list_index <- (1:length(var_list))
		#List of variable
		var_names <- unique(var_list)
		#Number of variable (minus cte)
		nterms <- length(var_names) - 1
		SCtot <- sum(diag(g_matrix))
		#complete_model <- nterms+1
		#If we permute dissmatrix, we can build all hat_matrix in advance
		if (permutation == "dissmatrix") {
			hat_matrix_list <- list()
			k_list <- list()
			#We build all "backward" hat_matrix based on QR decomposition
			for (var in 1:(nterms+1)) {
				pred <- predictor[, c(var_list_index[var_list!=var])]
				k_list[var] <- ncol(pred)
				qr_matrix <- qr(pred)
				q_matrix <- qr.Q(qr_matrix)
				hat_matrix_list[[var]] <- tcrossprod(q_matrix)
			}

			#Performing gc before all permutations
			gc()
			bts <- boot(predictor, internalbootmatrixregression, R=R, sim="permutation",
				stype="i", g_matrix=g_matrix, SCtot=SCtot, hat_matrix_list=hat_matrix_list,
				k_list=k_list, n=n, nterms=nterms)
		} else if (permutation=="model") {
			#If we permute the models, then hat_matrix is in function...
			gc()
			bts <- boot(predictor, internalbootmatrixregression2, R=R, sim="permutation",
				stype="i", g_matrix=g_matrix, SCtot=SCtot, n=n, nterms=nterms,
				var_list_index=var_list_index, var_list=var_list)
		} else {
			stop("Unknown permutation method")
		}
		#Computing Pvalue based on permutations tests
		p_value <- double()
		for (var in 1:(nterms+1)) {
			p_value[var] <- sum(bts$t[, var]>bts$t0[var])/bts$R
		}
		#Results
		mfac <- data.frame(Variable=var_complete_name,
				PseudoF=bts$t0[1:(nterms+1)],
				PseudoR2=bts$t0[(nterms+2):length(bts$t0)],
				p_value=p_value
			)
		#Var names
		dmf <- list(mfac = mfac, call = match.call(), perms=bts, perm_method=permutation)
		class(dmf) <- "dissmultifactor"
		return(dmf)

}
print.dissmultifactor <- function(x, ...) {
	print(x$mfac)
	invisible(x)
}
print.dissregression <- function(x, ...) {
	warning("dissreg function is deprecated. It has been renamed dissmfac.")
	print(x$mreg)
	invisible(x)
}
internalbootmatrixregression <- function(predictor, ind, g_matrix, SCtot,
																					hat_matrix_list, k_list, n, nterms) {
	g_perm <- g_matrix[ind, ind]
	R2_list <- numeric()
	F_list <- numeric()
	complete_model <- nterms+1
	m <- k_list[[complete_model]]
	SCexpC <- sum(hat_matrix_list[[complete_model]] * g_perm)
	SCresC <- SCtot-SCexpC
	R2_list[complete_model] <- (SCexpC/SCtot)
	F_list[complete_model] <- (SCexpC/(m-1)) / ((SCresC)/(n-m))
	if (nterms==1) {
		R2_list[1] <- R2_list[complete_model]
		F_list[1] <- F_list[complete_model]
		return(c(F_list, R2_list))
	}
	for (var in 1:(nterms)) {
		k <- k_list[[var]]
		#No need to transpose since centered_matrix is symmetric
		SCexp <- sum(hat_matrix_list[[var]] * g_perm)
		#traceHatCentered <- .Call("tmrpermtracematrix", hat_matrix, centered_matrix, 1:n)
		F_list[[var]] <- ((SCexpC-SCexp)/(m-k)) / ((SCresC)/(n-m-1))
		#calculate proportion of variance explained
		R2_list[[var]] <- (SCexpC-SCexp) / SCtot
		# Printing single regression analysis output
	}
	return(c(F_list, R2_list))
}
internalbootmatrixregression2 <- function(predictor, ind, g_matrix, SCtot, n,
																					nterms, var_list_index, var_list) {
	perm_predictor <- predictor[ind, ]
	complete_model <- nterms+1
	R2_list <- numeric()
	F_list <- numeric()
	m <- ncol(perm_predictor)
	qr_matrix <- qr(perm_predictor)
	q_matrix <- qr.Q(qr_matrix)
	hat_matrix <- tcrossprod(q_matrix)
	SCexpC <- sum(hat_matrix * g_matrix)
	SCresC <- SCtot-SCexpC
	R2_list[complete_model] <- (SCexpC/SCtot)
	F_list[complete_model] <- (SCexpC/(m-1)) / ((SCresC)/(n-m))
	if (nterms==1) {
		R2_list[1] <- R2_list[complete_model]
		F_list[1] <- F_list[complete_model]
		return(c(F_list, R2_list))
	}
	for (var in 1:nterms) {
		pred <- perm_predictor[, c(var_list_index[var_list!=var])]
		k <- ncol(pred)
		qr_matrix <- qr(pred)
		q_matrix <- qr.Q(qr_matrix)
		hat_matrix <- tcrossprod(q_matrix)
		SCexp <- sum(hat_matrix * g_matrix)
		F_list[[var]] <- ((SCexpC-SCexp)/(m-k)) / ((SCresC)/(n-m-1))
		#calculate proportion of variance explained
		R2_list[[var]] <- (SCexpC-SCexp) / SCtot
		# Printing single regression analysis output
	}
	return(c(F_list, R2_list))
}
