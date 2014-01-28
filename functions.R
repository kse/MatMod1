# Given a vector of values, calculates the fractile test of the data.
fractile <- function (x) {
	obs <- sort(unique(x));#{{{
	k <- c();
	p <- c();
	a <- c();

	n <- length(x);

	for(i in seq(1, length(obs))) {
		a[i] = length(which(x == obs[i]));

		if(i == 1) {
			k[i] = a[i];
		} else {
			k[i] = k[i-1] + a[i];
		}

		klow <- i - 1;
		if(klow == 0) {
			klow <- 1;
		}

		p[i] = 100* (sum(k[klow : i])/(2*n));
	}

	q <- qnorm(p*0.01);

	y <- data.frame(Observation=obs, 
					Number=a,
					Cumulative=k,
					Probability=p,
					Fractile=q
					);

	return(y);#}}}
}

# Given a vector of data points, calculates S, USS; SSD s2 and so on.
# Also gives confidence intervals
mm_analyzePoints <- function(x, f = NULL) {

	n <- length(x);#{{{
	df <- 0;

	if(!is.null(f)) {
		df <- f
	} else {
		df <- n - 1
	}

	S = sum(x);
	SSD = sum((x - mean(x))^2);
	s2 = (1/df)*SSD;
	USS = SSD + S^2/n;
	stderr = sqrt(s2/n);

	top <- (df*s2);

	mu_low = mean(x) - sqrt(s2/n)*qt(0.975, df);
	mu_high = mean(x) + sqrt(s2/n)*qt(0.975, df);
	
	sigsqr_low = top/qchisq(0.975, n-1)
	sigsqr_high = top/qchisq(0.025, n-1)

	y <- data.frame(n=n,
					  S=S,
					  xmean=S/n,
					  USS=USS,
					  S2_div_n = S^2/n,
					  SSD=SSD,
					  f=df,
					  s2=s2,
					  mu_low = mu_low,
					  mu_high = mu_high,
					  sigsqr_low = sigsqr_low,
					  sigsqr_high = sigsqr_high,
					  stderr= stderr
					  );
	return(y);#}}}
}

mm_analyzeSums <- function(USS = NULL, S = NULL, n = NULL, f = NULL) {
	SSD = USS - S^2/n;#{{{
	if(is.null(f)) {
		f = n - 1;
	}

	s2 = SSD/f;
	stderr = sqrt(s2/n);
	mean = S/n;

	top <- f*s2;

	mu_low = mean - sqrt(s2/n)*qt(0.975, f);
	mu_high = mean + sqrt(s2/n)*qt(0.975, f);
	
	sigsqr_low = top/qchisq(0.975, f)
	sigsqr_high = top/qchisq(0.025, f)

	y <- data.frame(n=n,
					  S=S,
					  xmean=S/n,
					  USS=USS,
					  S2_div_n=S^2/n,
					  SSD=SSD,
					  f=f,
					  s2=s2,
					  mu_low = mu_low,
					  mu_high = mu_high,
					  sigsqr_low = sigsqr_low,
					  sigsqr_high = sigsqr_high,
					  stderr= stderr
					  );
	return(y);#}}}
}

mm_addPointsAnalysis <- function(rows, x) {
	new <- mm_analyzePoints(x);#{{{

	d <- rbind(rows, new);
	return(d);#}}}
}

mm_addSumsAnalysis <- function(rows = NULL, USS = NULL, S = NULL, n = NULL) {
	new <- mm_analyzeSums(USS, S, n);#{{{

	d <- rbind(rows, new);
	return(d);#}}}
}

mm_estimateMean <- function(mean, data) {
	y <- (data$S/data$n - mean)/(sqrt(data$s2/data$n));#{{{

	pobs <- 2*(1 - pt(abs(y), df = data$n - 1))

	cat("t(x) =", y, ", Pobs=", pobs, "\n");#}}}
}

mm_estimateVariance <- function(var, data) {
	val <- ((data$n-1)*data$s2)/var;#{{{
	res = 0

	if(val <= (data$n - 1)) {
		res <- 2* pchisq(val, data$n-1)
	} else {
		res <- 2* (1 - pchisq(val, data$n-1))
	}

	cat("fs^2/sigma^2_0: =", val, ", Pobs =", res, "\n")#}}}
}

mm_2samples <- function(r, sigsqr = NULL, mu = NULL) {
	samples <- length(r[,1])#{{{

	if(samples != 2) {
		cat("Unexpected length of data, need length 2\n");
		return();
	}

	# Estimate if the variance is the same
	max <- if(r$s2[1] >  r$s2[2]) 1 else 2;
	min <- if(r$s2[1] <= r$s2[2]) 1 else 2;
	F <- r$s2[max]/r$s2[min]
	pF <- pf(F, r$f[max], r$f[min])
	pobs_SameVariance = 2*(1 - pF);
	cat ("------ Estimation of common variance: ------\n");
	cat ("F:", F, "\nF_F(f_num, f_denom)(F)", "\n");
	cat ("\t= F_F(", r$f[max], ",", r$f[min], ")(", F, ") =", pF, "\n");
	cat ("p_obs common variance:", pobs_SameVariance, "\n")
	s1sqr <- (r$f[1]*r$s2[1] + r$f[2]*r$s2[2])/(r$f[1] + r$f[2]);

	if(pobs_SameVariance > 0.05) {
		cat("It can be assumed that the variance is the same between the two",
			"samples\n");
		cat("New variance for both samples:", s1sqr, "\n");
		if(!is.null(sigsqr)) {
			f <- sum(r$n) - 2;
			cat("\tTest for sigma^2 =", sigsqr, "\n");
			val <- (f*s1sqr)/sigsqr;
			res = 0

			if(val <= f) {
				res <- 2* pchisq(val, f)
			} else {
				res <- 2 * (1 - pchisq(val, f))
			}

			cat("\tfs^2/sigma^2_0: =", val, ", Pobs =", res, "\n")
		}
	} else {
		cat("We cannot assume that the variance is the same between the two",
			"samples\n");
	}
	cat("\n");

	# Estimate if the mean is the same
	cat ("------ Estimation of common mean: ------\n");
	tx <- 0
	pobs_SameMean <- 0
	df <- 0 # fBar
	if(pobs_SameVariance <= 0.05) {
		tx <- (r$S[1]/r$n[1] - r$S[2]/r$n[2])/sqrt(r$s2[1]/r$n[1] + r$s2[2]/r$n[2])

		df <- (r$s2[1]/r$n[1] + r$s2[2]/r$n[2])^2
		df <- df/((r$s2[1]/r$n[1])^2/r$f[1] + (r$s2[2]/r$n[2])^2/r$f[2])
		pobs_SameMean <- 2*(1-pt(abs(tx), df))
		cat("fBar:", df, "\n")
	} else {
		tx <- (r$S[1]/r$n[1] - r$S[2]/r$n[2])/sqrt(s1sqr*(1/r$n[1] + 1/r$n[2]))
		pobs_SameMean <- 2*(1-pt(abs(tx), sum(r$f)))
	}

	cat("s1sqr:", s1sqr, "t(x):", tx, "\n");
	cat("p_obs for same mean:", pobs_SameMean, "\n");
	if(pobs_SameMean > 0.05) {
		cat("It can be assumed that the mean is the same for both\n")
	} else {
		cat("We cannot assume the mean is the same for both rows\n")
		if(pobs_SameVariance > 0.05) {
			confval <- sqrt(s1sqr/r$n[1])*qt(0.975, sum(r$n) - 2)
			cat("\tC95(mu1) = [", r$xmean[1] - confval, ",", r$xmean[1] + confval, "]\n")

			confval <- sqrt(s1sqr/r$n[2])*qt(0.975, sum(r$n) - 2)
			cat("\tC95(mu2) = [", r$xmean[2] - confval, ",", r$xmean[2] + confval, "]\n")
		}
	}

	# Confidence interval for the variance
	if(pobs_SameVariance > 0.05) {
		cat("\n------ C95(sigma^sqrt) ------\n")
		num <- (sum(r$f)*s1sqr)

		low <- num/qchisq(0.975, sum(r$f))
		high <- num/qchisq(0.025, sum(r$f))
		cat("[", low, ",", high, "]\n");
	}

	# Confidence interval for the mean
	x1mean <- r$S[1]/r$n[1];
	x2mean <- r$S[2]/r$n[2];
	if(pobs_SameVariance > 0.05) {
		val <- sqrt(s1sqr*(1/r$n[1] + 1/r$n[2])) * qt(0.975, r$n[1] + r$n[2] -
													  2);
		low  <- (x1mean - x2mean) - val
		high <- (x1mean - x2mean) + val

		cat("\n------ C95(mu1 - mu2) ------\n")
		cat("[", low, ",", high, "]\n");

		low  <- (x2mean - x1mean) - val
		high <- (x2mean - x1mean) + val
		cat("\n------ C95(mu2 - mu1) ------\n")
		cat("[", low, ",", high, "]\n");
	} else {
		val <- sqrt(r$s2[1]/r$n[1] + r$s2[2]/r$n[2]) * qt(0.975, df)
		low  <- (x1mean - x2mean) - val
		high <- (x1mean - x2mean) + val

		cat("\n------ C95(mu1 - mu2) ------\n")
		cat("[", low, ",", high, "]\n");

		low  <- (x2mean - x1mean) - val
		high <- (x2mean - x1mean) + val
		cat("\n------ C95(mu2 - mu1) ------\n")
		cat("[", low, ",", high, "]\n");
	}

	#}}}
}

mm_enumerateSums <- function(r) {
	j <- NULL;

	i <- 1
	while(i <= length(w[,1])) {
		j <- mm_addSumsAnalysis(j, r[i,1], r[i,2], r[i,3]);
		i <- i+1;
	}

	return(j);
}

# Takes a dataframe generated by mm_addSumsAnalysis or mm_analyzeSums/points
# and returns values for the end row, including confidence intervals.
mm_ksamples <- function(r) {
	#{{{
	k <- length(r[,1])

	if(k == 2) {
	}

	f1 <- sum(r$f)
	s1_sqr <- sum(r$SSD)/f1;

	lnq <- f1 * log(s1_sqr) - sum(r$f * log(r$s2));
	mmC <- 1 + (1/(3*(k - 1))) * (sum(1/(r$f)) - 1/f1)

	ba <- lnq/mmC;

	cat("---- Test for same Variance ----\n");
	cat("s1_sqr:", s1_sqr, "\n")
	cat("LNQ: ", lnq, "\nC:", mmC, "\nBartlett:", ba, "~~ ChiSq(", k-1, ")\n");
	pobs_SameVariance <- 1 - pchisq(ba, k-1)
	cat("pobs:", pobs_SameVariance, "\n")

	if(pobs_SameVariance > 0.05) {
		cat("It is safe to assume that the samples have the same variance\n\n")
	} else {
		cat("We cannot assume that the samples have the same variance\n\n");
	}

	SSD2 <- sum(r$S2_div_n) - sum(r$S)^2/sum(r$n);
	s2_sqr <- SSD2/(k-1)
	F <- s2_sqr/s1_sqr
	pobs_SameMean <- 1 - pf(F, k-1, f1)
	cat("---- Test for same mean ----\n");
	cat("SSD2:", SSD2, "\n")
	cat("s2_sqr:", s2_sqr, "\n");
	cat("F:", F, "~~ F(", k- 1, ",", f1, ")\n");
	cat("pobs:", pobs_SameMean, "\n");
	if(pobs_SameMean > 0.05) {
		cat("We can safely assume that the mean is the same\n\n");
	} else {
		cat("We cannot say that the samples have the same mean\n\n");
	}

	if(pobs_SameMean > 0.05 && pobs_SameVariance > 0.05) {
		cat("We now have that all the samples can be though of as one sample\n\n")
		j <- mm_analyzeSums(sum(r$USS), sum(r$S), sum(r$n), sum(r$n) - 1);
		return(j)
	}
	#}}}
}

mm_linreg <- function(x = NULL, t = NULL, n = NULL, Sx = NULL, St = NULL, 
					  USSx = NULL, USSt = NULL, SPxt = NULL, slope = NULL,
					  intercept = NULL) {
	#{{{

	#n;
	#Sx;
	#St;
	#USSx;
	#USSt;
	#SPxt;

	if(!is.null(x) && !is.null(t)) {
		if(!is.vector(x) || !is.vector(t)) {
			cat("X and T need to be vectors\n");
			return();
		}
		if(length(x) != length(t)) {
			print("ERROR: x and t are not equal length");
			return();
		}

		n <- length(x);
		Sx <- sum(x);
		St <- sum(t);
		USSx <- sum(x^2);
		USSt <- sum(t^2);
		SPxt <- sum(x*t);
	} else if(!is.null(n)
			  && !is.null(Sx)
			  && !is.null(St)
			  && !is.null(USSx)
			  && !is.null(USSt)
			  && !is.null(SPxt)) {
	} else {
		cat("You need to provide either x or t vectors, or the given values\n");
		return();
	}

	f <- n - 2;
	SSDx <- (USSx - Sx^2/n)
	SSDt <- (USSt - St^2/n)
	SPDxt <- SPxt - (Sx*St)/n;
	bhat <- SPDxt / SSDt;
	ahat <- 1/n * (Sx - bhat*St);
	SSD02 <- SSDx - SPDxt^2/SSDt;
	s02_square <- (1/f) * SSD02;

	cat(mm_space(19), "x", mm_space(30), "t\n");

	cat("         n    ", mm_space(15), sprintf("%-0 10f", n), "\n");
	cat("         S    ", sprintf("%-0 10f", Sx), mm_space(20), 
		sprintf("%-0 10f", St), "\n");
	cat("       USS    ", sprintf("%-0 10.4f", USSx), mm_space(20), 
		sprintf("%-0 10f", USSt), "\n")
	cat("        SP    ", mm_space(15), sprintf("%-0 10f", SPxt), "\n");
	cat("       SSD    ", sprintf("%-0 10f", SSDx), mm_space(20), 
		sprintf("%-0 10f", SSDt), "\n");
	cat("       SPD    ", mm_space(15), sprintf("%-0 10f", SPDxt), "\n");
	cat("      bhat    ", mm_space(15), sprintf("%-0 10f", bhat), "\n");
	cat("      ahat    ", mm_space(15), sprintf("%-0 10f", ahat), "\n");
	cat("     SSD02    ", mm_space(15), sprintf("%-0 10f", SSD02), "\n");
	cat("s02_square    ", mm_space(15), sprintf("%-0 10f", s02_square), 
		"\n\n");

	# Confidence interval for Alpha
	alphaconfval <- sqrt(s02_square * (1/n + (St/n)^2/SSDt)) * qt(0.975, n-2);
	cat("C95(alpha)   = [", ahat - alphaconfval, ",", ahat + alphaconfval, "]\n")

	betaconfval <- sqrt(s02_square/SSDt) * qt(0.975, n-2);
	cat("C95(beta)    = [", bhat - betaconfval, ",", bhat+betaconfval, "]\n")

	# Confidence interval for Sigma^2
	cat("C95(sigma^2) = [", s02_square/(qchisq(0.975, f)/f), ",", s02_square/(qchisq(0.025, f)/f), "]\n\n")

	if(!is.null(slope)) {
		cat("M2->M3\n");
		tofx <- (bhat - slope)/sqrt(s02_square/SSDt);
		cat("\tt(x) ", tofx, "~~ t(", f, ")\n");
		pobstx <- 2*(1 - pt(abs(tofx), f));
		cat("\tpobs(x): ", pobstx, "\n");

		if(pobstx > 0.05) {
			s03_square <- (1/(n-1))*(SSD02 + (bhat - slope)^2*SSDt);
			ahatM3 = Sx/n - slope*(St/n)
			cat("\ts03_square:", s03_square, "\n\tahatm3:", ahatM3, "\n");

			muconfval <- sqrt(s03_square * (1/n + (St/n)^2/SSDt)) * qt(0.975, n-1);
			cat("\tC95(alpha) = [", ahatM3 - muconfval, ",", ahatM3 + muconfval, "]\n\n")

			if(!is.null(intercept)) {
				cat("M3->M4\n");
				tofx <- (ahat - intercept)/sqrt(s03_square/n);
				cat("\tt(x) ", tofx, "~~t(", n-1, ")\n");
				pobstx <- 2*(1 - pt(abs(tofx), n - 1));
				cat("\tpobs(x): ", pobstx, "\n");

				if(pobstx > 0.05) {
					s04_square <- (1/n)*(USSx + n*intercept^2 + slope^2*USSt -
										 2*intercept*Sx - 2*slope*SPxt +
										 2*intercept*slope*St);
					cat("\ts04_square:", s04_square, "\n\n");
				} else {
					cat("\n")
				}
			}
		} else {
			cat("\n")
		}
	}

	if(!is.null(intercept)) {
		cat("M2->M3*\n");
		tval <- (ahat - intercept)/sqrt(s02_square*(1/n+((St/n)^2)/SSDt))
		cat("\tt(x) ", tval, "~~t(", n-2,")\n")
		pobs_int <- 2*(1-pt(abs(tval), n-2));
		cat("\tpobs(x): ", pobs_int, "\n")

		if(pobs_int > 0.05) {
			bhatM3 <- (SPxt - intercept*St)/USSt # ~~ N(beta, sigma^2/USSt)
			s03_2star <- (1/(n-1))*(USSx + n*intercept^2 - 2*intercept*Sx -
									bhatM3^2*USSt)
			cat("\tbhatM3*:", bhatM3, "\n")
			cat("\ts03*^2:", s03_2star, "\n")

			sigsqrconfval <- sqrt(s03_2star / SSDt) * qt(0.975, n-1);
			cat("\tsigsqrconf:", sigsqrconfval ,"\n")
			cat("\tC95(beta)  = [", bhatM3 - sigsqrconfval, ",", bhatM3 +
				sigsqrconfval, "]\n\n")

			if(!is.null(slope)) {
				cat("M3*->M4\n");
				tval <- (bhatM3 - slope)/sqrt(s03_2star/USSt)
				cat("\tt(x) ", tval, "~~t(", n-1,")\n")
				pobs_slope <- 2*(1-pt(abs(tval), n-1))
				cat("\tpobs(x): ", pobs_slope, "\n")

				if(pobs_slope > 0.05) {
					s04_square <- (1/n)*(USSx + n*intercept^2 + slope^2*USSt -
										 2*intercept*Sx - 2*slope*SPxt +
										 2*intercept*slope*St);
					cat("\ts04_square:", s04_square, "\n\n");
				} else {
					cat("\n")
				}
			}
		} else {
			cat("\n")
		}
	}

	return(data.frame(n, Sx, St, USSx, USSt, SPxt, SSDx, SSDt, SPDxt, bhat, 
					  ahat, SSD02, s02_square, f));
	#}}}
}

mm_linregAddRow <- function(a, b) {
	a <- rbind(a, b);
}

mm_linreg2samples <- function(r) {
	#{{{
	# Do a test for common variance, end in M1
	max <- 0;
	min <- 0;
	s01_sqr <- 0;

	if(r$s02_square[1] > r$s02_square[2]) {
		max <- 1;
		min <- 2;
	} else {
		max <- 2;
		min <- 1;
	}

	F <- r$s02_square[max]/r$s02_square[min];
	pobs_SameVariance <- 2*(1 - pf(F, r$f[max], r$f[min]));

	cat("------ Common Variance ------\n");
	cat("F:", F, "\n")
	cat("pobs_SameVariance:", pobs_SameVariance, "\n");
	if(pobs_SameVariance > 0.05) {
		s01_sqr <- (r$f[min]*r$s02_square[min] + r$f[max]*r$s02_square[max])/(r$f[min] + r$f[max])
		cat("We can assume that the variance is the same between the two linear regressions\n");
		cat("New variance estimate for M1:", s01_sqr, "\n\n");
	} else {
		cat("It is not likely that the linear regressions have the same variance. EOF\n");
		return();
	}

	# Do a test for common Alpha, and end in M2
	F <- 1

	# Do a test for common Beta, and end in M2*

	# Do a test for common Beta M2 -> M3

	# Do a test for common Alpha M2* -> M3
	#}}}
}

mm_space <- function(i) {
	return(paste(rep(" ", i), collapse=''));
}

mm_ftest <- function(v1, v2, f1, f2) {
	#{{{
	maxv <- 0
	maxf <- 0

	minv <- 0
	minf <- 0

	if(v1 > v2) {
		maxv <- v1
		maxf <- f1
		minv <- v2
		minf <- f2
	} else {
		maxv <- v2
		maxf <- f2
		minv <- v1
		minf <- f1
	}

	F <- maxv/minv

	cat("pobs:", 2*(1 - pf(F, maxf, minf)), ", F:", F, "\n")
	#}}}
}

mm_FromTotest <- function(SSDto, SSDfrom, fto, ffrom, varfrom) {
	#{{{
	F <- (SSDto - SSDfrom)/(fto - ffrom)/(varfrom);
	cat("\t  ", mm_format(SSDto), " - ", mm_format(SSDfrom),"\n")
	cat("\t", "  ----------------------------\n")
	cat("\t  ", mm_format(fto), " - ", mm_format(ffrom))
	cat(mm_space(8), "=  ", F, "~~ F(", fto-ffrom, ",", ffrom, ")\n")
	cat("\t", "--------------------------------\n")
	cat("\t", mm_space(11), varfrom, "\n")

	cat("\n\npobs(x): ", 1-pf(F, fto-ffrom, fto), "\n")
	#}}}
}

mm_format <- function(x) {
	return(format(x, width=10, justify="left"))
}

mm_binomConfidence <- function(x, n) {
	#{{{
	low <- (1/(n + 1.96^2))*((x - 0.5) + 1.96^2/2 - 1.96*sqrt(((x - 0.5)*(n - x + 0.5)) / n + 1.96^2/4) )
	high <- (1/(n + 1.96^2))*((x + 0.5) + 1.96^2/2 + 1.96*sqrt(((x + 0.5)*(n - x - 0.5)) / n + 1.96^2/4) )

	cat("C95(pi) (continuity approximation): [", low, ",", high, "]\n");

	low <- (1/(n + 1.96^2)) * (x + 1.96^2/2 - 1.96*sqrt((x*(n-x))/n + 1.96^2/4))
	high <- (1/(n + 1.96^2)) * (x + 1.96^2/2 + 1.96*sqrt((x*(n-x))/n + 1.96^2/4))

	cat("C95(pi): [", low, ",", high, "]\n");
	#}}}
}

mm_multinomial <- function(x) {
	b <- 3;

	dims <- dim(x);

	# Calculate expected values
	e <- matrix(nrow=dims[1], ncol=dims[2])
	for(i in 1:dims[1]) {
		for(j in 1:dims[2]) {
			e[i,j] = (sum(x[i,]) * sum(x[,j]))/sum(x)
		}
	}

	width1 <- max(nchar(unlist(dimnames(x)[1])));

	rsum   <- c() ## Calculate the total sum of rows
	for(i in 1:length(x[,1])) {
		rsum[i] <- sum(x[i,])
	}

	csum <- c()
	cWidth <- c()
	for(i in 1:length(x[1,])) {
		csum[i] <- sum(x[,i])
		cWidth[i] <- max(max(nchar(e[,i])), max(nchar(x[,i])))
	}

	print(cWidth)

	rowName <- unlist(dimnames(x)[1])
	colName <- unlist(dimnames(x)[2])

	# Draw the top row
	cat(sep='', rep(" ", width1 + b), "|")
	for(i in 1:dims[2]) {
			cat(sep='', format(colName[i], width=cWidth[j]+b, justify="centre"), "|")
	}
	cat("\n")

	cat(paste(rep("-", width1 + b), collapse=''), "+", sep='')
	for(j in 1:dims[2]) {
		cat( paste( rep("-", cWidth[j]+b), collapse=''), sep='', "+")
	}
	cat("\n")

	for(i in 1:dims[1]) {
		cat(format(rowName[i], width=width1+2, justify="right"), "|")
		for(j in 1:dims[2]) {
			cat(sep='', format(x[i,j], width=cWidth[j]+b), "|")
		}
		cat(sep='', "   ", format(rsum[i], justify="left"))
		cat("\n")

		cat(sep='', rep(" ", width1 + b), "|")
		for(j in 1:dims[2]) {
			cat(sep='', format(e[i,j], width=cWidth[j]+b), "|")
		}
		cat("\n")


		cat(paste(rep("-", width1 + b), collapse=''), "+", sep='')
		for(j in 1:dims[2]) {
			cat( paste( rep("-", cWidth[j]+b), collapse=''), sep='', "+")
		}

		cat("\n")
	}

	cat(sep='', rep(" ", width1 + b+1))
	for(i in 1:dims[2]) {
			cat(sep='', format(sum(x[,i]), width=cWidth[j]+b, justify="centre"), " ")
	}
		cat(sep='', "   ", format(sum(x), justify="left"))
	cat("\n\n")

	if(min(e) <= 5) {
		cat("WARNING, THERE IS AN EXPECTED VALUE LESS THAN 5\n\n")
	}

	lnq <- 2*( sum(x * log(x/e) ))
	f <- (dims[1] - 1)*(dims[2] - 1)
	cat("Performing Likelihood Ratio Chi-Square test:\n")
	cat(sep='', "-2lnQ(x): ", lnq, " ~~ Chisq(", f, ")\n")
	cat(sep='', "pobs: ", 1 - pchisq(lnq, f), "\n")
}

#x <- c(2.31, 2.33, 2.34, 2.41, 2.44, 2.45, 2.46, 2.47, 2.51, 2.54);
#x <- c(2.2,3.6,3.8,4.4,4.9,5.9,6.9,7.5,7.5,7.7,8.9,9.7);
#keys <- c(1215, 1181, 1195, 1189, 1185, 1206, 1192, 1203, 1213, 1156, 1164, 1177, 1167, 1206, 1193)

# Opgave på ugeseddel 5:
# x <- c(19.8, 22.8, 24.5, 27.3, 31.0, 35, 35.1, 37.1, 38.5, 39)
# t <- c(20, 22.5, 25, 28.5, 31, 33.5, 35.5, 37, 38, 40)
