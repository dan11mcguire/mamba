#loadtmb<-function(dir=getwd()){
#  compile(paste0(dir, "/tau2f.cpp"))
#  dyn.load(dynlib(paste0(dir, "/tau2f")))
#  data("res2ref_cpd")
#
#}


#M<-50000; alpha<-10
#lambda<-0.95; tau2<-2.5e-4; p<-0.01; resref=NULL;n=rep(50*10^3,10)
#generate_data_mamba(M,p,tau2,lambda,alpha,n=rep(50*10^3,10))

#' Simulate summary statistics according to the MAMBA model.
#' @param M number of SNPs for which to simulate summary stats.
#' @param p proportion of SNPs with real non-zero effects.
#' @param tau2 variance of real associated non-zero SNPs.
#' @param lambda proportion of NON-outlier studies for non-replicable SNPs. 
#' @param n vector of sample sizes from the contributing studies. 
#' @param alpha variance inflation factor of outlier studies.
#' @param resref an optional vector of residual variances to sample from (with replacement) 
#'   when generating the standard errors for the summary stats.  Default is NULL, in which case the residual variances calculated from a Cigarretes Per Day GWAS (Liu, Jiang 2019)  are used.
#' @return 
#'  A list containing 
#'  An Mxk matrix of effect size estimates \code{betajk}, 
#'  An Mxk matrix of effect size estimate  variances \code{sjk2}, 
#'  M-length vector inverse-variance weighted meta-analysis z-scores \code{meta.z}, 
#'  an M-length binary vector indicating real / non-real effect at each SNP \code{Rj},
#'  an M-length binary vector indicating true effect-size at each SNP \code{muj}.  
#'  Binary matrix \code{Ojk} indicating whether study was normal(Ojk=1) or outlier (Ojk=0).
#' @examples
#' generate_data_mamba(M=100, p=0.01, tau2=2.5e-4, n=rep(10^4, 10), lambda=0.975, alpha=15) 
#' @export


generate_data_mamba<-function(M=50*10^3, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10), 
			     lambda=0.975, alpha=5){

  if(missing(resref)){
	 data(res2ref_cpd, package="mamba",envir=environment())
	 resref<-sqrt(res2ref_cpd$res2);
   } 
  Rj<-rbinom(M, size=1, prob=p)
  muj<-Rj*rnorm(M, 0, sqrt(tau2))
  k<-length(n)
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  Ojk<-lapply(1:M, function(j){rbinom(k, size=1, prob=lambda)})
  betajk<-lapply(1:M, function(j){
	Rj[j] * rnorm(k, muj[j], sd=sqrt(sjk2[[j]])) + 
      (1-Rj[j]) * ((Ojk[[j]])*rnorm(k, 0, sqrt(sjk2[[j]])) + 
	(1-Ojk[[j]])*rnorm(k, 0, sqrt(alpha)*sqrt(sjk2[[j]])))})
  
  meta.z<-sapply(1:length(betajk), function(j){
    sum(betajk[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })

  betajk<-matrix(unlist(betajk), ncol=k, byrow=TRUE)
  sjk2<-matrix(unlist(sjk2), ncol=k, byrow=TRUE)
  Ojk<-matrix(unlist(Ojk), ncol=k, byrow=TRUE)
 
  return(list(betajk=betajk,
              sjk2=sjk2,
              meta.z=meta.z, 
              Rj=Rj, 
              Ojk=Ojk,
              muj=muj))
}

#' Simulate summary statistics according to the (inverse-variance weighted) fixed effects model. 
#' @param M number of SNPs for which to simulate summary stats.
#' @param p proportion of SNPs with real non-zero effects.
#' @param tau2 variance of real associated non-zero SNPs.
#' @param n vector of sample sizes from the contributing studies. 
#' @param resref an optional vector of residual variances to sample from (with replacement) 
#'   when generating the standard errors for the summary stats.  Default is NULL, in which case the residual variances calculated from a Cigarretes Per Day GWAS (Liu, Jiang 2019)  are used.
#' @return
#'  A list containing:
#'  An Mxk matrix of effect size estimates \code{betajk}, 
#'  An Mxk matrix of effect size estimate  variances \code{sjk2}, 
#'  M-length vector inverse-variance weighted meta-analysis z-scores \code{meta.z}, 
#'  an M-length binary vector indicating real / non-real effect at each SNP \code{Rj},
#'  an M-length binary vector indicating true effect-size at each SNP \code{muj}.  
#' @examples
#' generate_data_fe(M=100, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10)) 
#' @export


generate_data_fe<-function(M=50*10^3, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10)){
  
  if(missing(resref)){
	 data(res2ref_cpd, package="mamba",envir=environment())
	 resref<-sqrt(res2ref_cpd$res2);
   } 
  Rj<-rbinom(M, size=1, prob=p)
  muj<-Rj*rnorm(M, 0, sqrt(tau2))
  k<-length(n)
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  betajk<-lapply(1:M, function(j){
    rnorm(k, muj[j], sd=sqrt(sjk2[[j]]))
  })
  
  meta.z<-sapply(1:length(betajk), function(j){
    sum(betajk[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })

  betajk<-matrix(unlist(betajk), ncol=k, byrow=TRUE)
  sjk2<-matrix(unlist(sjk2), ncol=k, byrow=TRUE)
  
  return(list(betajk=betajk,
              sjk2=sjk2,
              meta.z=meta.z, 
              Rj=Rj, 
              muj=muj))
}

#' Simulate summary statistics according to the binary-effects (be) model.
#' @param M number of SNPs for which to simulate summary stats.
#' @param p proportion of SNPs with real non-zero effects.
#' @param tau2 variance of real associated non-zero SNPs.
#' @param n vector of sample sizes from the contributing studies. 
#' @param resref an optional vector of residual variances to sample from (with replacement) 
#'   when generating the standard errors for the summary stats.  Default is NULL, in which case the residual variances calculated from a Cigarretes Per Day GWAS (Liu, Jiang 2019)  are used.
#' @return
#'  A list containing:
#'  \itemize{
#'    \item An Mxk matrix of effect size estimates \code{betajk}, 
#'    \item An Mxk matrix of effect size estimate  variances \code{sjk2}, 
#'     \item M-length vector inverse-variance weighted meta-analysis z-scores \code{meta.z}, 
#'     \item an M-length binary vector indicating real / non-real effect at each SNP \code{Rj},
#'  \item an M-length binary vector indicating true effect-size at each SNP \code{muj}.  
#'  \item an M-length vector indicating the number of studies at each SNP with non-zero effects .
#'  \item an Mxk binary matrix indicating which studies at which SNP had non-zero effects.
#' }
#' @examples
#' generate_data_mamba(M=100, p=0.01, tau2=2.5e-4, n=rep(10^4, 10), lambda=0.975, alpha=15) 
#' @export

generate_data_be<-function(M=50*10^3, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10)){
  if(missing(resref)){
	 data(res2ref_cpd, package="mamba",envir=environment())
	 resref<-sqrt(res2ref_cpd$res2);
   } 
  Rj<-rbinom(M, size=1, prob=p)
  muj<-Rj*rnorm(M, 0, sqrt(tau2))

  k<-length(n)
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  Nej<-sapply(1:M, function(j){ Rj[j]*sample(1:k, 1) })
  bej<-lapply(1:M, function(j){
    be<-rep(0, k)
    if(Nej[j] > 0){
      be[sample(1:k, Nej[j])]<-1
    }
     return(be)
  })
  
  betajk<-lapply(1:M, function(j){bej[[j]] * rnorm(k, muj[j], sd=sqrt(sjk2[[j]])) + 
      (1-bej[[j]]) *rnorm(k, 0, sqrt(sjk2[[j]]))
  })
  
  meta.z<-sapply(1:length(betajk), function(j){
    sum(betajk[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })

  betajk<-matrix(unlist(betajk), ncol=k, byrow=TRUE)
  sjk2<-matrix(unlist(sjk2), ncol=k, byrow=TRUE)
  bej<-matrix(unlist(bej), ncol=k, byrow=TRUE)

  return(list(betajk=betajk,
              sjk2=sjk2,
              meta.z=meta.z, 
              Rj=Rj, 
              muj=muj, 
              Nej=Nej,
              bej=bej))
}

#' Simulate summary statistics according to random effects model.
#' @param M number of SNPs for which to simulate summary stats.
#' @param p proportion of SNPs with real non-zero effects.
#' @param tau2 variance of real associated non-zero SNPs.
#' @param n vector of sample sizes from the contributing studies. 
#' @param alpha variance inflation factor of outlier studies.
#' @param resref an optional vector of residual variances to sample from (with replacement) 
#'   when generating the standard errors for the summary stats.  Default is NULL, in which case the residual variances calculated from a Cigarretes Per Day GWAS (Liu, Jiang 2019)  are used.
#' @param I2 the I2 heterogeneity statistic for each SNP. 
#'  The variance of study-level effects around population level effect at each SNP 
#'  is specified given I2 level (between 0,1) and the simulated standard errors.   
#' @return
#'  A list containing:
#'  An Mxk matrix of effect size estimates \code{betajk}, 
#'  An Mxk matrix of effect size estimate  variances \code{sjk2}, 
#'  M-length vector inverse-variance weighted meta-analysis z-scores \code{meta.z}, 
#'  an M-length binary vector indicating real / non-real effect at each SNP \code{Rj},
#'  an M-length binary vector indicating true effect-size at each SNP \code{muj}.  
#'  an Mxk matrix of the true study-level effects \code{etajk}
#'  an M-length vector of the variance of the study-level effects around the SNP's population level effect \code{omega2j}
#' @examples
#'   generate_data_re1(M=100, p=0.01, tau2=2.5e-4, n=rep(10^4, 10), I2=0.2) 
#' @export

generate_data_re1<-function(M=50*10^3, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10),I2=0.3){
  
  if(missing(resref)){
	 data(res2ref_cpd, package="mamba",envir=environment())
	 resref<-sqrt(res2ref_cpd$res2);
   } 
  Rj<-rbinom(M, size=1, prob=p)
  muj<-Rj*rnorm(M, 0, sqrt(tau2))
  k<-length(n)
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  S2j<-sapply(1:M, function(j){
    ((k-1)*sum(1/sjk2[[j]]))/(sum(1/sjk2[[j]])^2 - sum(1/sjk2[[j]]^2))
  })
  omega2j<-sapply(1:M, function(j){
    I2/(1-I2) * S2j[j]
  })
  etajk<-lapply(1:M, function(j){
    rnorm(k, muj[j], sd=sqrt(omega2j[j])) 
  })
  
  betajk<-lapply(1:M, function(j){
    rnorm(k, etajk[[j]], sd=sqrt(sjk2[[j]]))
  })
  
  meta.z<-sapply(1:length(betajk), function(j){
    sum(betajk[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })


  betajk<-matrix(unlist(betajk), ncol=k, byrow=TRUE)
  sjk2<-matrix(unlist(sjk2), ncol=k, byrow=TRUE)
  etajk<-matrix(unlist(etajk), ncol=k, byrow=TRUE)

  return(list(betajk=betajk,
              sjk2=sjk2,
              meta.z=meta.z, 
              etajk=etajk,
              omega2j=omega2j,
              S2j=S2j,
              Rj=Rj, 
              muj=muj))
}

#' Simulate summary statistics according to the Han and Eskin's RE2 model.  (Heterogeneity only when the population mean effect is non-zero)
#' @param M number of SNPs for which to simulate summary stats.
#' @param p proportion of SNPs with real non-zero effects.
#' @param tau2 variance of real associated non-zero SNPs.
#' @param lambda proportion of NON-outlier studies for non-replicable SNPs. 
#' @param n vector of sample sizes from the contributing studies. 
#' @param alpha variance inflation factor of outlier studies.
#' @param resref an optional vector of residual variances to sample from (with replacement) 
#'   when generating the standard errors for the summary stats.  Default is NULL, in which case the residual variances calculated from a Cigarretes Per Day GWAS (Liu, Jiang 2019)  are used.
#' @return 
#'  A list containing:
#'  An Mxk matrix of effect size estimates \code{betajk}, 
#'  An Mxk matrix of effect size estimate  variances \code{sjk2}, 
#'  M-length vector inverse-variance weighted meta-analysis z-scores \code{meta.z}, 
#'  an M-length binary vector indicating real / non-real effect at each SNP \code{Rj},
#'  an M-length binary vector indicating true effect-size at each SNP \code{muj}.  
#'  an Mxk matrix of the true study-level effects \code{etajk}
#'  an M-length vector of the variance of the study-level effects around the SNP's population level effect \code{omega2j}
#' @examples
#' generate_data_re2(M=100, p=0.01, tau2=2.5e-4, n=rep(10^4, 10), I2=0.2) 
#' @export

generate_data_re2<-function(M=50*10^3, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10),I2=0.3){
  
  if(missing(resref)){
	 data(res2ref_cpd, package="mamba",envir=environment())
	 resref<-sqrt(res2ref_cpd$res2);
   } 
  Rj<-rbinom(M, size=1, prob=p)
  muj<-Rj*rnorm(M, 0, sqrt(tau2))
  k<-length(n)
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  S2j<-sapply(1:M, function(j){
    ((k-1)*sum(1/sjk2[[j]]))/(sum(1/sjk2[[j]])^2 - sum(1/sjk2[[j]]^2))
  })
  omega2j<-sapply(1:M, function(j){
    I2/(1-I2) * S2j[j]
  })
  etajk<-lapply(1:M, function(j){
   Rj[j] * rnorm(k, muj[j], sd=sqrt(omega2j[j])) 
  })
  
  betajk<-lapply(1:M, function(j){
    rnorm(k, etajk[[j]], sd=sqrt(sjk2[[j]]))
  })
  
  
  meta.z<-sapply(1:length(betajk), function(j){
    sum(betajk[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })


  betajk<-matrix(unlist(betajk), ncol=k, byrow=TRUE)
  sjk2<-matrix(unlist(sjk2), ncol=k, byrow=TRUE)
  etajk<-matrix(unlist(etajk), ncol=k, byrow=TRUE)
  
  return(list(betajk=betajk,
              sjk2=sjk2,
              meta.z=meta.z, 
              etajk=etajk,
              omega2j=omega2j,
              S2j=S2j,
              Rj=Rj, 
              muj=muj))
}

#d<-generate_data_mamba()
#betajk<-matrix(unlist(d$betajk),
#		byrow=TRUE,
#		nrow=length(d$betajk))
#sjk2<-matrix(unlist(d$sjk2),
#		byrow=TRUE,
#		nrow=length(d$betajk))
#parcores=1; 
#                   p=0.003;
#                   f=0.96;
#                   tau2=0.0002;
#                   alpha=3;
#                   #conv.eps=1e-3;
#                   rel.eps=1e-8;
#                   verbose=1;
#                   snpids=NA;
#                   maxIter=10^4L


#betajk<-betaij
#sjk2<-sij2
#
#parcores=1; 
#p=0.003;
#lambda=0.96;
#tau2=0.0002;
#alpha=3;
##conv.eps=1e-3;
#rel.eps=1e-8;
#verbose=1;
#snpids=NA;
#maxIter=10^4L


#' Fit the MAMBA model.
#' @param betajk Mxk matrix of effect size estimates, where row j corresponds to j-th SNP, and column k corresponds to k-th study .  Missing values can be represented with NA.
#' @param sjk2 Mxk matrix of effect size estimate variances, where row j corresponds to SNP j, and column k corresponds to study k.  Missing values can be represented with NA.
#' @param p initial value for EM algorithm, for the proportion of non-zero SNPs
#' @param lambda initial value for EM algorithm, for the proportion of non-replicable SNPs which are well behaved, or non-outliers.
#' @param tau2 initial value for EM algorithm, for the variance of replicable non-zero effect SNPs.
#' @param alpha initial value for EM algorithm, for the variance inflation of outlier summary statistics at non-replicable SNPs.
#' @param rel.eps threshold for when to end EM algorithm. rel.eps = (ll[i] - ll[i-1])/ll[i], where ll[i] is the log-likelihood at iteration i.
#' @param snpids an optional vector of SNP id names.  If not provided, the ID's will be 1:M, corresponding to the order of the matrix betajk.
#' @param maxIter maximum # of EM iterations.
#' @return
#' @examples 
#'   d<-generate_data_mamba()
#'   mod<-mamba(betajk=d$betajk, sjk2=d$sjk2)
#' @export 


mamba<-function(betajk, sjk2, 
                   parcores=1, 
                   p=0.003,
                   lambda=0.96,
                   tau2=0.0002,
                   alpha=3,
                   #conv.eps=1e-3,
                   rel.eps=1e-8,
                   verbose=1,
                   snpids=NA,
                   maxIter=10^4L){

  betajk<-as.matrix(betajk)
  sjk2<-as.matrix(sjk2)

## Replace any summary stats with standard error = 0 with missing
  zeroL<-mclapply(1:nrow(sjk2), function(j){
    which(sjk2[j,]==0)
  }, mc.cores = parcores)
 chk<- which(sapply(zeroL, length) > 0)
  if(length(chk) > 0){
    for(j in chk){
      sjk2[j,zeroL[[j]]]<-NA
      betajk[j,zeroL[[j]]]<-NA
    }
  }

## Replace any summary stats with standard error = inf or beta=inf with missing
 infL<-mclapply(1:nrow(sjk2), function(j){
   which(is.infinite(sjk2[j,]))
 }, mc.cores = parcores)
 chk<- which(sapply(infL, length) > 0)
 if(length(chk) > 0){
   for(i in chk){
     sjk2[j,infL[[j]]]<-NA
     betajk[j,infL[[j]]]<-NA
   }
 }
 infL<-mclapply(1:nrow(betajk), function(j){
   which(is.infinite(betajk[j,]))
 }, mc.cores = parcores)
 chk<- which(sapply(infL, length) > 0)
 if(length(chk) > 0){
   for(i in chk){
     sjk2[j,infL[[j]]]<-NA
     betajk[j,infL[[j]]]<-NA
   }
 }

### Identify indices of studies with non-missing summary stats at each SNP
  mis.inds<-mclapply(1:nrow(sjk2), function(j){
    which(!is.na(sjk2[j,]))
  }, mc.cores = parcores)
  
  
  b2s2<-sapply(1:nrow(betajk), function(j){
    sum(betajk[j,mis.inds[[j]]]^2/sjk2[j,mis.inds[[j]]])
  })
  bs22<-sapply(1:nrow(betajk), function(j){
    sum(betajk[j,mis.inds[[j]]]/sjk2[j,mis.inds[[j]]])^2
  })
  os22<-sapply(1:nrow(betajk), function(j){
    sum(1/sjk2[j,mis.inds[[j]]])
  })
  
  
  ll<-NA
  
  k_j<-sapply(mis.inds, length)
  
  MM<-nrow(betajk)
  
  mu.hat<-sapply(1:MM, function(j){
    sum(betajk[j,mis.inds[[j]]]/sjk2[j,mis.inds[[j]]])/(sum(1/sjk2[j,mis.inds[[j]]]) + 1/(tau2))
  })
  
  strt<-Sys.time()
  for(iter in 1:(maxIter)){
    
    deltajk<- mclapply(1:nrow(betajk), function(j){
	mamba:::deltis(betajk_j=betajk[j,], sjk2_j=sjk2[j,], alpha=alpha, lambda=lambda)
    })
    deltajk<-matrix(unlist(deltajk), ncol=MM, byrow=FALSE)
    
    missingdelta<-unlist(mclapply(1:MM, function(j){ 
			mean(is.na(deltajk[mis.inds[[j]],j]))
	      }, mc.cores = parcores))
    if(sum(missingdelta) > 0) {
      print("na's in deltajk")
      break
    }
    if(verbose > 0){
      print("deltajk calculated.")
    }
    if(iter==1){
      llbR1<-unlist(mclapply(1:nrow(betajk), function(j){
        mamba:::llbR1_j(betajk_j = betajk[j,mis.inds[[j]]], 
		sjk2_j=sjk2[j,mis.inds[[j]]], 
		tau2inv = 1/tau2)
      }, mc.cores=parcores)) - (k_j/2)*log(2*pi)
      
      llbR0<-unlist( mclapply(1:nrow(betajk), function(j){
        mamba:::llbR0_j(betajk_j = betajk[j,mis.inds[[j]]], 
		sjk2_j = sjk2[j,mis.inds[[j]]], 
		alpha = alpha, 
		lambda=lambda)
      }, mc.cores = parcores))
    }
    
    gammaj<-unlist(mclapply(1:MM, function(j){
      gam<- p*exp(llbR1[j])/(exp(mamba:::logsumexp(c(log(p) + llbR1[j], log(1-p) + llbR0[j]))))
      if(is.na(gam)){
        m<-which.max(c(log(1-p) + llbR0[j],
                       log(p) + llbR1[j]))                      
        gam<-(m-1)
      }
      gam
    }, mc.cores = parcores))
    if(verbose > 1){
      print(data.table(gammaj)[order(-gammaj)])
    }
    
    print("e step finished.")
    
    if(sum(is.na(gammaj)) > 0) break
    if(sum(1-gammaj)==0) break
    fn<-sum(unlist(mclapply(1:MM, function(j){
      (1-gammaj[j])*sum(deltajk[mis.inds[[j]],j])
    }, mc.cores = parcores)))
    fd<-sum((1-gammaj)*k_j)
    lambda<-fn/fd
    
    alpha<-sum((1-gammaj)*
                 unlist(mclapply(1:MM, function(j){
                   sum(
   (1-deltajk[mis.inds[[j]],j])*(betajk[j,mis.inds[[j]]]^2 / sjk2[j,mis.inds[[j]]])
   )
    }, mc.cores = parcores)))/
      sum((1-gammaj)*unlist(mclapply(1:MM, function(j){
        sum(1-deltajk[mis.inds[[j]],j])}, mc.cores = parcores)))
    
    nllk <- MakeADFun ( data = list ( b2s2=b2s2, bs22=bs22, os22=os22, gammaj=gammaj) , 
                        parameters = list ( logTau2=log(tau2)) , 
			DLL="tau2f",silent=TRUE)
    fit <- nlminb ( start = nllk $par , objective = nllk $fn , gradient = nllk $gr ,
                    lower =c(- Inf ,0) , upper =c( Inf , Inf ), control = list(trace=0))
    
    tau2<-exp(fit$par)
    tau2<-max(tau2, 1e-17)
    
    if(is.na(tau2)) break 
    
    mu.hat<-unlist(mclapply(1:MM, function(j){
      sum(betajk[j,mis.inds[[j]]]/sjk2[j,mis.inds[[j]]])/(sum(1/sjk2[j,mis.inds[[j]]]) + 1/tau2)
    }, mc.cores = parcores))
    
    p<-sum(gammaj)/MM
    if(verbose > 0){
      print("m step finished.")
      print(paste0("p=", round(p, 3), 
		   " lambda=", round(lambda, 3), 
		   " tau2=", round(tau2, 8), 
		   " alpha=", round(alpha, 3)))
    }
    
    llbR1<-unlist(mclapply(1:nrow(betajk), function(j){
	mamba:::llbR1_j(betajk_j = betajk[j,mis.inds[[j]]], 
	      sjk2_j=sjk2[j,mis.inds[[j]]], 
	      tau2inv = 1/tau2)
    }, mc.cores=parcores)) - (k_j/2)*log(2*pi)
    
    llbR0<-unlist( mclapply(1:nrow(betajk), function(j){
	    mamba:::llbR0_j(betajk_j = betajk[j,mis.inds[[j]]], 
	      sjk2_j = sjk2[j,mis.inds[[j]]], 
	      alpha = alpha, 
	      lambda=lambda)
    }, mc.cores = parcores))
    
    
    ll[iter]<-sum(unlist(mclapply(1:MM, function(j){
	mamba:::logsumexp(c(log(p) + llbR1[j], log(1-p) + llbR0[j]))
    }, mc.cores = parcores)))
    
    if(iter > 1){
      llMat<-rbindlist(list(llMat,
                            data.table(p=p,lambda=lambda,tau2=tau2,alpha=alpha,ll=ll[iter])))
    } else{
      llMat<-data.table(p=p,lambda=lambda,tau2=tau2,alpha=alpha,ll=ll[iter])
    } 
   
   # convergence criteria # 
    if(iter >=2){
      if(verbose > 0){
        print((ll[iter]-ll[iter-1])/ll[iter])
        print((ll[iter]-ll[iter-1]))
      }
      if(!is.na(ll[iter] - ll[iter-1])){
        if(((ll[iter] - ll[iter-1])/ll[iter] < rel.eps && ll[iter] - ll[iter-1] < 0.01) || 
		((ll[iter]-ll[iter-1]) < 0.001 && iter > 100) || 
		((ll[iter] - ll[iter-1])/ll[iter] < rel.eps && iter > 100)) break
      }
      #if(ll[iter] - ll[iter-1] < 0) break
    }
    if(verbose > 0){
      print(paste0("iteration ", iter, " complete."))
    } 
    #print(llMat[1:iter])
  }
  end<-Sys.time()
  se.hat<-unlist(mclapply(1:nrow(betajk), function(j){
    sqrt(1/(sum(1/sjk2[j,mis.inds[[j]]]) + 1/tau2))
  }, mc.cores = parcores))
  outliermat<-data.table::data.table(snp=1:ncol(deltajk),
                         t(deltajk))
  outliermat<-outliermat[,lapply(.SD, function(x) 1-x),.SDcols=paste0("V", 1:nrow(deltajk)),by=snp]
  colnames(outliermat)[2:ncol(outliermat)]<-paste0("Oj_",1:dim(deltajk)[1])
  outliermat[,ppr:=gammaj]
  #outlierprobs<-outliermat[,#(paste0("outlier_prob", 1:dim(deltajk)[1])):=
  #                         lapply(.SD, function(x){(1-gammaj)*(1-x)}),
  #                         .SDcols=paste0("deltaj_",1:dim(deltajk)[1]), by=.(snp)]
  #colnames(outlierprobs)[2:ncol(outlierprobs)]<-paste0("outlier_prob", 1:dim(deltajk)[1])
  # if(!is.null(colnames(betajk))){
  #   dnames<-gsub("beta.*\\_","",colnames(betajk))
  #   
  # }
  if(!is.na(snpids) && length(snpids)==nrow(outliermat)) {
    outliermat[,snp:=snpids]
  #  outlierprobs[,snp:=snpids]
  }
  return(list(ll=ll,
              p=p,
              ppr=gammaj,
              alpha=alpha,
              lambda=lambda,
              tau2=tau2, 
              mu.hat=mu.hat, 
              se.hat=se.hat,
              post.means=mu.hat*gammaj,
              outliermat=outliermat,
              #outlierprobs=outlierprobs, 
              time=end-strt,
              param_est_log=llMat[,`:=`(eps=ll-shift(ll,1),
                                rel.eps=((ll-shift(ll,1))/ll))]))
}

getpvals<-function(mod, 
		   s2, 
		   nModels, 
		   nullSNPsPerModel, 
		   numcores=parcores,
		   save_all=FALSE,
                   out){
  
  pdat<-nullMod<-list()
  for(j in 1:nModels){
    pdat[[j]]<-generate_data_S2(model=mod, sjk2=s2, Mnull=nullSNPsPerModel)
    nullMod[[j]]<-em_std_f(pdat[[j]]$betajk, pdat[[j]]$sjk2_sample, 
                           p=mod$p, lambda=mod$lambda, 
                           alpha=mod$alpha, tau2=mod$tau2,
                           parcores=numcores,verbose = TRUE)
    nullscoresj<-nullMod[[j]]$gammaj
    nullRij<-pdat[[j]]$Ri
    nullscoresj<-nullscoresj[nullRij==0] 
    fwrite(data.table(nullscoresj),
           file=paste0(out, "nullscores.",seed,".txt"),
           col.names = FALSE)
   
    system(paste0("cat ",
		  paste0(out, "nullscores.",seed,".txt")," >> ",
		  paste0(out, "nullscores.txt")))
    system(paste0(" echo $(wc -l ", 
		  paste0(out, "nullscores.txt) seed ",seed, " >> ",
			 paste0(out, "nullscores.log"))))
    
    system(paste0("rm ",paste0(out, "nullscores.",seed,".txt")))
    
    print(paste0("Model ",j, " of ", nModels, " complete.")) 
  }
  nullscores<-unlist(lapply(nullMod, "[[", "gammaj"))
  nullRi<-unlist(lapply(pdat, "[[", "Ri"))
  nullscores<-nullscores[nullRi==0]
  
  pvals<-sapply(mod$gammaj, function(score){
    mean(score < nullscores)
  })
  if(save_all){
    return(list(pvals=pvals,
                nullscores=nullscores,
                nullRi=nullRi,
                pdat=pdat,
                nullMod=nullMod))
  } else{
    
    return(list(nullscores=nullscores))
  }
}


generate_data_S2<-function(model, sjk2, Mnull=1000){
  
  zeroL<-mclapply(1:nrow(sjk2), function(j){
    which(sjk2[j,]==0)
  }, mc.cores = parcores)
  chk<- which(sapply(zeroL, length) > 0)
  if(length(chk) > 0){
    for(i in chk){
      sjk2[j,zeroL[[j]]]<-NA
    }
  }
  infL<-mclapply(1:nrow(sjk2), function(j){
    which(is.infinite(sjk2[j,]))
  }, mc.cores = parcores)
  chk<- which(sapply(infL, length) > 0)
  if(length(chk) > 0){
    for(i in chk){
      sjk2[j,infL[[j]]]<-NA
    }
  }
  p<-model$p
  f<-model$f
  tau2<-model$tau2
  alpha<-model$alpha
  M<-ceiling(Mnull/(1-p))
  sample_inds<-sample(nrow(sjk2),M,replace = TRUE)
  sjk2_sample<-sjk2[sample_inds,]
  mis.inds<-lapply(1:nrow(sjk2_sample), function(j){
    which(!is.na(sjk2_sample[j,]))
  })
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(tau2))
  k_j<-sapply(1:nrow(sjk2_sample), function(j){sum(!is.na(sjk2_sample[j,]))})
  sjk2<-lapply(1:M, function(j){sjk2_sample[j,]})
  Ojk<-lapply(1:M, function(j){
    vec<-rep(NA, dim(sjk2)[2])
    vec[mis.inds[[j]]]<-rbinom(k_j[j], size=1, prob=f)
    vec
  })
  
  betaj<-lapply(1:M, function(j){
    vec<-rep(NA, dim(sjk2)[2])
    vec[mis.inds[[j]]]<-
      Ri[j] * rnorm(k_j[j], mu[j], sd=sqrt(sjk2[[j]][mis.inds[[j]]])) + 
      (1-Ri[j]) * ((Ojk[[j]][mis.inds[[j]]])*rnorm(k_j[j], 0, sqrt(sjk2[[j]][mis.inds[[j]]])) + 
                     (1-Ojk[[j]][mis.inds[[j]]])*rnorm(k_j[j], 0, sqrt(alpha)*sqrt(sjk2[[j]][mis.inds[[j]]])))
    
    vec
  })
  
  meta.z<-sapply(1:length(betaj), function(j){
    sum(betaj[[j]][mis.inds[[j]]]/sjk2[[j]][mis.inds[[j]]]) / sqrt(sum(1/sjk2[[j]][mis.inds[[j]]]))
  })
  betajk<-matrix(unlist(betaj), ncol=dim(sjk2)[2], byrow=TRUE)
  return(list(betajk=betajk,
              sjk2_sample=sjk2_sample,
              sample_inds=sample_inds,
              meta.z=meta.z, 
              Rj=Rj, 
              Ojk=Ojk,
              muj=muj))
}
