loadtmb<-function(dir=getwd()){
  compile(paste0(dir, "/tau2f.cpp"))
  dyn.load(dynlib(paste0(dir, "/tau2f")))
  data("res2ref_cpd")

}


#M<-50000; alpha<-10
#lambda<-0.95; tau2<-2.5e-4; p<-0.01; resref=NULL;n=rep(50*10^3,10)
#generate_data_mamba(M,p,tau2,lambda,alpha,n=rep(50*10^3,10))
generateData_mamba<-function(M=50*10^3, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10), 
			     lambda=0.975, alpha=5){

  if(missing(resref)){
	 data(res2ref_cpd, package="emfuncs",envir=environment())
	 resref<-res2ref_cpd[,sqrt(res2)];
   } 
  Rj<-rbinom(M, size=1, prob=p)
  muj<-Rj*rnorm(M, 0, tau2)
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
  return(list(betajk=betajk,
              sjk2=sjk2,
              meta.z=meta.z, 
              Rj=Rj, 
              Ojk=Ojk,
              muj=muj))
}

generateData_fe<-function(M=50*10^3, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10)){
  
  if(missing(resref)){
	 data(res2ref_cpd, package="emfuncs",envir=environment())
	 resref<-res2ref_cpd[,sqrt(res2)];
   } 
  Rj<-rbinom(M, size=1, prob=p)
  muj<-Rj*rnorm(M, 0, tau2)
  k<-length(n)
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  Ojk<-lapply(1:M, function(j){rbinom(k, size=1, prob=1)})
  betajk<-lapply(1:M, function(j){
    rnorm(k, muj[j], sd=sqrt(sjk2[[j]]))
  })
  
  meta.z<-sapply(1:length(betajk), function(j){
    sum(betajk[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })
  return(list(betajk=betajk,
              sjk2=sjk2,
              meta.z=meta.z, 
              Rj=Rj, 
              Ojk=Ojk,
              muj=muj))
}


generateData_be<-function(M=50*10^3, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10)){
  if(missing(resref)){
	 data(res2ref_cpd, package="emfuncs",envir=environment())
	 resref<-res2ref_cpd[,sqrt(res2)];
   } 
  Rj<-rbinom(M, size=1, prob=p)
  muj<-Rj*rnorm(M, 0, tau2)
  k<-length(n)
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  Ojk<-lapply(1:M, function(j){rbinom(k, size=1, prob=1)})
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
  return(list(betajk=betajk,
              sjk2=sjk2,
              meta.z=meta.z, 
              Rj=Rj, 
              Ojk=Ojk,
              muj=muj, 
              Nej=Nej,
              bej=bej))
}


generateData_re1<-function(M=50*10^3, p=0.01, tau2=2.5e-4, resref=NULL, n=rep(50*10^3,10),I2=0.3){
  
  if(missing(resref)){
	 data(res2ref_cpd, package="emfuncs",envir=environment())
	 resref<-res2ref_cpd[,sqrt(res2)];
   } 
  Rj<-rbinom(M, size=1, prob=p)
  muj<-Rj*rnorm(M, 0, tau2)
  k<-length(n)
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  Ojk<-lapply(1:M, function(j){rbinom(k, size=1, prob=1)})
  S2j<-sapply(1:M, function(j){
    ((k-1)*sum(1/sjk2[[j]]))/(sum(1/sjk2[[j]])^2 - sum(1/sjk2[[j]]^2))
  })
  omega2j<-sapply(1:M, function(j){
    I2/(1-I2) * S2j[j]
  })
  etajk<-lapply(1:M, function(j){
    rnorm(k, muj[j], sd=sqrt(omega2j[j])) 
  })
  
  betaj<-lapply(1:M, function(j){
    rnorm(k, etajk[[j]], sd=sqrt(sjk2[[j]]))
  })
  
  meta.z<-sapply(1:length(betaj), function(j){
    sum(betaj[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })
  return(list(betajk=betajk,
              sjk2=sjk2,
              meta.z=meta.z, 
              etajk=etajk,
              omega2j=omega2j,
              S2j=S2j,
              Rj=Rj, 
              Ojk=Ojk,
              muj=muj))
}

### only heterogeneity under the alternative
generateData_re2<-function(M, p, h2g, f, alpha, I2=0.3,resref=res2ref_cpd[,sqrt(res2)], 
                                                     n=c(1192, 6613, 1829, 6181, 9607, 270820, 73331, 39480, 5459)){
  
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(h2g/(M*p)))
  k<-length(n)
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  Ojk<-lapply(1:M, function(j){rbinom(k, size=1, prob=1)})
  S2j<-sapply(1:M, function(j){
    ((k-1)*sum(1/sjk2[[j]]))/(sum(1/sjk2[[j]])^2 - sum(1/sjk2[[j]]^2))
  })
  omega2i<-sapply(1:M, function(j){
    I2/(1-I2) * S2j[j]
  })
  etajk<-lapply(1:M, function(j){
   Ri[j] * rnorm(k, mu[j], sd=sqrt(omega2i[j])) 
  })
  
  betaj<-lapply(1:M, function(j){
    rnorm(k, etajk[[j]], sd=sqrt(sjk2[[j]]))
  })
  
  
  meta.z<-sapply(1:length(betaj), function(j){
    sum(betaj[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })
  return(list(betaj=betaj,
              sjk2=sjk2,
              meta.z=meta.z, 
              etajk=etajk,
              omega2i=omega2i,
              S2j=S2j,
              Ri=Ri, 
              Ojk=Ojk,
              mui=mu))
}


generateData_t_re<-function(M, p, h2g, I2=0.3, resref=res2ref_cpd[,sqrt(res2)], f, 
                                              n=c(1192, 6613, 1829, 6181, 9607, 270820, 73331, 39480, 5459)){
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(h2g/(M*p)))
  k<-length(n)
  #maf<-runif(M, .05, 0.49)
  tau2<-h2g/(M*p)
  
  #ui<-rt(200, df=3)
  #ui<-var(ui/sqrt(3/var_ui_target))
  #nu<-re.p*tau2/(re.p*tau2 + re.p - 1 )
  #res.sd<-lapply(1:M, function(j){rlnorm(k, meanlog=-0.007891429, sdlog= 0.886734822)})
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  Ojk<-lapply(1:M, function(j){rbinom(k, size=1, prob=1)})
  
  S2j<-sapply(1:M, function(j){
    ((k-1)*sum(1/sjk2[[j]]))/(sum(1/sjk2[[j]])^2 - sum(1/sjk2[[j]]^2))
  })
  omega2i<-sapply(1:M, function(j){
    I2/(1-I2) * S2j[j]
  })
  etajk<-lapply(1:M, function(j){
    Ri[j] * rt(k, ncp=mu[j], df=2*(1-omega2i[j])/omega2i[j]) 
  })
  
  betaj<-lapply(1:M, function(j){
    rnorm(k, mean=etajk[[j]], sd=sqrt(sjk2[[j]]))
  })
  
  meta.z<-sapply(1:length(betaj), function(j){
    sum(betaj[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })
  return(list(betaj=betaj,
              sjk2=sjk2,
              meta.z=meta.z, 
              etajk=etajk,
              omega2i=omega2i,
              S2j=S2j,
              Ri=Ri, 
              Ojk=Ojk,
              mui=mu))
}


generateData_t_noisy<-function(M, p, h2g, re.p=0.3, resref=res2ref_cpd[,sqrt(res2)], f, 
                                              n=c(1192, 6613, 1829, 6181, 9607, 270820, 73331, 39480, 5459)){
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(h2g/(M*p)))
  k<-length(n)
  #maf<-runif(M, .05, 0.49)
  tau2<-h2g/(M*p)
  
  #ui<-rt(200, df=3)
  #ui<-var(ui/sqrt(3/var_ui_target))
  #nu<-re.p*tau2/(re.p*tau2 + re.p - 1 )
  #res.sd<-lapply(1:M, function(j){rlnorm(k, meanlog=-0.007891429, sdlog= 0.886734822)})
  res.sd<-lapply(1:M, function(j){ sample(resref, k, replace = TRUE) })
  sjk2<-lapply(1:M, function(j){res.sd[[j]]^2 / (n)})
  uij<-lapply(1:M, function(j){
    var_ui_target<-sjk2[[j]]*re.p/(1-re.p) 
    rt(k, df=3) / sqrt(3/var_ui_target)
  })
  Ojk<-lapply(1:M, function(j){rbinom(k, size=1, prob=1)})
  
  betaj<-lapply(1:M, function(j){
    etajk<-mu[j] + uij[[j]]
    rnorm(k, mean=etajk, sd=sqrt(sjk2[[j]]))
  })
  
  meta.z<-sapply(1:length(betaj), function(j){
    sum(betaj[[j]]/sjk2[[j]]) / sqrt(sum(1/sjk2[[j]]))
  })
  return(list(betaj=betaj,
              sjk2=sjk2,
              uij=uij,
              meta.z=meta.z, 
              Ojk=Ojk,
              Ri=Ri, 
              mui=mu))
}

em_std_f<-function(betajk, sij2, 
                   parcores=1, 
                   p=0.003,
                   f=0.96,
                   tau2=0.0002,
                   alpha=3,
                   #conv.eps=1e-3,
                   rel.eps=1e-8,
                   verbose=1,
                   snpids=NA,
                   maxIter=10^4L){
  betajk<-as.matrix(betajk)
  sij2<-as.matrix(sij2)
  
  zeroL<-mclapply(1:nrow(sij2), function(j){
    which(sij2[i,]==0)
  }, mc.cores = parcores)
 chk<- which(sapply(zeroL, length) > 0)
  if(length(chk) > 0){
    for(i in chk){
      sij2[i,zeroL[[j]]]<-NA
      betajk[i,zeroL[[j]]]<-NA
    }
  }
 infL<-mclapply(1:nrow(sij2), function(j){
   which(is.infinite(sij2[i,]))
 }, mc.cores = parcores)
 chk<- which(sapply(infL, length) > 0)
 if(length(chk) > 0){
   for(i in chk){
     sij2[i,infL[[j]]]<-NA
     betajk[i,infL[[j]]]<-NA
   }
 }
 infL<-mclapply(1:nrow(betajk), function(j){
   which(is.infinite(betajk[i,]))
 }, mc.cores = parcores)
 chk<- which(sapply(infL, length) > 0)
 if(length(chk) > 0){
   for(i in chk){
     sij2[i,infL[[j]]]<-NA
     betajk[i,infL[[j]]]<-NA
   }
 }
  mis.inds<-mclapply(1:nrow(sij2), function(j){
    which(!is.na(sij2[i,]))
  }, mc.cores = parcores)
  
  
  b2s2<-sapply(1:nrow(betajk), function(j){
    sum(betajk[i,mis.inds[[j]]]^2/sij2[i,mis.inds[[j]]])
  })
  bs22<-sapply(1:nrow(betajk), function(j){
    sum(betajk[i,mis.inds[[j]]]/sij2[i,mis.inds[[j]]])^2
  })
  os22<-sapply(1:nrow(betajk), function(j){
    sum(1/sij2[i,mis.inds[[j]]])
  })
  
  # tau2f2<-function(param, b2s2, bs22, os22, gammai=gammai){
  #   tau2<-exp(param)
  #   -1*  sum(  sapply(1:nrow(betajk), function(j){
  #     gammai[j]* (.5*log(1/tau2) -.5*log(os22[j] + 1/tau2) -
  #                   .5*(b2s2[j] - 
  #                         bs22[j]/(os22[j] + 1/tau2)))
  #   }))
  # }
  # 
  
  ll<-NA
  
  k_i<-sapply(mis.inds, length)
  
  MM<-nrow(betajk)
  
  mu.hat<-sapply(1:MM, function(j){
    sum(betajk[i,mis.inds[[j]]]/sij2[i,mis.inds[[j]]])/(sum(1/sij2[i,mis.inds[[j]]]) + 1/(tau2))
  })
  
  strt<-Sys.time()
  for(iter in 1:(maxIter)){
    
    deltaij<- mclapply(1:nrow(betajk), function(j){
      deltis(betajk_i=betajk[i,], sij2_i=sij2[i,], alpha=alpha, f=f)
    })
    deltaij<-matrix(unlist(deltaij), ncol=MM, byrow=FALSE)
    
    missingdelta<-unlist(mclapply(1:MM, function(j){ mean(is.na(deltaij[mis.inds[[j]],j]))}, mc.cores = parcores))
    if(sum(missingdelta) > 0) {
      print("na's in deltaij")
      break
    }
    if(verbose > 0){
      print("deltaij calculated.")
    }
    if(iter==1){
      llbR1<-unlist(mclapply(1:nrow(betajk), function(j){
        llbR1_i(betajk_i = betajk[i,mis.inds[[j]]], sij2_i=sij2[i,mis.inds[[j]]], tau2inv = 1/tau2)
      }, mc.cores=parcores)) - (k_i/2)*log(2*pi)
      
      llbR0<-unlist( mclapply(1:nrow(betajk), function(j){
        llbR0_i(betajk_i = betajk[i,mis.inds[[j]]], sij2_i = sij2[i,mis.inds[[j]]], alpha = alpha, f=f)
      }, mc.cores = parcores))
    }
    
    gammai<-unlist(mclapply(1:MM, function(j){
      gam<- p*exp(llbR1[j])/(exp(logsumexp(c(log(p) + llbR1[j], log(1-p) + llbR0[j]))))
      if(is.na(gam)){
        m<-which.max(c(log(1-p) + llbR0[j],
                       log(p) + llbR1[j]))                      
        gam<-(m-1)
      }
      gam
    }, mc.cores = parcores))
    if(verbose > 1){
      print(data.table(gammai)[order(-gammai)])
    }
    
    print("e step finished.")
    
    if(sum(is.na(gammai)) > 0) break
    if(sum(1-gammai)==0) break
    fn<-sum(unlist(mclapply(1:MM, function(j){
      (1-gammai[j])*sum(deltaij[mis.inds[[j]],i])
    }, mc.cores = parcores)))
    fd<-sum((1-gammai)*k_i)
    f<-fn/fd
    #f<-sum((1-gammai)*(colSums(deltaij, na.rm=TRUE))) / (sum((1-gammai)*k_i))   ## later, replace na.rm with index to sum
    
    alpha<-sum((1-gammai)*
                 unlist(mclapply(1:MM, function(j){
                   sum((1-deltaij[mis.inds[[j]],i])*(betajk[i,mis.inds[[j]]]^2 / sij2[i,mis.inds[[j]]]))}, mc.cores = parcores)))/
      sum((1-gammai)*unlist(mclapply(1:MM, function(j){
        sum(1-deltaij[mis.inds[[j]],i])}, mc.cores = parcores)))
    
    # system.time({
    # tau2Fit<-nlminb(start=log(tau2), objective = tau2f, betajk=betajk, sij2=sij2,  gammai=gammai, control=list(trace=2))
    
    nllk <- MakeADFun ( data = list ( b2s2=b2s2, bs22=bs22, os22=os22, gammai=gammai) , 
                        parameters = list ( logTau2=log(tau2)) , DLL="tau2f",silent=TRUE)
    fit <- nlminb ( start = nllk $par , objective = nllk $fn , gradient = nllk $gr ,
                    lower =c(- Inf ,0) , upper =c( Inf , Inf ), control = list(trace=0))
    tau2<-exp(fit$par)
    # tau2Fit2<-nlminb(start=log(tau2), objective = tau2f2, 
    #                 b2s2=b2s2, bs22=bs22, os22=os22, gammai=gammai, 
    #                 control=list(trace=2))
    # tau2<-exp(tau2Fit2$par)
    tau2<-max(tau2, 1e-17)
    
    if(is.na(tau2)) break 
    
    mu.hat<-unlist(mclapply(1:MM, function(j){
      sum(betajk[i,mis.inds[[j]]]/sij2[i,mis.inds[[j]]])/(sum(1/sij2[i,mis.inds[[j]]]) + 1/tau2)
    }, mc.cores = parcores))
    
    p<-sum(gammai)/MM
    if(verbose > 0){
      print("m step finished.")
      print(paste0("p=", round(p, 3), " f=", round(f, 3), " tau2=", round(tau2, 8), " alpha=", round(alpha, 3)))
    }
    
    llbR1<-unlist(mclapply(1:nrow(betajk), function(j){
      llbR1_i(betajk_i = betajk[i,mis.inds[[j]]], sij2_i=sij2[i,mis.inds[[j]]], tau2inv = 1/tau2)
    }, mc.cores=parcores)) - (k_i/2)*log(2*pi)
    
    llbR0<-unlist( mclapply(1:nrow(betajk), function(j){
      llbR0_i(betajk_i = betajk[i,mis.inds[[j]]], sij2_i = sij2[i,mis.inds[[j]]], alpha = alpha, f=f)
    }, mc.cores = parcores))
    
    
    ll[iter]<-sum(unlist(mclapply(1:MM, function(j){
      logsumexp(c(log(p) + llbR1[j], log(1-p) + llbR0[j]))
    }, mc.cores = parcores)))
    
    if(iter > 1){
      llMat<-rbindlist(list(llMat,
                            data.table(p=p,f=f,tau2=tau2,alpha=alpha,ll=ll[iter])))
    } else{
      llMat<-data.table(p=p,f=f,tau2=tau2,alpha=alpha,ll=ll[iter])
    } 
    
    if(iter >=2){
      if(verbose > 0){
        print((ll[iter]-ll[iter-1])/ll[iter])
        print((ll[iter]-ll[iter-1]))
      }
      if(!is.na(ll[iter] - ll[iter-1])){
        #if((ll[iter] - ll[iter-1]) < conv.eps) break
        if(((ll[iter] - ll[iter-1])/ll[iter] < rel.eps && ll[iter] - ll[iter-1] < 0.01) || ((ll[iter]-ll[iter-1]) < 0.001 && iter > 100) || ((ll[iter] - ll[iter-1])/ll[iter] < rel.eps && iter > 100)) break
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
    sqrt(1/(sum(1/sij2[i,mis.inds[[j]]]) + 1/tau2))
  }, mc.cores = parcores))
  outliermat<-data.table(snp=1:ncol(deltaij),
                         t(deltaij))
  colnames(outliermat)[2:ncol(outliermat)]<-paste0("deltai_",1:dim(deltaij)[1])
  outliermat[,gammai:=gammai]
  outlierprobs<-outliermat[,#(paste0("outlier_prob", 1:dim(deltaij)[1])):=
                           lapply(.SD, function(x){(1-gammai)*(1-x)}),
                           .SDcols=paste0("deltai_",1:dim(deltaij)[1]), by=.(snp)]
  colnames(outlierprobs)[2:ncol(outlierprobs)]<-paste0("outlier_prob", 1:dim(deltaij)[1])
  # if(!is.null(colnames(betajk))){
  #   dnames<-gsub("beta.*\\_","",colnames(betajk))
  #   
  # }
  if(!is.na(snpids) && length(snpids)==nrow(outliermat)) {
    outliermat[,snp:=snpids]
    outlierprobs[,snp:=snpids]
  }
  return(list(ll=ll,
              p=p,
              gammai=gammai,
              alpha=alpha,
              f=f,
              tau2=tau2, 
              mu.hat=mu.hat, 
              se.hat=se.hat,
              post.means=mu.hat*gammai,
              outliermat=outliermat,
              outlierprobs=outlierprobs, 
              time=end-strt,
              llMat=llMat[,`:=`(eps=ll-shift(ll,1),
                                rel.eps=(ll-shift(ll,1)/ll))]))
}

getpvals<-function(mod, s2, nModels, nullSNPsPerModel, numcores=parcores,save_all=FALSE,
                   out){
  
  pdat<-nullMod<-list()
  for(j in 1:nModels){
    pdat[[j]]<-generateDataS2(model=mod, sij2=s2, Mnull=nullSNPsPerModel)
    nullMod[[j]]<-em_std_f(pdat[[j]]$betajk, pdat[[j]]$sij2_sample, 
                           p=mod$p, f=mod$f, 
                           alpha=mod$alpha, tau2=mod$tau2,
                           parcores=numcores,verbose = TRUE)
    nullscoresj<-nullMod[[j]]$gammai
    nullRij<-pdat[[j]]$Ri
    nullscoresj<-nullscoresj[nullRij==0] 
    fwrite(data.table(nullscoresj),
           file=paste0(out, "nullscores.",seed,".txt"),
           col.names = FALSE)
    system(paste0("cat ",paste0(out, "nullscores.",seed,".txt")," >> ",paste0(out, "nullscores.txt")))
    system(paste0(" echo $(wc -l ", paste0(out, "nullscores.txt) seed ",seed, " >> ",paste0(out, "nullscores.log"))))
    system(paste0("rm ",paste0(out, "nullscores.",seed,".txt")))
    print(paste0("Model ",j, " of ", nModels, " complete.")) 
  }
  nullscores<-unlist(lapply(nullMod, "[[", "gammai"))
  nullRi<-unlist(lapply(pdat, "[[", "Ri"))
  nullscores<-nullscores[nullRi==0]
  
  pvals<-sapply(mod$gammai, function(score){
    mean(score < nullscores)
  })
  if(save_all){
    return(list(pvals=pvals,
                nullscores=nullscores,
                nullRi=nullRi,
                pdat=pdat,
                nullMod=nullMod))
  } else{
    
    return(list(nullscores=nullscores
                ))
  }
}


generateDataS2<-function(model, sij2, Mnull=1000){
  
  zeroL<-mclapply(1:nrow(sij2), function(j){
    which(sij2[i,]==0)
  }, mc.cores = parcores)
  chk<- which(sapply(zeroL, length) > 0)
  if(length(chk) > 0){
    for(i in chk){
      sij2[i,zeroL[[j]]]<-NA
    }
  }
  infL<-mclapply(1:nrow(sij2), function(j){
    which(is.infinite(sij2[i,]))
  }, mc.cores = parcores)
  chk<- which(sapply(infL, length) > 0)
  if(length(chk) > 0){
    for(i in chk){
      sij2[i,infL[[j]]]<-NA
    }
  }
  p<-model$p
  f<-model$f
  tau2<-model$tau2
  alpha<-model$alpha
  M<-ceiling(Mnull/(1-p))
  sample_inds<-sample(nrow(sij2),M,replace = TRUE)
  sij2_sample<-sij2[sample_inds,]
  mis.inds<-lapply(1:nrow(sij2_sample), function(j){
    which(!is.na(sij2_sample[i,]))
  })
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(tau2))
  k_i<-sapply(1:nrow(sij2_sample), function(j){sum(!is.na(sij2_sample[i,]))})
  sjk2<-lapply(1:M, function(j){sij2_sample[i,]})
  Ojk<-lapply(1:M, function(j){
    vec<-rep(NA, dim(sij2)[2])
    vec[mis.inds[[j]]]<-rbinom(k_i[j], size=1, prob=f)
    vec
  })
  
  betaj<-lapply(1:M, function(j){
    vec<-rep(NA, dim(sij2)[2])
    vec[mis.inds[[j]]]<-
      Ri[j] * rnorm(k_i[j], mu[j], sd=sqrt(sjk2[[j]][mis.inds[[j]]])) + 
      (1-Ri[j]) * ((Ojk[[j]][mis.inds[[j]]])*rnorm(k_i[j], 0, sqrt(sjk2[[j]][mis.inds[[j]]])) + 
                     (1-Ojk[[j]][mis.inds[[j]]])*rnorm(k_i[j], 0, sqrt(alpha)*sqrt(sjk2[[j]][mis.inds[[j]]])))
    
    vec
  })
  
  meta.z<-sapply(1:length(betaj), function(j){
    sum(betaj[[j]][mis.inds[[j]]]/sjk2[[j]][mis.inds[[j]]]) / sqrt(sum(1/sjk2[[j]][mis.inds[[j]]]))
  })
  betajk<-matrix(unlist(betaj), ncol=dim(sij2)[2], byrow=TRUE)
  return(list(betajk=betajk,
              sij2_sample=sij2_sample,
              sample_inds=sample_inds,
              meta.z=meta.z, 
              Ri=Ri, 
              Ojk=Ojk,
              mui=mu))
}
