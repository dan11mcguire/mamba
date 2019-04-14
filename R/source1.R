loadtmb<-function(dir=getwd()){
  compile(paste0(dir, "/tau2f.cpp"))
  dyn.load(dynlib(paste0(dir, "/tau2f")))
  data("res2ref_cpd")

}

generateData_standardized<-function(M, p, h2g, f, alpha, resref=res2ref_cpd[,sqrt(res2)], 
                                    n=c(1192, 6613, 1829, 6181, 9607, 270820, 73331, 39480, 5459)){
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(h2g/(M*p)))
  k<-length(n)
  #maf<-runif(M, .05, 0.49)
  #res.sd<-lapply(1:M, function(i){rlnorm(k, meanlog=-0.007891429, sdlog= 0.886734822)})
  res.sd<-lapply(1:M, function(i){ sample(resref, k, replace = TRUE) })
  sj2<-lapply(1:M, function(i){res.sd[[i]]^2 / (n)})
  f.ind<-lapply(1:M, function(i){rbinom(k, size=1, prob=f)})
  betaj<-lapply(1:M, function(i){Ri[i] * rnorm(k, mu[i], sd=sqrt(sj2[[i]])) + 
      (1-Ri[i]) * ((f.ind[[i]])*rnorm(k, 0, sqrt(sj2[[i]])) + (1-f.ind[[i]])*rnorm(k, 0, sqrt(alpha)*sqrt(sj2[[i]])))})
  
  meta.z<-sapply(1:length(betaj), function(i){
    sum(betaj[[i]]/sj2[[i]]) / sqrt(sum(1/sj2[[i]]))
  })
  return(list(betaj=betaj,
              sj2=sj2,
              meta.z=meta.z, 
              Ri=Ri, 
              f.ind=f.ind,
              mui=mu))
}

generateData_fe<-function(M, p, h2g, f, alpha, resref=res2ref_cpd[,sqrt(res2)], 
                                    n=c(1192, 6613, 1829, 6181, 9607, 270820, 73331, 39480, 5459)){
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(h2g/(M*p)))
  k<-length(n)
 
  res.sd<-lapply(1:M, function(i){ sample(resref, k, replace = TRUE) })
  sj2<-lapply(1:M, function(i){res.sd[[i]]^2 / (n)})
  f.ind<-lapply(1:M, function(i){rbinom(k, size=1, prob=1)})
  betaj<-lapply(1:M, function(i){
    rnorm(k, mu[i], sd=sqrt(sj2[[i]]))
    })
  
  meta.z<-sapply(1:length(betaj), function(i){
    sum(betaj[[i]]/sj2[[i]]) / sqrt(sum(1/sj2[[i]]))
  })
  return(list(betaj=betaj,
              sj2=sj2,
              meta.z=meta.z, 
              Ri=Ri, 
              f.ind=f.ind,
              mui=mu))
}


generateData_be<-function(M, p, h2g, f, alpha, resref=res2ref_cpd[,sqrt(res2)],
                                    n=c(1192, 6613, 1829, 6181, 9607, 270820, 73331, 39480, 5459)){
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(h2g/(M*p)))
  k<-length(n)
  res.sd<-lapply(1:M, function(i){ sample(resref, k, replace = TRUE) })
  sj2<-lapply(1:M, function(i){res.sd[[i]]^2 / (n)})
  f.ind<-lapply(1:M, function(i){rbinom(k, size=1, prob=1)})
  Nei<-sapply(1:M, function(i){ Ri[i]*sample(1:k, 1) })
  bei<-lapply(1:M, function(i){
    be<-rep(0, k)
    if(Nei[i] > 0){
      be[sample(1:k, Nei[i])]<-1
    }
     return(be)
  })
  
  betaj<-lapply(1:M, function(i){bei[[i]] * rnorm(k, mu[i], sd=sqrt(sj2[[i]])) + 
      (1-bei[[i]]) *rnorm(k, 0, sqrt(sj2[[i]]))
  })
  
  meta.z<-sapply(1:length(betaj), function(i){
    sum(betaj[[i]]/sj2[[i]]) / sqrt(sum(1/sj2[[i]]))
  })
  return(list(betaj=betaj,
              sj2=sj2,
              meta.z=meta.z, 
              Ri=Ri, 
              f.ind=f.ind,
              mui=mu, 
              Nei=Nei,
              bei=bei))
}


generateData_re1<-function(M, p, h2g, f, alpha, I2=0.3,resref=res2ref_cpd[,sqrt(res2)], 
                           n=c(1192, 6613, 1829, 6181, 9607, 270820, 73331, 39480, 5459)){
  
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(h2g/(M*p)))
  k<-length(n)
  res.sd<-lapply(1:M, function(i){ sample(resref, k, replace = TRUE) })
  sj2<-lapply(1:M, function(i){res.sd[[i]]^2 / (n)})
  f.ind<-lapply(1:M, function(i){rbinom(k, size=1, prob=1)})
  S2i<-sapply(1:M, function(i){
    ((k-1)*sum(1/sj2[[i]]))/(sum(1/sj2[[i]])^2 - sum(1/sj2[[i]]^2))
  })
  omega2i<-sapply(1:M, function(i){
    I2/(1-I2) * S2i[i]
  })
  etaij<-lapply(1:M, function(i){
    rnorm(k, mu[i], sd=sqrt(omega2i[i])) 
  })
  
  betaj<-lapply(1:M, function(i){
    rnorm(k, etaij[[i]], sd=sqrt(sj2[[i]]))
  })
  
  meta.z<-sapply(1:length(betaj), function(i){
    sum(betaj[[i]]/sj2[[i]]) / sqrt(sum(1/sj2[[i]]))
  })
  return(list(betaj=betaj,
              sj2=sj2,
              meta.z=meta.z, 
              etaij=etaij,
              omega2i=omega2i,
              S2i=S2i,
              Ri=Ri, 
              f.ind=f.ind,
              mui=mu))
}

### only heterogeneity under the alternative
generateData_re2<-function(M, p, h2g, f, alpha, I2=0.3,resref=res2ref_cpd[,sqrt(res2)], 
                                                     n=c(1192, 6613, 1829, 6181, 9607, 270820, 73331, 39480, 5459)){
  
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(h2g/(M*p)))
  k<-length(n)
  res.sd<-lapply(1:M, function(i){ sample(resref, k, replace = TRUE) })
  sj2<-lapply(1:M, function(i){res.sd[[i]]^2 / (n)})
  f.ind<-lapply(1:M, function(i){rbinom(k, size=1, prob=1)})
  S2i<-sapply(1:M, function(i){
    ((k-1)*sum(1/sj2[[i]]))/(sum(1/sj2[[i]])^2 - sum(1/sj2[[i]]^2))
  })
  omega2i<-sapply(1:M, function(i){
    I2/(1-I2) * S2i[i]
  })
  etaij<-lapply(1:M, function(i){
   Ri[i] * rnorm(k, mu[i], sd=sqrt(omega2i[i])) 
  })
  
  betaj<-lapply(1:M, function(i){
    rnorm(k, etaij[[i]], sd=sqrt(sj2[[i]]))
  })
  
  
  meta.z<-sapply(1:length(betaj), function(i){
    sum(betaj[[i]]/sj2[[i]]) / sqrt(sum(1/sj2[[i]]))
  })
  return(list(betaj=betaj,
              sj2=sj2,
              meta.z=meta.z, 
              etaij=etaij,
              omega2i=omega2i,
              S2i=S2i,
              Ri=Ri, 
              f.ind=f.ind,
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
  #res.sd<-lapply(1:M, function(i){rlnorm(k, meanlog=-0.007891429, sdlog= 0.886734822)})
  res.sd<-lapply(1:M, function(i){ sample(resref, k, replace = TRUE) })
  sj2<-lapply(1:M, function(i){res.sd[[i]]^2 / (n)})
  f.ind<-lapply(1:M, function(i){rbinom(k, size=1, prob=1)})
  
  S2i<-sapply(1:M, function(i){
    ((k-1)*sum(1/sj2[[i]]))/(sum(1/sj2[[i]])^2 - sum(1/sj2[[i]]^2))
  })
  omega2i<-sapply(1:M, function(i){
    I2/(1-I2) * S2i[i]
  })
  etaij<-lapply(1:M, function(i){
    Ri[i] * rt(k, ncp=mu[i], df=2*(1-omega2i[i])/omega2i[i]) 
  })
  
  betaj<-lapply(1:M, function(i){
    rnorm(k, mean=etaij[[i]], sd=sqrt(sj2[[i]]))
  })
  
  meta.z<-sapply(1:length(betaj), function(i){
    sum(betaj[[i]]/sj2[[i]]) / sqrt(sum(1/sj2[[i]]))
  })
  return(list(betaj=betaj,
              sj2=sj2,
              meta.z=meta.z, 
              etaij=etaij,
              omega2i=omega2i,
              S2i=S2i,
              Ri=Ri, 
              f.ind=f.ind,
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
  #res.sd<-lapply(1:M, function(i){rlnorm(k, meanlog=-0.007891429, sdlog= 0.886734822)})
  res.sd<-lapply(1:M, function(i){ sample(resref, k, replace = TRUE) })
  sj2<-lapply(1:M, function(i){res.sd[[i]]^2 / (n)})
  uij<-lapply(1:M, function(i){
    var_ui_target<-sj2[[i]]*re.p/(1-re.p) 
    rt(k, df=3) / sqrt(3/var_ui_target)
  })
  f.ind<-lapply(1:M, function(i){rbinom(k, size=1, prob=1)})
  
  betaj<-lapply(1:M, function(i){
    etaij<-mu[i] + uij[[i]]
    rnorm(k, mean=etaij, sd=sqrt(sj2[[i]]))
  })
  
  meta.z<-sapply(1:length(betaj), function(i){
    sum(betaj[[i]]/sj2[[i]]) / sqrt(sum(1/sj2[[i]]))
  })
  return(list(betaj=betaj,
              sj2=sj2,
              uij=uij,
              meta.z=meta.z, 
              f.ind=f.ind,
              Ri=Ri, 
              mui=mu))
}

em_std_f<-function(betaij, sij2, 
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
  betaij<-as.matrix(betaij)
  sij2<-as.matrix(sij2)
  
  zeroL<-mclapply(1:nrow(sij2), function(i){
    which(sij2[i,]==0)
  }, mc.cores = parcores)
 chk<- which(sapply(zeroL, length) > 0)
  if(length(chk) > 0){
    for(i in chk){
      sij2[i,zeroL[[i]]]<-NA
      betaij[i,zeroL[[i]]]<-NA
    }
  }
 infL<-mclapply(1:nrow(sij2), function(i){
   which(is.infinite(sij2[i,]))
 }, mc.cores = parcores)
 chk<- which(sapply(infL, length) > 0)
 if(length(chk) > 0){
   for(i in chk){
     sij2[i,infL[[i]]]<-NA
     betaij[i,infL[[i]]]<-NA
   }
 }
 infL<-mclapply(1:nrow(betaij), function(i){
   which(is.infinite(betaij[i,]))
 }, mc.cores = parcores)
 chk<- which(sapply(infL, length) > 0)
 if(length(chk) > 0){
   for(i in chk){
     sij2[i,infL[[i]]]<-NA
     betaij[i,infL[[i]]]<-NA
   }
 }
  mis.inds<-mclapply(1:nrow(sij2), function(i){
    which(!is.na(sij2[i,]))
  }, mc.cores = parcores)
  
  
  b2s2<-sapply(1:nrow(betaij), function(i){
    sum(betaij[i,mis.inds[[i]]]^2/sij2[i,mis.inds[[i]]])
  })
  bs22<-sapply(1:nrow(betaij), function(i){
    sum(betaij[i,mis.inds[[i]]]/sij2[i,mis.inds[[i]]])^2
  })
  os22<-sapply(1:nrow(betaij), function(i){
    sum(1/sij2[i,mis.inds[[i]]])
  })
  
  # tau2f2<-function(param, b2s2, bs22, os22, gammai=gammai){
  #   tau2<-exp(param)
  #   -1*  sum(  sapply(1:nrow(betaij), function(i){
  #     gammai[i]* (.5*log(1/tau2) -.5*log(os22[i] + 1/tau2) -
  #                   .5*(b2s2[i] - 
  #                         bs22[i]/(os22[i] + 1/tau2)))
  #   }))
  # }
  # 
  
  ll<-NA
  
  k_i<-sapply(mis.inds, length)
  
  MM<-nrow(betaij)
  
  mu.hat<-sapply(1:MM, function(i){
    sum(betaij[i,mis.inds[[i]]]/sij2[i,mis.inds[[i]]])/(sum(1/sij2[i,mis.inds[[i]]]) + 1/(tau2))
  })
  
  strt<-Sys.time()
  for(iter in 1:(maxIter)){
    
    deltaij<- mclapply(1:nrow(betaij), function(i){
      deltis(betaij_i=betaij[i,], sij2_i=sij2[i,], alpha=alpha, f=f)
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
      llbR1<-unlist(mclapply(1:nrow(betaij), function(i){
        llbR1_i(betaij_i = betaij[i,mis.inds[[i]]], sij2_i=sij2[i,mis.inds[[i]]], tau2inv = 1/tau2)
      }, mc.cores=parcores)) - (k_i/2)*log(2*pi)
      
      llbR0<-unlist( mclapply(1:nrow(betaij), function(i){
        llbR0_i(betaij_i = betaij[i,mis.inds[[i]]], sij2_i = sij2[i,mis.inds[[i]]], alpha = alpha, f=f)
      }, mc.cores = parcores))
    }
    
    gammai<-unlist(mclapply(1:MM, function(i){
      gam<- p*exp(llbR1[i])/(exp(logsumexp(c(log(p) + llbR1[i], log(1-p) + llbR0[i]))))
      if(is.na(gam)){
        m<-which.max(c(log(1-p) + llbR0[i],
                       log(p) + llbR1[i]))                      
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
    fn<-sum(unlist(mclapply(1:MM, function(i){
      (1-gammai[i])*sum(deltaij[mis.inds[[i]],i])
    }, mc.cores = parcores)))
    fd<-sum((1-gammai)*k_i)
    f<-fn/fd
    #f<-sum((1-gammai)*(colSums(deltaij, na.rm=TRUE))) / (sum((1-gammai)*k_i))   ## later, replace na.rm with index to sum
    
    alpha<-sum((1-gammai)*
                 unlist(mclapply(1:MM, function(i){
                   sum((1-deltaij[mis.inds[[i]],i])*(betaij[i,mis.inds[[i]]]^2 / sij2[i,mis.inds[[i]]]))}, mc.cores = parcores)))/
      sum((1-gammai)*unlist(mclapply(1:MM, function(i){
        sum(1-deltaij[mis.inds[[i]],i])}, mc.cores = parcores)))
    
    # system.time({
    # tau2Fit<-nlminb(start=log(tau2), objective = tau2f, betaij=betaij, sij2=sij2,  gammai=gammai, control=list(trace=2))
    
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
    
    mu.hat<-unlist(mclapply(1:MM, function(i){
      sum(betaij[i,mis.inds[[i]]]/sij2[i,mis.inds[[i]]])/(sum(1/sij2[i,mis.inds[[i]]]) + 1/tau2)
    }, mc.cores = parcores))
    
    p<-sum(gammai)/MM
    if(verbose > 0){
      print("m step finished.")
      print(paste0("p=", round(p, 3), " f=", round(f, 3), " tau2=", round(tau2, 8), " alpha=", round(alpha, 3)))
    }
    
    llbR1<-unlist(mclapply(1:nrow(betaij), function(i){
      llbR1_i(betaij_i = betaij[i,mis.inds[[i]]], sij2_i=sij2[i,mis.inds[[i]]], tau2inv = 1/tau2)
    }, mc.cores=parcores)) - (k_i/2)*log(2*pi)
    
    llbR0<-unlist( mclapply(1:nrow(betaij), function(i){
      llbR0_i(betaij_i = betaij[i,mis.inds[[i]]], sij2_i = sij2[i,mis.inds[[i]]], alpha = alpha, f=f)
    }, mc.cores = parcores))
    
    
    ll[iter]<-sum(unlist(mclapply(1:MM, function(i){
      logsumexp(c(log(p) + llbR1[i], log(1-p) + llbR0[i]))
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
  se.hat<-unlist(mclapply(1:nrow(betaij), function(i){
    sqrt(1/(sum(1/sij2[i,mis.inds[[i]]]) + 1/tau2))
  }, mc.cores = parcores))
  outliermat<-data.table(snp=1:ncol(deltaij),
                         t(deltaij))
  colnames(outliermat)[2:ncol(outliermat)]<-paste0("deltai_",1:dim(deltaij)[1])
  outliermat[,gammai:=gammai]
  outlierprobs<-outliermat[,#(paste0("outlier_prob", 1:dim(deltaij)[1])):=
                           lapply(.SD, function(x){(1-gammai)*(1-x)}),
                           .SDcols=paste0("deltai_",1:dim(deltaij)[1]), by=.(snp)]
  colnames(outlierprobs)[2:ncol(outlierprobs)]<-paste0("outlier_prob", 1:dim(deltaij)[1])
  # if(!is.null(colnames(betaij))){
  #   dnames<-gsub("beta.*\\_","",colnames(betaij))
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
    nullMod[[j]]<-em_std_f(pdat[[j]]$betaij, pdat[[j]]$sij2_sample, 
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
  
  zeroL<-mclapply(1:nrow(sij2), function(i){
    which(sij2[i,]==0)
  }, mc.cores = parcores)
  chk<- which(sapply(zeroL, length) > 0)
  if(length(chk) > 0){
    for(i in chk){
      sij2[i,zeroL[[i]]]<-NA
    }
  }
  infL<-mclapply(1:nrow(sij2), function(i){
    which(is.infinite(sij2[i,]))
  }, mc.cores = parcores)
  chk<- which(sapply(infL, length) > 0)
  if(length(chk) > 0){
    for(i in chk){
      sij2[i,infL[[i]]]<-NA
    }
  }
  p<-model$p
  f<-model$f
  tau2<-model$tau2
  alpha<-model$alpha
  M<-ceiling(Mnull/(1-p))
  sample_inds<-sample(nrow(sij2),M,replace = TRUE)
  sij2_sample<-sij2[sample_inds,]
  mis.inds<-lapply(1:nrow(sij2_sample), function(i){
    which(!is.na(sij2_sample[i,]))
  })
  Ri<-rbinom(M, size=1, prob=p)
  mu<-Ri*rnorm(M, 0, sqrt(tau2))
  k_i<-sapply(1:nrow(sij2_sample), function(i){sum(!is.na(sij2_sample[i,]))})
  sj2<-lapply(1:M, function(i){sij2_sample[i,]})
  f.ind<-lapply(1:M, function(i){
    vec<-rep(NA, dim(sij2)[2])
    vec[mis.inds[[i]]]<-rbinom(k_i[i], size=1, prob=f)
    vec
  })
  
  betaj<-lapply(1:M, function(i){
    vec<-rep(NA, dim(sij2)[2])
    vec[mis.inds[[i]]]<-
      Ri[i] * rnorm(k_i[i], mu[i], sd=sqrt(sj2[[i]][mis.inds[[i]]])) + 
      (1-Ri[i]) * ((f.ind[[i]][mis.inds[[i]]])*rnorm(k_i[i], 0, sqrt(sj2[[i]][mis.inds[[i]]])) + 
                     (1-f.ind[[i]][mis.inds[[i]]])*rnorm(k_i[i], 0, sqrt(alpha)*sqrt(sj2[[i]][mis.inds[[i]]])))
    
    vec
  })
  
  meta.z<-sapply(1:length(betaj), function(i){
    sum(betaj[[i]][mis.inds[[i]]]/sj2[[i]][mis.inds[[i]]]) / sqrt(sum(1/sj2[[i]][mis.inds[[i]]]))
  })
  betaij<-matrix(unlist(betaj), ncol=dim(sij2)[2], byrow=TRUE)
  return(list(betaij=betaij,
              sij2_sample=sij2_sample,
              sample_inds=sample_inds,
              meta.z=meta.z, 
              Ri=Ri, 
              f.ind=f.ind,
              mui=mu))
}
