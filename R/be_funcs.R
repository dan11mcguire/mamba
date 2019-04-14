M_shuffle<-function(t_config){
  bin<-c(1,0)
  k<-length(t_config)
  M_shuf<-sample(t_config)
  moves<-lapply(1:length(t_config), function(j){
    t<-t_config
    t[j]<-bin[bin!=t_config[j]]
    return(t)
  })
  moves[[k+1]]<-M_shuf
  return(moves[[sample(k+1, 1)]])
}

lPT<-function(TT=t_config, Alpha=1, Beta=1, k){
  lbeta(sum(TT) + Alpha, k-sum(TT) + Beta) - lbeta(Alpha, Beta)
}
lPX0<-function(i, t_config, betaij, sij2){
  sum(dnorm(betaij[i,t_config==0], 0, sqrt(sij2[i,t_config==0]), log=TRUE))
}

lPX1<-function(i, t_config, betaij, sij2, sigma2=0.2^2){
  n<-sum(t_config==1)
  xbar<-sum(betaij[i,t_config==1]/sij2[i,t_config==1]) / sum(1/sij2[i,t_config==1])
  vbar<-1/sum(1/sij2[i,t_config==1])
  nd<-dnorm(xbar, 0, sqrt(vbar + sigma2), log=TRUE)
  #n<-sum(t_config==1)
  lcbar<- .5*(-(n-1)*log(2*pi) +  
                sum(log(1/sij2[i,t_config==1])) - log(sum(1/sij2[i,t_config==1])) -
                (sum(betaij[i,t_config==1]^2 /sij2[i,t_config==1]) - 
                   (sum(betaij[i,t_config==1]/sij2[i,t_config==1])^2)/sum(1/sij2[i,t_config==1])))
  return(nd + lcbar)
}

exact_mi<-function(i, betaij, sij2, parcores=1){
  k<-dim(betaij)[2]
  t_configMat<-data.table(expand.grid(lapply(1:k, function(i)c(1,0))))
  colnames(t_configMat)<-paste0("t",1:k)
  t_configMat[,(paste0("t",1:k)):=lapply(.SD, as.integer)]
  gtvec<-unlist(mclapply(1:nrow(t_configMat), function(iter){
    t_config<-unlist(t_configMat[iter])
    lpx1<-lPX1(i=i, t_config = t_config, betaij=betaij, sij2=sij2)
    lpx0<-lPX0(i=i, t_config = t_config, betaij=betaij, sij2=sij2)  
    lpt<- lPT(TT=t_config, k=k)
    ll<-c(lpx1,lpx0,lpt)
    sum(ll[!is.na(ll)])
  }, mc.cores=parcores))
  t_configMat[,gt:=gtvec]
  mi<-sapply(1:k, function(j){
    exp(logsumexp(t_configMat[get(paste0("t",j))==1,gt])) / exp(logsumexp(t_configMat[,gt]))
  })
  return(list(mi=mi,
              t_configMat=t_configMat,
              metaz=sum(betaij[i,]/sij2[i,])/sqrt(sum(1/sij2[i,])),
              z=betaij[i,]/sqrt(sij2[i,])))
}


mh_mi<-function(i, ns, betaij, sij2, seed=sample(10^3,1)){
  set.seed(seed)
  k<-ncol(betaij)
  t_config<-  as.integer(sample(c(1, 0), prob = c(.5, .5), k, replace = TRUE))
  gtvec<-NA #rep(NA, ns)
  #t_configList<-list()
  t_configMat<-data.table(matrix(NA, ncol=k, nrow=ns))
  colnames(t_configMat)<-paste0("t",1:k)
  t_configMat[,(paste0("t",1:k)):=lapply(.SD, as.integer)]
  iter<-1
  while(iter < ns){
    lpx1<-lPX1(i=i, t_config = t_config, betaij=betaij, sij2=sij2)
    lpx0<-  lPX0(i=i, t_config = t_config, betaij=betaij, sij2=sij2)  
    lpt<-  lPT(TT=t_config,k=k)
    ll<-c(lpx1,lpx0,lpt)
    gtvec[iter]<-sum(ll[!is.na(ll)])
    t_configMat[iter,paste0("t",1:k):=as.list(t_config)]
    if(iter==1){
      t_config<-M_shuffle(t_config = t_config)
    } else {
      p<-min(1, exp(gtvec[iter] - gtvec[iter-1]))
      p<-max(p, 0)
      t_config<-list(t_config, M_shuffle(t_config=t_config))[[sample(c(2,1), size=1, prob=c(p, 1-p))]]
    }
    print(iter)
  }
  t_configMat[,gt:=gtvec]
  mi<-NA
  mi<-sapply(1:k, function(j){
    exp(logsumexp(t_configMat[get(paste0("t",j))==1,gt])) / exp(logsumexp(t_configMat[,gt]))
  })
  return(list(mi=mi,
              t_configMat=t_configMat,
              metaz=sum(betaij[i,]/sij2[i,])/sqrt(sum(1/sij2[i,])),
              z=betaij[i,]/sqrt(sij2[i,]), 
              seed=seed))
}


metasoft_format<-function(inds, betaij, sij2, parcores=1){
  mout<-  do.call(rbind, 
                  mclapply(inds, function(i){
                    do.call(cbind, 
                            lapply(1:ncol(betaij), function(j){
                              cbind(betaij[i,j], sqrt(sij2[i,j])) 
                            })
                    )
                  }, mc.cores = parcores))
  return(data.table(rs=paste0("rs",inds), mout))
}

# estimate fdr from Ri.hat:

myfdr<-function(gammai){
  rihat.seq<-seq(0+1e-10,1-1e-10,length.out = 5e3)
  fdr<-sapply(1:length(rihat.seq), function(u){
    data.table(rihat=myfit$gammai)[rihat >= rihat.seq[u]][,mean(1-rihat)]
  })
  fdr.ref<-data.table(fdr=list(fdr), ri=list(rihat.seq))
  fdr.ref
  gd<-data.table(gammai=myfit$gammai)
  gd[,fdr:=mapply(g=gammai, function(g){
    d<-abs(g-fdr.ref$ri[[1]])
    fdr.ref$fdr[[1]][which.min(d)]
  })][,fdr]
}


emp_fdr_control<-function(d, meth="PVALUE_RE2"){
  d[order(get(meth)),i:=1:nrow(m)]
  d[order(get(meth)),emp_Fdr:=cumsum(1-Ri)/i]
  d[order(get(meth))][emp_Fdr!=0][emp_Fdr < 0.05]
  d[Ri==1][,mean(emp_Fdr < 0.05)]
  d[Ri==0][,.(mean(c(m1, m2, m3, m4, m5, m6, m7, m8, m9) > 0.5)), by=1:nrow(d[Ri==0])][,summary(V1)]
  
}

