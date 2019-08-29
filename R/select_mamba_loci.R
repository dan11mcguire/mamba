#' Select loci to include in mamba() model, combining clumped and randomly pruned variants. 
#' 
#'  
#' Given 
#'  1. a set of clumped significant SNPs from meta-analysis which we want to assess with MAMBA, and 
#'  2. a random set of pruned markers in approximate LD, 
#' Identify pruned variants no closer than (bpdiff) away from clumped variants, 
#' to be used in fitting mamba, to be used in fitting mamba() model.  
#'  The inclusion of randomly pruned variants allows the non-replicable zero-effect cluster of the mixture model to be reliably estimated.
#' 
#'  
#' @param meta_file meta-analysis association file (with columns as outputted by mamba::make_mamba_summary_files()$assoc. i.e. plink assoc file format ) 
#' @param clump_file .clumped file output from plink identifying sentinel variants of interest for mamba model
#' @param prunedsnps_file a file with a set of SNP-ids (1 per line) which are approximately independent 
#'   according to a reference panel.  If NULL, a set of SNP id's (in chr:bp_ref_alt format) from HRC reference panel are used.  These were created using plink with parameters --indep-pairwise 500kb 1 0.1,  --maf 0.01  
#' @param bpdiff BP distance between a randomly pruned snp and a clumped snp of interest.
#' @param pthresh maximum p-value for index snps of interest.  If not provided, this is taken as maximum p-value of the clumped sentinel variants. 
#' @return a data.table with the inputted meta_file, subset to SNPs which are intended to be included in the mamba model.  
#'  The columns index_snp and prunesnp are indicators for whether the SNP belonged to either the clumped_file or is a randomly pruned marker.
#' 
#' @export
#'
#' @examples 
#'  \dontrun{
#'    ## see online vignette for detailed example. 
#'     
#' }





select_mamba_loci<-function(meta_file,clump_file,bpdiff=5*10^5L,prunedsnps_file=NULL,pthresh=NULL){

#meta_file="ivw_met_chr14.tsv"
#clump_file="example_chr14.clumped"
#if(missing(prunefile)) prunefile="prunesnps"
#prunefile="prunesnps"
#min_maf_clump<-0
#bpdiff=5*10^5L
#clump<-copy(as.data.table(clump_table))
#met<-copy(as.data.table(meta_table))


  strt<-Sys.time()
  clump<-fread(clump_file)
  clump<-clump[,c(1:11)][order(CHR,BP)]
  clump[,index_snp:=1]
  
  if(missing(prunedsnps_file)){
	 data(prunedsnps_HRC, package="mamba",envir=environment())
   } else {
     prunevars<-fread(prunedsnps_file)
     prunevars<-prunevars[,1]
     setnames(prunevars, old=names(prunevars), new="SNP")
   }

  met<-fread(meta_file)
  met[FRQ<=0.5,maf:=FRQ]
  met[FRQ>0.5,maf:=1-FRQ]

  if(missing(pthresh)) pthresh<-clump[,max(P)] + 1e-10

  metp<-met[P <= as.numeric(pthresh)]
  stp<-Sys.time()
  print("input read.  creating rema dataset..")
  print(stp-strt)
  
  strt<-Sys.time()
  metp<-merge(metp, clump[,.(CHR,BP, SNP,index_snp)], all=TRUE)
  metp[is.na(index_snp),index_snp:=0]
  #sig<-metp[maf < min_maf_clump | index_snp==1]
  sig<-metp[index_snp==1]
  sig<-sig[order(CHR,BP)]
  sig[,prunesnp:=0]

  prune<-merge(met, prunevars)
  prune[,index_snp:=0]
  prune[,prunesnp:=1]
  prune[,summary(P)]
  prune<-prune[!is.na(P)]
  setcolorder(prune, names(sig))
  
  modelsnps<-rbindlist(list(sig, prune))[order(CHR,BP)]
  rm1<-modelsnps[prunesnp==1 & P < as.numeric(pthresh),which=TRUE]
  if(length(rm1) > 0){
    modelsnps<-modelsnps[-rm1]
  }

  modelsnps[,.N,by=.(index_snp,prunesnp)]
  
  modelsnps[,ldif:=abs(BP-shift(BP,1, type="lag")),by=CHR]
  modelsnps[,lneighbor:=shift(as.numeric(prunesnp==0),1, type="lag"),by=CHR]
  rmneighbor<-length(modelsnps[lneighbor==1 & ldif < bpdiff & prunesnp==1, which=TRUE])
  if(rmneighbor > 0){
    while(rmneighbor > 0){
      rm2<-modelsnps[lneighbor==1 & ldif < bpdiff & prunesnp==1 , which=TRUE]
      modelsnps<-modelsnps[-rm2]
      modelsnps[,ldif:=abs(BP-shift(BP,1, type="lag")),by=CHR]
      modelsnps[,lneighbor:=shift(as.numeric(prunesnp==0),1, type="lag"),by=CHR]
      rmneighbor<-length(modelsnps[lneighbor==1 & ldif < bpdiff & prunesnp==1 ,which=TRUE])
      #print(rmneighbor)
    }
  }
  
  modelsnps[,rdif:=abs(BP-shift(BP,1, type="lead")),by=CHR]
  modelsnps[,rneighbor:=shift(as.numeric(prunesnp==0),1, type="lead"),by=CHR]
  rmneighbor<-length(modelsnps[rneighbor==1 & rdif < bpdiff & prunesnp==1,which=TRUE])
  if(rmneighbor > 0){
    while(rmneighbor > 0){
      rm2<-modelsnps[rneighbor==1 & rdif < bpdiff & prunesnp==1, which=TRUE]
      modelsnps<-modelsnps[-rm2]
      modelsnps[,rdif:=abs(BP-shift(BP,1, type="lead")),by=CHR]
      modelsnps[,rneighbor:=shift(as.numeric(prunesnp==0),1,type="lead"),by=CHR]
      rmneighbor<-length(modelsnps[rneighbor==1 & rdif < bpdiff & prunesnp==1,which=TRUE])
      #print(rmneighbor)
    }
  }
  modelsnps[,(c("ldif", "lneighbor", "rdif", "rneighbor")):=NULL]
  
  return(modelsnps)  
}



  
  
  
  
  
  

