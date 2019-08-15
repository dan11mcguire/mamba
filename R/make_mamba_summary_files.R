
#' Create input datasets for mamba() function and clumping
#'  
#' Take gwas summary stat files, 
#'   do ivw-meta-analysis, 
#'   create an assoc dataset which could be used for clumping association statistics in plink,
#'   create beta and s2 input datasets for mamba() model 
#'  
#' @param summary_files A list of summary statistic files with common column names
#' @param standardize_beta_and_s2_datasets in the outputted b and s2 data tables, 
#'   do you want the effect size estimates and variances standardized (i.e. beta = beta_k*sqrt(2*af_k*(1-af_k)))?  Default is TRUE.
#' @param beta_name name of effect size estimate column
#' @param s2_name name of effect size variance estimate column
#' @param ref_name name of reference allele column
#' @param alt_name name of alt allele column
#' @param chr_name name of chr column
#' @param n_name name of n column
#' @param af_name name of af column
#' @param bp_name name of bp column
#' @return a list of three data.tables (beta, s2, and assoc)
#' @export
#'
#' @examples 
#'  \dontrun{
#'     datadir<-system.file("extdata", package="mamba")
#'     datafiles<-paste0(datadir, "/", list.files(datadir))
#'     looky<-make_mamba_summary_files(summary_files=datafiles, 
#'                              standardize_beta_and_s2_datasets=TRUE,
#'                              beta_name="betaj", s2_name="sj2")
#'     looky$assoc[1:5]
#'     looky$b[1:5] 
#'     looky$s2[1:5]
#' }


make_mamba_summary_files<-function(summary_files, standardize_beta_and_s2_datasets=TRUE,
                                   beta_name="b", s2_name="s2", 
                                   n_name="n", af_name="af", 
                                   ref_name="ref", alt_name="alt", 
                                   chr_name="chr", bp_name="bp"){

#summary_files<-datafiles
#beta_name="betaj"; s2_name="sj2"; 
#n_name="n"; af_name="af"; 
#ref_name="ref"; alt_name="alt"; 
#chr_name="chr"; bp_name="bp"
#standardize_beta_and_s2_datasets=TRUE

needed<-c(beta_name, s2_name, n_name, af_name, ref_name, alt_name, chr_name, bp_name)
  
for(i in 1:length(summary_files)){
   d_tmp<-fread(summary_files[i])
   if(any(!needed %in% names(d_tmp))){
      stop(paste(needed[!needed %in% names(d_tmp)], " column not found in summary file ", i, " "))
   }
   d_tmp=d_tmp[,needed,with=F] 
   setnames(d_tmp, old=c(beta_name, s2_name, n_name, af_name, ref_name, alt_name), 
            new=c(paste0(c("b", "s2", "n", "af"), "_", i),"ref",  "alt" ))
    if(i == 1){
    d<-d_tmp
    }else{
    d<-merge(d, d_tmp, by=c("chr", "bp", "ref", "alt"), all=T)
    }
  }



beta<-melt(d[,c("chr", "bp", "ref", "alt",
                   grep("^b_",names(d),value=TRUE)),with=FALSE], 
             id.vars=c("chr", "bp", "ref", "alt"), 
             value.name = "b", 
             variable.name="cohort")[,cohort:=gsub("b_","",cohort)]

s2<-melt(d[,c("chr", "bp", "ref", "alt",
            grep("^s2_",names(d),value=TRUE)),with=FALSE], 
            id.vars=c("chr", "bp", "ref", "alt"), 
            value.name = "s2", 
            variable.name="cohort")[,cohort:=gsub("s2_","",cohort)]

bs<-merge(beta,s2,by=c("chr", "bp", "ref", "alt", "cohort"))[!is.na(b)][!is.na(s2)]

n<-melt(d[,c("chr", "bp", "ref", "alt",grep("n_",names(d),value=TRUE)),with=FALSE], 
           id.vars=c("chr", "bp", "ref", "alt" ), 
           value.name = "n", 
           variable.name="cohort")[,cohort:=gsub("n_","",cohort)]

af<-melt(d[,c("chr", "bp", "ref", "alt",grep("af_",names(d),value=TRUE)),with=FALSE], 
          id.vars=c("chr", "bp", "ref", "alt" ), 
          value.name = "af", 
          variable.name="cohort")[,cohort:=gsub("af_","",cohort)]

naf<-merge(n,af,by=c("chr", "bp", "ref", "alt","cohort"))[!is.na(n)][!is.na(af)]
bsnaf<-merge(bs, naf, by=c("chr", "bp", "ref", "alt","cohort"))


bs<-bs[,.(BETA=sum(b/s2)/sum(1/s2), SE=sqrt(1/sum(1/s2)),k=.N),by=.(chr,bp,ref,alt)] 
 
bs[,P:=2*pnorm(abs(BETA/SE),lower.tail=FALSE)]
bs[,z:=BETA/SE]
  
tot_n_af<-naf[,.(tot_ac=sum(2*n*af),tot_af=sum(2*n*af)/sum(2*n),
                 tot_n=sum(n)),by=.(chr,bp,ref,alt)] 

assoc<-merge(tot_n_af, bs, by=c("chr", "bp", "ref", "alt"))

assoc[,SNP:=paste0(chr,":",bp,"_",ref,"_",alt)]
assoc[,INF0:=1]
setcolorder(assoc, 
            c("chr", "bp", "SNP", "alt", "ref", 
              "tot_af", "INF0", "BETA", "SE", "P", "k", "z", 
              "tot_n", "tot_ac"))
setnames(assoc, old=c("chr", "bp", "alt", "ref", "tot_af"),
         new=c("CHR", "BP", "A1", "A2", "FRQ"))


  if(standardize_beta_and_s2_datasets){
    bsnaf[,b:=b*sqrt(2*af*(1-af))]
    bsnaf[,s2:=s2*2*af*(1-af)]
  }

keys<-c("chr", "bp","snp","ref", "alt")

bsnaf[,snp:=paste0(chr, ":", bp, "_", ref, "_", alt)]

b_out<-dcast(bsnaf,chr + bp + snp + ref + alt ~cohort, value.var="b")
cohorts<-setdiff(names(b_out), keys)
setcolorder(b_out, c(keys,paste0(sort(as.numeric(cohorts)))))
setnames(b_out, old=setdiff(names(b_out), keys), new=paste0("b_", setdiff(names(b_out), keys)))

s2_out<-dcast(bsnaf, chr + bp + snp + ref + alt ~cohort, value.var="s2")
cohorts<-setdiff(names(s2_out), keys)
setcolorder(s2_out, c(keys,paste0(sort(as.numeric(cohorts)))))
setnames(s2_out, old=setdiff(names(s2_out), keys), new=paste0("s2_", setdiff(names(s2_out), keys)))


  return(list(
    assoc=assoc,
    beta=b_out,
    s2=s2_out
  ))


}

