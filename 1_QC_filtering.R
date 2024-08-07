### 2) QC filtering 
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

## having obtained the top SNPs associated with EXP (genome-wide signficance), and ran PheWAS to get their effects on all UKBB traits we have (filtered for small sample size), we proceed to quality control.
## Data inputs needed: unstanderdised association matrix between IV and traits, SE matrix of association matrix, t-stat matrix of association matrix, P-value matrix of association matrix, trait info (trait, description, sample size), file path of traots for duplicate trait filtering

# variable set-up
EXP_pheno = "21001"

res_dir = "/data"
#--- EDIT - File missing
hail_gcorr_dir = paste0(res_dir,"/both_sexes_rsid_corrmat.csv")
#--- end edit
fpaths_fil_dir = paste0(res_dir,"/fpaths_fil_nfil.txt")
#setwd(res_dir)
exp_gcorr_thresh = 0.75

# read in data
print(getwd())
unstdBeta_df = as.matrix(data.table::fread(paste0(res_dir,"/unstdBeta_df.csv")), rownames=1)
unstdSE_df = as.matrix(data.table::fread(paste0(res_dir,"/unstdSE_df.csv")), rownames=1)
tstat_df = as.matrix(data.table::fread(paste0(res_dir,"/tstat_df.csv")), rownames=1)
pval_df = as.matrix(data.table::fread(paste0(res_dir,"/pval_df.csv")), rownames=1)
trait_info = data.table::fread(paste0(res_dir,"/trait_info_nfil.csv"))

# filter out traits that are NA (albumin and other measurements)
na_traits = which(is.na(trait_info$description)==T)
# not all are also NA in beta file
na_effect = which(is.na(unstdBeta_df)==T, arr.ind=T)
na_both = intersect(na_traits,unique(as.numeric(na_effect[,2])));
if(length(na_both)>0){
  unstdBeta_df = unstdBeta_df[,-na_both]
  unstdSE_df = unstdSE_df[,-na_both]
  tstat_df = tstat_df[,-na_both]
  pval_df = pval_df[,-na_both]
  trait_info = trait_info[-na_both,]
}

# keep only V2 of replicated traits - future update: those with larger N
v2_traits = which(duplicated(trait_info$phenotype, fromLast = F )==T)
fpaths_fil = fread(fpaths_fil_dir, header=F)
if(length(na_both)>0){fpaths_fil = fpaths_fil[-na_both,]} #continuity from above 
v2_traits_2 = which(grepl("v2", fpaths_fil$V1)==T)
# check
print(paste0("Duplicates in fpaths_fil and trait_info match: ",all(v2_traits==v2_traits_2))) #all v2 indices are accounted for, remove indices -1
v1_traits = v2_traits-1
unstdBeta_df = unstdBeta_df[,-v1_traits]
unstdSE_df = unstdSE_df[,-v1_traits]
tstat_df = tstat_df[,-v1_traits]
pval_df = pval_df[,-v1_traits]
trait_info = trait_info[-v1_traits,]
print(paste0("Dimension of unstd SNP-trait matrix: ", dim(unstdBeta_df)[1],", ",dim(unstdBeta_df)[2]))
print(paste0("Number of traits: ", dim(trait_info)[1]))

# check for duplicate exposure trait, use one with larger N
dup_exp = which(trait_info$description == trait_info$description[which(grepl(paste0("^",EXP_pheno), trait_info$phenotype)==TRUE)])
if(length(dup_exp)>1){
  dup_remove = dup_exp[!dup_exp %in% dup_exp[which.max(trait_info$n_eff[dup_exp])]]
  #dup_remove = dup_exp[!dup_exp %in% dup_exp[which.max(trait_info$n_complete_samples[dup_exp])]]
  unstdBeta_df = unstdBeta_df[,-dup_remove]
  unstdSE_df = unstdSE_df[,-dup_remove]
  tstat_df = tstat_df[,-dup_remove]
  pval_df = pval_df[,-dup_remove]
  trait_info = trait_info[-dup_remove,]
  print(paste0(length(dup_remove), " exposure-duplicate trait(s) removed."))
}

# removing traits with gcorr > 0.75 with EXP trait
#--- EDIT - This section has been edited to work with the download file as no file 
# available for working example
hail_df_wide = fread(hail_gcorr_dir, header = TRUE)
n <- dim(hail_df_wide)[2]
hail_df <- hail_df_wide %>% melt(
  na.rm = FALSE,
  measure.vars = names(hail_df_wide)[2:n],
  id.vars = c("ph1"),
  variable.name = , "ph2",
  variable.factor = FALSE,
  value.name = "rg"
)
#--- EDIT - string extraction and matching found no results. 
# - Currently sets ID as the numeric part of the phenotype and ph as the full phenotype
hail_df$ID1 = stringr::str_extract(hail_df$`ph1`, "[0-9]+")
hail_df$ID2 = stringr::str_extract(hail_df$`ph2`, "[0-9]+")
#---- EDIT - phenotype is 21001_irnt not 21001 so nothing is found with an == search.
# - replaced `Phenotype 1` with `ID1`
df_gcorr_EXP = hail_df[which(hail_df$ID1 == EXP_pheno | hail_df$ID2 == EXP_pheno), ]
#--- end edit
high_gcorr = which(abs(df_gcorr_EXP$rg)>=exp_gcorr_thresh)
df_gcorr_EXP_high = df_gcorr_EXP[high_gcorr,]
phenotype_remove = unique(c(df_gcorr_EXP_high$ph1,df_gcorr_EXP_high$ph2))
#--- EDIT - replace description with phonetype
trait_remove1 = trait_info$phenotype[which(trait_info$phenotype %in% phenotype_remove)]
#--- end edit
trait_remove1 = trait_remove1[-which(grepl(paste0("^",EXP_pheno), trait_remove1))] #remove exposure from traits to remove
print(paste0("Number of traits removed due to correlation ",length(trait_remove1)))


if(length(trait_remove1)>0){
  unstdBeta_df = unstdBeta_df[,-which(colnames(unstdBeta_df) %in% trait_remove1)] 
  unstdSE_df = unstdSE_df[,-which(colnames(unstdSE_df) %in% trait_remove1)]
  tstat_df = tstat_df[,-which(colnames(tstat_df) %in% trait_remove1)]
  pval_df = pval_df[,-which(colnames(pval_df) %in% trait_remove1)]
  trait_info = trait_info[-which(trait_info$phenotype %in% trait_remove1),]
  print(paste0(length(trait_remove1), " exposure-correlation high trait(s) removed."))
}

# getting the standardised effects (tstat/sqrt(N))  
exp_col = which(grepl(paste0("^",EXP_pheno), colnames(tstat_df))) 
stdBeta_df = as.matrix(tstat_df) %*% diag(1/(sqrt(trait_info$n_eff)))   #mat %*% diag(1 / dev)
colnames(stdBeta_df) = colnames(unstdBeta_df)
stdBeta_df_noEXP = stdBeta_df[,-exp_col] #removing exposure
stdBeta_EXP = stdBeta_df[,exp_col]
stdSE_noEXP = 1/(sqrt(trait_info$n_eff[-exp_col]))
stdSE_EXP = 1/(sqrt(trait_info$n_eff[exp_col]))

trait_info_noEXP = trait_info[-which(grepl(paste0("^",EXP_pheno), trait_info$phenotype)) ,]

print(paste0("Dimension of std SNP-trait matrix: ", dim(stdBeta_df)[1],", ",dim(stdBeta_df)[2]))
print(paste0("Number of traits: ", dim(trait_info)[1]))

# filtering for SNPs more associated with other traits using t-stat
tstat_dif = matrix(NA, nrow = nrow(stdBeta_df_noEXP), ncol = ncol(stdBeta_df_noEXP))
for(i in 1:nrow(stdBeta_df_noEXP)){ #1:i SNPs
  for(j in 1:ncol(stdBeta_df_noEXP)){ #1:j traits
    tstat_dif[i,j] = (abs(stdBeta_EXP[i]) - abs(stdBeta_df_noEXP[i,j])) / 
      sqrt(stdSE_EXP^2 + stdSE_noEXP[j]^2)
  }
}

tstat_pval_df = pnorm(tstat_dif)  # larger than on EXP
tstat_thresh = 0.05/dim(trait_info_noEXP)[1] # accuracy: divided by nSNPs* nIndp.traits
#--- EDIT - rename tstat_dif to tstat_dif_df, separating the instances of matrix and dataframe so the code is more backawards compatible.
tstat_dif_df = as.data.frame(tstat_dif)
colnames(tstat_dif_df) = colnames(stdBeta_df_noEXP)
rownames(tstat_dif_df) = rownames(stdBeta_df_noEXP)
sig_ind = which(tstat_pval_df<tstat_thresh, arr.ind = T)
sus_SNPs=rownames(tstat_dif_df)[unique(sig_ind[,1])] 
#--- end edit
sus_traits = unique(trait_info_noEXP$description[(sig_ind[,2])+1]) 
print(paste0("SNP filtering revealed ", length(sus_SNPs)," SNPs more strongly associated to traits other than exposure trait ", EXP_pheno,"."))

sus_SNP_ind = which(rownames(stdBeta_df_noEXP) %in% sus_SNPs) #same as sig_ind[1,]
save(list=c("stdBeta_df","stdBeta_df_noEXP","stdBeta_EXP","stdSE_EXP","stdSE_noEXP","sus_SNP_ind","sus_SNPs","sus_traits","trait_info", "pval_df","trait_info_noEXP", "tstat_df","tstat_dif","tstat_pval_df","tstat_thresh","unstdBeta_df","unstdSE_df", "hail_df"), file = paste0(res_dir,"/QCdata_",EXP_pheno,".Rdata"))


#--- EDIT - print the matrix dimensions that will be clustered
stdBeta_df_noEXP_nosus = stdBeta_df_noEXP[-sus_SNP_ind,]
stdBeta_df_nosus = stdBeta_df[-sus_SNP_ind,]
print(paste0("Dimension of std SNP-trait matrix: ", dim(stdBeta_df_noEXP[-sus_SNP_ind,])[1],", ",dim(stdBeta_df_noEXP[-sus_SNP_ind,])[2]))
print(paste0("Number of traits: ", dim(trait_info)[1]))
write.csv(stdBeta_df_noEXP_nosus, paste0(res_dir,"/eff_dfs_pre_transform.csv"))
write.csv(stdBeta_EXP, paste0(res_dir,"/stdBeta_EXP.csv"))
#--- end edit
