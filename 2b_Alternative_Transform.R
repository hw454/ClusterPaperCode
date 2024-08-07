library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(TwoSampleMR)

## after quality control od SNPs and traits in the PheWAS matrix, we proceed with clustering, enrichment analysis and MR est8imation per cluster given a certain outcome
## Data inputs needed: .RData file from previous step, pre-prepared exposure SNP effect file (genome-wide significant and clumped) and the corresponding effects in the outcome trait

# variable set-up
EXP_pheno= "21001"
OUT_pheno = "845"
res_dir = "./data"
alg1 = "Lloyd"
alg1_lab <- gsub('[[:punct:] ]+','', alg1)
which_transform1 = "diff_transform"
#setwd(res_dir)

alt_transform <- function(exp_dat, eff_df_in){
  eff_df <- data.frame(eff_df_in)
  neg_nums <- which(exp_dat<0)# Get the rows with negatives in exp_dat
  neg_names_dat <- names(neg_nums)
  neg_names <- rownames(exp_dat)[neg_nums]
  eff_df[neg_names,] <- -eff_df[neg_names,]
  eff_df_out <- t(scale(t(eff_df)))
  exp_dat_out <- exp_dat
  exp_dat_out[neg_names_dat]<- -exp_dat[neg_names_dat]
  exp_dat_out <- scale(exp_dat_out)
  out_list <- list("eff_df"= eff_df_out, "exp_dat" = exp_dat_out)
  return(out_list)
}


#load QC-filtered data
load(paste0(res_dir,"/QCdata_",EXP_pheno,".Rdata"))

# k-means clustering data set up
alt_out_list <- alt_transform(stdBeta_EXP,stdBeta_df_noEXP_nosus)
eff_dfs_alt <- alt_out_list$eff_df
exp_dat_alt <- alt_out_list$exp_dat
colnames(eff_dfs_alt) <- colnames(stdBeta_df_noEXP_nosus)
write.csv(eff_dfs_alt, paste0(res_dir,"/",which_transform1,"_eff_dfs.csv"))
write.csv(exp_dat_alt, paste0(res_dir,"/",which_transform1,"_stdBeta_EXP.csv"))

# Fitting K-Means clustering Model 
kmeansIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(data.frame(AIC = D + 2*m*k,
                    BIC = D + log(n)*m*k))
}

cluster_list = list()
IC_df = data.frame(matrix(data=NA, nrow = 49, ncol=1))
colnames(IC_df) = c("AIC")
IC_df$nCluster = 2:50

for(i in 2:50){
  set.seed(240) # setting seed
  kmeans.re <- kmeans(eff_dfs_alt, 
                      centers = i, 
                      nstart = 50, 
                      iter.max = 300,
                      algorithm = alg1)
  cluster_list[[length(cluster_list)+1]] = kmeans.re
  IC = kmeansIC(kmeans.re)
  IC_df[(i-1),1] = IC$AIC
}

# cluster number identification for each observation
kmeans.minAIC = cluster_list[[which(IC_df$AIC==min(IC_df$AIC))]]
nClust.AIC = max(kmeans.minAIC$cluster)

# plot AIC
pdf(file=paste0(res_dir,"/",which_transform1,"_AIC-ncluster_",EXP_pheno,".pdf"), width=7, height = 5)
plot(IC_df$nCluster,IC_df$AIC, xlab="nCluster", ylab="AIC")
abline(v=nClust.AIC,col="red")
dev.off()

# assigning SNPs to clusters
AICclusters_rsid2 = list()
for(i in 1:nClust.AIC){
  AICclusters_rsid2[[i]] = names(kmeans.minAIC$cluster)[which(kmeans.minAIC$cluster==i)]
}
print(paste0("Number of SNPs in each AIC grouped clusters: ", paste0(lengths(AICclusters_rsid2), collapse = ", ")))
AICclusters_rsid_df2 = t(plyr::ldply(AICclusters_rsid2, rbind))
write.csv(AICclusters_rsid_df2, paste0(res_dir,"/",which_transform1,"_",alg1_lab, "_AICclusters_rsid_",EXP_pheno,".csv"), row.names = FALSE)

