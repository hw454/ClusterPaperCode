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
alg0 = "Hartigan-Wong"
alg0_lab <- gsub('[[:punct:] ]+','', alg0)
which_transform0 = ""

# cochran's Qtest >>>>
q.meta.test = function(b_x, se_x){
  se_x[which(se_x==0)]=1
  w0 = se_x^-2
  w = w0/sum(w0)
  
  meta.bet_x = sum(b_x*w)
  meta.se_x = (sum(w^2*se_x^2))^0.5
  Q = sum( (b_x-meta.bet_x)^2 * w0 )
  df = length(b_x)-1
  Q.pval = 1 - pchisq(Q,df)
  return(list("Q"=Q, "Q.pval"=Q.pval))
}

#load QC-filtered data
load(paste0(res_dir,"/QCdata_",EXP_pheno,".Rdata"))

# k-means clustering data set up
stdBeta_df_noEXP_nosus = stdBeta_df_noEXP[-sus_SNP_ind,]
eff_df = abs(stdBeta_df_noEXP_nosus)
eff_dfs = t(scale(t(eff_df)))
# -- edit - save the data used for clustering
write.csv(eff_dfs, paste0(res_dir,"/",which_transform0,"eff_dfs.csv"))
# -- end edit

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
  kmeans.re <- kmeans(eff_dfs, 
                      centers = i, 
                      nstart = 50, 
                      iter.max = 300,
                      algorithm = alg0) # -- edit - added the algorithm parameter
  cluster_list[[length(cluster_list)+1]] = kmeans.re
  IC = kmeansIC(kmeans.re)
  IC_df[(i-1),1] = IC$AIC
}

# cluster number identification for each observation
kmeans.minAIC = cluster_list[[which(IC_df$AIC==min(IC_df$AIC))]]
nClust.AIC = max(kmeans.minAIC$cluster)

# plot AIC
pdf(file=paste0(res_dir,"/AIC-ncluster_",EXP_pheno,".pdf"), width=7, height = 5)
plot(IC_df$nCluster,IC_df$AIC, xlab="nCluster", ylab="AIC")
abline(v=nClust.AIC,col="red")
dev.off()

# assigning SNPs to clusters
AICclusters_rsid = list()
for(i in 1:nClust.AIC){
  AICclusters_rsid[[i]] = names(kmeans.minAIC$cluster)[which(kmeans.minAIC$cluster==i)]
}
print(paste0("Number of SNPs in each AIC grouped clusters: ", paste0(lengths(AICclusters_rsid), collapse = ", ")))
AICclusters_rsid_df = t(plyr::ldply(AICclusters_rsid, rbind))
write.csv(AICclusters_rsid_df, paste0(res_dir,"/",which_transform0,"_",alg0_lab,"AICclusters_rsid_",EXP_pheno,".csv"), row.names = FALSE)
