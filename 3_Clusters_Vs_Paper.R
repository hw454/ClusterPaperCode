#---------------------- FUNCTIONS

long_form_cluster_table <- function(df){
  clust_nums = 1:dim(df)[2]
  if (is.matrix(df)){df <- data.frame(df)}
  df_long <- df %>% reshape2::melt(
    measure.vars = names(df),
    variable.name = "Cluster",
    value.name = "rsid"
  )
  df_long <- df_long %>% drop_na("rsid")
  #tibble::column_to_rownames(df_long, var = "rsid")
  #print(df_long)
  return(df_long)
}

compare_df1_to_df2 <- function(df1, df2){
  # Since the clusters may have none integer create a mapping between the terms and their position. 
  cnums1 = unique(df1["Cluster"])
  cnums2 = unique(df2["Cluster"])
  rownames(cnums1) <- 1:nrow(cnums1)
  rownames(cnums2) <- 1:nrow(cnums2)
  res_mat = matrix(0,nrow(cnums1), nrow(cnums2))
  snp_list = as.vector(df1["rsid"])
  #print(snp_list)
  for (i in 1:nrow(df1)){
    snp <- df1[i, "rsid"]
    #print(snp)
    ind1 = which(df1["rsid"]==snp, arr.ind = T)[1]
    ind2 = which(df2["rsid"]==snp, arr.ind = T)[1]
    cn1 = as.character(df1[ind1,"Cluster"])
    cn2 = as.character(df2[ind2,"Cluster"])
    i = which(cn1 == cnums1["Cluster"])
    j = which(cn2 == cnums2["Cluster"])
    res_mat[i,j] = res_mat[i,j] + 1
  }
  res_mat_per <- t(apply(res_mat,1, function(x) x/sum(x)))
return(res_mat_per)
}
# ------------------

clusters_paper_df <- fread(paste0(res_dir,"/ClusterMembership_BMI.csv"),
                           header = TRUE)

# - Get both sets of clusters into long-form
AICclusters_rsid_df_long <- long_form_cluster_table(AICclusters_rsid_df)
df2 <- AICclusters_rsid_df_long 
df <- long_form_cluster_table(AICclusters_rsid_df2)
#clusters_paper_long <- long_form_cluster_table(clusters_paper_df)
#df <- clusters_paper_long
# - Find overlap between clusters
comp_df <- compare_df1_to_df2(df, df2)

# -- PLOT
# - plot the comparison table
## convert to tibble, add row identifier, and shape "long"
dat2 <-
  comp_df %>%
  as_tibble() %>%
  tibble::rownames_to_column("Alg1") %>%
  pivot_longer(-Alg1, names_to = "Alg0", values_to = "value") %>%
  mutate(
    Alg1 = factor(Alg1, levels = 1:10),
    Alg0 = factor(gsub("V", "", Alg0), levels = 1:10)
  )


alg0 = "paper"
p<- ggplot(dat2, aes(Alg0, Alg1)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "red")
p + ggtitle(paste(which_transform,alg0,"compared to",alg1))
p + xlab(paste(which_transform0,"_",alg0)) +ylab(which_transform1,"_",alg1)
plot_dir <- "./plots/"
dir.create(file.path(plot_dir), showWarnings = FALSE)
ggsave(paste0(plot_dir,which_transform,"_",alg0, "membership_compare_",alg1,"_R.png"))

