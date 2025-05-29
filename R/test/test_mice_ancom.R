## Testing the functions to Manually pool using Rubin's rules

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("ANCOMBC")
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
# BiocManager::install("phyloseq")

library(ANCOMBC)
library(phyloseq)
library(data.table)
library(tidyverse)


data(QMP, package = "ANCOMBC")
  set.seed(123)
  n = 150
  d = ncol(QMP)
  diff_prop = 0.1
  lfc_cont = 1
  lfc_cat2_vs_1 = -2
  lfc_cat3_vs_1 = 1
  
  # Generate the true abundances
  abn_data = sim_plnm(abn_table = QMP, taxa_are_rows = FALSE, prv_cut = 0.05, 
                      n = n, lib_mean = 1e8, disp = 0.5)
  log_abn_data = log(abn_data + 1e-5)
  rownames(log_abn_data) = paste0("T", seq_len(d))
  colnames(log_abn_data) = paste0("S", seq_len(n))
  
  # Generate the sample and feature meta data
  # Sampling fractions are set to differ by batches
  smd = data.frame(samp_frac = log(c(runif(n/3, min = 1e-4, max = 1e-3),
                                     runif(n/3, min = 1e-3, max = 1e-2),
                                     runif(n/3, min = 1e-2, max = 1e-1))),
                   cont_cov = rnorm(n),
                   cat_cov = as.factor(rep(seq_len(3), each = n/3))
                   )
  rownames(smd) = paste0("S", seq_len(n))
  
  fmd = data.frame(taxon = paste0("T", seq_len(d)),
                   seq_eff = log(runif(d, min = 0.1, max = 1)),
                   lfc_cont = sample(c(0, lfc_cont), 
                                     size = d,
                                     replace = TRUE,
                                     prob = c(1 - diff_prop, diff_prop)),
                   lfc_cat2_vs_1 = sample(c(0, lfc_cat2_vs_1), 
                                          size = d,
                                          replace = TRUE,
                                          prob = c(1 - diff_prop, diff_prop)),
                   lfc_cat3_vs_1 = sample(c(0, lfc_cat3_vs_1), 
                                          size = d,
                                          replace = TRUE,
                                          prob = c(1 - diff_prop, diff_prop))) %>%
    mutate(lfc_cat3_vs_2 = lfc_cat3_vs_1 - lfc_cat2_vs_1)
  
  # Add effect sizes of covariates to the true abundances
  smd_dmy = model.matrix(~ 0 + cont_cov + cat_cov, data = smd)
  log_abn_data = log_abn_data + outer(fmd$lfc_cont, smd_dmy[, "cont_cov"] )
  log_abn_data = log_abn_data + outer(fmd$lfc_cat2_vs_1, smd_dmy[, "cat_cov2"])
  log_abn_data = log_abn_data + outer(fmd$lfc_cat3_vs_1, smd_dmy[, "cat_cov3"])
  
  # Add sample- and taxon-specific biases
  log_otu_data = t(t(log_abn_data) + smd$samp_frac)
  log_otu_data = log_otu_data + fmd$seq_eff
  otu_data = round(exp(log_otu_data))
  
  
smd


ps<-phyloseq(otu_table(otu_data, taxa_are_rows = T),
             sample_data(smd))

output = ancombc2(data = ps, 
                  fix_formula = "cont_cov + cat_cov", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = FALSE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "cat_cov", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = F, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)

res_prim = output$res
res_pair = output$res_pair

res_prim = output$res %>%
  mutate_if(is.numeric, function(x) round(x, 2))
res_prim %>% select( taxon, contains("W_"), contains("p_")) 

## spike in missing:

## randomly setting 2 continuous variables to NA
smd[sample(150,2),2]<-NA

## randomly setting 2 categorical variables to NA
smd[sample(150,2),3]<-NA


## run MICE
library(mice)

smd<-rownames_to_column(smd)

smd$samp_frac<-as.numeric(smd$samp_frac)
smd$cont_cov<-as.numeric(smd$cont_cov)
smd$cat_cov<-factor(smd$cat_cov)

imp_prelim<-mice(smd, maxit = 5, seed = 123, m = 5)

imp_prelim$predictorMatrix[1,1]<-0


imp<-mice(smd, maxit = 5, seed = 123, m = 10,predictorMatrix = imp_prelim$predictorMatrix)

plot(imp)
imp$loggedEvents


impl<-complete(imp, action = "long")

split_imp<-split(impl, impl$.imp)

set.seed(123)

# Make the sample ID col a rowname and coerce back into phyloseq objects. 

split_imp_sd<-lapply(split_imp,`rownames<-`, split_imp[[1]]$rowname)
split_imp_sd<-lapply(split_imp_sd, sample_data)


# coerce to phyloseq

ps_imp<-lapply(split_imp_sd,phyloseq, otu_table(otu_data,taxa_are_rows = T))

output_imp <-lapply(ps_imp, ancombc2,
                  fix_formula = "cont_cov + cat_cov", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = F,
                  prv_cut = 0.1, lib_cut = 0, s0_perc = 0,
                  group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = F, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)


output_imp_res<-lapply(output_imp, function(x){x$res})

long_imp_res<-rbindlist(output_imp_res, idcol = T)

source("../_pool_manual.R")

## test continuous variable:
pooled_res<-pool_manual(result = long_imp_res, name_col = "taxon",imput_col = ".id",beta_col = "lfc_cont_cov", 
            sample_size = 150, se_col = "se_cont_cov", p_col = "p_cont_cov",
            alpha = 0.05, p_adjust = "holm")

res_merge<-merge(res_prim, pooled_res, by.x = "taxon", by.y = "name")
res_merge[,c("taxon","lfc_cont_cov", "q_cont_cov","pool.beta","pool.pvalue.fwe")]

# test categorical variable:

pooled_res2<-pool_manual(result = long_imp_res, name_col = "taxon",imput_col = ".id",beta_col = "lfc_cat_cov2", 
                        sample_size = 150, se_col = "se_cat_cov2", p_col = "p_cat_cov2",
                        alpha = 0.05, p_adjust = "holm")

res_merge<-merge(res_prim, pooled_res2, by.x = "taxon", by.y = "name")

res_merge[,c("taxon","lfc_cat_cov2", "q_cat_cov2","pool.beta","pool.pvalue.fwe")] %>%
  mutate(across(lfc_cat_cov2:pool.pvalue.fwe, round, 2))
