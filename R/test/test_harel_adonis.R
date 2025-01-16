
library(data.table)
library(tidyverse)

source("../_harel_adonis.R")

compare_digits<-function(x){
  x<-format(x, digits = 5,scientific = 5)
  x<-trimws(x)
  return(x)
}

setUp<- function(){
  self<-list()
  # The correlation testing and z-testing
  self[["fisher_r2"]] <- data.frame("r2" =0.25)
  self[["num_obs"]]  <- data.frame("n" =12)
  self["fisher_q"] <- 0.5493061443340549
  self["fisher_var"] <- 0.11111111111111
  
  imputed <- fread(
    paste0("../data", 
                 '/Example_imputed_dataset.tsv'),
    drop=1
  )
  self[["impute_res"]] <- imputed
  return(self)
}

test_pool_adonis_perr<-function(self){
  try(stop('The p_method is Jesper.\nAcceptable values are "mean", "median", and "max".\nPlease supply a correct option for p_method.'),silent = T)
  known <-geterrmessage
  try(pool_adonis( self[["impute_res"]], p_method='Jesper'),silent = T)
  res<-geterrmessage
  identical(known, res)
}

test_pool_adonis<-function(self){
  # known is same as Justines, but I add model = 1 to everything
  known <- t(data.frame(
    R1 = c(1.00000000e+00, 4.99810567e-01, 7.84331842e-03,
                    1.29042888e-03, 4.42434149e-02, 2.80000000e-03,
                    3.78857946e-04, 3.78889830e-04,1),
    R2=c(1.00000000e+00, 4.99635004e-01, 2.57162720e-02,
         1.36805651e-03, 7.77982566e-02, 1.00000000e-03,
         7.29961663e-04, 7.30079986e-04,1),
    R3=c(2.00000000e+00, 1.19852576e+00, 1.25530568e-02,
         1.50620749e-04, 5.42624016e-02, 3.60000000e-03,
         1.22833636e-03, 1.22867124e-03,1),
    R4 =c(3.00000000e+00, 1.99937671e+00, 1.19279284e-02,
                    2.27267347e-04, 5.29836414e-02, 3.65400000e-01,
                    3.11624567e-04, 3.11646141e-04,1),
    R5 =c(3.00000000e+00, 1.99892407e+00, 1.91005582e-02,
                    2.05947215e-04, 6.65030355e-02, 1.30000000e-03,
                    5.37899838e-04, 5.37964100e-04,1),
    R6 = c(2.39000000e+02, 2.36471917e+02, 9.22849967e-01,
           9.02063747e-01, 9.39371134e-01,NaN,
           2.20509005e-03, 2.20616820e-03,1)
    ))
  
  colnames(known)<-c('dof_obs', 'dof_br', 'r2', 'r2_ci_lo', 'r2_ci_hi','p', 'fmi', 'rvi', 'model')
  rownames(known)<-NULL
  known<-as.data.frame(known)
  
  known$name <-c("age", "ethnicity","bmi_class","alcohol_use","birth_country","Residual")
  known$var_order <-c("1", "2","3","4","5","6")
  known<-known[,c("name", "var_order",'dof_obs', 'dof_br', 'r2', 'r2_ci_lo', 'r2_ci_hi','p', 'fmi', 'rvi', 'model')]
  
  test <-pool_adonis(self[["impute_res"]], 
                   name_col ='name', var_order = 'var_order')
  ## are the first 7 digits the same?
  test_s<-sapply(test, compare_digits) # correct up to 7 significant figs?
  known_s<-sapply(as.data.frame(known), compare_digits) # correct up to 7 significant figs?
 
  identical(known_s, test_s)
}

test_fisher_z<-function(self){
  
  known_q<-sapply(data.frame("r2" = self[["fisher_r2"]], 
                    "ra" = 0.5, 
                    "q" = self[["fisher_q"]]), compare_digits) 
  known_v<-sapply(data.frame("n" = self[["num_obs"]], 
                      "v" = self[["fisher_var"]]), compare_digits) 
  
  z_out<- fisher_z(self[["fisher_r2"]], self[["num_obs"]], "r2")
  z_out_s<-lapply(z_out, function(x){sapply(x, compare_digits)}) 
  
  print(identical(known_q, z_out_s[["q"]]))
  print(identical(known_v, z_out_s[["v"]]))
}


test_release_z<-function(self){
  test_r2 <- release_z(self[["fisher_q"]])
  
  test_r2_s <- compare_digits(test_r2)
  known_s<-deframe(self[["fisher_r2"]]) %>% compare_digits
  
  identical(known_s, test_r2_s)
}


test_pool_variance_simp<-function(self){
    var <- t(data.frame(
        'A'=rep(0.00404858, 10)
      ))
    param <- t(data.frame(
      'A'=c(0.08973622, 0.08796234, 0.08785097, 0.09045728,
                      0.09035581, 0.08890988, 0.08730415, 0.08727498,
                      0.08937271, 0.08872715)
    ))
    coef_pool<-apply(param, 1,pool_est)
    
    # known result: 
    k_var_w = 0.00404858
    k_var_b = 1.39492858e-06
    k_var_f = 1.39492858e-07
    k_var_t = 0.00405012
    
    k_var_w_s<-compare_digits(k_var_w)
    k_var_b_s<-compare_digits(k_var_b)
    k_var_f_s<-compare_digits(k_var_f)
    k_var_t_s<-compare_digits(k_var_t)

    var <- pool_variance_simp(param, coef_pool, var, 10)
    var_s<-lapply(var, as.numeric)
    
    var_s<-lapply(var_s, compare_digits)

    print(identical(k_var_w_s, var_s[["var_w"]]))
    print(identical(k_var_b_s, var_s[["var_b"]]))
    print(identical(k_var_f_s, var_s[["var_f"]]))
    print(identical(k_var_t_s, var_s[["var_t"]]))
    
}

self<-setUp()

test_pool_adonis_perr(self)

test_pool_adonis(self)

test_fisher_z(self)

test_release_z(self)

test_pool_variance_simp(self)


