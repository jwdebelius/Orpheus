## Pool Adonis-2 results using Fisher Z-transformation and Rubin's rules
library(data.table)
library(tidyverse)

pool_adonis<- function(adonis, name_col='name', imput_col='imputation',var_order='var_order',model_level=NULL,
                       r2_col='r2', p_method='mean', df_col='df',p_col = 'pval',alpha = 0.05){
  # """
  #   Pools Adonis imputed MICE results using Harel's R2 pooling
  # 
  #   Parameters
  #   ----------
  #   adonis_res : DataFrame
  #       A dataframe with a multi level index including the adonis results.
  #       The frame should be indexed by a multiple level index which is 
  #       based on the `index_col` plus the `impute_level`, where the 
  #       variable name is in the `name_level` and the imputation identifier 
  #       is in the `impute_level` level.
  #       The dataframe should also include columns summarizing the 
  #       correlation coeffecient (`r2_col`), p-value (`p_col`), and 
  #       degrees of freedom (`df_col`).
  #   ci_alpha: float, (0, 1)
  #       The critical value to calculate the confidence interval for R2
  #       values
  #   p_method: {'mean', 'median', 'max'}
  #       The method for pooling the p-value across imputations. Options
  #       are mean, median, and max.
  #   name_col: str
  #       The level name for the index portion identifying each unique
  #       variable. This must be part of the index_levels list.
  #   var_order: str
  #       And identifier for a model or equation where the number of
  #       observations is contingent on the model. Relevant only when
  #       `adonis` contains observations from multiple models with different
  #       sample sizes. If this is included, it must be part of the 
  #       index_levels list.
  # 
  #   impute_col: str
  #       The name of a column containing a dummy variable identifying 
  #       each unique imputation. 
  # 
  #   model_level: str
  #       And identifier for a model or equation where the number of
  #       observations is contingent on the model. Relevant only when
  #       `adonis` contains observations from multiple models with different
  #       sample sizes. If this is included, it must be part of the 
  #       index_levels list.
  # 
  #   df_col, r2_col, p_col: str
  #       Column names in adonis_res file that describes the variance 
  #       explained (r2_col), permutative p-value (p_col), and the degrees
  #       of freedom (df_col)
  # 
  #   Returns
  #   -------
  #   DataFrame
  #   
  # 
  #   References
  #   ----------
  #   [1] Harel (2009). "The estimation of R2 and adjusted R2 in incomplete 
  # data sets using multiple imputation". Journal of Applied Statitics. 
  #       36:1109. doi: 10.1080/02664760802553000
  #   """
  
  
  adonis<-adonis %>% rename(
    imputation =.data[[imput_col]], 
    name=.data[[name_col]], 
    var_order=.data[[var_order]], 
    r2=.data[[r2_col]],
    pval=.data[[p_col]],
    df=.data[[df_col]]
  )
  
  if(!is.null(model_level)){
    adonis<-adonis %>% 
      rename(model_level=.data[[model_level]]) %>%
      mutate(var_order= paste(var_order, model_level,sep="_"))
  }
  
  adonis <- data.table(adonis)
  
  num_obs <- get_num_obs(adonis)
  
  adonis_nototal<- adonis %>% filter(.data[[name_col]]  != "Total")
  
  # Checks the p method
  `%nin%`<-Negate(`%in%`)
  
  p_opts <- c('mean', 'max', 'median')
  if(p_method %nin% p_opts){
    stop(paste0('The p_method is ',p_method,'.\nAcceptable values are "mean", "median", and "max".\nPlease supply a correct option for p_method.'))
  }else{
    pfunc <- get(p_method)
    p_pool<-adonis_nototal %>% group_by(name, var_order) %>% summarise(p = pfunc(pval, na.rm = T))
  }
  
  
  long_dat <- fisher_z(adonis_nototal, num_obs, r2_col)
  
  # Pivots the variable to longer so we have each column as an imputed
  # transformed value
  
  
  q_wide <- long_dat[["q"]] %>% pivot_wider(id_cols=c('name', 'var_order'),
                                            names_from ='imputation',
                                            values_from=q)
  v_wide <- long_dat[["v"]] %>% pivot_wider(id_cols  =c('name', 'var_order'),
                                            names_from ='imputation',
                                            values_from =v)
  
  # Calculates the pooled estimate
  `%nin%` <- Negate(`%in%`)
  
  m<- sum(names(q_wide) %nin% c('name', 'var_order'))
  
  # Calculates the pooled q estimate
  q_pool<-apply(q_wide[,names(q_wide) %nin% c('name', 'var_order')],1, pool_est)
  
  
  
  # Calculates the pooled variance estimate
  q_wide_nokey<-q_wide[,names(q_wide) %nin% c('name', 'var_order')]
  v_wide_nokey<-v_wide[,names(v_wide) %nin% c('name', 'var_order')]
  
  v_wide_nokey
  
  variances <- invisible(pool_variance_simp(q_wide_nokey, q_pool, v_wide_nokey, m))
  
  bse_t <-sqrt(variances$var_t)
  
  
  # Gets the fmi and rvi
  missing_info<-calc_fmi_rvi(variances$var_b, variances$var_f, variances$var_t, variances$var_w, m)
  
  
  # Calculates the CI on the q-value
  q_width <- qnorm(1 - alpha / 2) * bse_t
  
  q_ci_lo <-  q_pool - q_width
  q_ci_hi <-  q_pool + q_width
  
  # Extracts the degrees of freedom for each model
  
  dof <- adonis_nototal %>% 
    filter(imputation == 1) %>% 
    group_by(name, var_order) %>%
    summarise(dof_obs=first(df), .groups = "drop")
  
  
  dof<- merge(data.frame(name = q_wide$name, var_order = q_wide$var_order), dof, by = c('name', 'var_order'), sort = FALSE)
  
  
  dof_br <- pool_dof(missing_info[[2]], missing_info[[1]], m, dof$dof_obs)
  
  
  ## i am getting CIs that do not contain the estimate for low R^2
  
  dof_p<-merge(dof, p_pool, by = c('name', 'var_order'))
  
  # Summarizes the results
  r2_summary <- data.frame('name' = q_wide[, c("name")],
                           'var_order' = q_wide[, c("var_order")],
                           'r2' = release_z(q_pool),
                           'r2_ci_lo' =  release_z(q_ci_lo),
                           'r2_ci_hi' =  release_z(q_ci_hi),
                           'fmi'=  missing_info[[1]],
                           'rvi'= missing_info[[2]],
                           'dof_br'=dof_br
                           
  )
  
  
  r2_summary_m<-merge( r2_summary, dof_p, by = c("name","var_order"))
  
  r2_summary_m<-r2_summary_m %>% select(name, var_order, dof_obs, dof_br, r2, r2_ci_lo, r2_ci_hi, p, fmi, rvi) %>%
    mutate(model = as.numeric(ifelse(str_detect(pattern = "_", var_order), str_split_i(pattern = "_",var_order,2), 1))) %>%
    arrange(model, var_order) %>%
    data.table()
  
  ## if Lower CI > R^2, set LCI to 0 
  r2_summary_m[r2_ci_lo>r2, r2_ci_lo:=0]
  
  r2_summary_m<-as.data.frame(r2_summary_m)
  
  return(r2_summary_m)
  
}

fisher_z<- function(adonis_nototal,  num_obs, r2_col){
  # """
  # Performs a fisher's Z transform for R2 correlation data
  # 
  # Parameters
  # ----------
  # r2: pd.Series
  #     The correlation coeffecients for the analysis
  # num_obs: int
  #     The number of observations in the orginal data set
  # 
  # Returns
  # pd.Series
  #     The z-transformed correlation coeffecient
  # pd.Series
  #     The standard deviation in the z-transformed correlation coeffecient
  # 
  # References
  # ----------
  # [1] Harel (2009). "The estimation of R2 and adjusted R2 in incomplete 
  #     data sets using multiple imputation". Journal of Applied Statitics. 
  #     36:1109. doi: 10.1080/02664760802553000
  # [2] Fisher Z transform. (2 June 2024) In Wikipedia.
  #     https://en.wikipedia.org/wiki/Fisher_transformation 
  # 
  # Also See
  # --------
  # _throwback_z
  # """
  
  adonis_nototal <-adonis_nototal %>% mutate(ra = sqrt(.data[[r2_col]]))
  adonis_q <- adonis_nototal %>% mutate(q = 0.5 * log((1 + ra) / (1 - ra)))
  
  adonis_v<-num_obs %>% mutate(v = (1 / sqrt(n- 3))^2)
  
  
  q_v<-list()
  q_v[['q']]<-adonis_q
  q_v[['v']]<-adonis_v
  
  return(q_v)
}




pool_est<- function(params){
  # """
  # Pools the parameter estimate
  # 
  # References
  # ----------
  # eq. 7.12 ref [1]_
  # pg 9 ref [3]_
  # """
  pooled<-mean(params)
  return(pooled)
}

pool_variance_simp<- function(params, coef_pool,variance, m){
      
    # """
    # Pools the variance estimates, not variance/covariance
    # """
    var_w <-apply(as.matrix(variance), 1, mean)
    
    diffs<-apply(params, MARGIN = 2, function(x){x - coef_pool})
    
    
    if("matrix" %in% class(diffs)){
      vcov_b <-(diffs %*% t(diffs)) / (m - 1)
      var_b <- diag(vcov_b)
    }else{
      var_b <- sum((diffs * t(diffs))) / (m - 1) 
    }
    
    var_f <- var_b / m
    var_t <- var_b + var_w + var_f
    
    results<-list()

    results[["var_w"]]<-var_w
    results[["var_b"]]<-var_b
    results[["var_f"]]<-var_f
    results[["var_t"]]<-var_t

    return(results)
}

calc_fmi_rvi<-function(var_b, var_f, var_t, var_w, m){
  # """
  # Calculates the fraction of missing data and relative increase in variance
  # 
  # References
  # ----------
  # [1]
  # [2]
  # """
  fmi <- (var_b + var_f) / var_t
  df_r <- (m - 1) * (1 + 1 / (fmi^2))
  rvi <- ((var_b + var_f) / var_t) * (df_r + 1) / (df_r + 3) + 2 / (df_r + 3)
  
  return(list(fmi, rvi, df_r))
}

pool_dof<-function(riv, fmi, m, dof_com){
  # """
  # Calculates the pooled degrees of freedom
  # 
  # References
  # ----------
  # [1] https://www.daviddisabato.com/blog/2021/2/13/analyzing-and-pooling-results-from-multiply-imputed-data
  # """
  # dof_com <- dof_res * (dof_res + 1)/(dof_res + 3)*(1 - riv)
  
  v <- (m - 1) / fmi^2
  vhat <- (1 - fmi) * ((dof_com + 1) / (dof_com + 3)) * dof_com
  dof_br <- (v^(-1) + vhat^(-1))^(-1)
  
  
  return(dof_br)
  
}


release_z<-function(q){
  # """
  # Undoes the Fisher's Z transformation
  # 
  # Parameters
  # ----------
  # q: pd.Series
  #     The fisher z-transforemd correlation coeffecient
  # 
  # Returns
  # -------
  # pd.Series
  #     The untransfromed correlation coeffecient values    
  # 
  # References
  # ----------
  # [1] Harel (2009). "The estimation of R2 and adjusted R2 in incomplete 
  #     data sets using multiple imputation". Journal of Applied Statitics. 
  #     36:1109. doi: 10.1080/02664760802553000
  # [2] Fisher Z transform. (2 June 2024) In Wikipedia.
  #     https://en.wikipedia.org/wiki/Fisher_transformation 
  # 
  # Also See
  # --------
  # _fisher_z
  # """
  ra <- tanh(q)
  r2 <- ra^2
  
  return(r2)
  
}


get_num_obs<- function(adonis){
  # """
  # Gets the number of observations
  # """
  df_<-adonis[name == 'Total', ]
  
  df_$n<-df_$df+1
  
  ## get rid of the unnessecary cols
  
  
  keys<-c("name","var_order","imputation")
  
  ## merge the index rows with n variables
  
  adonis_m<-merge(adonis[,..keys], df_,allow.cartesian=TRUE)
  adonis_m<-adonis_m %>% distinct()
  
}

