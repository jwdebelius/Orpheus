## Pool results using Rubin's rules
## 5-23-25
## Meredith Palmore

library(data.table)
library(tidyverse)

pool_manual<- function(result, name_col='name', imput_col='imputation',
                      sample_size, beta_col='beta', se_col = 'se',p_col = 'pval',
                      alpha = 0.05, p_adjust = 'holm'){
  # """
  #   Pools imputed MICE results using Rubin's rules
  # 
  #   Parameters
  #   ----------
  #   result : DataFrame
  #       A dataframe with a multi level index including the results.
  #       The frame should be indexed by a multiple level index which is 
  #       based on the `impute_level`, where the variable name is in the 
  #       `name_level` and the imputation identifier is in the 
  #       `impute_level` level. The dataframe should also include columns 
  #       summarizing the effect estimate (`beta`), standard error (`se_col`), and p-value (`p_col`).
  #   name_col: str
  #       The feature/taxon name column.
  #   impute_col: str
  #       The name of a column containing a dummy variable identifying 
  #       each unique imputation. 
  #   sample_size: num
  #       number of participants or samples included in the model.
  #   p_col: num
  #       The raw p-value before adjustment. 
  #   beta_col, se_col, p_col: str
  #       Column names in result file that describes the effect estimate
  #       (beta_col), standard error of the effect (se_col), 
  #       parametric p-value (p_col), and the degrees of freedom (df_col)
  #   p_adjust: str 
  #      p_adj_method character. method to adjust p-values. Default is "holm".
  #      ' Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  #      ' "fdr", "none". See \code{?stats::p.adjust} for more details.
  # 
  #   Returns
  #   -------
  #   data.frame 
  #   
  # 
  
  result<-result %>% rename(
    imputation =.data[[imput_col]],
    name=.data[[name_col]],
    beta=.data[[beta_col]],
    se = .data[[se_col]],
    pval=.data[[p_col]]
  )
  
  result <- data.table(result)
  
  # model parameters
  k<-result %>% filter(imputation==1) %>% nrow()
  
  # number of imputations
  m<-max(as.numeric(result$imputation))
  
  # calculate the pooled effect estimate, within and between variance
  pooled<-result%>% group_by(name) %>%
    summarise(pool.beta = sum(beta)/m, 
              Var_W = sum(se**2/m),
              Var_B = sum((beta - pool.beta)**2)/(m-1)
    ) %>% 
    ungroup()
  
  pooled_output = pooled %>%
    mutate(
      Var_F = (Var_B/m),
      # total variance
      Var_T = Var_W + Var_B + Var_F,
      # pooled standard error
      pool.SE = sqrt(Var_T),
      # pooled wald stat
      pool.Wald = pool.beta/pool.SE,
      # lambda value which is necessary to calculate the degrees of freedom
      lambda = (Var_B + (Var_B/m))/Var_T,
      # old degrees of freedom
      df_old = (m-1)/(lambda**2),
      # observation degrees of freedom
      df_obs = ((sample_size-k + 1)/(sample_size-k + 3))*(sample_size-k)*(1-lambda),
      #adjusted degrees of freedom
      df_adj = (df_old*df_obs)/(df_old+df_obs),
      # pooled p-value
      pool.pvalue = 2*pt(q=abs(pool.Wald),df=df_adj, lower.tail = F),
      # family-wise error adjusted pooled p-value
      pool.pvalue.fwe = p.adjust(pool.pvalue, method = p_adjust),
      # pooled lower CI
      pool.lwr = pool.beta  - (qt(p=alpha/2, df=df_adj,lower.tail=F)*pool.SE),
      # pooled upper CI
      pool.upr = pool.beta + (qt(p=alpha/2, df=df_adj,lower.tail=F)*pool.SE)
    )
  
  
  # Gets the fmi and rvi
  missing_info<-calc_fmi_rvi(pooled_output$Var_B, pooled_output$Var_F, 
                             pooled_output$Var_T, pooled_output$Var_W, 
                             m)
  
  pooled_output_fmi<-data.frame(cbind(pooled_output,missing_info))
  
  
  summary_tab<-pooled_output_fmi %>%select(name, pool.beta, pool.SE, pool.Wald, 
                                           pool.pvalue,pool.lwr, pool.upr, 
                                           pool.pvalue.fwe, df_adj, fmi, rvi) %>% 
    arrange(pool.pvalue.fwe) %>%
    rename_with(~paste0(., ".",beta_col), -c("name")) %>%
    rename(name_col=name)
    
  
  return(summary_tab)
  
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
  
  tab<-cbind(fmi, rvi, df_r)
  
  return(tab)
}
