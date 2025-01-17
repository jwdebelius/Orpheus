library(data.table)

var_mod <- data.frame(
  'r2'=rep(0.0017775014, 10),
  'name'=rep("sex",10),
  'df'=rep(1,10),
  'imp'=1:10,
  'pval'=c(0.901,0.902,0.875,0.891,0.886,0.904,0.895,0.884,0.895,0.88),
  "var_order"=rep(1,10)
)

var_resid <- data.frame(
  'r2'=rep(0.998222499, 10),
  'name'=rep("Residual",10),
  'df'=rep(257,10),
  'imp'=1:10,
  'pval'=rep(NA,10),
  "var_order"=rep(2,10)
)

var_total<- data.frame(
  'r2'=rep(1, 10),
  'name'=rep("Total",10),
  'df'=rep(258,10),
  'imp'=1:10,
  'pval'=rep(NA,10),
  "var_order"=rep(NA,10)
)

adonis_res<-data.table(rbind(var_mod,var_resid, var_total ))


source("../_harel_adonis.R")

lci_notincl<-pool_adonis(adonis_res, imput_col = "imp")


print(lci_notincl)

lci_notincl %>%
  filter(name == "sex") %>%
  ggplot(data=., aes(y=factor(name), x = as.numeric(r2), xmin=as.numeric(r2_ci_lo),xmax = as.numeric(r2_ci_hi))) +
    geom_point(position = position_dodge(width=0.5)) +
    geom_errorbar(position=position_dodge(width = 0.1), width = 0.3) +
    scale_y_discrete(limits = rev)+
    labs(x = expression(italic(R)^2), y = "") +
    geom_vline(xintercept = 0, color = 'black', linetype='dashed', alpha = 0.5, linewidth=0.5)+
    theme_classic() +
    guides(color=guide_legend(reverse = T)) 
    