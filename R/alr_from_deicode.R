## Author: Meredith Palmore
## Date: 5/13/25
## Purpose: Calculate additive log ratios for presence of features, ranking features along PC coordinates


library(data.table)
library(dplyr)

## Check that the sample names of the coordinate and feature tables match
## and that the Feature ID is in the first column in both dataset

check_tables<-function(features, coords, featurcol= "FeatureID"){
  if(names(features)[1]!=featurcol){
    stop("Feature IDs should be the first column in feature table")
  }
  if(names(coords)[1]!=featurcol){
    stop("Feature IDs should be the first column in coords table")
  }
  if(!identical(features[,featurcol],coords[,featurcol])){
    stop("Order of Feature IDs in feature table should match order of Feature IDs in coordinates")
  }
}



## add PC#

## Sample | PC1 | PC2 

sort_pc_asc<-function(coords, pc){
  # Sort along PC in an ascending fashion
  coords<-data.table(coords)
  
  #sort
  setorderv(coords, pc, 1) # ascending 
  return(coords)
}


sort_pc_desc<-function(coords, pc){
  # Sort along PC in an descending fashion
  coords<-data.table(coords)
  
  # sort
  setorderv(coords, pc, -1) # descending 
  return(coords)
}

# Convert to presence/absence 

pres_abs<-function(features, thresh){
  # = 1 if feature is >=thresh
  features_pa<-sapply(features[,-c(1)], function(x){ifelse(x>=thresh, 1,0)})

  return(features_pa)
}


cummax_feat<-function(features, coords_sorted, prop_samp, thresh){
  # features_pa is the presence-abscence table
  # coords_sorted is the coordinates sorted by PCs
  
  # Extract feature order
  order<-unlist(coords_sorted[,1])
  
  # arrange presence absence and features dataset by order in the sorted coordinates file
  feat_sort<-features %>%
    dplyr::arrange(factor(FeatureID, 
                          levels = order))
  
  #convert to presence absence
  pres_abs_sort<- pres_abs(features = feat_sort, thresh = thresh)
  
  # print(pres_abs_sort)
  
  # column-wise (Sample) cumulative max
  cummax_feat_samp<-apply(pres_abs_sort, 2, cummax)
  
  # Take sample mean across features
  row_index<-which(apply(cummax_feat_samp,1, mean, na.rm = T)>=prop_samp) %>% min
  
  # Take all of the rows up to (and including) where the first row’s minimum num 
  # of OTUs present equals 1
  selected_feats<-feat_sort[1:row_index,]
  
  ## return list of feature IDs instead
  return(selected_feats)
}


alr<-function(features,coords, pc, featurcol="FeatureID", thresh=1, prop_samp=1){
  # Main function, sourcing the others above
  
  # 'coords' is a df or tibble, rows are features, columns are PCs, first column should be the feature name
  # 'features' is a a df or tibble, rows are features, columns are samples, first column should be the feature name
  # 'pc' is the column name of the PC of interest
  # 'featurcol' is the column name of the Feature IDs in features and coordinates table
  # 'thresh' is the minimum number of features to be counted as "present". Default is 1. 
  # 'prop_samp' is the proportion of samples that must have a feature present for the feature to contribute to the numerator, ranging from 0 to 1. Default is 1 (all samples). 
  #  We assume at least two features in the numerator and denominator
  
  if(thresh<=0){stop("Thresold for number of features to be counted as 'present' must be >0.")}
  
  check_tables(features,coords, featurcol)
  
  coords_asc<-sort_pc_asc(coords, pc=pc)
  coords_desc<-sort_pc_desc(coords, pc=pc) 
  
  
  ## Ascending PC1 
  asc_feats<-cummax_feat(features = features, coords_sorted = coords_asc, prop_samp = prop_samp, thresh = thresh)
  
  ## Descending PC1
  desc_feats<-cummax_feat(features = features, coords_sorted = coords_desc, prop_samp = prop_samp, thresh = thresh)
  
  ## warning if features are in num AND denom 
  mutual_features<-asc_feats[,featurcol][asc_feats[,featurcol]  %in% desc_feats[,featurcol]]
  
  if(length(mutual_features)!=0){
    warning("Features contributing to both the numerator and denominator:\n",paste0(mutual_features,"\n"))
  }
  
  # For each sample ; sum of PC1asc and PC1desc 
  ## actual feature table to compute ALRs -- subset to selected ascending and descending features
  
  asc_feats_n<-sapply(asc_feats[,-c(1)], as.numeric)
  desc_feats_n<-sapply(desc_feats[,-c(1)], as.numeric)
  
  # Can assume is matrix because never going to have just one feature
  
  sampsums_pc_a<-apply(asc_feats_n,2, sum)

  sampsums_pc_d<-apply(desc_feats_n,2, sum) 
  
  # take the ratio of the sum of features along the PCs
  ar<-(sampsums_pc_a + 1)/(sampsums_pc_d + 1) # add 1 as a buffer when the sum of features = 0
  
  # Log2 of the ratio
  log2alr<-log(ar,
              base = 2)
  
  return(log2alr)
}
