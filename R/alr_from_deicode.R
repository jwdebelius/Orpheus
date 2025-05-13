## Author: Meredith Palmore
## Date: 4/30/25
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
  # = 1 if feature x>0
  features_pa<-sapply(features[,-c(1)], function(x){ifelse(x>thresh, 1,0)})
  
  # re-add the features column
  features_pa<-as.data.frame(cbind(features[,1], features_pa))
  names(features_pa)[1]<-"FeatureID"
  
  return(features_pa)
}


cumsum_feat<-function(features, feature_pa, coords_sorted){
  # features_pa is the presence-abscence table
  # coords_sorted is the coordinates sorted by PCs
  
  # Extract feature order
  order<-unlist(coords_sorted[,1])
  
  # arrange presence absence and features dataset by order in the sorted coordinates file
  pres_abs_sort<-feature_pa %>%
    dplyr::arrange(factor(FeatureID, 
                          levels = order))
  
  feat_sort<-features %>%
    dplyr::arrange(factor(FeatureID, 
                          levels = order))
  
  # column-wise (Sample) cumulative sum
  cumsum_feat_samp<-apply(pres_abs_sort[,-c(1)], 2, cumsum)
  
  # Take min across row (samples)
  row_index<-which(apply(cumsum_feat_samp,1, min)==1) %>% min
  
  # Take all of the rows up to (and including) where the first row’s minimum num 
  # of OTUs present equals 1
  selected_feats<-feat_sort[1:row_index,]
  
  ## return list of feature IDs instead
  return(selected_feats)
}


alr<-function(features,coords, pc, featurcol="FeatureID", thresh=0){
  # Main function, sourcing the others above
  
  # 'coords' is a df or tibble, rows are features, columns are PCs, first column should be the feature name
  # 'features' is a a df or tibble, rows are features, columns are samples, first column should be the feature name
  # 'pc' is the column name of the PC of interest
  # 'featurcol' is the column name of the Feature IDs in features and coordinates table
  # 'thresh' is the minimum count of features required to be "present". Default is 0. 
  #  We assume at least two features in the numerator and denominator
  
  check_tables(features,coords, featurcol)
  
  coords_asc<-sort_pc_asc(coords, pc=pc)
  coords_desc<-sort_pc_desc(coords, pc=pc) 
  
  # convert to presence/absence
  
  pres_abs_feat<-pres_abs(features, thresh)
  
  ## Ascending PC1 
  asc_feats<-cumsum_feat(features = features, feature_pa = pres_abs_feat, coords_sorted = coords_asc)
  
  ## Descending PC1
  desc_feats<-cumsum_feat(features = features, feature_pa = pres_abs_feat, coords_sorted = coords_desc)
  
  # For each sample ; sum of PC1asc and PC1desc 
  ## actual feature table to compute ALRs -- subset to selected ascsenidng and descending features
  
  asc_feats_n<-sapply(asc_feats[,-c(1)], as.numeric)
  desc_feats_n<-sapply(desc_feats[,-c(1)], as.numeric)
  
  # Can assume is matrix because never going to have just one feature
  # If matrix, (more than one OTU selected) take the column-wise sum (sum OTUs present per sample). Else, keep the vector
  
  sampsums_pc_a<-apply(asc_feats_n,2, sum)

  sampsums_pc_d<-apply(desc_feats_n,2, sum)
  
  # take the ratio of the sum of features along the PCs
  ar<-(sampsums_pc_a/sampsums_pc_d)
  
  # Log2 of the ratio
  log2alr<-log(ar,
              base = 2)
  
  return(log2alr)
}


