#' Merge_PVR_Traits
#'
#' Merge eigenvectors to traits dataframe 
#'
#' @param my_tree phylogenetic tree
#' @param my_df traits data frame
#' @param ID Column name containing species names
#' @export

Merge_PVR_Traits <- function(my_tree,my_df,ID){
  x <- PVR::PVRdecomp(my_tree)
  pvrs <- x@Eigen$vectors
  
  rownames(pvrs) <- my_tree$tip.label
  rownames(my_df) <- my_df[,ID]
  
  pvr_df <- merge(my_df,pvrs,by = "row.names")
  rownames(pvr_df) <- pvr_df$Row.names
  pvr_df <- pvr_df[, -1]  
  
  return(pvr_df)
}