#' Union_EV
#'
#' Get the union eigenvectors sets of X and Y
#'
#' @param my_df traits data frame
#' @param my_tree phylogenetic tree
#' @param ID Column name containing species names
#' @param Y Column name of dependent variable
#' @param X Column name of independent variable
#' @param method Eigenvector selection method, "ESRRV" or "mMorI"
#' @export

Union_EV <- function(my_df,my_tree,ID,Y,X,method){

  x <- PVR::PVRdecomp(my_tree)
  pvrs <- x@Eigen$vectors

  rownames(pvrs) <- my_tree$tip.label
  rownames(my_df) <- my_df[,ID]

  pvr_df <- merge(my_df,pvrs,by = "row.names")
  rownames(pvr_df) <- pvr_df$Row.names
  pvr_df <- pvr_df[, -1]

  colnames(pvrs) <- paste0("c",c(1:dim(pvrs)[2]))
  c_num <- length(my_tree$tip.label)-1

  if (method == "ESRRV"){
    pvr_name_Y <- Select_EV_ESRRV(Y,pvr_df,my_tree)
    pvr_name_X <- Select_EV_ESRRV(X,pvr_df,my_tree)
    pvr_XY_union <- union(pvr_name_X,pvr_name_Y)
  }else if (method == "mMorI"){
    pvr_name_Y <- Select_EV_MoranI(x,trait = pvr_df[,Y],phy = my_tree)
    pvr_name_X <- Select_EV_MoranI(x,trait = pvr_df[,X],phy = my_tree)
    pvr_XY_union <- union(pvr_name_X,pvr_name_Y)
  }

  return(pvr_XY_union)

}
