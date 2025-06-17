#' PVR_Robust
#'
#' Apply PVR using union eigenvector sets and robust esimators 
#'
#' @param my_df traits data frame
#' @param my_tree phylogenetic tree
#' @param ID Column name containing species names
#' @param Y Column name of dependent variable
#' @param X Column name of independent variable
#' @param method Eigenvector selection method, "ESRRV" or "mMorI", default is "mMorI"
#' @param estimator Estimator used in regression, "L1","L2","S","M" or "MM", default is "L2"
#' @export

PVR_Robust <- function(my_df,my_tree,ID,Y,X,method="mMorI",estimator="L2"){
  
  x <- PVR::PVRdecomp(my_tree)
  pvrs <- x@Eigen$vectors
  
  rownames(pvrs) <- my_tree$tip.label
  rownames(my_df) <- my_df[,ID]
  
  pvr_df <- merge(my_df,pvrs,by = "row.names")
  rownames(pvr_df) <- pvr_df$Row.names
  pvr_df <- pvr_df[, -1]  
  
  if (method == "ESRRV"){
    pvr_name_Y <- Select_EV_ESRRV(Y,pvr_df,my_tree)
    pvr_name_X <- Select_EV_ESRRV(X,pvr_df,my_tree)
    pvr_XY_union <- union(pvr_name_X,pvr_name_Y)
  }else if (method == "mMorI"){
    pvr_name_Y <- Select_EV_MoranI(x,trait = pvr_df[,Y],phy = my_tree)
    pvr_name_X <- Select_EV_MoranI(x,trait = pvr_df[,X],phy = my_tree)
    pvr_XY_union <- union(pvr_name_X,pvr_name_Y)
  }
  
  pvr_XY_npar <- length(pvr_XY_union)
  
  xy_pvr_formula <- NULL
  
  if (pvr_XY_npar == 0){
    xy_pvr_formula <- as.formula(paste0(Y,"~",X))
    
  }else{
    consider_xy_pvr <- paste(sort(pvr_XY_union),collapse="+")
    xy_pvr_formula <- as.formula(paste0(Y,"~",X,"+",consider_xy_pvr))
  } 
  
  if (estimator == "L2"){
    pvr_xy_model <- lm(xy_pvr_formula,data=pvr_df)
  } else if (estimator == "L1"){
    pvr_xy_model <- lad(xy_pvr_formula,data=pvr_df, method = "EM")
  } else if (estimator == "S") {
    pvr_xy_model <- lmRob(xy_pvr_formula,data=pvr_df,estim="Initial")
  } else if (estimator == "M") {
    pvr_xy_model <- rlm(xy_pvr_formula, data=pvr_df, method = "M", scale.est = "Huber", k2 = 1.345)
  } else if (estimator == "MM") {
    pvr_xy_model <- lmRob(xy_pvr_formula,data=pvr_df)
  }
  
  return(pvr_xy_model)
  
}