#' Select_EV_ESRRV
#'
#' Select eigenvectors using ESRRV
#'
#' @param pvr_df merged data frame,output of Merge_PVR_Traits
#' @param Y the traits used to select eigenvectors
#' @param my_tree phylogenetic tree

Select_EV_ESRRV <- function(Y,pvr_df,my_tree){
  c_num <- length(my_tree$tip.label)-1
  EV_vec <- NULL
  ####### calculate corelate PVRs #########
  for (c in c(1:c_num)){
    col_name <- paste("c",c,sep="")
    norm_p_y <- shapiro.test(pvr_df[,Y])$p.value
    norm_p_EV <- shapiro.test(pvr_df[,col_name])$p.value
    if (norm_p_y>=0.05 & norm_p_EV>=0.05){
      fit <- cor.test(pvr_df[,Y],pvr_df[,col_name],method="pearson")
    }else{
      fit <- cor.test(pvr_df[,Y],pvr_df[,col_name],method="spearman")
    }
    if (fit$p.value<0.05){
      EV_vec <- c(EV_vec,col_name)
    }
  }
  return(EV_vec)
}
