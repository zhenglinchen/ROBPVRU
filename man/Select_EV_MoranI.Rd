% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Select_EV_MoranI.R
\name{Select_EV_MoranI}
\alias{Select_EV_MoranI}
\title{Select_EV_MoranI}
\usage{
Select_EV_MoranI(
  x,
  phy = my_tree,
  trait = NULL,
  sig = TRUE,
  sig.t = 0.05,
  MI.t = 0.05
)
}
\arguments{
\item{pvr_df}{merged data frame,output of Merge_PVR_Traits}

\item{my_tree}{phylogenetic tree}
}
\description{
Select eigenvectors using mMorI method
}
