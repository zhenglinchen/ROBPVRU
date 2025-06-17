
> ROBust Phylogentic eigenVector Regression with Union eigenvector sets

---

## 🔍 Features

- `PVR_Robust()`: Conduct PVR with specified estimator and eigenvector selection method
- `Merge_PVR_Traits()`: Merge phylogenetic eigenvectors into traits dataframe
- `Union_EV()`: Get the union set of phylogenetic eigenvectors selected from X and Y


---

## ⚡ Quick Start

### Installation

```r
# Development version
remotes::install_github("zhenglinchen/ROBPVRU")
library(ROBPVRU)

```
### Usage
```r
# View example data
View(my_df)
View(my_tree)

# Merge phylogenetic eigenvectors into traits dataframe
Merge_PVR_Traits(my_tree,my_df,"species")

# Get the union set of phylogenetic eigenvectors selected from X and Y
Union_EV(my_df,my_tree,ID = "species",Y = "FG",X="SO",method = "ESRRV")
Union_EV(my_df,my_tree,ID = "species",Y = "FG",X="SO",method = "mMorI")
```

Do PVR robust regression with different eigenvector selection methods
```r
PVR_Robust(my_df,my_tree,"species","FG","SO",method = "mMorI",estimator = "L1")
PVR_Robust(my_df,my_tree,"species","FG","SO",method = "mMorI",estimator = "L2")
PVR_Robust(my_df,my_tree,"species","FG","SO",method = "mMorI",estimator = "S")
PVR_Robust(my_df,my_tree,"species","FG","SO",method = "mMorI",estimator = "M")
PVR_Robust(my_df,my_tree,"species","FG","SO",method = "mMorI",estimator = "MM")

PVR_Robust(my_df,my_tree,"species","FG","SO",method = "ESRRV",estimator = "L1")
PVR_Robust(my_df,my_tree,"species","FG","SO",method = "ESRRV",estimator = "L2")
PVR_Robust(my_df,my_tree,"species","FG","SO",method = "ESRRV",estimator = "S")
PVR_Robust(my_df,my_tree,"species","FG","SO",method = "ESRRV",estimator = "M")
PVR_Robust(my_df,my_tree,"species","FG","SO",method = "ESRRV",estimator = "MM")
```
