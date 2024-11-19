# veneeR

## This package was developed to denoise stem data and develop veneer rolls, optimised for maximal wood usage, in the stem until crown base height.

## 

## Get it running

Please install the package from github and load the required packages.

``` r
install.packages('devtools')
devtools::install_github('https://github.com/elenaWWD/veneeR')
library(veneeR)
```

## Introduction

The package is written to process the 3D point cloud .las file in a step-wise manner. It starts with reading in the .las file and the mandatory information on crown base height.

``` r
#data path: 
folder_out = "veneer_results"

# create direction, where the output is saved at:
dir.create(folder_out)
setwd(folder_out)

# required input data:

#1. read in manually derived single tree information (DBH, tree height, CBH, coordinates) of you point cloud:
laser_segmented_tree_info_test = cbind(read.csv(paste0(path_in, folder_out,"/laser_segmented_tree_info_test.csv"))) 

#set the coordinate of your field to remove high values:
x_move = mean(laser_segmented_tree_info_test$X)
y_move = mean(laser_segmented_tree_info_test$Y)

#2. get all the tree IDs out:
all_TID =  laser_segmented_tree_info_test$TreeID

#test tree ID: 
TID = 1
cbh_tree_loop = laser_segmented_tree_info_test[which(laser_segmented_tree_info_test$TreeID == TID),]

#create data.frame:
final_tree_information = data.frame(
  TreeID = TID, 
  radius_at_bh = cbh_tree_loop$Radius, 
  tree_height= cbh_tree_loop$h_man, 
  mean_x = cbh_tree_loop$X, 
  mean_y = cbh_tree_loop$Y, 
  cbh = cbh_tree_loop$cbh)

```
