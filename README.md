# veneeR

## This package was developed to denoise stem data and develop veneer rolls, optimised for maximal veneer roll usage, in the stem until crown base height.

## 

## Get it running

Please install the package from github and load the required packages.

``` r
install.packages('devtools')
devtools::install_github('https://github.com/elenaWWD/veneeR')
library(veneeR)
```

## Data preparation

The package is written to process the 3D point cloud .las file in a step-wise manner. It starts with reading in the .las file and the mandatory information on crown base height.

``` r

#create a input folder:
path_in = NA

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

## Example of the workflow

``` r
#read in example point cloud: ########################################################################################### 
tls_loop = lidR::readLAS(paste0(path_in, "1.las"))
tls_all  = tls_loop [which(tls_loop@data$TreeID == TID),]

#try it: ###############################################################################################################
tree_f   = veneer_noise (tree_las = tls_all)
tree_f   = veneer_incli(tree_f = tree_f, plot=T)
radii    = veneer_ransac(tree_f, steps = 0.05, species = "FaSy")
tree_f   = veneer_outlier(radii=radii, tree_f = tree_f, steps=0.05)
tree_new = veneer_card(tree_f = tree_f, radii = radii)
taper_vol_list = taper_volume(radii = radii)
tree_new_spline = veneer_spline(tree_new)

#the final veneer function: ##############################################################################################

<<<<<<< HEAD
tls_all_withcrown = lidR::readLAS(paste0(path_in, "test_1_with_crown.las"))
=======
tls_all_withcrown = lidR::readLAS(paste0( "path_in, "test_1_with_crown.las"))
>>>>>>> 1b4999472905fad4b870888f7848021f2cefb924
veneer_rolls = analyse_veneer_potential(pred_all = tree_new)

```
