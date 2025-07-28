# veneeR

<img align="right" width="259" height="300" alt="veneer" src="https://github.com/user-attachments/assets/a681c027-a940-46fe-bb00-af83a00982c5" />
</br>
</p>


## This package was developed to denoise stem data and develop veneer rolls, optimised for maximal veneer roll usage, in the stem until crown base height.

Have a look at our article for the explanation of the workflow and the veneer estimation results of the pilot study: Larysch, E., Frey, J., Schindler, Z. et al. Quantifying and mapping the ready-to-use veneer volume of European beech trees based on terrestrial laser scanning data. Eur J Forest Res (2025). https://doi.org/10.1007/s10342-025-01796-z


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

#fill in the final tree information based on the ransac results: 

#get the radius at breast height:
final_tree_information$radius_ransac_at_bh  = round(mean(radii$R [which(radii$Z > 1 & radii$Z < 1.5)], na.rm=T),4)
# get the radius at crown base height:
final_tree_information$radius_ransac_at_cbh = round(mean(radii$R [which(radii$Z > 0.9*final_tree_information$cbh)], na.rm=T),4)
# get the mitten diameter: 
final_tree_information$mitten_diameter = mean(radii$R [which (radii$Z > 0.9* ((max(radii$Z)- min(radii$Z))/2) & radii$Z< 1.1* ((max(radii$Z)- min(radii$Z))/2))], na.rm=T)
# get (if possible based on the tree height) diameter along the stem height:
if(final_tree_information$tree_height>6)  final_tree_information$radius_ransac_at_5  = round(mean(radii$R [which(radii$Z > 4.5 & radii$Z < 5.5)], na.rm=T),4)
if(final_tree_information$tree_height>10) final_tree_information$radius_ransac_at_10 = round(mean(radii$R [which(radii$Z > 9.5 & radii$Z < 10.5)], na.rm=T),4)
if(final_tree_information$tree_height>15) final_tree_information$radius_ransac_at_15 = round(mean(radii$R [which(radii$Z > 14.5 & radii$Z < 15.5)], na.rm=T),4)
if(final_tree_information$tree_height>20) final_tree_information$radius_ransac_at_20 = round(mean(radii$R [which(radii$Z > 19.5 & radii$Z < 20.5)], na.rm=T),4)

#move to next tree and safe tree dimensional information, if tree is too thin or stem length lower the minimum required stem length for Pollmeier: 
if(is.na(final_tree_information$radius_ransac_at_bh) ){ 
rm(tree_f, tls_all) 
write.csv(final_tree_information, "final_tree_information.csv")
gc() #make space for new calc
}
if(is.na(final_tree_information$radius_ransac_at_bh) ) next
if(c(final_tree_information$cbh - 0.3)  <1.8288 ) {
  rm(tree_f, tls_all)
  write.csv(final_tree_information, "final_tree_information.csv")
  gc()
}
if(c(final_tree_information$cbh - 0.3)  <1.8288 ) next #move to the next tree but keep in information table.




tree_f   = veneer_outlier(radii=radii, tree_f = tree_f, steps=0.05)
tree_new = veneer_card(tree_f = tree_f, radii = radii)
taper_vol_list = taper_volume(radii = radii)
tree_new_spline = veneer_spline(tree_new)

#the final veneer function: ##############################################################################################

tls_all_withcrown = lidR::readLAS(paste0(path_in, "test_1_with_crown.las"))

veneer_rolls = analyse_veneer_potential(pred_all = tree_new_spline, final_tree_information = final_tree_information, TID=TID)

```
