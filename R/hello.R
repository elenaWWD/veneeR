##### Functions ####

#get the local curvature points: 
inflect <- function(x, threshold = 1){
  up   <- sapply(1:threshold, function(n) c(x[-(seq(n))], rep(NA, n)))
  down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), x[-seq(length(x), length(x) - abs(n) + 1)]))
  a    <- cbind(x,up,down)
  list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
}

#not in:
'%!in%' <- function(x,y)!('%in%'(x,y))

#calculate angle:
angle2 <- function(M,N){ atan2(N[2],N[1]) - atan2(M[2],M[1]) }

#change angle to direction:
d2c.2 <- function(x) {
  upper <- seq(from = 11.25, by = 22.5, length.out = 17)
  card1 <- c('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW', 'N')
  ifelse(x>360 | x<0,NA,card1[findInterval(x,upper,rightmost.closed = T)+1])
}

# #color palette:
col2mtl <- function(col) { paste(paste(c("Ka","Kd","Ks"),paste(col2rgb(col)/255, collapse = " ")) , collapse = "\n")}
# pal=viridis::magma(16,begin=.1)
# pal_shuffle = sample(pal)
# for (c in 1:13) {
#   writeLines(paste0("newmtl c",c))
#   writeLines(col2mtl(pal_shuffle[c]))
# }

#export mesh:
write_single_obj <- function(pol, file, mtl_name = "Material1", mtllib = "colors.mtl"){
  # generate header information
  obj_head <- c("# File created with R write_obj", paste("#", Sys.Date()),paste("mtllib", mtllib), paste("usemtl", mtl_name), 'o object1')
  #check if its a list or just one polygon
  if(is.list(pol)){
    l <- 0 # index for vertex numbers already used
    lines <- character() # obj text lines
    p_i <- 1 # polygon index
    p = pol
    #g <- paste0("o polygon", p_i) # create a group for evry polygon
    v = apply(round(p$vb,5), 2, function(x) paste("v", paste(x[1:3], collapse = " "))) # add the coordinates round(p$vb,5)
    f = apply(p$ib+l, 2, function(x) paste("f", paste(x[1:4], collapse = " "))) # add the triangle coordinate order
    f_rev = apply(p$ib+l, 2, function(x) paste("f", paste(x[4:1], collapse = " "))) # add the reverse order for backsides
    l <- l+ncol(p$vb) # increase the vertex index
    lines <- c(lines, c(v,f,f_rev)) # add the lines
    p_i <- p_i + 1 # increase the polygon index
    
    writeLines(c(obj_head,lines), file) # write to file
  } else {
    v = apply(pol$vb, 2, function(x) paste("v", paste(x[1:3], collapse = " ")))
    f = apply(pol$it, 2, function(x) paste("f", paste(x[1:3], collapse = " ")))
    writeLines(c(obj_head,v,f), file)
  }
  
}

#export multiple mesh:
write_obj <- function(pol, file, mtl_name = "Material1"){
  # generate header information
  obj_head <- c("# File created with R write_obj", paste("#", Sys.Date()), paste("usemtl", mtl_name), paste("usemtl", mtl_name), 'o object1')
  #check if its a list or just one polygon
  l <- 0 # index for vertex numbers already used
  lines <- character() # obj text lines
  p_i <- 1 # polygon index
  for(pi in 1:length(pol)){
    
    p = pol[[pi]]
    #g <- paste0("o polygon", p_i) # create a group for evry polygon
    v = apply(round(p$vb,5), 2, function(x) paste("v", paste(x[1:3], collapse = " "))) # add the coordinates round(p$vb,5)
    f = apply(p$ib+l, 2, function(x) paste("f", paste(x[1:4], collapse = " "))) # add the triangle coordinate order
    f_rev = apply(p$ib+l, 2, function(x) paste("f", paste(x[4:1], collapse = " "))) # add the reverse order for backsides
    l <- l+ncol(p$vb) # increase the vertex index
    lines <- c(lines, c(v,f,f_rev)) # add the lines
    p_i <- p_i + 1 # increase the polygon index
  }
  
  writeLines(c(obj_head,lines), file) # write to file
  
  
}

# Helper function to remove outliers based on boxplot stats
remove_outliers <- function(disc_df) {
  # Outlier removal for X and Y
  out_x <- boxplot.stats(disc_df$X)$out
  out_x_ind <- which(disc_df$X %in% out_x)
  disc_o <- if (length(out_x_ind) != 0) disc_df[-out_x_ind, ] else disc_df
  
  out_y <- boxplot.stats(disc_o$Y)$out
  out_y_ind <- which(disc_o$Y %in% out_y)
  if (length(out_y_ind) != 0) disc_o <- disc_o[-out_y_ind, ]
  
  return(disc_o)
}

# Helper function to update the final tree information
update_tree_info <- function(radii, species) {
  final_tree_information$radius_ransac_at_bh <- round(mean(radii$R[radii$Z > 1 & radii$Z < 1.5], na.rm = TRUE), 4)
  final_tree_information$radius_ransac_at_cbh <- round(mean(radii$R[radii$Z > 0.9 * final_tree_information$cbh], na.rm = TRUE), 4)
  
  # Mitten diameter calculation
  final_tree_information$mitten_diameter <- mean(radii$R[radii$Z > 0.9 * ((max(radii$Z) - min(radii$Z)) / 2) & radii$Z < 1.1 * ((max(radii$Z) - min(radii$Z)) / 2)], na.rm = TRUE)
  
  # Radius at various stem heights
  if (final_tree_information$tree_height > 6) final_tree_information$radius_ransac_at_5 <- round(mean(radii$R[radii$Z > 4.5 & radii$Z < 5.5], na.rm = TRUE), 4)
  if (final_tree_information$tree_height > 10) final_tree_information$radius_ransac_at_10 <- round(mean(radii$R[radii$Z > 9.5 & radii$Z < 10.5], na.rm = TRUE), 4)
  if (final_tree_information$tree_height > 15) final_tree_information$radius_ransac_at_15 <- round(mean(radii$R[radii$Z > 14.5 & radii$Z < 15.5], na.rm = TRUE), 4)
  if (final_tree_information$tree_height > 20) final_tree_information$radius_ransac_at_20 <- round(mean(radii$R[radii$Z > 19.5 & radii$Z < 20.5], na.rm = TRUE), 4)
  
  # Clean up tree if not valid
  if (is.na(final_tree_information$radius_ransac_at_bh)) {
    rm(tree_f, tls_all)
    write.csv(final_tree_information, "final_tree_information.csv")
    gc()  # Clean up memory
  }
  if (is.na(final_tree_information$radius_ransac_at_bh)) next
  if ((final_tree_information$cbh - 0.3) < 1.8288) {
    rm(tree_f, tls_all)
    write.csv(final_tree_information, "final_tree_information.csv")
    gc()
  }
}

