
##### 1 noise removal: ####

#' Remove Noise and Outliers from Tree LiDAR Data
#'
#' This function processes LiDAR data from a tree scan by filtering out low-verticality points, 
#' high-planarity points, and neighbor noise, and by removing outliers in the X, Y, and Z dimensions.
#' It also resets the Z coordinates using the original Z reference data.
#'
#' @param tree_las A LAS object representing a tree scan. The object must include `Verticality`, `Planarity`, and `Zref` attributes.
#' @return A LAS object with noise and outlier points removed.
#' @examples
#' # Assuming 'tree_las' is a valid LAS object with Verticality, Planarity, and Zref attributes:
#' clean_tree <- veneer_noise(tree_las)
#' @export
#' 
#' @import lidR 
#' @import CspStandSegmentation 
#' @import tidyr

#' 
veneer_noise <- function(tree_las) {
  
  invisible(lapply(c('tidyr','CspStandSegmentation', 'lidR'), require, character.only = TRUE))
  
  # Add geometry metrics to the LAS data if not already present
  tls_all <- CspStandSegmentation::add_geometry(tree_las)
  
  # Apply filters: keep only points with high verticality and planarity, and remove classified noise
  tls_all <- tls_all %>%
    filter_poi(Verticality > 0.85) %>%
    filter_poi(Planarity > 0.4) %>%
    classify_noise(ivf(0.1, 6)) %>%
    filter_poi(Classification != LASNOISE)
  
  # Restore Z to its original value and normalize to start from zero
  tls_all@data$Z <- tls_all@data$Zref - min(tls_all@data$Zref)
  
  # Identify and filter outliers in the X and Y dimensions
  x_outliers <- tls_all@data$X %in% boxplot(tls_all@data$X, plot = FALSE)$out
  y_outliers <- tls_all@data$Y %in% boxplot(tls_all@data$Y, plot = FALSE)$out
  tls_all <- tls_all[!(x_outliers | y_outliers), ]
  
  # Optionally filter by CBH if final_tree_information exists
  if (exists("final_tree_information") && "cbh" %in% names(final_tree_information)) {
    tls_all <- tls_all[which(tls_all@data$Z < final_tree_information$cbh), ]
  }
  
  return(tls_all)
}


##### 2 correct stem inclination: ####

#' Correct Tree Inclination Using Stem Axis
#'
#' This function estimates and corrects the inclination of a tree point cloud by aligning it to a new coordinate system 
#' based on the tree's stem axis. It also provides a 2D plot comparison of original and corrected inclination.
#'
#' @param tree_f A LAS object containing tree scan data with `X`, `Y`, and `Z` coordinates.
#' @param plot Logical, if `TRUE`, generates and saves a 2D plot comparison of the original and corrected inclination.
#' @return A LAS object with inclination-corrected coordinates (`X_cor`, `Y_cor`, `Z_cor`).
#' @examples
#' # Assuming 'tree_f' is a valid LAS object:
#' corrected_tree <- veneer_incli(tree_f, plot = TRUE)
#' @export
#' @import ggplot2
#' @import viridis


veneer_incli <- function(tree_f, plot = TRUE) {
  invisible(lapply(c('ggplot2','viridis', 'ggpubr'), require, character.only = TRUE))
  
  # Extract and normalize point coordinates
  points <- as.matrix(tree_f@data[, c("X", "Y", "Z")])
  points <- sweep(points, 2, apply(points, 2, min))
  
  # Estimate stem direction using bottom and top sections of the tree
  tree_height <- max(points[,3])
  bottom_points <- points[points[,3] < 0.1 * tree_height, ]
  top_points <- points[points[,3] > 0.8 * tree_height & points[,3] < tree_height, ]
  
  # Calculate stem axis vector
  stem_axis <- colMeans(top_points) - colMeans(bottom_points)
  stem_axis <- stem_axis / sqrt(sum(stem_axis^2))
  
  # Define orthogonal vectors for the rotation matrix
  vec_orth <- c(0, -stem_axis[3], stem_axis[2])
  vec_norm <- c(
    vec_orth[2] * stem_axis[3] - vec_orth[3] * stem_axis[2],
    vec_orth[3] * stem_axis[1] - vec_orth[1] * stem_axis[3],
    vec_orth[1] * stem_axis[2] - vec_orth[2] * stem_axis[1]
  )
  
  # Normalize auxiliary vectors
  vec_orth <- vec_orth / sqrt(sum(vec_orth^2))
  vec_norm <- vec_norm / sqrt(sum(vec_norm^2))
  
  # Arrange rotation matrix
  rotation_matrix <- rbind(vec_norm, vec_orth, stem_axis)
  
  # Rotate each point in the dataset
  rotated_points <- t(rotation_matrix %*% t(points))
  
  # Generate and save plots if requested
  if (plot) {
    plot_data <- data.frame(
      X = points[,1], Y = points[,2], Z = points[,3],
      X_rot = rotated_points[,1] - mean(tree_f@data$X), 
      Y_rot = rotated_points[,2] - mean(tree_f@data$Y), 
      Z_rot = rotated_points[,3]
    )
    
    p1 <- ggplot(plot_data, aes(x = X, y = Y, color = Z)) +
      geom_point(pch = ".") +
      labs(title = "Uncorrected Inclination", x = "X (m)", y = "Y (m)", color = "Z (m)") +
      theme_bw() + scale_color_viridis(option = "A")
    
    p2 <- ggplot(plot_data, aes(x = X_rot, y = Y_rot, color = Z_rot)) +
      geom_point(pch = ".") +
      labs(title = "Corrected Inclination", x = "X (m)", y = "Y (m)", color = "Z (m)") +
      theme_bw() + scale_color_viridis(option = "A")
    
    ggpubr::ggarrange(p1, p2, labels = "AUTO", common.legend = TRUE, legend = "bottom", nrow = 1)
    ggsave("Inclination_correction_comparison.png", width = 8, height = 4)
  }
  
  # Update tree data with corrected coordinates
  tree_f@data$X_cor <- rotated_points[, 1]
  tree_f@data$Y_cor <- rotated_points[, 2]
  tree_f@data$Z_cor <- rotated_points[, 3] - min(rotated_points[, 3])
  
  # Save rotation matrix and offset details
  offset_min <- c("x" = min(tree_f@data$X), "y" = min(tree_f@data$Y))
  offset_z <- c("z" = min(rotated_points[, 3]))
  saveRDS(list(offset_min, offset_z, rotation_matrix), file = "rotation_stuff.rds")
  
  return(tree_f)
}

##### 3 fit multiple circles along the stem height ####
#' Fit RANSAC Circles to Tree Point Cloud Data and Calculate Radii Metrics
#'
#' This function uses RANSAC to fit circles to stem slices and calculates several radius-related metrics, 
#' such as the convex hull area, circle area, and bark thickness, at different stem heights.
#'
#' @param tree_f A LAS object containing tree scan data, with `Z_cor`, `X_cor`, and `Y_cor` coordinates.
#' @param steps Numeric, the step size for slicing the stem along the Z-axis (e.g., 0.05 for 5 cm steps).
#' @param species Character, species name used to calculate bark thickness (e.g., "FaSy").
#' @return A data.frame with calculated radii and other metrics at different stem heights.
#' @examples
#' # Assuming 'tree_f' is a valid LAS object:
#' radii_metrics <- veneer_ransac(tree_f, steps = 0.05, species = "FaSy")
#' @export
#' @import TreeLS
#' @import sp

veneer_ransac <- function(tree_f, steps, species) {
  invisible(lapply(c('TreeLS','sp', 'ggplot2'), require, character.only = TRUE))
  
  # Initialize empty data frame to store results
  radii <- data.frame(
    Z = numeric(), R = numeric(), X = numeric(), Y = numeric(),
    chull.area = numeric(), circle.area = numeric(),
    proportion_circle_hull = numeric(), raw.x = numeric(), 
    raw.y = numeric(), bark = numeric()
  )
  
  # Counter for results
  i <- 1
  
  # Loop through stem heights with the specified step size
  for (z in seq(min(tree_f@data$Z_cor) + 0.01, max(tree_f@data$Z_cor), steps)) {
    
    # Extract slice of points around the current height z
    disc <- filter_poi(tree_f, Z_cor > (z - steps / 2), Z_cor < (z + steps / 2))
    
    # Continue only if enough points are available for circle fitting
    if (length(disc@data$X_cor) > 5) {
      
      disc@data$X <- disc@data$X_cor
      disc@data$Y <- disc@data$Y_cor
      disc@data$Z <- disc@data$Z_cor
      
      # Fit circle to the current slice using RANSAC
      circ <- TreeLS::shapeFit(disc, shape = "circle", algorithm = "ransac")
      
      if (!is.na(circ[1, 3])) {
        
        # Store the circle's radius and other calculated metrics
        radii[i, ] <- c(
          z, circ[1, 3], circ[1, 1], circ[1, 2], NA, NA, NA,
          mean(disc@data$X), mean(disc@data$Y), NA
        )
        
        # Calculate bark thickness for the given species (if applicable)
        if (species == "FaSy") {
          radii$bark[i] <- 2.61029 + 0.28522 * radii$R[i] * 200
        }
        
        # Perform outlier detection using boxplot
        disc_df <- data.frame(disc@data[, 1:3])
        disc_df$znew <- z
        
        # Remove X and Y outliers using boxplot.stats
        disc_o <- remove_outliers(disc_df)
        
        # Calculate convex hull
        chull_pts <- chull(disc_o[, c(1, 2, 4)])
        chull_pts <- c(chull_pts, chull_pts[1])  # Close the polygon
        disc_hull <- disc_o[chull_pts, ]
        chull.poly <- Polygon(disc_hull[, c(1, 2)], hole = FALSE)
        
        # Store convex hull area and circle area metrics
        radii[i, 5] <- chull.poly@area
        radii[i, 6] <- pi * radii$R[i]^2
        radii[i, 7] <- radii[i, 5] / radii[i, 6]  # Proportion of hull area to circle area
        
        # Move to the next iteration
        i <- i + 1
      }
    }
  }
  
  # Update final tree information (e.g., radius at various heights)
  update_tree_info(radii, species)
  
  return(radii)
}



##### 3.1 outlier detection based on splines for axis along tree height:##### 
#' Detect and Remove Outliers from Tree Point Cloud Data
#'
#' This function detects outliers based on distance from a fitted spline to the points in a tree's point cloud data.
#'
#' @param radii A data.frame with calculated radii, X, Y, Z values, and other metrics.
#' @param tree_f A LAS object containing tree scan data, with coordinates (X_cor, Y_cor, Z_cor).
#' @param steps Numeric, the step size for slicing the stem along the Z-axis.
#' @return A LAS object with outliers removed based on calculated distance deviations.
#' @examples
#' # Assuming 'tree_f' is a valid LAS object and 'radii' is computed:
#' tree_f_cleaned <- veneer_outlier(radii, tree_f, steps = 0.05)
#' @export
#' @import ggpubr

veneer_outlier <- function(radii, tree_f, steps) {
  
  # Smooth splines for X and Y coordinates based on Z
  new_x <- smooth.spline(x = radii$Z, y = radii$X, nknots = 10)
  new_y <- smooth.spline(x = radii$Z, y = radii$Y, nknots = 10)
  
  radii$X_spline <- predict(new_x, newdata = radii$Z)$y
  radii$Y_spline <- predict(new_y, newdata = radii$Z)$y
  
  # Initialize columns in tree_f data
  tree_f@data$dist <- NA
  tree_f@data$diff_dist <- NA
  tree_f@data$diff_dist_raus <- FALSE
  
  # Loop over unique Z values in the radii dataset to calculate distances
  for (radii_z in unique(radii$Z)) {
    rz = radii [which(radii$Z == radii_z),]
    cc_x = rz$X 
    cc_y = rz$Y
    dist=c()
    
    tree_f_sub = tree_f@data[which(tree_f@data$Z_cor>= radii_z & tree_f@data$Z_cor < radii_z+ steps ),]
    for (i in 1:nrow(tree_f_sub))  dist[i] =  dist(rbind(cbind(cc_x, cc_y), cbind(tree_f_sub$X_cor[i], tree_f_sub$Y_cor[i])))
    
    tree_f@data$dist      [which(tree_f@data$Z_cor>= radii_z & tree_f@data$Z_cor < radii_z+ steps )] = dist
    tree_f@data$diff_dist [which(tree_f@data$Z_cor>= radii_z & tree_f@data$Z_cor < radii_z+ steps )] = dist - rz$R
    
    
  }
  
  
  # Mark outliers where the difference in distance exceeds a threshold
  tree_f@data$diff_dist_raus <- tree_f@data$diff_dist > 0.1
  
  
  # Plot outlier detection results
  a <- ggplot(tree_f@data) + 
    geom_point(aes(x = X_cor, y = Y_cor, col = diff_dist_raus), alpha = 0.3) + 
    theme_minimal()
  
  b <- ggplot(tree_f@data) + 
    geom_point(aes(x = X_cor, y = Z_cor, col = diff_dist_raus), alpha = 0.3) + 
    theme_minimal()
  
  c <- ggplot(tree_f@data) + 
    geom_point(aes(x = Y_cor, y = Z_cor, col = diff_dist_raus), alpha = 0.3) + 
    theme_minimal()
  
  ggarrange(a, b, c, ncol = 3, common.legend = TRUE, legend = "bottom")
  
  # Save the plot
  ggsave("fine_outlier_detection.png", width = 10, height = 4)
  
  # Remove outliers from the tree_f data
  tree_f = tree_f[which(tree_f@data$diff_dist_raus==F),]
  
  return(tree_f)
}


##### 4 cardinal directions: ####

#' Compute Cardinal Directions for Tree Points with Optional Radii Recalculation
#'
#' This function calculates the cardinal directions (e.g., N, NE, E, etc.) for each point along a tree's trunk 
#' based on 3D coordinates and radii information. If the tree scan is incomplete, additional points are generated 
#' using circular fit radii. The function can generate both 2D and 3D plots of the calculated cardinal directions.
#'
#' @param tree_f A `LAS` object or similar, containing 3D coordinates (`X_cor`, `Y_cor`, and `Z_cor`) of tree points.
#' @param radii A data frame containing radii information for calculating tree circular fits. Must contain columns `Z`, `X`, and `Y`.
#' @param steps Numeric, the height interval step (in meters) for computing directions along the tree (default is `0.05`).
#' @param plot2d Logical, if `TRUE`, generates a 2D plot of the cardinal directions (default is `TRUE`).
#' @param plot3d Logical, if `TRUE`, generates a 3D plot of the cardinal directions (default is `TRUE`).
#' @return Data frame with updated `car_ang` (angle) and `car_dir` (direction) columns for each point.
#'         If necessary, points are interpolated to fill gaps in tree data.
#'
#' @details This function checks whether radii need recalculating based on the scan completeness.
#'          It assigns a cardinal direction to each point based on its position relative to the tree center at each height layer.
#' 
#' @importFrom ggplot2 ggplot aes geom_point ggtitle theme_bw
#' @importFrom plotly plot_ly add_trace layout
#' @importFrom plyr rbind.fill
#' @importFrom conicfit calculateCircle
#' @importFrom viridis magma
#' @importFrom htmlwidgets saveWidget

#' @export
#' @examples
#' \dontrun{
#' # Example usage:
#' tree_data <- readLAS("path/to/tree_points.las")
#' radii_data <- data.frame(Z = seq(0, 10, 0.5), X = runif(20), Y = runif(20))
#' veneer_card(tree_f = tree_data, radii = radii_data, plot2d = TRUE, plot3d = TRUE)
#' }
#' 

veneer_card <- function(tree_f, radii, steps = 0.05, plot_it_2d = TRUE, plot_it_3d = TRUE) {
  invisible(lapply(c('plotly','sp', 'ggplot2'), require, character.only = TRUE))
  
  # Initialize columns for cardinal angles and directions
  tree_f@data$car_ang <- 0
  tree_f@data$car_dir <- "P"
  
  # Check if the tree is fully scanned; if not, create additional points from radii circle fits
  if (length(which(table(tree_f@data$car_dir) < 2000)) > 3) {
    circle_tree <- c()
    
    for (i in 1:nrow(radii)) {
      # Add conditions for appending radii points to create missing tree points
      if (i > 2 && (radii$R[i] > radii$R[i - 1] || 
                    (is.na(radii$R[i - 1]) && radii$R[i] > 1.05 * radii$R[i - 2]) || 
                    radii$R[i] > 1.2 * mean(radii$R[radii$Z > 1 & radii$Z < 1.3]))) {
        next
      }
      
      # Generate additional points using radii circles
      new_points <- cbind(radii$Z[i], conicfit::calculateCircle(radii$X[i], radii$Y[i], radii$R[i], steps = 100))
      colnames(new_points) <- c("Z_cor", "X_cor", "Y_cor")
      circle_tree <- rbind(circle_tree, new_points)
    }
    
    # Append new points to the original tree data
    circle_tree <- as.data.frame(circle_tree)
    tree_f@data <- rbind.fill(tree_f@data, circle_tree)
  }
  
  # Calculate cardinal angles and directions for the combined data
  for (radii_z in unique(radii$Z)) {
    rz <- radii[which(radii$Z == radii_z), ]
    cc_x <- rz$X
    cc_y <- rz$Y
    tree_f_sub <- tree_f@data[which(tree_f@data$Z_cor >= radii_z & tree_f@data$Z_cor < radii_z + steps), ]
    
    # Calculate angles relative to North
    cardinal_angle <- apply(matrix(1:nrow(tree_f_sub)), 1, function(i) {
      stempoint <- tree_f_sub[i, ]
      branchv <- as.numeric(c((stempoint$X_cor - cc_x), (stempoint$Y_cor - cc_y)))
      northv <- as.numeric(c(cc_x - cc_x, (cc_y + 100) - cc_y))
      angle <- 180 * (angle2(branchv, northv)) / pi
      ifelse(angle < 0, 360 - abs(angle), angle)
    })
    
    # Map angles to cardinal directions
    cardinal_dir <- d2c.2(cardinal_angle)
    
    # Update data frame with calculated angles and directions
    tree_f@data$car_ang[which(tree_f@data$Z_cor >= radii_z & tree_f@data$Z_cor < radii_z + steps)] <- cardinal_angle
    tree_f@data$car_dir[which(tree_f@data$Z_cor >= radii_z & tree_f@data$Z_cor < radii_z + steps)] <- cardinal_dir
  }
  
  # Filter out rows without assigned directions
  tree_f <- tree_f[(tree_f@data$car_dir != "P"), ]
  tree_f@data$car_dir <- factor(tree_f@data$car_dir, levels = c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"))
  
  # Optional 2D plotting of cardinal directions
  if (plot_it_2d) {
    ggplot(tree_f@data, aes(x = X_cor, y = Y_cor, col = car_dir)) + 
      geom_point() + 
      ggtitle("Cardinal Directions") + 
      scale_colour_viridis_d(option = "A") + 
      theme_bw()
    ggsave("Cardinal_direction.pdf", width = 10, height = 5)
  }
  
  # Optional 3D plotting of cardinal directions
  if (plot_it_3d) {
    p <-  plotly::plot_ly(colors = setNames(viridis::magma(30), unique(tree_f@data$car_dir)))  %>%
      plotly::add_trace(
        x = tree_f@data$X_cor - mean(tree_f@data$X_cor),
        y = tree_f@data$Y_cor - mean(tree_f@data$Y_cor),
        z = tree_f@data$Z_cor,
        type = "scatter3d",
        mode = "markers",
        color = tree_f@data$car_dir
      ) %>%
      plotly::layout(
        scene = list(
          xaxis = list(title = "X (m)"),
          zaxis = list(title = "Tree height (m)"),
          yaxis = list(title = 'Y (m)'),
          aspectmode = "manual",
          aspectratio = list(z = 1, x = 0.2, y = 0.2)))

    htmlwidgets::saveWidget(plotly::partial_bundle(p), file = "Tree_Cardinal_Direction.html", selfcontained = TRUE)
  }
  
  return(tree_f@data)
}

##### 5 fit splines along vertical partly point clouds: ####
#' @title Veneer Spline Fitting and Visualization
#' @description This function fits splines to tree coordinate data (X, Y) at different tree heights (Z) 
#'              based on cardinal directions, visualizes the fitted splines in a 3D plot, and returns the predicted coordinates.
#' @param tree_new A data frame with columns 'X_cor', 'Y_cor', 'Z_cor' for coordinates and 'car_dir' for cardinal direction.
#' @return A data frame containing the predicted X, Y coordinates at different Z heights for each cardinal direction.
#' @import plotly
#' @import magma
#' @import smooth
#' @import htmlwidgets
#' @examples
#' # Example usage
#' # tree_data <- data.frame(X_cor = ..., Y_cor = ..., Z_cor = ..., car_dir = ...)
#' # result <- veneer_spline(tree_data)

veneer_spline = function(tree_new) {
  
  # Initialize an empty list to store the predicted spline values
  pred_all = c()
  
  # Generate a color palette based on the number of unique cardinal directions (car_dir)
  pal = magma(length(unique(tree_new$car_dir)))
  
  # Center the coordinates by subtracting the mean of X and Y
  tree_new$X_cor = tree_new$X_cor - mean(tree_new$X_cor)
  tree_new$Y_cor = tree_new$Y_cor - mean(tree_new$Y_cor)
  
  # Loop through each unique cardinal direction (car_dir)
  for (cd in unique(tree_new$car_dir)) {
    
    # Skip directions with fewer than 11 data points (necessary for spline fitting)
    if(nrow(tree_new[which(tree_new$car_dir == cd),]) < 11) {
      message(paste("Skipping direction:", cd, "due to insufficient data"))
      next
    }
    
    # Fit a spline to the X coordinates based on Z coordinates for this direction
    s_X <- smooth.spline(tree_new$Z_cor[which(tree_new$car_dir == cd)], tree_new$X_cor[which(tree_new$car_dir == cd)], df = 15)
    
    # Fit a spline to the Y coordinates based on Z coordinates for this direction
    s_Y <- smooth.spline(tree_new$Z_cor[which(tree_new$car_dir == cd)], tree_new$Y_cor[which(tree_new$car_dir == cd)], df = 15)
    
    # Generate a sequence of Z values to predict X and Y coordinates
    z_here = seq(min(tree_new$Z_cor[which(tree_new$car_dir == cd)]), 
                 max(tree_new$Z_cor[which(tree_new$car_dir == cd)]), 
                 by = 0.01)
    
    # Predict the X and Y coordinates for the generated Z values using the splines
    modelled_tree_x = predict(s_X, z_here)$y
    modelled_tree_y = predict(s_Y, z_here)$y
    
    # Store the predicted coordinates along with the direction and height (z)
    pred_1 = data.frame("cd" = cd, "z" = z_here, "x" = modelled_tree_x, "y" = modelled_tree_y)
    
    # Append the predictions for this direction to the overall list
    pred_all = rbind(pred_all, pred_1)
  }
  
  # Create a 3D plot using Plotly
  p = plotly::plot_ly(colors = setNames(viridis::magma(30), unique(tree_new$car_dir))) %>%
    plotly::layout(scene = list(
      xaxis = list(title = 'X (m)'), 
      zaxis = list(title = "Tree height (m)"),
      yaxis = list(title = 'Y (m)'),
      camera = list(eye = list(x = 1.25, y = 1.25, z = 1.25), center = list(x = 0, y = 0, z = 0))
    )) %>%
    plotly::layout(scene = list(aspectmode = "manual", aspectratio = list(z = 1, x = 0.2, y = 0.2)))
  
  # Loop through each unique direction in the predictions and add traces to the plot
  for (cd in unique(pred_all$cd)) {
    p <- p %>% plotly::add_trace(
      x = pred_all$x[which(pred_all$cd == cd)], 
      y = pred_all$y[which(pred_all$cd == cd)], 
      z = pred_all$z[which(pred_all$cd == cd)], 
      type = "scatter3d", 
      mode = "markers", 
      color = cd, 
      hoverinfo = 'skip', 
      opacity = 1, 
      size = 1
    )
  }
  
  # Display the 3D plot
  p
  
  # Save the plot as an interactive HTML file with the current site and TID for naming
  htmlwidgets::saveWidget(p, paste( TID, "_splines_3d_size1.html"), selfcontained = F, libdir = "lib")
  
  # Return the predicted values for all cardinal directions
  return(pred_all)
}

##### 6 taper and volume: #####

#' Estimate Taper and Volume of a Tree
#'
#' This function estimates the taper and volume of a tree based on RANSAC-fitted radii.
#'
#' @param radii A data frame containing tree radii, convex hull areas, circle areas, bark thickness, and other measurements.
#' @return A list with taper information and volume calculations using various methods.
#' @export


taper_volume <- function(radii) {
  # Validate input
  if (!is.data.frame(radii)) stop("The input 'radii' must be a data frame.")
  if (!all(c("Z", "R", "chull.area", "circle.area", "bark") %in% colnames(radii))) {
    stop("The input 'radii' must contain columns: Z, R, chull.area, circle.area, and bark.")
  }
  
  # Compute taper
  lm_t2 <- lm(R * 200 ~ Z, data = radii[which(radii$Z > 1 & radii$Z < 0.98 * max(radii$Z)), ])
  D1 <- mean(radii$R[which(radii$Z > 0.9 & radii$Z < 1.1)] * 200)  # Diameter at 1m
  D2 <- mean(radii$R[which(radii$Z > 0.9 * max(radii$Z) & radii$Z < 0.98 * max(radii$Z))] * 200)  # Diameter at stem top
  L <- max(radii$Z) - 1  # Length of the stem
  taper_mean <- (D2 - D1) / L  # Taper value (cm/m)
  
  taper_method1 <- list(
    mean_taper = round(taper_mean, 2),
    lm_taper = round(lm_t2$coefficients[2], 2)
  )
  
  # Adjust convex hull area where the hull is larger than the circle area
  radii$chull.area[which(radii$proportion_circle_hull > 2)] <- radii$circle.area[which(radii$proportion_circle_hull > 2)]
  
  # Volume calculations
  radii$chull.volume            <- radii$chull.area * 0.05  # Convex hull volume (per 5cm slice)
  radii$circle.volume.nobark    <- (radii$R - (radii$bark / 2000))^2 * pi * 0.05  # Circle fit volume without bark
  radii$circle.volume           <- radii$circle.area * 0.05
  radii$circle.vol.bark.diff    <- (radii$circle.volume - radii$circle.volume.nobark) / radii$circle.volume
  radii$chull.volume.nobark     <- radii$chull.volume - radii$chull.volume * radii$circle.vol.bark.diff
  
  # Summarize volumes
  fm_polygon <- sum(radii$chull.volume.nobark, na.rm = TRUE)
  fm_circle  <- sum(radii$circle.volume.nobark, na.rm = TRUE)
  
  # Conventional volume calculation
  mid_diameter <- mean(radii$R[which(radii$Z > 0.45 * max(radii$Z) & radii$Z < 0.55 * max(radii$Z))]) * 2
  mid_bark     <- mean(radii$bark[which(radii$Z > 0.45 * max(radii$Z) & radii$Z < 0.55 * max(radii$Z))]) / 1000
  fm_tree      <- ((mid_diameter - mid_bark)^2) * pi / 4 * max(radii$Z)
  
  # Output
  return(list(
    taper = taper_method1,
    volumes = list(
      fm_polygon = round(fm_polygon, 4),
      fm_circle = round(fm_circle, 4),
      fm_conventional = round(fm_tree, 4)
    ),
    radii = radii  # Return updated radii with volume calculations
  ))
}

##### 7 the important function: #####

#' Analyse Veneer Potential and Cylinder Fitting in Tree Stems
#'
#' This function analyzes point cloud data from TLS of tree stems to identify
#' optimal veneer log combinations by simulating cut lengths, estimating radii,
#' and calculating veneer vs. industrial wood volumes.
#'
#' @param pred_all A data.frame of TLS spline points with at least `x`, `y`, `z` coordinates.
#' @param tls_all_withcrown A LAS object containing TLS points with tree crowns.
#' @param final_tree_information A data.frame with existing metadata on the tree.
#' @param TID Tree ID to filter TLS data.
#' @param x_move, y_move Numeric values to reposition TLS data.
#' @param site Character string for site name (used in file naming).
#' @param stump Numeric height of the stump (default: 0.3 m).
#' @param length1, length2 Numeric lengths (in meters) for possible veneer log cuts.
#' @param pal Optional vector of ggplot2 color codes.
#' @param fm_tree,fm_circle Optional volume estimates for comparison.
#' @param taper_mean, lm_t2 Optional taper model results.
#'
#' @return A list containing best combination information and final_tree_information.
#' @export
#' @importFrom dplyr group_by summarize right_join reframe
#' @importFrom geometry convhulln
#' @importFrom conicfit CircleFitByKasa calculateCircle
#' @importFrom rgeos gCentroid
#' @importFrom sf st_point st_polygon st_covered_by
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom grDevices ggsave
#' @importFrom graphics points
#' @import ggplot2
#' 
analyse_veneer_potential <- function(pred_all, tls_all_withcrown, final_tree_information,
                                     TID, x_move, y_move, 
                                     stump = 0.3, length1 = 1.8288, length2 = 2.4384,
                                     pal = magma(16), fm_tree = NA, fm_circle = NA,
                                     taper_mean = NA, lm_t2 = list(coefficients = NA)) {

  #Remaining roll has a diameter of 95mm / 9.5 cm
  #Median diameter (German: Mittendurchmesser) > 30 cm 
  length1 = 1.8288 #(m) length 1 6ft or 1.85m
  length2 = 2.4384 #(m) length 2 8ft or 2.45m
  max_length = max(tree_new$Z_cor, na.rm=T) #max. Stem length (m)
  stump = 0.3
  
  if (c(max(pred_all$z)-min(pred_all$z)-.3) < length1){
    write.table(final_tree_information, "final_tree_information.csv", col.names=T, row.names=F)
    tls_with_crown = as.data.frame(tls_all_withcrown@data [which(tls_all_withcrown@data$TreeID == TID),])
    tls_with_crown$X = tls_with_crown$X - x_move
    tls_with_crown$Y = tls_with_crown$Y - y_move
    writeLAS(LAS(tls_with_crown), paste0( "cylinder_mesh/", TID, "_tls_with_crown_grey.las"))
  }
  
  if (c(max(pred_all$z)-min(pred_all$z)-.3) < length1) next
  
  if (final_tree_information$mitten_diameter*2 < .3){
    write.table(final_tree_information, "final_tree_information.csv", col.names=T, row.names=F)
    tls_with_crown = as.data.frame(tls_all_withcrown@data [which(tls_all_withcrown@data$TreeID == TID),])
    tls_with_crown$X = tls_with_crown$X - x_move
    tls_with_crown$Y = tls_with_crown$Y - y_move
    writeLAS(LAS(tls_with_crown), paste0("cylinder_mesh/", TID, "_tls_with_crown_grey.las"))
  }
  
  if (final_tree_information$mitten_diameter*2 < .3) next
  
  ##### Prepare the data to find centroid along the stem####
  
  splines   = pred_all
  splines$Z = round(splines$z, digits=3)
  sp0       = splines [order(splines$Z),]
  #begin with the maximum radius at the bottom
  #choose the minimum radius of all cardinal directions at certain height
  
  veneer_all = c()
  for (h in seq(0,max_length,by=0.01 )){
    
    #get centroid of the splines:
    t_l = h
    t_h = h +0.04 #4cm bereich
    sp0_z1_base = sp0 [which (sp0$Z< t_h & sp0$Z> t_l),]
    if(nrow(sp0_z1_base) == 0) next
    sp0_z1_base_mean = sp0_z1_base %>% group_by (cd) %>% summarize (mx = mean(x), my = mean(y))
    sp0_z1_base_mean = sp0_z1_base_mean [order(sp0_z1_base_mean$mx),]
    chull.poly = Polygon(sp0_z1_base_mean[c(2,3)] [c(chull(sp0_z1_base_mean$mx, sp0_z1_base_mean$my), chull(sp0_z1_base_mean$mx, sp0_z1_base_mean$my)[1]),], hole=F)
    Poly.2 <- SpatialPolygons(list(Polygons(list(chull.poly), ID = 1)))
    trueCentroids = rgeos::gCentroid(Poly.2)
    sp0_z1_base_mean$center_x = trueCentroids@coords[1]
    sp0_z1_base_mean$center_y = trueCentroids@coords[2]
  
    sp0_z1_base_mean$dist = 0
    for (j in 1:nrow(sp0_z1_base_mean)) {
      sp0_z1_base_mean$dist[j] = dist(rbind(c(trueCentroids@coords[1], trueCentroids@coords[2]), c(sp0_z1_base_mean[j,2],  sp0_z1_base_mean[j,3])))
    }
    
    
    veneer_disk = data.frame(z=h, x=sp0_z1_base_mean$center_x[1], y=sp0_z1_base_mean$center_y[1], r=min(sp0_z1_base_mean$dist)  )
    veneer_all = rbind(veneer_all, veneer_disk)
    
  }
  
  #veneer to use: 
  veneer_0.3 = veneer_all [which(veneer_all$z >stump),]
  x <- veneer_0.3$x
  y <- veneer_0.3$y
  z <- veneer_0.3$z
  radi <- veneer_0.3$r
  
  ##### include in the calculation of the best combination also the radii of the potential cylinders.####
  generate_combinations <- function(length1, length2, max_length, current_sequence, remaining_length, combinations) {
    if (remaining_length < length1) {
      combinations <- append(combinations, list(list(sequence = current_sequence, remaining_length = remaining_length)))
      return(combinations)
    }
    
    if (length(current_sequence) + length1 <= max_length && remaining_length - length1 >= 0) {
      combinations <- generate_combinations(length1, length2, max_length, c(current_sequence, length1), remaining_length - length1, combinations)
    }
    
    if (length(current_sequence) + length2 <= max_length && remaining_length - length2 >= 0) {
      combinations <- generate_combinations(length1, length2, max_length, c(current_sequence, length2), remaining_length - length2, combinations)
    }
    
    return(combinations)
  }
  get_all_combinations  <- function(length1, length2, max_length) {
    remaining_length <- max_length  - stump
    combinations <- list()
    combinations <- generate_combinations(length1, length2, max_length, numeric(0), remaining_length, combinations)
    return(combinations)
  }
  all_combinations      <- get_all_combinations(length1, length2, max_length)
  
  # Iterate through all possible combinations:
  results = data.frame(combination=0, theor_vol = 0, remaining_length=0)
  
  for (i in 1:length(all_combinations)) { 
    # i=1
    s = all_combinations[[i]]$sequence
    z_base <- stump
    combination_volume <- 0
    
    # Iterate through each cylinder in the current combination
    for (j in 1:length(s)) {
      # j=1
      # Determine the length and top position of the cylinder
      L <- s[j]
      z_top <- z_base + L
      
      # Calculate the minimum radius and theoretical volume of the cylinder
      mr <- min(veneer_0.3$r[which(veneer_0.3$z > z_base & veneer_0.3$z < z_top)])
      tv <- pi * mr^2 * L
      
      # Add the volume to the combination volume
      combination_volume <- combination_volume + tv
      
      # Update the base position for the next cylinder
      z_base <- z_top
      
      
    }
    
    # Add the results to the table
    results[i,] <-  c(combination = i,  theor_vol =  combination_volume, remaining_length = all_combinations[[i]]$remaining_length)
    
  }
  
  best_vol = results$combination[which(results$theor_vol > 0.999* round(max(results$theor_vol),  digits=3))]
  schlechtestes_ergebnis = min(results$theor_vol)
  bestes_ergebnis = max(results$theor_vol)
  results$relative_loss =  (bestes_ergebnis-results$theor_vol) / bestes_ergebnis
  best_vol= results$combination[which(results$relative_loss < 0.05)] #ergebnisse welche einen theoretischem volumen verlust unter 5% vom besten ergebnis erzielen.
  if (length(best_vol)>20) best_vol= best_vol[1:20]
  # Create a parallel backend with 24 cores for parallel processing
  
  begin_fe = Sys.time()
  cl <- makeCluster(15)  # Use 24 cores for parallel processing
  registerDoParallel(cl)
  
  # Define the point_inside_cylinder function
  point_inside_cylinder <- function(point, center, radius) {
    rad = dist(rbind(point, center))
    inside = rad < radius
    return(inside)
  }
  
  
  #Define optim function
  optim_fn <- function(cent, circ){
    cent <- c(X,Y)
    dists <- sqrt((circ[,1]-cent[1])^2+(circ[,2]-cent[2])^2)
    chullc <- conicfit::calculateCircle(cent[1],cent[2],true_radius) 
    chullc <- rbind(chullc, chullc[1,])
    check1 = as.logical(sf::st_covered_by(sf::st_point(c(cent[1],cent[2])), 
                                          sf::st_polygon(list(chullc))))
    
    if(cent[1] < max(circ[,1]) & cent[1] > min(circ[,1]) & cent[2] <  max(circ[,2]) & cent[2] > min(circ[,2]) & check1 == F) { #bounding box
      
      radius <- -min(dists)} else {
        radius <- min(dists)
      }
    
    
    return(radius)
  }
  
  #Define optim function
  optim_fn <- function(cent, circ){
    
    
    dists <- sqrt((circ[,1]-cent[1])^2+(circ[,2]-cent[2])^2)
    
    
    check1 = as.logical(sf::st_covered_by(sf::st_point(c(cent[1],cent[2])), 
                                          sf::st_polygon(list(chullc))))
    #print(check1)
    #  if(cent[1] < max(circ[,1]) & cent[1] > min(circ[,1]) & cent[2] <  max(circ[,2]) & cent[2] > min(circ[,2]) & check1 == F) { #bounding box
    if(!is.na(check1) & check1) { #bounding box
      
      radius <- -min(dists)} else {
        radius <-  min(dists)
      }
    
    
    return(radius)
  }
  
  # Initialize a container for the results
  ci = list()
  
  # Parallelize the for loop using foreach
  ci <- foreach(b = best_vol, .combine = rbind, .packages = c("ggplot2","geometry", "plyr", "dplyr","CVD", "pracma", "ggforce")) %dopar% {
    # Rest of your code as it is
    #  b = best_vol[1]
    veneer_0.3 = veneer_all [which(veneer_all$z >stump),]
    cutting_lengths = as.numeric(all_combinations[[b]]$sequence)
    cutting_heights = c(stump, stump + cumsum(cutting_lengths)) #m
    restroll = as.numeric(all_combinations[[b]]$remaining_length) #m
    
    veneer_0.3$group <- cut(veneer_0.3$z, breaks = c(-Inf, cutting_heights, Inf),  right = FALSE)
    
    table_group = table(veneer_0.3$group)
    
    vr = veneer_0.3 %>% group_by (group) %>% summarize (min_r = min(r, na.rm=T))
    
    veneer_0.3 = right_join(veneer_0.3, vr, by = "group")
    
    cylinder_info = data.frame( veneer_0.3 %>% 
                                  group_by ( group) %>% 
                                  reframe (min_r = mean(min_r), 
                                           base_x = first(x), base_y=first(y), base_z = first(z),
                                           top_x = last(x), top_y= last (y), top_z = last(z), 
                                           min_x =mean(x[which(r==min_r)]), min_y = mean(y[which(r==min_r)]), min_z = mean(z[which(r==min_r)])))
    
    if( nrow(cylinder_info) == length(cutting_heights) )  {
      cylinder_info$cutting_heights =   cutting_heights 
      cylinder_info$height = c(cutting_lengths, restroll) 
    }
    
    if( nrow(cylinder_info) != length(cutting_heights) )  {
      
      cylinder_info$cutting_heights =   cutting_heights[table_group>1]
      
      
      cylinder_info$height =  cutting_lengths
    }
    
    cylinder_info = cylinder_info [which(cylinder_info$height %in% c(length1, length2)),]
    # Plot 3D the cutting heights (minimum restrolle) along the stem: ####
    x <- mean(pred_all$x)  # x-coordinate of plane position
    y <- mean(pred_all$y) # y-coordinate of plane position
    z <- cutting_heights  # z-coordinates of plane centers
    
    # do the cylinder fit for all cylinders per tree####
    
    #prepare for the following parameters needed to be found:
    cylinder_info$true_radius = cylinder_info$true_centroid_x = cylinder_info$true_centroid_y = 0 #based on the true center of the point cloud
    cylinder_info$optimal_radius = cylinder_info$optimal_centroid_x = cylinder_info$optimal_centroid_y = 0 #optimized moved center to increase the radius
    
    cylinder_info = cylinder_info[which(cylinder_info$height > 1.8),]
    cylinder_info$optimal_bark = 0
    cylinder_info$optimal_radiuswithoutbark = 0
    cylinder_info$pointcloud_vol = 0 
    cylinder_info$combination = b
    cylinder_info$c = 1:nrow(cylinder_info)
    
    # Parallelize the for loop using foreach
    
    for (c in 1:nrow(cylinder_info)) {
      
      point_cloud <- cbind(splines[which(splines$z> cylinder_info$base_z [c] & splines$z < cylinder_info$top_z[c]), c(3, 4, 2)])
      point_cloud_2d <- point_cloud[, c(1, 2)]
      height <- cylinder_info$height[c]   # the height of the inscribed cylinder
      true_centroid <- colMeans(point_cloud_2d)   # Find the centroid of the point cloud
      point_cloud_distance <- apply(point_cloud_2d, 1, function(x) dist(rbind(x, true_centroid)))
      radius <- max(point_cloud_distance)   # Calculate the radius based on the candidate centroid
      max_radius_center = radius
      # Check if there are points inside the initial cylinder:
      points_within_cylinder <- apply(point_cloud_2d, 1, point_inside_cylinder, center = true_centroid, radius = radius)
      # Adjust the cylinder radius to avoid inlying points
      reduction_factor <- 0.999  # Starting reduction factor
      reduction_step <- 0.001
      
      while (any(points_within_cylinder) && reduction_factor >= 0.2) {
        radius <- radius * reduction_factor
        points_within_cylinder <- sqrt((point_cloud_2d[, 1] - true_centroid[1])^2 + (point_cloud_2d[, 2] - true_centroid[2])^2) <= radius
        reduction_factor <- reduction_factor - reduction_step  # Decrease the reduction factor step by step
      }
      
      true_radius = radius
      cylinder_info$true_radius[c] = true_radius
      cylinder_info$true_centroid_x[c]= true_centroid[1]
      cylinder_info$true_centroid_y[c]= true_centroid[2]
      hull = geometry::convhulln(point_cloud, options="FA")
      cylinder_info$pointcloud_vol [c] = hull$vol
      
      
      bounding_circle <- conicfit::CircleFitByKasa(as.matrix(point_cloud_2d))
      chullc <- conicfit::calculateCircle(bounding_circle[1],bounding_circle[2],bounding_circle[3]) 
      chullc <- rbind(chullc, chullc[1,])
      
      
      #optimization of the centroid to potentially increase the radius:####
      opt <- optim(c(X = true_centroid[1], Y = true_centroid[2]), optim_fn, circ = as.matrix(point_cloud_2d))
      
      cylinder_info$optimal_radius [c]    = abs(opt$value)
      cylinder_info$optimal_centroid_x [c]= opt$par[1]
      cylinder_info$optimal_centroid_y [c]= opt$par[2]
      cylinder_info$optimal_bark [c]      = c(2.61029  + 0.28522* max(abs(opt$value)*200))/2 #radius in cm to obtain oneside bark thickness in mm
      cylinder_info$optimal_radiuswithoutbark [c]  = cylinder_info$optimal_radius [c] -   c(cylinder_info$optimal_bark [c]/1000)  #radius in cm to obtain oneside bark thickness in mm
      
      # ggplot() + geom_circle(aes(x0=cylinder_info[c,]$optimal_centroid_x, y0=cylinder_info[c,]$optimal_centroid_y, r=cylinder_info[c,]$optimal_radius), col="blue") + geom_point(aes(x=point_cloud_2d$x, y=point_cloud_2d$y)) + 
      #   geom_circle(aes(x0=cylinder_info[c,]$true_centroid_x, y0=cylinder_info[c,]$true_centroid_y, r=cylinder_info[c,]$true_radius)) +   geom_circle(color="red", aes(x0=bounding_circle[1], y0=bounding_circle[2], r=bounding_circle[3]))
      
      
    }
    
    cylinder_info$veneer_volume = cylinder_info$optimal_radius^2*pi*(cylinder_info$height)
    cylinder_info$veneer_volumewithout_bark = cylinder_info$optimal_radiuswithoutbark^2*pi*(cylinder_info$height)
    
    return(cylinder_info)
  }
  
  # Stop the parallel backend
  stopCluster(cl)
  
  #calculate final information on volumes: 
  point_cloud <- splines
  if(nrow(point_cloud[which(point_cloud$z <stump),c(3,4,2)])>0) {
    ass = ashape3d(as.matrix(data.frame(point_cloud[which(point_cloud$z <stump),c(3,4,2)])), pert=T, alpha=.5)
    ci$stump_volume = volume_ashape3d(ass) 
  } else {ci$stump_volume = 0 }
  
  ci$ven_vol = pi*ci$height*ci$optimal_radius^2
  ci$restrolle_vol = pi*ci$height*c(0.095/2)^2 #restrolle hat einen durchmesser von 95mm / 9.5 cm
  ci$long_term_c_vol = ci$ven_vol-ci$restrolle_vol
  ci$radius_improvement = 100*c(ci$optimal_radius-ci$true_radius)/ci$true_radius
  ci$sum_ven_vol = 0
  ci$sum_longterm_c_vol = 0
  ci$stem_volume = 0
  
  for (i in unique(ci$combination)){
    
    cci = ci [which(ci$combination==i),]
    ci$sum_ven_vol [which(ci$combination==i)] = sum(cci$ven_vol)
    ci$sum_longterm_c_vol [which(ci$combination==i)] = sum(cci$long_term_c_vol)
    ci$stem_volume  [which(ci$combination==i)] = sum(cci$pointcloud_vol) 
    
  }
  
  ci$share_veneer_rolls = ci$sum_ven_vol/ ci$stem_volume #give out a volume and the relative share of veneer wood / industry wood / bark 
  ci$share_longterm_c = ci$sum_longterm_c_vol/ ci$stem_volume #give out a volume and the relative share of veneer wood / industry wood / bark 
  ci$restrolle_vol = pi*ci$height*c(0.095/2)^2 #restrolle hat einen durchmesser von 95mm / 9.5 cm
  
  best_combination = ci[which(ci$share_longterm_c == max(ci$share_longterm_c)),]
  best_combination$remaining_stem_disk_length = all_combinations[[unique(best_combination$combination)]]$remaining_length
  write.table(best_combination, paste0(site, "_", TID, "_bestcombi_veneer.csv"), col.names=T, row.names=F, dec=".", sep=";")
  
  ggplot(data=best_combination) + geom_col( aes(x=base_z, y=100*c(optimal_radius-true_radius)/true_radius), fill=pal[3]) +  ylab("Improvement Radius (%)") + xlab("Cylinder base height (m)") + theme_minimal()+ ggtitle("Relative change in radius after improvement") 
  ggplot(data=best_combination) + geom_col( aes(x=c, y=100*c(optimal_radius-true_radius)/true_radius), fill=pal[3]) +  ylab("Improvement Roll Radius (%)") + xlab("Cylinder number") + theme_minimal()
  ggsave("best_combination.pdf", width=3, height=3)
  ggsave("best_combination.png", width=3, height=3)
  
  worst_combination = ci[which(ci$share_longterm_c == min(ci$share_longterm_c)),]
  worst_combination$remaining_stem_disk_length = all_combinations[[unique(worst_combination$combination)]]$remaining_length
  write.table(worst_combination, paste0(site, "_", TID, "_worstcombi_veneer.csv"), col.names=T, row.names=F, dec=".", sep=";")
  
  
  ##### get final information per tree together ####
  
  #final tree information
  final_tree_information$DBH_spline = 2*mean(veneer_all$r [which(veneer_all$z >1.2 & veneer_all$z < 1.4)])
  if(max(veneer_all$z) >5)final_tree_information$D5_spline = 2*mean(veneer_all$r [which(veneer_all$z >4.9 & veneer_all$z < 5.1)])
  if(max(veneer_all$z) >10) final_tree_information$D10_spline = 2*mean(veneer_all$r [which(veneer_all$z >9.8 & veneer_all$z < 10.2)])
  
  
  # volume
  final_tree_information$stem_volume = best_combination$stem_volume[1]
  final_tree_information$longterm_c_volume = best_combination$sum_longterm_c_vol[1]
  final_tree_information$industry_volume = best_combination$stem_volume[1]-best_combination$sum_longterm_c_vol[1]
  final_tree_information$stump_vol = best_combination$stump_volume[1]
  # optim improvements:
  final_tree_information$combination_improvement_percent = 100* unique(best_combination$share_longterm_c) - unique(worst_combination$share_longterm_c) # Anteil am Gesamtvolumen (beste kombination) - anteil am gesamtvolumen (schlechteste combi)
  final_tree_information$radius_improvement_percent  = mean(best_combination$radius_improvement) #reine verbesserung innerhalb der besten combination 
  final_tree_information$stem_volume_conventional = fm_tree
  final_tree_information$stem_volume_circleransac = fm_circle
  # taper
  final_tree_information$taper_conventional = taper_mean
  final_tree_information$taper_lm_splines   = lm_t2$coefficients[2]
  final_tree_information$restrolle_vol = sum(best_combination$restrolle_vol)
  write.table(final_tree_information, "final_tree_information.csv", col.names=T, row.names=F)
  
  
  ggplot(data=final_tree_information) + 
    geom_col( aes(x="Stem", y=stem_volume), fill=pal[8], col= pal[8], width=1) +
    geom_col( aes(x="1: Veneer", y=longterm_c_volume ), fill=pal[5]) +
    geom_col( aes(x="2: Non-veneer incl. Restroll", y=industry_volume), fill=pal[12])+ 
    geom_col( aes(x="2: Non-veneer incl. Restroll", y=restrolle_vol), fill=pal[14])+ 
    ylab("Volume (m³)") + xlab("Type") + theme_minimal() 
  ggsave("Share_veneer_industry.pdf", width=5, height=4)
  ggsave("Share_veneer_industry.png", width=5, height=4)
}
