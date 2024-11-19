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


