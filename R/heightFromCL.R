# Metadata

heightFromCL = function(coastlines, DEM){
  # function that ...
  n = length(coastlines)
  result = data.frame(matrix(rep(NA, n*2), nrow = n))
  colnames(result) = c('Mean','STDEV')
  for (i in 1:n){
    result$Mean[i]  = extract(DEM, coastlines[i,], fun = mean, na.rm = TRUE)
    result$STDEV[i] = extract(DEM, coastlines[i,], fun = sd, na.rm = TRUE)  
  }
  return(result)
}