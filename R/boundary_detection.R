# Metadata

boundary_detection = function(input_stack){
  #function that....
  result = list()
  for (i in 1:nlayers(input_stack)){
    result[[i]] = rasterToContour(input_stack[[i]], nlevels = 1)
  }
  return(result)
}