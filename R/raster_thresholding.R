# Metadata

raster_thresholding = function(input_stack){
  #function that....
  result = input_stack
  for (i in 1:nlayers(input_stack)){
    raster = result[[i]]
    threshold = autoThreshold(raster, 0.25)[3]
    raster[raster >= threshold] = 1
    raster[raster <  threshold] = 0
    result[[i]] = raster
  }
  return(result)
}