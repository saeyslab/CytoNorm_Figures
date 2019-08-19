gaussNorm <- function (flowset, channel.names, max.lms = 2, base.lms = NULL, 
                       peak.density.thr = 0.05, peak.distance.thr = 0.05, 
                       debug = FALSE, fname = ""){
  remb.flowset = flowStats:::remBoundary(flowset, channel.names)
  expr.list = c(1:length(flowset))
  if (length(max.lms) == 1) {
    max.lms = rep(max.lms, times = length(channel.names))
  }
  if (length(max.lms) != length(channel.names)) {
    cat("Error: length of max.lms and channel.names doesn't match\n")
    return(NULL)
  }
  names(max.lms) = channel.names
  lms = flowStats:::extract.landmarks(remb.flowset$res, expr.list, channel.names, 
                          max.lms, peak.density.thr, peak.distance.thr)
  if (is.null(base.lms)) 
    base.lms = flowStats:::extract.base.landmarks(lms$filter, channel.names, 
                                      max.lms)
  matched.lms = flowStats:::match.all.lms(lms, base.lms, channel.names, 
                              max.lms)
  #confidence = compute.confidence(matched.lms, base.lms)
  confidence = NULL
  cat("\nAdjusting the distance between landmarks\n")
  newRange = matrix(ncol = length(channel.names), nrow = 2)
  colnames(newRange) = channel.names
  newRange[, channel.names] <- c(Inf, -Inf)
  for (i in expr.list) {
    cat(".")
    if (fname != "") 
      file.name = paste(fname, as.character(i), sep = ".")
    else file.name = ""
    exprs(remb.flowset$res[[i]]) = flowStats:::normalize.one.expr(exprs(remb.flowset$res[[i]]), 
                                                      base.lms, lms$filter[, i], lms$original[, i], matched.lms[, 
                                                                                                                i], channel.names, max.lms, file.name, debug)
    for (p in channel.names) {
      newRange[1, p] <- min(newRange[1, p], min(exprs(remb.flowset$res[[i]])[, 
                                                                             p], na.rm = TRUE))
      newRange[2, p] <- max(newRange[2, p], max(exprs(remb.flowset$res[[i]])[, 
                                                                             p], na.rm = TRUE))
    }
  }
  cat("\n")
  restoreBoundary(flowset, remb.flowset, channel.names)
  for (i in expr.list) {
    for (p in channel.names) {
      ip <- match(p, pData(parameters(remb.flowset$res[[i]]))$name)
      tmp <- parameters(remb.flowset$res[[i]])
      oldRanges <- unlist(range(flowset[[i]])[, p])
      pData(tmp)[ip, c("minRange", "maxRange")] <- c(min(oldRanges[1], 
                                                         newRange[1, p]), max(oldRanges[2], newRange[2, 
                                                                                                     p]))
      remb.flowset$res[[i]]@parameters <- tmp
    }
  }
  list(flowset = remb.flowset$res, confidence = confidence)
}


restoreBoundary <- function (org.flowset, remb.flowset, channel.names) 
{
  for (i in 1:length(org.flowset)) {
    ids <- is.na(exprs(remb.flowset$res[[i]]))
    exprs(remb.flowset$res[[i]])[ids] <- exprs(org.flowset[[i]])[ids]
    
    # for (p in channel.names) {
    #   exprs(remb.flowset$res[[i]])[remb.flowset$index[[i]], p] = 
    #     exprs(org.flowset[[i]])[remb.flowset$index[[i]], p]
    # }
  }
}
