#################### SANEM ##########################

na.lomf <<- function(object, na.rm = F) {
     numCols <- ""
     for (c in 1:ncol(object)){
          numCols[c] <- is.numeric(object[,c])
     }
     
     na.lomf.0 <- function(object) {
          idx <- which(!is.na(object))
          if (is.na(object[1])) idx <- c(1, idx)
          rep.int(object[idx], diff(c(idx, length(object) + 1)))
     }    
     dimLen <- length(dim(object))
     object <- if (dimLen == 0) na.lomf.0(object) else apply(object, dimLen, na.lomf.0)
     if (na.rm) na.trim(object, sides = "left", is.na = "all") else object
     
     object <- as.data.frame(object)
     for (c in 1:ncol(object)){
          if (numCols[c]){
               object[,c] <- as.numeric(object[,c])    
          }
     }
     return(object)
}