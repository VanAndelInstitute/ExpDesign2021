#' list memory-resident objects by footprint
#' 
#' @param ... optional arguments to pass on to .ls.objects()
#' @param n   optional integer specifying how many to show (default is 10)
#' 
#' @return  a pretty-printed list of objects and their sizes 
#' 
#' @export
#' 
lsos <- function(n=10, ...) {
  .ls.objects(order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

.ls.objects <- function(pos=1, pattern="", order.by, decreasing=F, head=F, n=5){
  napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos=pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by)) 
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head == TRUE) out <- head(out, n)
  return(out)
}
