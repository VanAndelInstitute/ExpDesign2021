#' A version of identical() that accepts arbitrarily many inputs.  
#' (It turns out that Reduce(identical, ...) doesn't do what you might like.)
#' 
#' @param ...   (list) The things to test for identicality.
#' @return      (boolean) Whether they are all identical.
#' 
#' @examples
#' 
#' identical(1:3*1, c(1, 2, 3), seq(3, 9, 3) / 3)           # TRUE
#' identical(1:3, c(1, 2, 3), seq(1, 3))                    # FALSE
#' identical(c(1, 2, 3), c(2, 4, 6) / 2, c(3, 6, 9) / 3)    # TRUE 
#' identical(1:3, c(1, 2, 3), seq_len(3))                   # FALSE 
#' identical(1:3, seq_len(3), seq_along(1:3))               # TRUE 
#'
#' @export
#' 
identical <- function(...) {
  if (length(list(...)) < 3) {
    return(base::identical(list(...)[[1]], list(...)[[2]]))
  } else { 
    eq1 <- function(x) base::identical(x, list(...)[[1]])
    return(all(unlist(lapply(list(...), eq1)) == TRUE) )
  }
}
