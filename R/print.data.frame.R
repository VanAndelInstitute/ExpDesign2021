#' pretty-print data.frames like IRanges or Julia do 
#' 
#' @param df  a data.frame
#' 
#' @import IRanges
#' 
#' @export 
#'
print.data.frame <- function(df) {
  if (ncol(df) > 0 && require("IRanges")) {
    prev.max.print <- getOption("max.print")
    on.exit(options(max.print=prev.max.print))
    options(max.print=ncol(df) * 20)
    x <- capture.output(print(as(df, "DataFrame")))
    cat(sub("DataFrame", "data frame", x[[1]]), x[-1], sep="\n")
  } else {
    base::print.data.frame(df)
  } 
}
