# 
# h i s t o r i z e . R  Historize function
# 


Historize <- function(fun, len = 0, ...) {
    stopifnot(is.function(fun))
    if (len != floor(len) || len < 0)
        stop("Argument 'len' must be a non-negative integer.")
    local({
        H <- numeric()
        myFun <- function(x, ...) {
            if(missing(x)) {
                return(H)
            } else if (is.null(x)) {
                H <<- numeric()
                return(invisible(NULL))
            } else {
                y <- fun(x, ...)
                if (len == 0) {
                    H <<- c(H, y)
                } else {
                    if (length(x) != len)
                        stop("Incorrect parameter length.")
                    H <<- rbind(H, c(x, y))
                }
                return(y)
            }
        }
        return(myFun)
    })
}
