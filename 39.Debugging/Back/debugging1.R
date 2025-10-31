# https://www.youtube.com/watch?v=QQxTf3o07NU
fun_b <- function(y){
    y * y 
}


fun_a <- function(x){
    x_squared <- fun_b(x)
    result <- x * x_squared
    browser()
    for(value in 1:10){
        result <- value * result
    }
    return(result)
}

z=3
debug(fun_a)
#fun_a(3)

fun_a(z)
