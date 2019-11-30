# Jake Yeung
# MathFunctions.R
# Math funtions here  
# 2016-10-21

# http://www.win-vector.com/blog/2012/03/modeling-trick-the-signed-pseudo-logarithm/
pseudoLog10 <- function(x) { asinh(x/2)/log(10) }

