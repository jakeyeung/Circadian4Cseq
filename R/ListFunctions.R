# appending to list is super slow
# thank you http://stackoverflow.com/questions/2436688/append-an-object-to-a-list-in-r-in-amortized-constant-time
# Jan Kanis:
# use this as follows
# > l <- expandingList()
# > l$add("hello")
# > l$add("world")
# > l$add(101)
# > l$as.list()
# [[1]]
# [1] "hello"
# 
# [[2]]
# [1] "world"
# 
# [[3]]
# [1] 101

expandingList <- function(capacity = 10) {
  buffer <- vector('list', capacity)
  length <- 0
  
  methods <- list()
  
  methods$double.size <- function() {
    buffer <<- c(buffer, vector('list', capacity))
    capacity <<- capacity * 2
  }
  
  methods$add <- function(val) {
    if(length == capacity) {
      methods$double.size()
    }
    
    length <<- length + 1
    buffer[[length]] <<- val
  }
  
  methods$as.list <- function() {
    b <- buffer[0:length]
    return(b)
  }
  
  methods
}

linkedList <- function() {
  head <- list(0)
  length <- 0
  
  methods <- list()
  
  methods$add <- function(val) {
    length <<- length + 1
    head <<- list(head, val)
  }
  
  methods$as.list <- function() {
    b <- vector('list', length)
    h <- head
    for(i in length:1) {
      b[[i]] <- head[[2]]
      head <- head[[1]]
    }
    return(b)
  }
  methods
}