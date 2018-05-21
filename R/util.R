# Copyright 2018 Kliebenstein Lab - UC Davis Plant Sciences. All rights reserved.
#
# You can redistribute it and/or modify it under the terms of the
# GNU General Public License Version 2. You should have received a copy of the
# GNU General Public License Version 2 along with Octopus project.
# If not, you can get one at https://github.com/WeiZhang317/octopus
#


.info.instance <- new.env(parent=emptyenv())
assign("info.OFF", FALSE, envir=.info.instance)

#' @export
info.enable <- function(enable = TRUE) {
  if (enable) {
    assign("info.OFF", FALSE, envir=.info.instance)
  } else{
    assign("info.OFF", TRUE, envir=.info.instance)
  }
}

#' Print info with timestamp
#'
#' Print debug message, progress reports... with timestamp.
#' \itemize{
#'   \item Default on
#'   \item Turn off by \code{info.enable(FALSE)}
#'   \item Turn on by \code{info.enable(TRUE)}
#' }
#' @param prefix prefix of message.
#'
#' @return object of class \code{info}
#' @export
#' @examples --------------------
#'
#' tutorial("info") # print tutorial to console
#'
create.info <- function(prefix = "") {
  structure(prefix, class = "info")
}

#' @export
print.info <- function(x, ...) {
  if (!.info.instance$info.OFF) {
    cat(paste0(Sys.time()), "[", x, "\t] INFO:" , paste(...) , "\n")
  }
}

#' @export
tutorial.info <- function(topic=NULL) {
  tutorial.info.do <- function(){
    info.self <- create.info() #create object of info class
    print(info.self,"print message without prefix.")

    info.self <- create.info("tutorial.info") #create object of info class with "tutorial.info" as prefix
    print(info.self,"print message with prefix, use identifier for prefix usually.")

    # turn it off
    info.enable(FALSE) #turn information printing off all together
    print(info.self,"will print nothing from now on.")
    print(info.self,"conform print nothing.")

    # turn it on
    info.enable(TRUE) #turn information printing on all together
    print(info.self,"start printing again.")
  }
  print(tutorial.info.do)
}

#' @export
tutorial.create.info <- tutorial.info

#' Load package, install automatically if missing.
#'
#' Load \code{package}, install highest version
#' from \href{https://www.bioconductor.org/install/}{bioconductor}
#' or \href{http://cran.rstudio.com/}{rstudio} if missing.
#' @param package the name of a package
#'
#' @return invisible \code{TRUE}
#' @export
#'
#' @examples --------------------
#'
#' tutorial("auto.library") # print tutorial to console
#'
auto.library <- function(package) {
  info.self <- create.info("library.auto")
  package <- as.character(substitute(package))
  if (package %in% rownames(installed.packages()) == FALSE) {
    print(info.self, "missing library:", package)
    print(info.self, "Trying to install it from bioconductor")
    source("https://www.bioconductor.org/biocLite.R")
    biocLite(package)
  }
  if (package %in% rownames(installed.packages()) == FALSE) {
    print(info.self, "Trying to install it from cran.rstudio.com")
    install.packages(package, repos = "http://cran.rstudio.com/")
  }
  if (package %in% rownames(installed.packages()) == FALSE) {
    print(info.self, "failed to install ", package)
  } else{
    library(package, character.only = TRUE)
  }
  invisible(TRUE)
}

#' @export
tutorial.auto.library <- function(topic=NULL) {
  tutorial.auto.library.do <- function(){
    auto.library(stringr) #use package name
    auto.library("stringr") #use package name string
    auto.library("package_no_exist") # Print warning message when package is not found
  }
  print(tutorial.auto.library.do)
}


#' Slice list into sublists by n
#'
#' @param x list
#' @param n desired length of sublist
#'
#' @return A list of sublists
#' @export
#'
#' @examples --------------------
#'
#' tutorial("slice") # print tutorial to console
#'
slice <- function(x,n){
  len <- length(x)
  if (is.null(len) || len< n) {
    return(x)
  }
  ret <- vector(mode = "list", length = ceiling(len/n))
  for(i in 1:floor(len/n)){
    b <- (i - 1) * n + 1
    e <- i * n
    ret[i] <- list(x[b:e])
  }
  if(i != ceiling(len / n)) {
    b <- i * n + 1
    e <- len
    ret[i+1] <- list(x[b:e])
  }
  ret
}

#' @export
tutorial.slice <- function(topic=NULL) {
  tutorial.slice.do <- function(){
    x <- 5:505
    slice(x,100) # 5 sublists
    slice(x,99)  # 6 sublists
    slice(x,101) # 4 sublists
  }
  print(tutorial.slice.do)
}

create.counter <- function(total = 0, name = "") {
  if (total < 0) {
    total <- 0
  }
  structure(
    list(
      total = total,
      name = name,
      count = 0,
      time.cur = Sys.time() -10,
      time.start = Sys.time()
    ),
    class = "counter"
  )
}

"+.counter" <- function(x, n) {
  x$count <- x$count + n
  x
}

print.counter <- function(x, prefix = "") {

  if (.info.instance$info.OFF) {
    return(invisible(NULL))
  }

  if (as.integer(difftime(Sys.time(), x$time.cur, units = "secs")) < 5) {
    return(invisible(NULL))
  }

  x$time.cur <- Sys.time()
  time.passed <- as.integer(difftime(Sys.time(), x$time.start, units = "secs"))
  cps <- as.integer(x$count) / time.passed
  cpn <- paste0(round(cps), "/s")
  cpm <- (60 * as.integer(x$count)) / time.passed
  if (cps < 1) {
    cpn <- paste0(round(cpm) , "/m")
  }
  cph <- (3600 * as.integer(x$count)) / time.passed
  if (cps < 1) {
    cpn <- paste0(round(cph) , "/h")
  }
  if (x$total < x$count) {
    x$total <- x$count
  }
  eta <- Sys.time() + (x$total - x$count) / cps
  print(
    paste0(Sys.time(),x$name," : ",
      prefix, " [" , x$count, " of ",x$total,"]",
      " ETA: ", eta ,
      ", speed: ",cpn ,
      ", time cost: ", time.passed," seconds."
    )
  )
}

