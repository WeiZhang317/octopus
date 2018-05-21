# Copyright 2018 Kliebenstein Lab - UC Davis Plant Sciences. All rights reserved.
#
# You can redistribute it and/or modify it under the terms of the
# GNU General Public License Version 2. You should have received a copy of the
# GNU General Public License Version 2 along with Octopus project.
# If not, you can get one at https://github.com/WeiZhang317/octopus
#


#' Print Tutorial
#'
#' Print tutorial for function/class and return tutorial code as a function. It is a generic function.
#'
#' @param topic tutorial topic, it is usually a function name.
#'
#' @return tutorial code as a function
#' @export
#'
#' @examples --------------------
#'
#' tutorial("tutorial") # print tutorial to console
#'
tutorial <- function(topic) {
  UseMethod("tutorial")
}

#' @export
tutorial.default <- function(topic) {
  if(!is.null(topic)){
    warning("No tutorial found for function/class ",class(topic))
  }
  invisible(NULL)
}

#' @export
tutorial.character <- function(topic) {
  tutorial(structure("", class = topic))
}

#' @export
tutorial.function <- function(topic) {
  tutorial(as.character(substitute(topic)))
}

#' @export
tutorial.tutorial <- function(topic=NULL) {
  tutorial.tutorial.do <- function(){
    # access tutorial code for class info
    tutorial("info")          # use class name string
    tutorial(create.info())   # use object

    tutorial(create.info)     # use function name
    tutorial("create.info")   # use function name string

    tutorial_func <- tutorial("info")
    tutorial_func() # run tutorial code, suggest read it first
  }
  print(tutorial.tutorial.do)
}


