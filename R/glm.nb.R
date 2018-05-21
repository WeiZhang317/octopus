# Copyright 2018 Kliebenstein Lab - UC Davis Plant Sciences. All rights reserved.
#
# You can redistribute it and/or modify it under the terms of the
# GNU General Public License Version 2. You should have received a copy of the
# GNU General Public License Version 2 along with Octopus project.
# If not, you can get one at https://github.com/WeiZhang317/octopus
#


#' Fit a Negative Binomial Generalized Linear Model
#'
#' Fit glm.nb , add anova to result, add lsm to result.\cr
#' Before fitting, \code{data} and \code{factors} will be merged by row.names, mismatching rows will be removed.\cr
#' If formula start with ~, will fit model with each data column of data as respond y.
#'
#' @param formula formula
#' @param data normalized reads, expecting gene as column, sample as row
#' @param factors factors, expecting factor name as column, sample as row
#' @param specs4lsmeans specs for lsmeans
#'
#' @return model if \code{formula} has respond y; flattened anova table, lsm table and se table if formula starts with ~.
#' @export
#'
#' @examples --------------------
octopus.glm.nb <- function(formula,data,factors,specs4lsmeans=NULL) {
  if(is.null(data) || length(data) == 0 ){ stop("data is empty.") }
  if(is.null(factors) || length(factors) == 0 ){ stop("factors is empty.") }
  if(is.null(formula) ){ stop("formula is empty.") }
  info.self <- create.info("octopus.glm.nb")
  do.glm.nb <- function(...){
    result <- tryCatch({ do.call(glm.nb,list(...)) }, error=function(e) { print(info.self,"glm.nb",e) ; NULL }) # glm.nb
    if(is.null(result)) { return(NULL) } # Skip

    # Pull p.values/variances/DFs from the model. Can't pull variances/deviances with Type II SS, so stuck with Type I
    result$anova <- tryCatch({ anova(result) }, error=function(e) { print(info.self,"anova",e) ; NULL }) # anova
    if(is.null(result$anova)) { return(NULL) } # Skip

    if(!is.null(specs4lsmeans)) { # lsmeans
      result$lsmeans <- tryCatch({lsmeans(result,specs4lsmeans)}, error=function(e) { print(info.self,"lsmeans",e) ; NULL })
    }
    result
  }

  genes <- colnames(data)
  data <- merge(factors,data, by = "row.names")
  row.names(data) <- data[,1]
  data <- data[,-1]

  formula_len <- length(as.character(formula(formula)))
  if(formula_len==3){
    return(do.glm.nb(formula,data))
  }

  formula <- as.character(formula(formula))
  formula <- formula(paste("data[,gene] ~",formula[2])) # build formula

  ret <- NULL
  sp <- create.counter(length(genes)," glm.nb ")
  for(gene in genes){
    sp <- sp + 1
    print(sp)
    result <- do.glm.nb(formula,factors)
    anova.flat <- data.frame(p=c(result$anova[,1], result$anova[,2], result$anova[,3], result$anova[,4], result$anova[,5])
                             ,row.names = apply(expand.grid(rownames(result$anova), colnames(result$anova)),1,paste,collapse="_"))
    if(!is.null(specs4lsmeans)) { lsm <- summary(result$lsmeans) }
    if(is.null(ret)){
      ret = list(anova=NULL,lsm=NULL,se=NULL)
      # initialize memory in advance to avoid data copying
      ret$anova <- data.frame(matrix(nrow = dim(anova.flat)[1] ,ncol = length(genes) ,dimnames=list(row.names(anova.flat),genes)))
      if(!is.null(specs4lsmeans)) {
        ret$lsm <- matrix(nrow = dim(lsm)[1] ,ncol = length(genes) ,dimnames=list(row.names(lsm),genes))
        ret$lsm <- cbind(lsm[1:(length(names(lsm))-5)],ret$lsm) # attach factor colmuns
        ret$se <- matrix(nrow = dim(lsm)[1] ,ncol = length(genes),dimnames=list(row.names(lsm),genes))
        ret$se <- cbind(lsm[1:(length(names(lsm))-5)],ret$se)  # attach factor colmuns
      }
    }
    ret$anova[gene] <- anova.flat[1]
    if(!is.null(specs4lsmeans)) {
      ret$lsm[gene] <- lsm$lsmean
      ret$se[gene] <- lsm$SE
    }
  }
  ret
}


#' Normalize
#'
#' NA will be treated as 0
#'
#' @param x number matrix, expecting gene as row, sample as column
#' @param dist_dir
#'
#' @return number matrix
#' @export
#'
#' @examples tutorial("octopus.glm.nb") # to print tutorial
#' @seealso \link[octopus]{octopus.glm.nb}
octopus.normalize <- function(x,dist_dir="results/"){
  #Normalize using edgeR
  info.self <- create.info("octopus.normalize")
  x[is.na(x)] <- 0
  x <- DGEList(x) ###Transform to DGEobject
  print(info.self,"TMM Normalization ...")
  x <- calcNormFactors(x, method ="TMM")
  x <- estimateCommonDisp(x, verbose=TRUE)
  print(info.self,"write out TMM Normalization info to norm.reads.info.csv")
  dir.create(dist_dir, showWarnings = FALSE, recursive = TRUE )
  write.csv(x$samples,file = paste0(dist_dir,"norm.reads.info.csv"))
  as.data.frame(x$pseudo.counts)
}

#' @export
tutorial.octopus.glm.nb <- function(topic=NULL) {
  tutorial.octopus.glm.nb.do <- function(){
    data(sample_keys) # facotrs
    for(cur_col in colnames(sample_keys)){ sample_keys[,cur_col] <- as.factor(sample_keys[,cur_col]) }

    data(sample_reads) # reads
    # reads normalization
    norm.reads <- octopus.normalize(sample_reads)

    print(row.names(sample_keys))
    print(row.names(norm.reads))
    # transpose it to meet input requirements of founction octopus.glm.nb()
    norm.reads <- data.frame(t(norm.reads))
    print(row.names(norm.reads))

    # fit model with first gene column, column name is Bcin01g00040.1
    result <- octopus.glm.nb(Bcin01g00040.1 ~ Experiment + Experiment/GrowingFlat + Experiment/GrowingFlat/AgarFlat + Isolate + HostGenotype + HostGenotype*Isolate
                             ,norm.reads
                             ,sample_keys)
    print(result$anova)

    # calculate LSMean
    result$lsmeans <- lsmeans(result,~ Isolate | HostGenotype)
    summary(result$lsmeans)

    # fit model and pull LSMean in one call
    result <- octopus.glm.nb(Bcin01g00040.1 ~ Experiment + Experiment/GrowingFlat + Experiment/GrowingFlat/AgarFlat + Isolate + HostGenotype + HostGenotype*Isolate
                             ,norm.reads
                             ,sample_keys
                             ,specs4lsmeans = ~ Isolate | HostGenotype)
    print(result$anova)
    summary(result$lsmeans)

    # iterate through all data
    result <- octopus.glm.nb( ~ Experiment + Experiment/GrowingFlat + Experiment/GrowingFlat/AgarFlat + Isolate + HostGenotype + HostGenotype*Isolate
                              ,norm.reads
                              ,sample_keys
                              ,specs4lsmeans = ~ Isolate | HostGenotype)
    print(result$anova)
    print(result$lsm)
    print(result$se)

  }

  print(tutorial.octopus.glm.nb.do)

}

