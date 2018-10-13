# Copyright 2018 Kliebenstein Lab - UC Davis Plant Sciences. All rights reserved.
#
# You can redistribute it and/or modify it under the terms of the
# GNU General Public License Version 2. You should have received a copy of the
# GNU General Public License Version 2 along with Octopus project.
# If not, you can get one at https://github.com/WeiZhang317/octopus
#


#' Extract Short Reads
#'
#' Extract short reads from \code{seq_files} using \href{http://bowtie-bio.sourceforge.net/manual.shtml}{bowtie}.\cr
#' Output 2 files, \code{dist_dir}/orig.reads.csv and \code{dist_dir}/orig.reads.info.csv \cr
#' Remove \code{dist_dir}/reusing/ folder if don't want reuse data from that folder \cr
#' Setup ./tmp/ as an RAM disk folder will avoid lots of disk IOs, speed things up and protect you SSD.\cr
#' Example for linux :
#' \itemize{
#'   \item cd to curent folder.
#'   \item create ram disk in console : mount -t tmpfs -o size=4g tmpfs ./tmp/
#' }
#' Example for mac :
#' \itemize{
#'   \item cd to curent folder.
#'   \item mkdir -p tmp
#'   \item sudo mount -t tmpfs -o size=4096M tmpfs ./tmp/
#' }
#' For windows, there are a number of RAM disk softerwares you can use.
#'
#' @param seq_files  Sequencing files , accepts .fastq or .gz format for files.\cr
#' If \code{seq_files} is a \code{List} or \code{Vector}, index of \code{seq_files} is assumed to be sampe names of sequencing files.\cr
#' If \code{seq_files} is a \code{data.frame}, row.names is assumed to be sampe names of sequencing files, the first column is assumed to be sequencing files, when \code{type} is \code{paired} the second column is assumed to be the second mate pair sequences.
#' @param references A comma-separated list of FASTA files containing the reference sequences to be aligned to
#' @param type Could be one of c("single", "paired", "crossbow").\cr
#' If single, the input sequences are interpreted as single reads.\cr
#' If paired, they are supposed to be mate pair reads.\cr
#' If crossbow, they are considered to be Crossbow-style reads.
#' @param ... Additional arguments to be passed on to the binaries. See ... of \link[Rbowtie]{bowtie}
#' @param dist_dir folder for result file orig.reads.csv and orig.reads.info.csv
#' @export
#'
#' @examples --------------------
octopus.short_reads <- function(seq_files,references,...,type="single",dist_dir="results/"){
  # check input data
  if (is.data.frame(seq_files)){
    sample_names <- row.names(seq_files)
    if(type=="paired") {
      sample_seq_files <- seq_files[,c(1,2)]
    }else{
      sample_seq_files <- seq_files[,1]
    }
  }else if(is.vector(seq_files) || is.list(seq_files) ){
    sample_seq_files <- seq_files
    sample_names <- 1:length(seq_files)
    if(type=="paired") {
      stop("")
    }
  }else{
    print("seq_files must be one of vector/list/dta.frame type.")
    return(NULL)
  }

  # create directory for results
  short_reads_path <- paste0(dist_dir,"reusing/short_reads/")
  dir.create(short_reads_path, showWarnings = FALSE, recursive = TRUE )
  dir.create("tmp", showWarnings = FALSE, recursive = TRUE )

  refs_sign_str <- NULL
  for (str in references) {
    if(is.null(refs_sign_str)){
      refs_sign_str <- readChar(str,10000000)
    }else{
      refs_sign_str <- c(refs_sign_str,readChar(str,10000000))
    }
  }

  # detect changing of input data
  signs <- list(refs=sha1(refs_sign_str),seq_files=sha1(seq_files),ext_args=sha1(list(...)),type=type) # TODO list(2) -> list(...)
  signs_file <- paste0(short_reads_path,"sha1_signatures")
  if(!file.exists(signs_file)){ write.csv(signs,signs_file,row.names = FALSE) }
  signs_old <- read.csv(signs_file)
  is_mismatch = signs$refs != signs_old$refs || signs$seq_files != signs_old$seq_files || signs$ext_args != signs_old$ext_args || signs$type != signs_old$type
  if(is_mismatch){
    warning("parameter has been changed...
            remove folder ",short_reads_path," and restart if changes can affect results.
            remove ",signs_file," and restart to stop this message if changes won't affect results.")
  }
  mismatch <- function(){
    if(!is_mismatch){ return() }
    print("parameter has been changed...")
    print(paste("remove folder ",short_reads_path," and restart if changes can affect results."))
    print(paste("remove file ",signs_file," and restart to stop this message if changes won't affect results."))
  }
  mismatch()

  # build index
  refs <- octopus.bowtie_build(references,dist_dir)
  # alignment
  for(i in 1:length(sample_names)){
    sample_name <- sample_names[i]
    if(type=="paired") {
      sample_seq_file <- sample_seq_files[i,]
    } else {
      sample_seq_file <- sample_seq_files[i]
    }

    print(paste("process sample",sample_name,", seq_file",sample_seq_file))
    reads_file <- paste0(short_reads_path,sample_name,".reads.csv")

    if(file.exists(reads_file)) { # try reuse
      print(paste("reusing sample",sample_name,", file ",reads_file))
      mismatch()
    }else{ #align
      ret <- octopus.bowtie_align(sample_seq_file,refs,type=type,...)
      colnames(ret) <- c("gene",sample_name) #set column name
      print(paste("saving into file ",reads_file))
      write.csv(ret,file=reads_file,row.names=FALSE)
    }
  }
  # merge
  octopus.short_reads.merge(short_reads_path,dist_dir)
}

octopus.short_reads.merge <- function(src_dir,dist_dir="results/"){

  print(paste("start merging reads files in folder ",src_dir))
  orig.reads = NULL

  files <- list.files(src_dir, pattern="\\.reads\\.csv$", recursive = TRUE , all.files = FALSE)

  for(f in files) {
    f <- paste0(src_dir,"/",f)
    ret <- read.csv(f,check.names = FALSE)
    colnames(ret)

    if (is.null(orig.reads)) {
      orig.reads <- ret
    }else{
      orig.reads <- merge(orig.reads,ret,all=T,by="gene")
    }
    names <- names(orig.reads)
    print(paste("merge" , f , " into orig.reads. names [",names[1],"...",names[length(names)],"]"))
  }

  if(is.null(orig.reads )){ octopus.info("orig.reads is empty ...") ; return(FALSE) }
  ret_map <- orig.reads[orig.reads$gene!="*",]

  print(paste0("saving results to ",dist_dir,"orig.reads.csv"))
  dir.create(dist_dir, showWarnings = FALSE, recursive = TRUE )
  write.csv(ret_map,file=paste0(dist_dir,"orig.reads.csv"),row.names=FALSE)

  print(paste("summarize mapped and unmapped reads ..."))

  unmapped <- orig.reads[orig.reads$gene=="*",-1]
  unmapped <- t(unmapped)
  colnames(unmapped) <- "unmapped"

  mapped <- apply(ret_map[-1],2,sum,na.rm=T)
  mapped <- data.frame(mapped)

  mapped_percent <- round(mapped/(mapped+unmapped)*100,1)
  mapped_percent <- cbind(rownames(mapped_percent),mapped_percent[,1])
  colnames(mapped_percent) <- c("SampleName","Qulity")

  orig.reads.info <- data.frame(mapped_percent,mapped,unmapped)
  print(paste0("saving qulity info to ",dist_dir,"orig.reads.info.csv"))
  write.csv(orig.reads.info,file=paste0(dist_dir,"orig.reads.info.csv"),row.names=FALSE)

  print(paste("merging done!!!"))
}

octopus.bowtie_build <- function(references,dist_dir="results/"){
  # sha1 signature
  # build index
  refs_sign_str <- NULL
  for (str in references) {
    if(is.null(refs_sign_str)){
      refs_sign_str <- readChar(str,10000000)
    }else{
      refs_sign_str <- c(refs_sign_str,readChar(str,10000000))
    }
  }
  refs_sign <- sha1(refs_sign_str)
  refs <- paste0(dist_dir,"/reusing/refs/",refs_sign)
  if(dir.exists(refs)){
    print(paste("reusing bowtie index ",refs," for references ",references))
  }else{
    print(paste("build bowtie index for references ",references))
    refs_tmp <- paste0(refs,"_tmp")
    bowtie_build(references=references, outdir=refs_tmp, force=TRUE)
    file.rename(refs_tmp,refs)
  }

  # use ram disk for bowtie index
  refs_ram <- "tmp/refs"
  dir.create(refs_ram, showWarnings = FALSE, recursive = TRUE )
  file.copy(paste0(refs,"/",list.files(refs)),paste0(refs_ram,"/",list.files(refs)),overwrite = TRUE)
  paste0(refs_ram,"/index")
}


octopus.bowtie_align <- function(sample_seq_file,refs,...,type="single"){

  # unzip
  tmp_files <- NULL
  if(type=="paired") {
    # sample_seq_file <- sample_seq_files[1,]
    sample_seq_file_pair <- as.character(sample_seq_file[2])
    if(file_ext(sample_seq_file_pair) == "gz"){# unzip
      fastq_file_pair <- "tmp/orig_pair.fastq"
      gunzip(sample_seq_file_pair,fastq_file_pair,overwrite=TRUE,remove=FALSE)
      tmp_files <- c(tmp_files,fastq_file_pair)
    }else{
      if(file_ext(sample_seq_file_pair) != "fastq"){ octopus.info("expecting .gz/.fastq file.") ; return(FALSE) }
      fastq_file_pair <- sample_seq_file_pair
    }
    sample_seq_file <- as.character(sample_seq_file[1])
  }

  if(file_ext(sample_seq_file) == "gz"){# unzip
    fastq_file <- "tmp/orig.fastq"
    gunzip(sample_seq_file,fastq_file,overwrite=TRUE,remove=FALSE)
    tmp_files <- c(tmp_files,fastq_file)
  }else{
    if(file_ext(sample_seq_file) != "fastq"){ octopus.info("expecting .gz/.fastq file.") ; return(FALSE) }
    fastq_file <- sample_seq_file
  }


  #do alignment
  sam_file <- "tmp/orig.sam"
  tmp_files <- c(tmp_files,sam_file)

  # S=TRUE # SAM format ,
  # q=TRUE # query input files are FASTQ .fq/.fastq
  ext_args <- list(...)
  if(type=="paired") {
    print(paste("fist mate pair sequences:",sample_seq_file))
    print(paste("second mate pair sequences:",sample_seq_file_pair))
    ext_args <- append(list(sequences=list(fastq_file,fastq_file_pair),index=refs,outfile=sam_file, type=type,force=TRUE,S=TRUE,q=TRUE),ext_args) # ,y=TRUE
  }else{
    ext_args <- append(list(sequences=fastq_file,index=refs,outfile=sam_file, type=type,force=TRUE,S=TRUE,q=TRUE),ext_args)
  }
  do.call(bowtie,ext_args)

  #read the sam file.
  ret <- scan(sam_file,what=list(NULL,NULL,""),comment.char="@",sep="\t",flush=T)[[3]] # We only care about the third column.also discard the header info (rows starting with "@")  -- Mike Covington - Julin Maloofs lab
  file.remove(tmp_files) # clean up

  if(length(ret)==0){ print(paste("empty results ...")) ; return(ret) }
  ret <- as.data.frame(table(ret))
  colnames(ret) <- c("gene","reads") #set column name
  ret
}

tutorial.octopus.short_reads <- function(topic=NULL){
  tutorial.octopus.short_reads.do <- function(){

    # single
    references <- "seq_data/cdna/Arabidopsis_thaliana.TAIR10.25.cdna.all.fa"
    seq_files <- data.frame(seq_file=c("seq_data/1_AACGTGAT_L003_R1_001.fastq.gz","seq_data/1_AACGTGAT_L007_R1_001.fastq","seq_data/3_AACGTGAT_L003_R1_001.fastq.gz")
                            ,sample_name=c("sample1","sample2","sample3")
                            ,stringsAsFactors = FALSE)
    row.names(seq_files) <- seq_files$sample_name

    octopus.short_reads(seq_files,references
                        ,p=3 # number of alignment threads to launch
                        ,`phred33-quals`=TRUE # input quals are Phred+33
                        ,t=TRUE # print wall-clock time taken by search phases
                        ,quiet=TRUE # print nothing but the alignments
                        ,trim5=10 # trim <int> bases from 5' (left) end of reads
    )

    # paired

    references <- "seq_data/cdna/Arabidopsis_thaliana.TAIR10.25.cdna.all.fa"
    seq_files <- data.frame(seq_file=c("seq_data/A9_S1_L001_R1_001.fastq.gz","seq_data/xxx_R1_001.fastq.gz","seq_data/xxxx_R1_001.fastq.gz")
                            ,seq_file_pair=c("seq_data/A9_S1_L001_R2_001.fastq.gz","seq_data/xxx_R2_001.fastq.gz","seq_data/xxxx_R2_001.fastq.gz")
                            ,sample_name=c("sample1","sample2","sample3")
                            ,stringsAsFactors = FALSE
                            )
    row.names(seq_files) <- seq_files$sample_name

    octopus.short_reads(seq_files,references
                        ,p=3 # number of alignment threads to launch
                        ,`phred33-quals`=TRUE # input quals are Phred+33
                        ,t=TRUE # print wall-clock time taken by search phases
                        ,quiet=TRUE # print nothing but the alignments
                        ,trim5=10 # trim <int> bases from 5' (left) end of reads
                        # ,y=TRUE # more sensitive but much slower, see http://bowtie-bio.sourceforge.net/manual.shtml#bowtie-options-y
                        ,type="paired"
                        ,dist_dir="results_paired/"
    )

    # multiple referencing files
    references  <- c("seq_data/ref_sequences/Botrytisfusarivirus1.txt","seq_data/ref_sequences/BotrytisHypovirus1.txt")
    seq_files <- data.frame(seq_file=c("seq_data/1_AACGTGAT_L003_R1_001.fastq.gz","seq_data/1_AACGTGAT_L007_R1_001.fastq","seq_data/3_AACGTGAT_L003_R1_001.fastq.gz")
                            ,sample_name=c("sample1","sample2","sample3")
                            ,stringsAsFactors = FALSE)
    row.names(seq_files) <- seq_files$sample_name

    octopus.short_reads(seq_files,references
                        ,p=3 # number of alignment threads to launch
                        ,`phred33-quals`=TRUE # input quals are Phred+33
                        ,t=TRUE # print wall-clock time taken by search phases
                        ,quiet=TRUE # print nothing but the alignments
                        ,trim5=10 # trim <int> bases from 5' (left) end of reads
    )


  }

  print(tutorial.octopus.short_reads.do)

}


