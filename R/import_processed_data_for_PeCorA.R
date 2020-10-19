#' @title import_processed_data_for_PeCorA
#' @description reads peptide.txt MaxQuant output example file inside the package and process them for PeCorA
#' @param x peptide.txt MaxQuant output file
#' @author Maria Dermit <maria.dermit@qmul.ac.uk>
#' @return dataframe ready for PeCorA analysis
#' @details output files containing columns: Peptide.Modified.Sequence", Protein, variable,value,Condition,BioReplicate
#' @examples
#' \dontrun{
#' if(interactive()){
#' t<-import_processed_data_for_PeCorA("peptides.txt")
#'  }
#' }
#' @rdname import_processed_data_for_PeCorA
#' @importFrom reshape melt
#'@export
import_processed_data_for_PeCorA <- function(x){
  mypath <- system.file(
    "extdata",
    x,
    package = "PeCorA"
  )

  p <- read.delim(
    mypath,
    stringsAsFactors = FALSE,
  )
  grep(colnames(p), pattern="LFQ")
  colnames(p)[1:50]
  idcol = c("Sequence", "Leading.razor.protein")

  pmelt <-reshape::melt(p,
               id.vars=idcol,
               measure.vars=grep(colnames(p), pattern="LFQ"))
  head(pmelt)

  levels(as.factor(pmelt$variable))

  pmelt$Condition = rep(0, nrow(pmelt))
  pmelt$Condition[grep(pmelt$variable, pattern="_COVID19")]<-"COVID"
  head(pmelt)
  pmelt$Condition[grep(pmelt$variable, pattern="NON.COVID19")]<-"NONCO"
  pmelt$Condition[grep(pmelt$variable, pattern="control")]<-"CNTRL"



  levels(as.factor(pmelt$Condition))

  splt_last = function(x){tail(unlist(strsplit(x, "_")), n=1)}
  charnames<-as.character(pmelt$variable)
  charnames[1]
  splt_last(charnames[1])
  pmelt$BioReplicate<-lapply(FUN=splt_last, charnames)
  head(pmelt)
  ?write.table
  df <- apply(pmelt,2,as.character)
  head(df)
  df[,1]
  colnames(df)
  colnames(df)[1]<-"Peptide.Modified.Sequence"
  colnames(df)[2]<-"Protein"
  return(df)
}

