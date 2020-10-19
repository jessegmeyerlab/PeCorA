#' @title import_processed_data_for_PeCorA_LFQ
#' @description reads peptide.txt MaxQuant output example file inside the package and process them for PeCorA.
#' @param LFQfile peptide.txt MaxQuant output file
#' @param condition1 experimental condition 1
#' @param condition2 experimental condition 2
#' @param condition3 experimental condition 3
#' @author Maria Dermit <maria.dermit@qmul.ac.uk>
#' @return dataframe ready for PeCorA analysis
#' @details output files containing columns: Peptide.Modified.Sequence", Protein, variable,value,Condition,BioReplicate
#' @examples
#' \dontrun{
#' if(interactive()){
#' t <- import_processed_data_for_PeCorA_LFQ( LFQfile = "peptides.txt",condition1="CONTROL",condition2="_COVID",condition3="NON.COVID")
#'  }
#' }
#' @rdname import_processed_data_for_PeCorA_LFQ
#' @importFrom reshape melt
#'@export
import_processed_data_for_PeCorA_LFQ <- function(LFQfile,condition1,condition2,condition3){
  mypath <- system.file(
    "extdata",
    LFQfile,
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
  #pmelt$Condition[grep(pmelt$variable, pattern="_COVID19")]<-"COVID"
  #pmelt$Condition[grep(pmelt$variable, pattern="NON.COVID19")]<-"NONCO"
  #pmelt$Condition[grep(pmelt$variable, pattern="control")]<-"CNTRL"

  pmelt$Condition[grep(pmelt$variable, pattern=condition1)]<-gsub("[^0-9A-Za-z///' ]","" , condition1 ,ignore.case = TRUE)
  pmelt$Condition[grep(pmelt$variable, pattern=condition2)]<-gsub("[^0-9A-Za-z///' ]","" , condition2 ,ignore.case = TRUE)
  pmelt$Condition[grep(pmelt$variable, pattern=condition3)]<-gsub("[^0-9A-Za-z///' ]","" , condition3 ,ignore.case = TRUE)



  levels(as.factor(pmelt$Condition))

  splt_last = function(x){tail(unlist(strsplit(x, "_")), n=1)}
  charnames<-as.character(pmelt$variable)
  charnames[1]
  splt_last(charnames[1])
  pmelt$BioReplicate<-lapply(FUN=splt_last, charnames)
  head(pmelt)
  df <- apply(pmelt,2,as.character)
  head(df)
  df[,1]
  colnames(df)
  colnames(df)[1]<-"Peptide.Modified.Sequence"
  colnames(df)[2]<-"Protein"
  df<-as.data.frame(df)
  return(df)
}

