#' @title import_processed_data
#' @description reads files already processed for PeCorA inside the package
#' @param x csv file
#' @author Maria Dermit <maria.dermit@qmul.ac.uk>
#' @return dataframe ready for PeCorA analysis
#' @details file containing columns: Peptide, Protein, Peptide.Modified.Sequence,Begin.Pos,End.Pos,Condition,BioReplicate,Normalized.Area
#' @examples
#' \dontrun{
#' if(interactive()){
#'  t<-import_processed_data("PeCorA_noZ.csv")
#'  }
#' }
#' @rdname import_processed_data
#' @importFrom utils read.csv
#'@export
import_processed_data <- function(x){
  mypath <- system.file(
    "extdata",
    x,
    package = "PeCorA"
  )

  file <- read.csv(
    mypath,
    stringsAsFactors = FALSE,
    encoding = "UTF-8"
  )
  return(file)
}

