#' @title PeCorA_preprocessing
#' @description steps to pre-process data for PeCorA
#' @param t data in PeCorA format
#' @param area_column_name column number with peak areas
#' @param threshold_to_filter threshould value of peak areas to filter
#' @param control_name control reference
#' @author Maria Dermit <maria.dermit@qmul.ac.uk>
#' @return dataframe output of pre-process data ready for PeCorA analalysis
#' @details Peak areas are log transformed and scaled to center.After global scaling, each peptide scaling was performed to center each peptide relative to the mean of the control group’s peak area.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  scaled_peptides <- PeCorA_preprocessing(t,
#'                                          area_column_name=8,
#'                                          threshold_to_filter=100,
#'                                          control_name="cntrl")
#'  }
#' }
#' @rdname PeCorA_preprocessing
#' @importFrom standardize scale_by
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats p.adjust
#' @export
PeCorA_preprocessing <- function (t,area_column_name,threshold_to_filter,control_name){
  # make new column that combines petide and charge
  t["modpep_z"] = paste(t$Peptide.Modified.Sequence, "all", sep = "_")
  # make column that combines condition and replicate
  t["condition_rep"] = paste(t$Condition, t$BioReplicate, sep = "_")
  n_reps_total <- length(unique(t$condition_rep))
  ## filter the data to remove NA, <100 area, and then reps without all measures
  if(suppressWarnings(length(which(is.na(as.numeric(t[,area_column_name]))==TRUE))>0)){
    t <- t[-suppressWarnings(which(is.na(as.numeric(t[,area_column_name]))==TRUE)),]
  }
  t <- t[-which(as.numeric(t[,area_column_name])<=threshold_to_filter),]
  t <- t[t$modpep_z %in% names(which(table(t$modpep_z)==n_reps_total)),]


  ## scale the data and then plot
  t$ms1log2<-log2(as.numeric(t[,area_column_name]))
  t$ms1scaled <- standardize::scale_by(ms1log2 ~ Condition*BioReplicate, t)

  ########### data now scaled to 0 center ###################
  ######I am going to do a subset #######
  ###  next subtract the mean value of control from each peptide
  ### takes 3 minutes for 31k peptides, 10k/min
  ### change to vectors
  #control_name="cntrl"
  #control_name="Condition1"
  #control_name="NONCO"

  ms1adj <- rep(0, nrow(t))
  i=1
  t_cntrl<-t[t$Condition==control_name,]
  ms1scaled_cntrl<-t_cntrl[,"ms1scaled"]
  ms1scaled_full<-t[,"ms1scaled"]
  t0<-Sys.time()
  print("scaling peptides to control == 0")
  pb <- txtProgressBar(min = 0, max = length(unique(t$modpep_z)), style = 3)
  for(x in unique(t$modpep_z) ){
    ms1adj[which(t$modpep_z==x)] <- suppressWarnings(ms1scaled_full[which(t$modpep_z==x)] - mean(ms1scaled_cntrl[which(t_cntrl$modpep_z==x)]))  # subtract the ave of control
    setTxtProgressBar(pb, i)
    i=i+1
  }
  print(" ")
  print(paste("peptide scaling completed in", round(Sys.time()-t0,2), "minutes", sep=" "))
  t$ms1adj<-ms1adj
  t
}
