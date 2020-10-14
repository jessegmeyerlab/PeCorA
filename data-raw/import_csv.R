## code to prepare `import_csv` dataset goes here
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

#t<-import_processed_data("PeCorA_noZ.csv")
usethis::use_data(import_processed_data, overwrite = TRUE)
