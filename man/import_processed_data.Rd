% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_processed_data.R
\name{import_processed_data}
\alias{import_processed_data}
\title{import_processed_data}
\usage{
import_processed_data(x)
}
\arguments{
\item{x}{csv file}
}
\value{
dataframe ready for PeCorA analysis
}
\description{
reads files already processed for PeCorA inside the package
}
\details{
file containing columns: Peptide, Protein, Peptide.Modified.Sequence,Begin.Pos,End.Pos,Condition,BioReplicate,Normalized.Area
}
\examples{
\dontrun{
if(interactive()){
 t<-import_processed_data("PeCorA_noZ.csv")
 }
}
}
\author{
Maria Dermit \href{mailto:maria.dermit@qmul.ac.uk}{maria.dermit@qmul.ac.uk}
}
