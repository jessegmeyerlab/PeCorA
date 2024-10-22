% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PeCorA_plotting.R
\name{PeCorA_plotting}
\alias{PeCorA_plotting}
\title{PeCorA_plotting}
\usage{
PeCorA_plotting(w, u, v)
}
\arguments{
\item{w}{disagree peptides}

\item{u}{selection of disagree peptides}

\item{v}{scaled peptides}
}
\value{
PeCorA plot
}
\description{
plot PeCorA output
}
\details{
This function generates boxplots for a given protein representing the peptides that are statistically different from the quantities of all the other peptides in green and all other peptides in grey.
}
\examples{
\dontrun{
if(interactive()){
 PeCorA_plotting_plot<-PeCorA_plotting(disagree_peptides,disagree_peptides[1,],scaled_peptides)
 }
}
}
\author{
Maria Dermit \href{mailto:maria.dermit@qmul.ac.uk}{maria.dermit@qmul.ac.uk}
}
