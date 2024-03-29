% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_MC_Report.R
\name{get_MC_Report}
\alias{get_MC_Report}
\title{Report about Missed Cleavages}
\usage{
get_MC_Report(input_list, metric = c("absolute", "percentage"))
}
\arguments{
\item{input_list}{A list with data frames and respective missed cleavage information.}

\item{metric}{\code{"absolute"} for absolute numbers or \code{"percentage"} for displaying percentages. Default is absolute.}
}
\value{
This function returns a list. For each analysis a respective data frame including information of missed cleavages is stored in the generated list.
\itemize{
\item Analysis - analysis name.
\item Missed.Cleavage - categorical entry with number of missed cleavages.
\item mc_count - number of missed cleavages per categorical missed cleavage entry - absolute or in percentage.
}
}
\description{
Generates report with information about number of missed cleavages
}
\details{
For each submitted data a report is generated with information about the number of missed cleavages.
}
\examples{
# Load libraries
library(tibble)
library(magrittr)
library(stringr)

# Example data
data <- list(
Spectronaut = list(
  filename = "C",
  software = "Spectronaut",
  data = list(
     "Spectronaut" = tibble::tibble(
        Stripped.Sequence_mpwR = c("A", "B", "C", "D", "E"),
        Missed.Cleavage_mpwR = c(0, 1, 1, 2, 2)
     )
  )
)
)

# Result
output <- get_MC_Report(
  input_list = data,
  metric = "absolute"
)
}
\author{
Oliver Kardell
}
