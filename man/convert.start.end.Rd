% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sizing.R
\name{convert.start.end}
\alias{convert.start.end}
\title{Convert peaks start and end points to base-pair sizing units.}
\usage{
convert.start.end(peaks, data)
}
\arguments{
\item{peaks}{List of peak tables that you want to convert the start and end values to basepairs}

\item{data}{List of data tables genereated by read.data and modified by apply.size.}
}
\value{
The input peak table is returned with added columns for the base pair sized start and end values.
}
\description{
This is a simple function to convert the start and end values of a peak table to base-pair units.
}
\seealso{
\code{\link{get.models, precompute.avgs}}
}
