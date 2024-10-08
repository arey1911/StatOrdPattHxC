% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MovingAverage.R
\name{mov.av}
\alias{mov.av}
\title{Moving average}
\usage{
mov.av(x, order)
}
\arguments{
\item{x}{time series}

\item{order}{window size (default: 1, no filter)}
}
\value{
the filtered time series, smaller that the input (only values are returned)
}
\description{
Computes the moving average of a time series over a specified window size
}
\details{
Used only for the examples
}
\examples{
set.seed(1234567890, kind="Mersenne-Twister")
x <- rnorm(1000)
x.ma <- mov.av(x, order=11)

}
