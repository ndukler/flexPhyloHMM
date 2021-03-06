% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maskNoiseCalls.R
\name{maskNoiseCalls}
\alias{maskNoiseCalls}
\title{Heuristically find noisy state calls}
\usage{
maskNoiseCalls(thmm, noise.calls, noise.threshold = 16)
}
\arguments{
\item{thmm}{the tree hmm the the masking will be performed on}

\item{noise.calls}{A list of bins that are noisy and the corresponding state}

\item{noise.threshold}{The minimum value of sf.sum required to not be masked}
}
\value{
a list of transition matricies, one entry for each specified branch length
}
\description{
Heuristically find noisy state calls
}
