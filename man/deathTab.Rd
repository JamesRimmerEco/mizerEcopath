% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deathTab.R
\name{deathTab}
\alias{deathTab}
\alias{deathTabUI}
\title{Serve tab with death plots}
\usage{
deathTab(input, output, session, params, logs, diet = NULL, ...)

deathTabUI(...)
}
\arguments{
\item{input}{Reactive holding the inputs}

\item{output}{Reactive holding the outputs}

\item{session}{Shiny session}

\item{params}{Reactive value holding updated MizerParams object}

\item{logs}{Environment holding the log of steady states.}

\item{diet}{A diet matrix to match the diet of the model to.}

\item{...}{Unused}
}
\description{
Serve tab with death plots
}
