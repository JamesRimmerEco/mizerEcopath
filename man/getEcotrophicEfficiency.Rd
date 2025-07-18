% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.R
\name{getEcotrophicEfficiency}
\alias{getEcotrophicEfficiency}
\title{Get Ecotrophic Efficiency for each species}
\usage{
getEcotrophicEfficiency(params)
}
\arguments{
\item{params}{A MizerParams object}
}
\value{
A named vector of ecotrophic efficiency for each species
}
\description{
For each species, returns the ecotrophic efficiency, the proportion of
production that is not lost to mortality external to the model:
\deqn{EE_i = 1-\frac{M0B_i}{P_i}}
where \eqn{M0B_i} is the rate of biomass loss due to external mortality
and \eqn{P_i} is the rate of production.
}
\examples{
getEcotrophicEfficiency(NS_params)
}
\seealso{
Other rate functions: 
\code{\link{getConsumption}()},
\code{\link{getDietMatrix}()},
\code{\link{getGonadicProduction}()},
\code{\link{getM0B}()},
\code{\link{getM2B}()},
\code{\link{getMetabolicRespiration}()},
\code{\link{getOffspringProduction}()},
\code{\link{getProduction}()},
\code{\link{getReproductiveLoss}()},
\code{\link{getRespiration}()},
\code{\link{getSomaticProduction}()},
\code{\link{getTotalProduction}()},
\code{\link{getUnassimilated}()},
\code{\link{getZB}()}
}
\concept{rate functions}
