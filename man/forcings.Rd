% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aadoc.R, R/RcppExports.R
\name{forcings}
\alias{forcings}
\alias{RF_TOTAL}
\alias{RF_T_ALBEDO}
\alias{RF_CO2}
\alias{RF_N2O}
\alias{RF_H2O}
\alias{RF_O3}
\alias{RF_BC}
\alias{RF_OC}
\alias{RF_SO2D}
\alias{RF_SO2I}
\alias{RF_SO2}
\alias{RF_VOL}
\alias{RF_CH4}
\title{Identifiers for capabilities in the Hector forcing component}
\usage{
RF_TOTAL()

RF_T_ALBEDO()

RF_CO2()

RF_N2O()

RF_H2O()

RF_O3()

RF_BC()

RF_OC()

RF_SO2D()

RF_SO2I()

RF_SO2()

RF_VOL()

RF_CH4()
}
\description{
These identifiers specify forcing values that can be provided by hector via
the forcing component.  All of the values corresponding to these identifiers
are read-only (\emph{i.e.}, they can only appear in \code{\link{GETDATA}}
messages.)
}
\section{Functions}{
\itemize{
\item \code{RF_TOTAL}: Total radiative forcing

\item \code{RF_T_ALBEDO}: Radiative forcing due to albedo

\item \code{RF_CO2}: Radiative forcing due to CO2

\item \code{RF_N2O}: Radiative forcing due to N2O

\item \code{RF_H2O}: Radiative forcing due to water vapor

\item \code{RF_O3}: Radiative forcing due to ozone

\item \code{RF_BC}: Radiative forcing due to black carbon

\item \code{RF_OC}: Radiative forcing due to organic carbon

\item \code{RF_SO2D}: Direct contribution of SO2 to radiative forcing

\item \code{RF_SO2I}: Indirect contribution of SO2 to radiative forcing

\item \code{RF_SO2}: Total contribution (direct and indirect) of SO2 to radiative forcing

\item \code{RF_VOL}: Radiative forcing due to volcanic activity

\item \code{RF_CH4}: Radiative forcing due to methane
}}

\section{Note}{

Because these identifiers are provided as \code{#define} macros in the hector code,
these identifiers are provided in the R interface as function.  Therefore,
these objects must be called to use them; \emph{e.g.}, \code{GETDATA()}
instead of the more natural looking \code{GETDATA}.
}

\seealso{
\link{haloforcings} for forcings from halocarbons.

Other capability identifiers: \code{\link{carboncycle}},
  \code{\link{concentrations}}, \code{\link{emissions}},
  \code{\link{haloemiss}}, \code{\link{haloforcings}},
  \code{\link{methane}}, \code{\link{ocean}},
  \code{\link{parameters}}, \code{\link{so2}},
  \code{\link{temperature}}
}
\concept{capability identifiers}
