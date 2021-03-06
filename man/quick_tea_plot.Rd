% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TeaSetClass.R
\name{quick_tea_plot}
\alias{quick_tea_plot}
\title{quick_tea_plot}
\usage{
quick_tea_plot(data, dataLabels = NA, main = "", col = "gradient",
  pch = 17, cex = 1.5, ...)
}
\arguments{
\item{data}{the data with which to create the ternary plot}

\item{dataLabels}{include if you want to plot text instead of symbols.  Defaults to NA}

\item{main}{the title for your plot}

\item{col}{takes normal colors plus "contrast", "gradient", "gradient90", "gradient180" and "gradient270"}

\item{pch}{as in plot and points functions.  Determines symbol of points plotted.  Defautls to 17 (triangle)}

\item{cex}{as in plot, points and text functions, increases size of these.  Defaults to 1.5 here.}

\item{...}{any additional variables that may work wtih TeaSet's tea_plot() method or the points() function}
}
\value{
NONE
}
\description{
creates a tea_plot given ternary data
}
\note{
essentially accesses tea_plot without having to create a TeaSet object oneself
}
