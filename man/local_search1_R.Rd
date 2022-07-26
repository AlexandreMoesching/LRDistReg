% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_search_functions.R
\name{local_search1_R}
\alias{local_search1_R}
\title{Local search (row), R version}
\usage{
local_search1_R(theta, l, m, n, mM, lL, PP, w, w_ul)
}
\arguments{
\item{theta}{Log-parameter}

\item{l}{Number of unique covariates}

\item{m}{Number of unique responses}

\item{n}{Sample size}

\item{mM}{(m_j,M_j) index pairs}

\item{lL}{(l_k,L_k) index pairs}

\item{PP}{Reduced index space}

\item{w}{Sample weights}

\item{w_ul}{Cumulative row-sums of sample weights}
}
\value{
New proposal Psi and step-size delta
}
\description{
Local search (row), R version
}
