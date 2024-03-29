% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusters.R
\docType{methods}
\name{clusters}
\alias{clusters}
\alias{clusters,Rscexv-method}
\title{description of function clusters}
\usage{
clusters(
  dataObj,
  clusterby = "raw",
  useGrouping = NULL,
  mds.proj = NULL,
  groups.n = 3,
  ctype = "hierarchical clust",
  onwhat = "Expression",
  cmethod = "ward.D"
)
}
\arguments{
\item{dataObj}{TEXT MISSING}

\item{clusterby}{TEXT MISSING default="raw"}

\item{useGrouping}{do not calculate a new grouping - use this column in the samples table (default=NULL)}

\item{groups.n}{TEXT MISSING default= 3}

\item{ctype}{TEXT MISSING default='hierarchical clust'}

\item{onwhat}{TEXT MISSING default="Expression"}

\item{cmethod}{TEXT MISSING default='ward.D'}
}
\description{
create the grouping based on the MDS or raw data.
}
