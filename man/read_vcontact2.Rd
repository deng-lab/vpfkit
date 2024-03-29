% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_readvpf.R
\name{read_vcontact2}
\alias{read_vcontact2}
\title{Import vConTACT2 results}
\usage{
read_vcontact2(fin, assembler_label = "_NODE_", version = "v1")
}
\arguments{
\item{fin}{File `genome_by_genome_overview.csv` created by vConTACT2}

\item{assembler_label}{"_NODE_" for SPAdes}

\item{version}{version of the tool used to annotate taxonomy}
}
\value{
list of results
}
\description{
Import vConTACT2 results
}
