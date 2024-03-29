% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datPivot.R
\name{datPivot}
\alias{datPivot}
\title{Data - Generates a pivot table.}
\usage{
datPivot(
  x,
  pvar,
  xvar,
  yvar,
  pfun = sum,
  xfilter = NULL,
  NAto0 = TRUE,
  dropNAxvar = TRUE,
  dropNAyvar = TRUE,
  pvar.round = 6,
  returnDT = FALSE,
  savedata = FALSE,
  savedata_opts = NULL,
  gui = FALSE
)
}
\arguments{
\item{x}{Dataframe. Table with pivot variables.}

\item{pvar}{String. The name of the variable for pivot table values.}

\item{xvar}{String. The name of the variable for rows.}

\item{yvar}{String. The name of the variable for columns.}

\item{pfun}{Function. The name of the function to use for pivot table values
(ex. sum, mean, max).}

\item{xfilter}{String. A filter to subset the datatable table x before
pivoting (ex. "STATUSCD == 1").}

\item{NAto0}{Logical. If TRUE, converts NA values to 0.}

\item{dropNAxvar}{Logical. If TRUE, removes columns that are NA.}

\item{dropNAyvar}{Logical. If TRUE, removes rows that have NA values.}

\item{pvar.round}{Integer. Number to round pvar values to.}

\item{returnDT}{Logical. If TRUE, returns a datatable.}

\item{savedata}{Logical. If TRUE, writes output data to outfolder.}

\item{savedata_opts}{List. See help(savedata_options()) for a list
of options. Only used when savedata = TRUE.  If out_layer = NULL, 
default = 'datpivot'.}

\item{gui}{Logical. If TRUE, pop-up windows will appear for user-interface.}
}
\value{
\item{ptab}{ Matrix. The pivot table. }
}
\description{
Generates a pivot table of values by x row and y column.
}
\examples{
# Pivot WYcond table
datPivot(x = FIESTA::WYcond,
         pvar = "CONDPROP_UNADJ",
         xvar = "FORTYPCD",
         yvar = "STDSZCD")

# Pivot WYtree table
datPivot(x = FIESTA::WYtree,
         pvar = "TPA_UNADJ",
         xvar = "SPCD",
         yvar = "STATUSCD",
         pfun = mean,
         NAto0 = TRUE)
}
\author{
Tracey S. Frescino
}
\keyword{data}
