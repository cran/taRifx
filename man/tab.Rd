\name{tab}
\alias{tab}
\title{Table function which lists NA entries by default
This is a simple wrapper to change defaults from the base R table()}
\usage{
  tab(..., exclude = NULL,
    useNA = c("no", "ifany", "always"), deparse.level = 1)
}
\arguments{
  \item{\dots}{one or more objects which can be interpreted
  as factors (including character strings), or a list (or
  data frame) whose components can be so interpreted. (For
  as.table and as.data.frame, arguments passed to specific
  methods.)}

  \item{exclude}{levels to remove for all factors in ....
  If set to NULL, it implies useNA = "always". See
  'Details' for its interpretation for non-factor
  arguments.}

  \item{useNA}{whether to include NA values in the table.
  See 'Details'.}

  \item{deparse.level}{controls how the default dnn is
  constructed. See 'Details'.}
}
\value{
  tab() returns a contingency table, an object of class
  "table", an array of integer values
}
\description{
  Table function which lists NA entries by default This is
  a simple wrapper to change defaults from the base R
  table()
}
\seealso{
  table
}

