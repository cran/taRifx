\name{distinct}
\alias{distinct}
\title{
Returns number of distinct observations in each column of a data frame or in a vector
}
\description{
Returns number of distinct observations in each column of a data frame or in a vector
}
\usage{
distinct(input, na.rm = TRUE)
}
\arguments{
  \item{input}{
data.frame or vector
}
  \item{na.rm}{
remove nas or not
}
}
\value{
Num of distinct obs
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (input, NA.rm = TRUE) 
{
    if (NA.rm != TRUE) {
        exclude = c()
    }
    else {
        exclude = c(NA, NaN)
    }
    return(switch(class(input), data.frame = unlist(lapply(input, 
        function(x) length(table(x, exclude = exclude)))), numeric = length(table(input, 
        exclude = exclude)), integer = length(table(input, exclude = exclude)), 
        character = length(table(input, exclude = exclude))))
  }
}
