\name{categorize}
\alias{categorize}
\title{
Categorize a vector based on a data.frame with two columns, the low and high end points of each category
}
\description{
Categorize a vector based on a data.frame with two columns, the low and high end points of each category.  
}
\usage{
categorize(vec, cutpoints.df, match.min = TRUE, names = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
vector to categorize
}
  \item{cutpoints.df}{
quantile_cutpoints will create a data.frame of the proper format here
}
  \item{match.min}{
%%     ~~Describe \code{match.min} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
}
\value{
Categorized values
}
\seealso{
\code{\link{quantile_cutpoints}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vec, cutpoints.df, match.min = TRUE, names = TRUE) 
{
    cat.one <- function(x, cutpoints.df, names) {
        if (length(x) != 1 & class(x) != "numeric") {
            stop("x must be a single number.")
        }
        if (match.min) {
            cutpoints.df[1, 1] <- cutpoints.df[1, 1] - 1e-05
        }
        selector <- (x > cutpoints.df[, 1]) & (x <= cutpoints.df[, 
            2])
        if ("TRUE" \%in\% names(table(selector))) {
            if (table(selector)[["TRUE"]] != 1) {
                warning(x, "matched more than one category.")
                return(NA)
            }
        }
        else {
            warning(x, "matched zero categories")
            return(NA)
        }
        if (names) {
            return(rownames(cutpoints.df)[selector])
        }
        else {
            return(seq(nrow(cutpoints.df))[selector])
        }
    }
    sapply(vec, cat.one, cutpoints.df = cutpoints.df, names = names)
  }
}
