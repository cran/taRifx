\name{daysofweek}
\alias{daysofweek}
\title{
 Return a vector of the days of the week, in order
}
\description{
 Return a vector of the days of the week, in order
}
\usage{
daysofweek(start.day = "Monday")
}
\arguments{
  \item{start.day}{
	Day of the week to begin the week with (as a text item)
}
}
\value{
Character vector of length 7
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (start.day = "Monday") 
{
    wkdays <- c("Monday", "Tuesday", "Wednesday", "Thursday", 
        "Friday", "Saturday", "Sunday")
    wkdays <- rep(wkdays, 2)
    selector <- unlist(lapply(wkdays, function(x) x == start.day))
    selector.numeric <- seq(length(wkdays))[selector]
    selector.numeric[2] <- selector.numeric[2] - 1
    return(wkdays[selector.numeric[1]:selector.numeric[2]])
  }
}

