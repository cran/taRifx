\name{sides}
\alias{sides}
\alias{sides.default}
\alias{sides.formula}
\title{Figure out how many "sides" a formula has
See also SimonO101's answer at http://stackoverflow.com/a/16376939/636656}
\usage{
  sides(x, ...)

  \method{sides}{default} (x, ...)

  \method{sides}{formula} (x, ...)
}
\arguments{
  \item{x}{The object to calculate the sidedness of}

  \item{\dots}{Other items to pass along}
}
\value{
  An integer of the number of sides
}
\description{
  Figure out how many "sides" a formula has See also
  SimonO101's answer at
  http://stackoverflow.com/a/16376939/636656
}
\examples{
test <- list( ~ a + b, a ~ b + c, b + c ~ a, ~ a ~ b, a ~ b ~ c, a~b+c|d~c~d~e~f~g )
sapply(test,sides)
}

