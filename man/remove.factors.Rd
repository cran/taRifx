\name{remove.factors}
\alias{remove.factors}
\title{
Converts all factors in a data.frame to character.
}
\description{
Converts all factors in a data.frame to character.
}
\usage{
remove.factors(df)
}
\arguments{
  \item{df}{
  	  A data.frame
  }
}
\value{
data.frame
}
\examples{
my.test.df <- data.frame(grp=rep(c("A","B"),10),data=runif(20))
remove.factors(my.test.df)
}
