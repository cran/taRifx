\name{xtable.lme}
\alias{xtable.lme}
\alias{xtable.summary.lme}
\title{
Add in methods to handle LME objects in xtable
}
\description{
Add in methods to handle LME objects in xtable
}
\usage{
xtable.lme(x, caption = NULL, label = NULL, align = NULL, digits = NULL, display = NULL, beta.names = NULL, ...)
xtable.summary.lme(x, caption = NULL, label = NULL, align = NULL, digits = NULL, display = NULL, beta.names = NULL, ...)
}
\arguments{
  \item{x}{
	Model object
}
  \item{caption}{
	Caption for table
}
  \item{label}{
	See ?xtable
}
  \item{align}{
See ?xtable
}
  \item{digits}{
See ?xtable
}
  \item{display}{
See ?xtable
  }
  \item{beta.names}{
See ?xtable
  }
  \item{\dots}{
  	  Arguments to pass to xtable
  }
}
\value{
xtable object
}
\seealso{
\code{\link[xtable]{xtable}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
    display = NULL, beta.names = NULL, ...) 
{
    return(xtable.summary.lme(summary(x), caption = caption, 
        label = label, align = align, digits = digits, display = display, 
        beta.names = beta.names))
  }
}
