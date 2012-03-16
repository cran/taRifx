\name{title.page.new}
\alias{title.page.new}
\title{
Plot a title page containing the given text.  Good for breaking up sections of plot PDFs.
}
\description{
Plot a title page containing the given text.  Good for breaking up sections of plot PDFs.
}
\usage{
title.page.new(title.text = "")
}
\arguments{
  \item{title.text}{
Text to plot on its own page
}
}
\value{
	Plot
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (title.text = "") 
{
    plot.new()
    text(0.5, 0.7, title.text)
  }
}
