\name{EDAgraph}
\alias{EDAgraph}
\alias{EDAgraph,EDADiscrete,EDADiscrete-method}
\alias{EDAgraph,EDADiscrete,EDAContinuous-method}
\alias{EDAgraph,EDAContinuous,EDADiscrete-method}
\alias{EDAgraph,EDAContinuous,EDAContinuous-method}
\alias{EDAgraph,NULL,ANY-method}

\title{
Exploratory Data Graphs
}
\description{
Plots how a response variable varies with an explanatory variable.
}
\usage{
EDAgraph(v, vname, response, responseName, usePlotly = TRUE, ...)
\S4method{EDAgraph}{EDADiscrete,EDAContinuous}(v, vname, response, responseName, usePlotly = TRUE, ...)
\S4method{EDAgraph}{EDADiscrete,EDADiscrete}(v, vname, response, responseName, usePlotly = TRUE, ...)
\S4method{EDAgraph}{EDAContinuous,EDAContinuous}(v, vname, response, responseName, usePlotly = TRUE, ...)
\S4method{EDAgraph}{EDAContinuous,EDADiscrete}(v, vname, response, responseName, usePlotly = TRUE, ...)
\S4method{EDAgraph}{NULL,ANY}(v, vname, response, responseName, usePlotly = TRUE, ...)
}

\arguments{
  \item{v}{
The explanatory variable. If it's not \code{\link{numeric}}, \code{\link{integer}} or a date/time class, it will be assumed to be a factor.
}
  \item{vname}{
Explanatory variable name (used in annotating the plot). Defaults to \code{substitute(v)}, with some trimming in case that is extremely long.
}
  \item{response}{
Response variable. Same class rules as for \code{v}.
}
  \item{responseName}{
Response name. Same rules as for \code{vname}.
}
  \item{usePlotly}{
Logical. If TRUE (the default), use the plotly library which results in much prettier plots but which render in HTML and potentially take more memory to render. If FALSE, default R plotting is used.
}
  \item{\dots}{
Other arguments and graphical parameters to pass to \code{\link{plot}}. Not used for plotly plots.
}
}
\details{
Missing values in \code{v} will always be treated separately, whether it's continuous or a factor.

For continuous \code{v}, if there are fewer than 20 unique values, each is assigned its own bin. Else, the breakpoints are chosen as the unique odd quattuordeciles of \code{v}. This should ensure where possible that there are similar amounts of data in each bin.

The vertical lines are the mean +/- two SDs of the response in that bin.

If there is only one datum in a bin, that point will be coloured red for a continuous response.

\code{EDAgraph} is a (S4) generic function, dispatching different methods depending on whether \code{v} and \code{response} are continuous or discrete. This is determined with the hidden class unions \code{EarnixContinuous} and \code{EarnixDiscrete}.
}
\value{
NULL. The plot is drawn as a side-effect.
}

\author{
James Lawrence
}
\note{
This function supersedes the function \code{\link{EarnixModelAnalysisPlot}}, which may become defunct in future versions of the \code{rsai} package.
}

\seealso{
\code{link{boxplot}}.
}
\examples{
## continuous response
x <- rnorm(1000)
y <- x^2 + rnorm(1000)
EDAgraph(x,"x",y,"y")
## discrete response
y2 <- as.factor(round(y/(abs(x)+1)))
EDAgraph(x,"x",y2,"y2")
}

\keyword{ EDA }