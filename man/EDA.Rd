\name{EDA}
\alias{EDA}

\title{
Exploratory Data Analysis
}
\description{
Produces summaries of a dataset and, optionally, a function to clean similar data and/or a html report giving one-way plots of each variable and their effect on the response.

This function has many arguments, but only the first is absolutely mandatory. It is designed as far as possible to "guess right", and then give you options to override that guess if it is wrong.
}
\usage{
EDA(
	data,
	response,
	responseClass = c("numeric", "factor"),
	correlations = ncol(data) <= 250L,
	displayCorrelations = 5L,
	ignoreColumns = character(0),
	ignoreColumnSpec = c("name", "partial", "regex"),
	returnData = FALSE,
	returnFunction = TRUE,
	returnRmd = TRUE,
	returnCorrelations = FALSE,
	outfile = paste0("EDA", strtrim(make.names(substitute(data)), 29), ".Rmd"),
	RmdHeader,
	splitRmd = 100L,
	response1 = FALSE,
	NAIsALevel = FALSE,
	fixDollars = TRUE,
	fixCommas = TRUE,
	fixRanges = TRUE,
	fixCorrelated = correlations,
	fixUnorderedNoise = TRUE,
	excludeDollars = character(0),
	excludeCommas = character(0),
	excludeRanges = character(0),
	excludeCorrelated=character(0),
	currencySymbols = c("$"),
	maxFactorLevels = Inf,
	renderQuietly = TRUE,
	usePlotly = TRUE,
	showHistogram = !usePlotly
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A \code{\link{data.frame}} (or something that inherits from that, like a \code{data.table}).
}
  \item{response}{
Which variable is the response variable. Can either be a vector (the response itself) or a scalar (the name or column index of the response).
}
  \item{responseClass}{
The class of the response. Normally auto-detection works fine, but you can over-ride this if desired. Must be either "numeric" or "factor".
}
  \item{correlations}{
Logical. Should correlations between the columns of \code{data} be computed? The default is yes unless there are more than 250 columns (where it would take a long time; the complexity is quadratic in the number of columns)

This argument must be \code{TRUE} if you want to specify \code{returnCorrelations} or \code{fixCorrelated}, and will be set to TRUE (with a warning) if \code{correlations=FALSE}.
}
  \item{displayCorrelations}{
Logical. How many correlations should be displayed in the markdown. A smaller number will be used if fewer columns actually exist in the data. Default 5.
}
  \item{ignoreColumns}{
Which columns you want to exclude from the routine. See argument \code{ignoreColumnSpec} below.
}
  \item{ignoreColumnSpec}{
Either:
	- "name", whence \code{ignoreColumns} is a vector of names of columns to exclude from analysis.
	- "partial", whence any column name partially matched by anything in \code{ignoreColumns} is ignored.
	- "regex", whence any column name matching \code{ignoreColumns} as a regular expression is ignored.
}
  \item{returnData}{
logical. If TRUE, returns a cleaned data set.
}
  \item{returnFunction}{
logical. If TRUE, returns a function that can clean the specified data and data in a similar format.
}
  \item{returnRmd}{
logical. If TRUE, returns a HTML document and Rmarkdown summarising the data. Requires \code{knitr} and \code{rmarkdown}, and probably you need to have pandoc installed as well.
}
  \item{returnCorrelations}{
logical. If TRUE, the output will include the covariance matrices (up to three of them) between the columns of \code{data}.
}
  \item{outfile}{
If \code{returnRmd=TRUE}, the filename to write the markdown to. Automatically generated if not specified.
}
  \item{RmdHeader}{
A header to use in the R markdown. One is generated automatically if this is not set.
}
  \item{splitRmd}{
Numeric. If the data has more than this many columns, markdown output will be split into several files with this many columns in each one. Default 100.
}
  \item{response1}{
If \code{response} is of length one, it is assumed to be a column name or index. If that isn't the case, set \code{response1=TRUE}.
}
  \item{NAIsALevel}{
Logical. Should the response variable have NA as an explicit level?
}
  \item{fixDollars}{
Logical. Should data with format such as "$1000" be converted to numeric?
}
  \item{fixCommas}{
Logical. Should character data which appears to be numeric with thousands commas be converted to numeric?
}
  \item{fixRanges}{
Logical. Should character data which appears to be a range (e.g. "15 to 20") be converted to numeric? If \code{TRUE}, the midpoint of the range is taken, unless it is a semi-infinite range (e.g. "up to 100"), whence a value close to the finite endpoint is chosen.
}
  \item{fixCorrelated}{
Logical or numeric. Any numeric column \code{x1} with a correlation of more than \code{fixCorrelated} with another numeric column \code{x2} will be removed from the data, unless either (a) it is more correlated with \code{response} than \code{x2} was, whence \code{x2} will be removed instead; or (b) either \code{x1} or \code{x2} is the response. There is a small numerical tolerance of 1E-6 implemented to allow for numerical error.

Using 1 or TRUE will remove only numeric columns that are effectively duplicates of each other.

The default is to use 1 if \code{correlations==TRUE}, else zero. A nonzero value implies \code{correlations} must be TRUE.

This may be extended to include factor columns in future.
}
  \item{fixUnorderedNoise}{
Logical. If TRUE, any column which is categorical but never \code{duplicated} is removed, since it would be useless for modelling and may take up quite a lot of memory.
}
  \item{excludeDollars}{
Vector of column names to be ignored by the check in \code{fixDollars} (see above).
}
  \item{excludeCommas}{
Vector of column names to be ignored by the check in \code{fixCommas} (see above).
}
  \item{excludeRanges}{
Vector of column names to be ignored by the check in \code{fixRanges} (see above).
}
  \item{excludeCorrelated}{
Vector of column names to be ignored when calculating correlations (and when checking to see if columns can be removed for correlation).
}
  \item{currencySymbols}{
Character vector of currency symbols used by \code{fixDollars}. Following the UK's decision to leave the EU, the Euro and pound symbols have been removed from this list. Depending on the result of the US elections, the dollar symbol may soon follow.
}
  \item{maxFactorLevels}{
Numeric. Maximum number of permissible levels for a categorical variable. The default is no limit, but that can lead to large outputs and memory usage if there are factors with many levels.
}
  \item{renderQuietly}{
Logical. If TRUE, the output from \code{render} will be suppressed. Bear in mind that rendering will typically take longer than the analysis itself, since that is where all the graphs are created.
}
  \item{usePlotly}{
Logical. If TRUE, plotly plots are used (which look much prettier). If FALSE, default R graphics are used instead. These will be quicker to render, but take up more disk space. Conversely the plotly files tend to be smaller for large datasets, but take longer to render if there are many columns.
}
  \item{showHistogram}{
Logical. If TRUE, show a histogram of each column as well as the reponse plot. Defaults to TRUE when default R plotting is used, but FALSE for plotly plotting (since a histogram is included in those plots anyway)
}
}
\details{
Some of the "guessing" is detailed here.

The response column, if not specified, is first taken to be any column with a name "response" (ignoring case). If there are more than one such column, an error is returned. If there are no such columns, the last column is taken to be the response.

Any column which is always missing or always the same value will be ignored, and stripped out of the data to save space in the \code{$clean} element of the output.

A "range" is currently defined as anything which regex-matches either "^[[(](-?[0-9.e]+),(-?[0-9.e]+)[])]$" or "^(-?[0-9.e]+|up|-inf) *(to|-+|or below|or above) *(-?[0-9.e]+|inf)?$".
}
\section{Note}{
Since plotly have completely changed their API as of October 2016, this requires plotly version 4.5.2 or later. If you do not wish to upgrade plotly, then you will have to use an earlier version of the \code{rsai} package.
}
\value{
A list, which depending on the options \code{return*} will have some or all of the components:
- \code{data}: a cleaned data set
- \code{clean}: a function that would produce \code{.$data} if applied to the input, and is suitable for other datasets of similar format
- \code{outfile}: the rmarkdown file
- \code{corNN}: correlation matrix between numeric columns
- \code{corNF}: conditional entropy ratio matrix between numeric and factor columns
- \code{corFF}: cross-entropy marix between factor columns (note, not symmetric)
}

\author{
James Lawrence
}

\seealso{
\code{\link{summary}}. \code{\link{EDAgraph}}
}
\examples{
## simplest example, categorical response
## the last column is the response, so
## the default guess works fine
EDA(iris)

## an example where we specify the response column
EDA(mtcars,response="mpg")
}

\keyword{ EDA }