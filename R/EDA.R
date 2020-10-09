## stacked barplots in plotly. 
## these are absolutly icky, but since plotly has
## moved to NSE I see no other way!
stackedPlotlyBarPlot <- function(M,names,x.centres,fit.length,axisX,axisY1,axisY2,title,layout=TRUE){
	p <- plotly::plot_ly(M)
	p <- p %>% plotly::add_lines(x=~x.centres,y=~fit.length,name="num.obs",yaxis="y2")
	for(nm in names){eval(parse(text=paste0(
		"p <- p %>% plotly::add_trace(x=~x.centres,y=~",
		nm,
		",type=\"bar\",name=\"",
		nm,
		"\")"
	)))}
	if(layout) p %>% plotly::layout(
		barmode="stack",
		xaxis=axisX,
		yaxis=axisY1,
		yaxis2=axisY2,
		title=title
	) else p
}

setClassUnion("EDAContinuous",c("numeric","POSIXt","Date")) ## difftime isn't a class??
setClassUnion("EDADiscrete",c("character","factor","logical"))

## entropy functions

entropy <- function(s){ ## s should come from summary.factor
	s <- s[s>0]
	if(!length(s)) return(0)
	n <- sum(s)
	sum(s/n*log(s/n))
}

semiEntropy <- function(x,y){ ## x is numeric, y is factor
	sy <- summary.factor(y)
	sx <- sd(x,na.rm=TRUE)
	sb <- by(x,y,sd,na.rm=TRUE)
	sy <- sy[!is.na(sb)]
	sb <- sb[!is.na(sb)]
	sy <- sy / sum(sy)
	1-sum(sb*sy)/sx
} ## close to 1 means x|y is near-deterministic

crossEntropyRatio <- function(x,y){
	sx <- summary.factor(x,maxsum=Inf)
	sy <- summary.factor(y,maxsum=Inf)
	sb <- by(x,y,summary.factor,maxsum=Inf)
	ent.x <- entropy(sx)
	ent.y <- entropy(sy)
	sy <- sy / sum(sy)
	ent.xy <- sapply(sb,entropy)
	return(sum(sy*ent.xy)/ent.y)
} ## close to 1 means x|y is near-deterministic

## functions to produce EDA graphs for individual variables

EDAgraphCC <- function(v,vname,response,responseName,usePlotly=TRUE,...){
## cts input, cts response
    graph.title <- paste0(responseName, " vs. ", vname)
    discrete <- FALSE
    if(length(unique(v)) <= 20){
        return(EDAgraphDC(addNA(v,ifany=TRUE),vname,response,responseName,usePlotly))
    }
    x.cuts <- unique(quantile(v, seq(0.025, 0.975, length.out = 19), na.rm = TRUE))
    if (length(x.cuts) < 2L && any(v != x.cuts[1],na.rm=TRUE)) 
        x.cuts <- unique(c(x.cuts, quantile(v[v != x.cuts[1]], 
            seq(0.025, 0.975, length.out = 19), na.rm = TRUE)))
    if(max(x.cuts) == min(x.cuts) && length(x.cuts == 1L)) x.cuts <- c(x.cuts[1]-0.5,x.cuts[1]+0.5)
    if(min(v,na.rm=TRUE) == x.cuts[1]) x.cuts <- c(min(v,na.rm=TRUE)-(range(x.cuts)[2]-range(x.cuts)[1])/20,x.cuts)
    x.centres <- sapply(1:(length(x.cuts) - 1), function(i) (x.cuts[i + 
        1] + x.cuts[i])/2)
    v[!is.na(v) && v==min(x.cuts)] <- v[!is.na(v) && v==min(x.cuts)]*(1+.Machine$double.eps)
    v[!is.na(v) && v==max(x.cuts)] <- v[!is.na(v) && v==max(x.cuts)]*(1-.Machine$double.eps)
    x2 <- cut(as.numeric(v), unique(x.cuts),right=TRUE)
    fit.means <- by(response, x2, mean, na.rm = TRUE)
    fit.upper <- by(response, x2, function(x) {
        mean(x, na.rm = TRUE) + 2 * sd(x, na.rm = TRUE)/sqrt(length(x))
    })
    fit.lower <- by(response, x2, function(x) {
        mean(x, na.rm = TRUE) - 2 * sd(x, na.rm = TRUE)/sqrt(length(x))
    })
    fit.length <- by(response,x2,length)
    if(anyNA(v)){
        x.centres <- c(x.centres,max(x.centres) + (max(x.centres) - min(x.centres))*0.075)
        fit.lower <- c(fit.lower,mean(response[is.na(v)]) - 2*sd(response[is.na(v)]))
        fit.upper <- c(fit.upper,mean(response[is.na(v)]) + 2*sd(response[is.na(v)]))
        fit.means <- c(fit.means,mean(response[is.na(v)]))
    }
    if(usePlotly){
	axisY1 <- list(overlaying = "y2", side = "left", title=responseName)
  	axisY2 <- list(side = "right", showgrid=FALSE, title="num observations",rangemode="tozero",fixedrange=TRUE,range=c(0,2*max(fit.length)))
	axisX <- list(title=vname, showline=TRUE)
	plotly::plot_ly(x=~x.centres,y=~fit.means,error_y=list(array=(fit.upper-fit.lower)/2),name=~responseName,yaxis="y1") %>%
	plotly::add_lines(x=~x.centres,y=~fit.means,error_y=list(array=(fit.upper-fit.lower)/2),name=~responseName,yaxis="y1") %>%
	plotly::add_bars(x=~x.centres,y=~fit.length,name="Num. Obs",error_y=NULL,yaxis="y2") %>%
	plotly::layout(title=graph.title, xaxis=axisX, yaxis=axisY1, yaxis2=axisY2)
    } else {
	plot(x.centres, fit.means, type = "o", pch = 1, 
     		xaxp = c(1, nlevels(x2), 1), xaxt = "n", xlab = vname, 
     		ylab = responseName, main = graph.title, col=2-is.finite(fit.upper),...)
	axisLabels <- signif(x.centres, 3)
	if(anyNA(v)){
     		axisLabels <- c(axisLabels[-length(axisLabels)],"NA")
	}
	axis(1, labels = axisLabels, at = x.centres)
	segments(x.centres, y0 = fit.upper, y1 = fit.lower)
    }
}

EDAgraphDC <- function(v,vname,response,responseName,usePlotly=TRUE,...){
## dsc input, cts response
    if (!inherits(v, "factor")) 
        v <- as.factor(v)
    graph.title <- paste0(responseName, " vs. ", vname)
    fit.means <- by(response, v, mean, na.rm = TRUE)
    fit.upper <- by(response, v, function(x) {
        mean(x, na.rm = TRUE) + 2 * sd(x, na.rm = TRUE)/sqrt(length(x))
    })
    fit.lower <- by(response, v, function(x) {
        mean(x, na.rm = TRUE) - 2 * sd(x, na.rm = TRUE)/sqrt(length(x))
    })
    fit.length <- by(response,v,length)
    if(usePlotly){
	x.centres <- seq_along(fit.means)
	axisY1 <- list(overlaying = "y2", side = "left", title=responseName)
  	axisY2 <- list(side = "right", showgrid=FALSE, title="num observations",rangemode="tozero",fixedrange=TRUE,range=c(0,2*max(fit.length)))
	axisX <- list(title=vname, showline=TRUE,type="category",ticktext=levels(v),tickmode="array",tickvals=x.centres)
	plotly::plot_ly(x=~x.centres,y=~fit.means,error_y=list(array=(fit.upper-fit.lower)/2),name=~responseName,yaxis="y1") %>%
	plotly::add_lines(x=~x.centres,y=~fit.means,error_y=list(array=(fit.upper-fit.lower)/2),name=~responseName,yaxis="y1") %>%
	plotly::add_bars(x=~x.centres,y=~fit.length,name="Num. Obs",error_y=NULL,yaxis="y2") %>%
	plotly::layout(title=graph.title, xaxis=axisX, yaxis=axisY1, yaxis2=axisY2)
    } else {
      plot(fit.means, type = "o", pch = 1, col = 1, xaxp = c(1, 
          nlevels(v), 1), xaxt = "n", xlab = vname, ylab = responseName, 
          main = graph.title, ...)
      axis(1, labels = levels(v), at = 1:nlevels(v), xaxp = c(1, 
          nlevels(v), 1))
      segments(1:length(fit.lower), y0 = fit.upper, y1 = fit.lower, col=2-is.finite(fit.upper))
    }
}

EDAgraphCD <- function(v,vname,response,responseName,usePlotly=TRUE,...){
## cts input, dsc response
    response <- as.factor(response)
    graph.title <- paste0(responseName, " vs. ", vname)
    if(length(unique(v)) <= 20){
        return(EDAgraphDD(addNA(v,ifany=TRUE),vname,response,responseName))
    }
    x.cuts <- unique(quantile(v, seq(0.025, 0.975, length.out = 19), 
        na.rm = TRUE))
    if (length(x.cuts) < 2L && any(v != x.cuts[1])) 
        x.cuts <- unique(c(x.cuts, quantile(v[v != x.cuts[1]], 
            seq(0.025, 0.975, length.out = 19), na.rm = TRUE)))
    if(max(x.cuts) == min(x.cuts) && length(x.cuts == 1L)) x.cuts <- c(x.cuts[1]-0.5,x.cuts[1]+0.5)
    if(min(v,na.rm=TRUE) == x.cuts[1]) x.cuts <- c(min(v,na.rm=TRUE)-(range(x.cuts)[2]-range(x.cuts)[1])/20,x.cuts)
    x.centres <- sapply(1:(length(x.cuts) - 1), function(i) (x.cuts[i + 
        1] + x.cuts[i])/2)
    v[!is.na(v) && v==min(x.cuts)] <- v[!is.na(v) && v==min(x.cuts)]*(1+.Machine$double.eps)
    v[!is.na(v) && v==max(x.cuts)] <- v[!is.na(v) && v==max(x.cuts)]*(1-.Machine$double.eps)
    x2 <- cut(as.numeric(v), unique(x.cuts),right=TRUE)
    fit.length <- by(response,x2,length)
    fs <- by(response,x2,summary,maxsum=Inf)
    summary0 <- rep(0,nlevels(response))
    names(summary0) <- levels(response)
    for(i in which(sapply(fs,is.null))) fs[[i]] <- summary0
    fs2 <- rbind(0,sapply(fs,function(x) if(is.null(x)) rep(0,nlevels(response)) else cumsum(x)/sum(x)))
    fs3 <- sapply(fs,function(x) x / sum(x))
    xind <- apply(fs2,2,function(v) all(v==0))
    fs2 <- fs2[,!xind,drop=FALSE]
    fs3 <- data.frame(n=c(unlist(fs3)),xc=x.centres[rep(seq_along(x.centres),each=nlevels(response))],yt=as.factor(levels(response)[rep(seq_len(nlevels(response)),length(x.centres))]),xn=as.factor(rep(seq_along(x.centres),each=nlevels(response))))
    class(fs3) <- c("grouped_df","tbl_df","tbl","data.frame") ## dplyr not necessary
	attr(fs3,"vars") <- list(as.name(substitute(xc))) 
    x.centres <- x.centres[!xind]
    col <- rainbow(nlevels(response))
    if(anyNA(v)){
        x.centres <- c(x.centres,max(x.centres) + (max(x.centres) - min(x.centres))*0.075)
        s1 <- summary(response[is.na(v)],maxsum=Inf)
        fs2 <- cbind(fs2,cumsum(c(0,s1))/sum(s1))
    }
    xl <- signif(x.centres, 3)
    if(anyNA(v)){
         xl <- c(xl[-length(xl)],"NA")
    }
    if(usePlotly){
	n2 <- fs3$n
	dim(n2) <- c(length(n2)/length(fit.length),length(fit.length))
	rownames(n2) <- if(is.factor(response)) levels(response) else unique(response)
	n2 <- as.data.frame(t(n2))
## remove missing entries
	x.centres <- x.centres[t1 <- !is.na(x.centres)]
	fit.length <- fit.length[t1]
	n2 <- n2[t1,]
	axisY1 <- list(side = "left", title=responseName,range=c(0,1),fixedrange=TRUE)
  	axisY2 <- list(overlaying = "y", 
		side = "right", showgrid=FALSE, title="num observations",rangemode="tozero",fixedrange=TRUE,range=c(0,2*max(fit.length,na.rm=TRUE)))
	axisX <- list(title=vname, showline=TRUE)
#	plotly::plot_ly(fs3,x=~xc,y=~n,name=responseName,yaxis="y1",type="bar",color=~yt) %>%
#	plotly::add_bars(fs3,x=~xc,y=~n,name=responseName,yaxis="y1",type="bar",color=~yt) %>%
#	plotly::add_lines(x=~x.centres,y=~fit.length,name="Num. Obs",error_y=NULL,yaxis="y2") %>%
	stackedPlotlyBarPlot(n2,names=colnames(n2),x.centres=x.centres,fit.length=fit.length,axisX,axisY1,axisY2,title=graph.title,layout=TRUE)## %>% ## yaxis="y1"
#	plotly::add_lines(x=~x.centres,y=~fit.length,name="Num. Obs",error_y=NULL,yaxis="y2") %>%
#	plotly::layout(barmode="stack",xaxis=axisX, yaxis=axisY1, yaxis2=axisY2,title=graph.title)
#	p <- plotly::plot_ly(fs3,x=~xc,y=~n,name=responseName,yaxis="y1",type="bar",color=~yt)
#	layout(title=graph.title, xaxis=axisX, yaxis=axisY1, yaxis2=axisY2,barmode="stack")
    } else {
      plot(NA,xaxt="n",xlim=range(x.centres)+c(0,(range(x.centres)[2]-range(x.centres)[1])*0.2),ylim=c(0,1),xlab=vname,ylab=responseName,main=graph.title)
      axis(1, labels = xl, at = x.centres)
      for(i in seq_len(nrow(fs2)-1)){
          polygon(
              c(x.centres,rev(x.centres)),
              c(fs2[i,,drop=FALSE],rev(fs2[i+1,,drop=FALSE])),
              col=col[i]
          )
      }
      legend(
          "right",
          levels(response),
          fill=col
      )
    }
}

EDAgraphDD <- function(v,vname,response,responseName,usePlotly=TRUE,...){
## dsc input, dsc response
    response <- as.factor(response)
    v <- as.factor(v)
    graph.title <- paste0(responseName, " vs. ", vname)
    fs <- by(response,v,summary,maxsum=Inf)
    fs2 <- sapply(fs,function(x) if(is.null(x)) rep(0,nlevels(response)) else x/sum(x))
    x.centres <- seq_len(nlevels(v))
    fit.length <- by(response,v,length)
    fs3 <- sapply(fs,function(x) x / sum(x))
    fs3 <- data.frame(n=c(fs3),xc=x.centres[rep(seq_along(x.centres),each=nlevels(response))],yt=as.factor(levels(response)[rep(seq_len(nlevels(response)),length(x.centres))]),xn=as.factor(rep(seq_along(x.centres),each=nlevels(response))))
    class(fs3) <- c("grouped_df","tbl_df","tbl","data.frame") ## dplyr not necessary
	attr(fs3,"vars") <- list(as.name(substitute(xc))) 
    xind <- apply(fs2,2,function(v) all(v==0))
    fs2 <- fs2[,!xind,drop=FALSE]
    fs3 <- fs3[fs3$xn %in% which(!xind),,drop=FALSE]
    x.centres <- x.centres[!xind]
    if(usePlotly){
	n2 <- fs3$n
	dim(n2) <- c(length(n2)/length(fit.length),length(fit.length))
	rownames(n2) <- if(is.factor(response)) levels(response) else unique(response)
	n2 <- as.data.frame(t(n2))
	axisY1 <- list(side = "left", title=responseName,range=c(0,1),fixedrange=TRUE)
  	axisY2 <- list(overlaying = "y", side = "right", showgrid=FALSE, title="num observations",rangemode="tozero",fixedrange=TRUE,range=c(0,2*max(fit.length)))
	axisX <- list(title=vname, showline=TRUE,type="category",ticktext=levels(v),tickmode="array",tickvals=x.centres)
	stackedPlotlyBarPlot(n2,names=colnames(n2),x.centres=x.centres,fit.length=fit.length,axisX,axisY1,axisY2,title=graph.title,layout=TRUE)
#	plotly::plot_ly(fs3,x=~xc,y=~n,name=responseName,yaxis="y1",type="bar",barmode="stack",color=~yt) %>%
#	plotly::add_bars(fs3,x=~xc,y=~n,name=responseName,yaxis="y1",type="bar",barmode="stack",color=~yt) %>%
#	plotly::add_lines(x=~x.centres,y=~fit.length,name="Num. Obs",error_y=NULL,yaxis="y2") %>%
#	plotly::layout(title=graph.title, xaxis=axisX, yaxis=axisY1, yaxis2=axisY2,barmode="stack")
    } else {
	    col <- rainbow(nlevels(response))
	    barplot(
	    fs2,
     	    col=col,
     	    xlab=vname,
     	    ylab=responseName,
          main=graph.title,
          xlim=c(0,25)
      )
      legend(
          "right",
          levels(response),
          fill=col
      )
    }
}

EDAgraphNoGraph <- function(v,response,...){
	print("This data was removed from the analysis, because it was highly correlated with another column that explained the response better. See the log for details")
}

EDAgraph <- function(v,vname,response,responseName,usePlotly=TRUE,...){}
setGeneric("EDAgraph",signature=c("v","response"))
setMethod("EDAgraph",signature(v="EDAContinuous",response="EDAContinuous"),EDAgraphCC)
setMethod("EDAgraph",signature(v="EDADiscrete",response="EDAContinuous"),EDAgraphDC)
setMethod("EDAgraph",signature(v="EDAContinuous",response="EDADiscrete"),EDAgraphCD)
setMethod("EDAgraph",signature(v="EDADiscrete",response="EDADiscrete"),EDAgraphDD)
setMethod("EDAgraph",signature(v="NULL",response="ANY"),EDAgraphNoGraph)

## main function

EDA <- function(
	data,
	response,
	responseClass=c("numeric","factor"),
	correlations=ncol(data) <= 250L,
	displayCorrelations=5L,
	ignoreColumns=character(0),
	ignoreColumnSpec=c("name","partial","regex"),
	returnData=FALSE,
	returnFunction=TRUE,
	returnRmd=TRUE,
	returnCorrelations=FALSE,
	outfile=paste0("EDA",strtrim(make.names(substitute(data)),29),".Rmd"),
	RmdHeader,
	splitRmd=100L,
	response1=FALSE,
#	loss=c("MSE","MAE","AUC","log","mlog"),
	NAIsALevel=FALSE,
	fixDollars=TRUE,
	fixCommas=TRUE,
	fixRanges=TRUE,
	fixCorrelated=correlations,
	fixUnorderedNoise=TRUE,
#	findNA=TRUE,
	excludeDollars=character(0),
	excludeCommas=character(0),
	excludeRanges=character(0),
	excludeCorrelated=character(0),
#	excludeFindNA=character(0),
	currencySymbols=c("$"),
	maxFactorLevels=Inf,
	renderQuietly=TRUE,
	usePlotly=TRUE,
	showHistogram=!usePlotly
	){
	## Rmd header
	if(returnCorrelations && !correlations){
		warning("'returnCorrelations' is TRUE but 'correlations' is FALSE. Can't return something I haven't calculated!\nSetting 'correlations' to TRUE...\n")
		correlations <- TRUE
	}
	if(fixCorrelated && !correlations){
		warning("'fixCorrelated' is nonzero but 'correlations' is FALSE. Can't fix something I haven't calculated!\nSetting 'correlations' to TRUE...\n")
	}
	if(returnRmd){
		rmdPart <- 1L
		dataName <- strtrim(make.names(substitute(data)),29)
		requiresSplit <- ncol(data) > splitRmd
#		require(knitr)
#		require(rmarkdown)
		## open outfile for writing
		if(requiresSplit){
			outfileRoot <- outfile
			outfile <- gsub("[Rr]md$",paste0(rmdPart,".Rmd"),outfile)
		}
		filesToRender <- outfile
		f <- file(outfile,"wt")
		on.exit(close(f))
		## define a function to write to rmd
		## this makes it clearer what is being written to where
		## cat writes to the console, catf to the Rmd file
		catf <- function(...)cat(...,file=f)
		## Rmd header
		if(missing(RmdHeader)){
			if(requiresSplit){
				catf("---\ntitle: \"EDA of",dataName,"part",rmdPart,"\"\nauthor: \"",as.list(Sys.info())$nodename,"\"\noutput:\n  html_document:\n    theme: spacelab\n---\n")
			} else {
				catf("---\ntitle: \"EDA of",dataName,"\"\nauthor: \"",as.list(Sys.info())$nodename,"\"\noutput:\n  html_document:\n    theme: spacelab\n---\n")
			}
		} else {
			catf(RmdHeader)
		}
	}
	if(fixDollars){
		dollarRegex <- paste0(c("^[",currencySymbols,",0-9]+$"),collapse="")
		dollarRegexRepl <- paste0("[",currencySymbols,",]",collapse="")
	}
	if(fixRanges){ ## use tolower
		rangeRegex1 <- "^[[(](-?[0-9.e]+),(-?[0-9.e]+)[])]$"
		rangeRegex2 <- "^(-?[0-9.e]+|up|-inf) *(to|-+|or below|or above) *(-?[0-9.e]+|inf)?$"
	}
	testRowsNum <- min(10L,nrow(data))
	if(returnFunction){
		returnF <- function(x){} ## put things into it with substitute()
		lf <- 2L
	}
	## warn if data isn't a data frame
	## things which inherit from data.frame
	## (like data.table) are ok
	## be sure to use portable syntax, such as data[[name]] rather than data[,name]
	if(!is.data.frame(data)){
		warning("'data' isn't a data.frame, Coercing to same.")
	}
	data <- as.data.frame(data)
	## first, calculate the response
	if(missing(response)){
		cat("'response' is missing. Attempting to work it out from the data...\n")
		if(any(t1 <- toupper(colnames(data)) == "RESPONSE")){
			if(sum(t1) > 1) stop("multiple columns called 'response'! Don't know which to choose...")
			cat("using column",colnames(data)[t1],"\n")
			response <- data[[colnames(data)[t1]]]
			responseName <- colnames(data)[t1]
		} else {
			cat("no column called \"response\", using the last column, ",colnames(data)[length(colnames(data))],"\n")
			response <- data[[colnames(data)[length(colnames(data))]]]
			responseName <- colnames(data)[length(colnames(data))]
		}
	}
	if(length(response) == 1L && !response1){
		cat("'response' has length one. Assuming it is the name or index of a column (specify response1=TRUE to override this behaviour)\n")
		if(is.character(response)) responseName <- response
		if(is.numeric(response)) responseName <- colnames(data)[response]
		response <- data[[response]]
	}
	## now, analyse the response
	cat("'response' has class",class(response),"\n")
	if(missing(responseClass)){
	## coerce it to either numeric or factor
		switch(class(response)[1],
			"numeric"=,
			"difftime"=,
			"POSIXct"=,
			"POSIXlt"=,
			"POSIXt"={
				cat("coercing 'response' to numeric\n")
				response <- as.numeric(response)
				if(returnFunction){
					body(returnF)[[lf]] <- substitute(x[[response]] <- as.numeric(x[[response]]),list(response=responseName[1]))
					lf <- lf + 1L
				}
			},
			{
				cat("coercing 'response' to a factor\n")
				response <- as.factor(response)
				if(returnFunction){
					body(returnF)[[lf]] <- substitute(x[[response]] <- as.factor(x[[response]]),list(response=responseName[1]))
					lf <- lf + 1L
				}
			}
		)
	} else {
		responseClass <- match.arg(responseClass)
		switch(responseClass,
			"numeric"={
				cat("coercing 'response' to numeric\n")
				response <- as.numeric(response)
				if(returnFunction){
					body(returnF)[[lf]] <- substitute(x[[response]] <- as.numeric(x[[response]]),list(response=responseName[1]))
					lf <- lf + 1L
				}
			},
			"factor"={
				cat("coercing 'response' to a factor\n")
				response <- as.factor(response)
				if(returnFunction){
					body(returnF) <- substitute(x[[response]] <- as.factor(x[[response]]),list(response=responseName[1]))
					lf <- lf + 1L
				}
			}
		)
	}
	if(is.na(NAIsALevel)) stop("'NAIsALevel' can't be NA. Were you trying to be funny?")
	if(NAIsALevel && is.factor(response)) response <- addNA(response)
	## see if response has any NA in it
	if(anyNA(response) && !NAIsALevel) warning("'response' has some NA in it. These rows will probably be ignored when modelling. If 'response' is a factor, then you can specify NAIsALevel=TRUE to have a separate level for NA.")
	## now, we can actually look at the response
	if(returnRmd){
		catf("## Response summary\n")
		catf("Response name: ",responseName,"\n\n")
		catf("Response class: ",class(response),"\n\n")
		if(is.numeric(response)){
			catf("```{r response-summary}\n summary(response)\n```\n")
			catf("```{r response-hist}\n hist(response)\n```\n")
		} else if(is.factor(response)){
			catf("Num unique: ",nlevels(response),"\n\n")
			catf("```{r response-summary}\n summary(response,maxsum=100)\n```\n")
		}
	}
	## analyse each explanatory column
	modelVars <- character(0)
	ignoreColumnSpec <- match.arg(ignoreColumnSpec)
	for(v in colnames(data)){
		z <- data[[v]]
		vnna <- sum(is.na(z))
		vl <- length(z)
		vind <- match(v,colnames(data))
		if(returnRmd && requiresSplit && vind %% splitRmd == 0L && vind != ncol(data)){
			cat("starting new file...\n")
			close(f)
			rmdPart <- rmdPart + 1L
			## open outfile for writing
			outfile <- gsub("[Rr]md$",paste0(rmdPart,".Rmd"),outfileRoot)
			filesToRender <- c(filesToRender,outfile)
			f <- file(outfile,"wt")
			on.exit(close(f))
			## define a function to write to rmd
			## this makes it clearer what is being written to where
			## cat writes to the console, catf to the Rmd file
			catf <- function(...)cat(...,file=f)
			## Rmd header
			if(missing(RmdHeader)){
				if(requiresSplit){
					catf("---\ntitle: \"EDA of",dataName,"part",rmdPart,"\"\nauthor: \"",as.list(Sys.info())$nodename,"\"\noutput:\n  html_document:\n    theme: spacelab\n---\n")
				} else {
					catf("---\ntitle: \"EDA of",dataName,"\"\nauthor: \"",as.list(Sys.info())$nodename,"\"\noutput:\n  html_document:\n    theme: spacelab\n---\n")
				}
			} else {
				catf(RmdHeader)
			}
		}
		cat("analysing variable",v," ",vind,"/",ncol(data),"\n")
		if(v == responseName){
			cat("this is the response variable. Skipping...\n")
			next
		}
		if(
			ignoreColumnSpec=="name" && v %in% ignoreColumns || 
			ignoreColumnSpec=="partial" && pmatch(v,ignoreColumns,nomatch=0L) ||
			ignoreColumnSpec=="regex" && any(sapply(ignoreColumns,grepl,v))
		){
			cat("skipping variable",v,"\n")
			next
		}
		## univariate fix routines
		if(inherits(z,"character") || inherits(z,"factor")){
			if(fixDollars && !(v %in% excludeDollars) && all(grepl(dollarRegex,z[1:testRowsNum])) && all(grepl(dollarRegex,z))){
				cat(v,"looks like a numeric specified with currency symbols. Converting to numeric...\n")
				data[[v]] <- as.numeric(gsub(dollarRegexRepl,"",z))
			}
			if(fixCommas && !(v %in% excludeCommas) && all(grepl("^[0-9,]+$",data[[v]][1:testRowsNum])) && all(grepl("^[0-9,]+$",data[[v]]))){
				cat(v,"looks like a numeric with commas in it. Converting to numeric...\n")
				data[[v]] <- as.numeric(gsub(",","",z))
			}
			if(fixRanges && !(v %in% excludeRanges) && all(grepl(rangeRegex1,tolower(z[1:testRowsNum]))) && all(grepl(rangeRegex1,tolower(z)))){
				cat(v,"looks like a range variable with [a,b] format. Converting to numeric...\n")
				vl <- as.numeric(gsub(rangeRegex1,"\\1",z))
				vu <- as.numeric(gsub(rangeRegex1,"\\2",z))
				vl[!is.finite(vl)] <- min(vu,vl,na.rm=TRUE)-1
				vu[!is.finite(vu)] <- max(vu,vl,na.rm=TRUE)+1
				data[[v]] <- 0.5*(vl + vu)
				rm(vl)
				rm(vu)
			}
			if(fixRanges && !(v %in% excludeRanges) && all(grepl(rangeRegex2,tolower(z[1:testRowsNum]))) && all(grepl(rangeRegex2,tolower(z)))){
				cat(v,"looks like a range variable with \"a to b\" format. Converting to numeric...\n")
				vl <- as.numeric(gsub(rangeRegex1,"\\1",z))
				vu <- as.numeric(gsub(rangeRegex1,"\\3",z))
				vl[!is.finite(vl)] <- min(vu,vl,na.rm=TRUE)-1
				vu[!is.finite(vu)] <- max(vu,vl,na.rm=TRUE)+1
				data[[v]] <- 0.5*(vl + vu)
				rm(vl)
				rm(vu)
			}
		}
		if(class(data[[v]]) == "numeric" || class(data[[v]]) == "integer"){

		}
		## end univariate fix routines
		if(returnRmd) catf("\n## ",v,"(column",vind,"of",ncol(data),")\n\n")
		if(vnna == vl){
			cat("Variable is nothing but NA:",v,"\n")
			if(returnRmd) catf("This variable is always NA.\n")
			if(returnFunction){
				body(returnF)[[lf]] <- substitute(x[[v]] <- NULL,list(v=v))
				lf <- lf + 1L
			}
			data[[v]] <- NULL
			next
		}
		if(vnna == 0 && all(data[[v]]==data[[v]][1])){
			cat("Variable is always the same:",v,"=",data[[v]][1],"\n")
			if(returnRmd) catf("This variable is always",data[[v]][1],".\n")
			if(returnFunction){
				body(returnF)[[lf]] <- substitute(x[[v]] <- NULL,list(v=v))
				lf <- lf + 1L
			}
			data[[v]] <- NULL
			next
		}
		if(fixUnorderedNoise && is.factor(data[[v]]) && !is.ordered(data[[v]]) && nlevels(data[[v]]) >= length(data[[v]])){
			cat("Variable is an unordered factor that is never the same.\n")
			if(returnRmd) catf("This variable is an unordered factor but always different.\n")
			if(returnFunction){
				body(returnF)[[lf]] <- substitute(x[[v]] <- NULL,list(v=v))
				lf <- lf + 1L
			}
			data[[v]] <- NULL
			next
		}
		if(returnRmd){
			catf("- class: ",class(data[[v]]),"\n")
			if(vnna) catf("- missing:",vnna,"out of",vl,"\n")
			catf("\n### summary\n")
			if(is.numeric(data[[v]])){
		## todo: something to look for values like -999 which are probably NA
		## todo: analyse columns in parallel or multi-threaded
				catf("```{r ",v,"-summary}\n summary(data[[\"",v,"\"]])\n```\n",sep="")
				if(showHistogram) catf("```{r ",v,"-hist}\n hist(data[[\"",v,"\"]])\n```\n",sep="")
			} else {
				catf("Num unique: ",nlevels(as.factor(data[[v]])),"\n",sep="")
				catf("```{r ",v,"-summary}\n summary.factor(data[[\"",v,"\"]],maxsum=100)\n```\n",sep="")
			}
			if(correlations){
				catf("```{r ",v,"-corr}\nmostCorrelatedNumeric(\"",v,"\",",displayCorrelations,")\nmostCorrelatedFactor(\"",v,"\",",displayCorrelations,")\n```\n",sep="")
			}
			catf("\n### response\n")
			catf("```{r ",v,"-graph,fig.width=10}\nEDAgraph(data[[\"",v,"\"]],\"",v,"\",response,responseName,usePlotly=",usePlotly,")\n```\n",sep="")
		}
	}
	## correlations
	## we do this after the univariate stuff so we are correlating cleaned columns, not uncleaned ones
	if(correlations){
		cat("calculating correlations...\n")
		numericCols <- sapply(data,is.numeric)
		factorCols <- !numericCols ## in case something slipped through that was neither numeric nor factor
		numericCols[match(names(numericCols),excludeCorrelated)] <- FALSE
		factorCols[match(names(factorCols),excludeCorrelated)] <- FALSE
		cat("numeric-numeric\n")
		if(any(numericCols)){
			NNcorMatrix <- abs(cor(as.data.frame(data)[,numericCols],use="pairwise.complete.obs"))
			diag(NNcorMatrix) <- 0 ## so that variable drops to the bottom of the ranking
		}
		## todo: replace this with stuff from the entropy package
		cat("numeric-factor\n")
		if(any(numericCols) && any(factorCols)){
			NFcorMatrix <- outer(which(numericCols),which(factorCols),Vectorize(function(i,j){
				semiEntropy(data[,i],data[,j])
			}))
		}
		cat("factor-factor\n")
		if(any(factorCols)){
			FFcorMatrix <- outer(which(factorCols),which(factorCols),Vectorize(function(i,j){
				if(i==j) return(0) else crossEntropyRatio(data[,i],data[,j])
			}))
		}
		numCompareNumeric <- min(displayCorrelations,max(0L,sum(numericCols)-1L))
		numCompareFactor <- min(displayCorrelations,max(0L,sum(factorCols)-1L))
		mostCorrelatedNumeric <- function(z,n){
			n <- min(n,numCompareNumeric)
			if(n <= 0) return(NULL)
			if(z %in% names(numericCols[numericCols])){
				mrow <- match(z,rownames(NNcorMatrix))
				ri <- order(NNcorMatrix[mrow,],decreasing=TRUE)[1:min(n,ncol(NNcorMatrix))]
				return(t(NNcorMatrix[mrow,ri,drop=FALSE]))
			} else {
				mrow <- match(z,colnames(NFcorMatrix))
				ri <- order(NFcorMatrix[,mrow],decreasing=TRUE)[1:min(n,nrow(NNcorMatrix))]
				return(NFcorMatrix[ri,mrow,drop=FALSE])
			}
		}
		mostCorrelatedFactor <- function(z,n){
			n <- min(n,numCompareFactor)
			if(n <= 0) return(NULL)
			if(z %in% names(numericCols[numericCols])){
				mrow <- match(z,rownames(NFcorMatrix))
				ri <- order(NFcorMatrix[mrow,],decreasing=TRUE)[1:min(n,ncol(NFcorMatrix))]
				return(t(NFcorMatrix[mrow,ri,drop=FALSE]))
			} else {
				mrow <- match(z,colnames(FFcorMatrix))
				ri <- order(FFcorMatrix[,mrow],decreasing=TRUE)[1:min(n,nrow(FFcorMatrix))]
				return(FFcorMatrix[ri,mrow,drop=FALSE])
			}
		}
		## fix correlated
		if(fixCorrelated && numCompareNumeric > 0L) for(z in names(numericCols)){
			if(!(z %in% names(data))) next ## it already got deleted in an earlier step
			t1 <- mostCorrelatedNumeric(z,1)
			altVar <- rownames(t1)[1]
			if(altVar == responseName) next ## don't remove something because it is correlated to the response!
			if(z == responseName) next ## also, don't remove the response... it causes problems
			if(z == altVar) next ## a variable shouldn't be removed just because it is correlated with itself
			## decide whether to remove z or altVar; remove the one less correlated to response
			if(is.factor(response)){
				zcor <- NFcorMatrix[match(z,rownames(NFcorMatrix)),match(responseName,colnames(NFcorMatrix))]
				avcor <- NFcorMatrix[match(altVar,rownames(NFcorMatrix)),match(responseName,colnames(NFcorMatrix))]
			} else {
				zcor <- NNcorMatrix[match(z,rownames(NNcorMatrix)),match(responseName,colnames(NNcorMatrix))]
				avcor <- NNcorMatrix[match(altVar,rownames(NNcorMatrix)),match(responseName,colnames(NNcorMatrix))]
			}
			if(t1+1E-6 >= fixCorrelated){
				if(isTRUE(avcor > zcor)){
					outVar <- z
					keepVar <- altVar
				} else {
					outVar <- altVar
					keepVar <- z
				}
				cat("removing",outVar,"because it exceeds correlation criteria against",keepVar,"\n")
				data[[outVar]] <- NULL
				if(returnFunction){
					body(returnF)[[lf]] <- substitute(x[[z]] <- NULL,list(z=outVar))
					lf <- lf + 1L
				}
				if(nrow(NNcorMatrix)){
					NNr <- match(outVar,rownames(NNcorMatrix))
					NNc <- match(outVar,colnames(NNcorMatrix))
					NNcorMatrix <- NNcorMatrix[-NNr,-NNc]
				}
				if(nrow(NFcorMatrix)){
					NFr <- match(outVar,rownames(NFcorMatrix))
					NFcorMatrix <- NFcorMatrix[-NFr,]
				}
			}
		}
		## reset number of comparisons; we might have lost columns
		numCompareNumeric <- min(displayCorrelations,max(0L,sum(numericCols)-1L))
		numCompareFactor <- min(displayCorrelations,max(0L,sum(factorCols)-1L))
	}
	ro <- list()
	if(returnRmd){
		for(i in filesToRender){
			cat("rendering file",i,"\n")
			rmarkdown::render(i,rmarkdown::html_document(),quiet=renderQuietly)
		}
		ro$outfile <- filesToRender
	}
	if(returnFunction){
		body(returnF)[[lf]] <- substitute(return(x))
		ro$clean <- returnF
	}
	if(returnData) ro$data <- data ## has already effectively been through clean()
	if(returnCorrelations){
		ro$corNN <- NNcorMatrix
		ro$corNF <- NFcorMatrix
		ro$corFF <- FFcorMatrix
	}
	return(ro)
}
## return list
## [[1]] Rmd files
## [[2]] modified data
## [[3]] function(data) to do the same to other sets
## [[4]]-[[6]] correlation matrices (numeric-numeric, numeric-factor, factor-factor)
## as a side effect, the R markdown and html
