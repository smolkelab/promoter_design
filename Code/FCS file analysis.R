# Functions for analyzing FCS files
# load flowframe
# subset flowframe
# extract values as a table
# plotting functions e.g. basicPlot

# Basic workflow:
# Load flowframe with load.fromfilename
# If needed, use test.DAPI to choose cutoff values for DAPI filtering
# use filter.DAPI to enforce a DAPI cutoff
# use extract.data to convert the gated flowframe object to a numerical matrix
# use basicPlot to examine the actual data

# Possible improvements: find an easier way of doing gating, generalize gating to include FSC/SSC.

library(flowCore)
library(flowViz)
library(data.table)
library(ggplot2)
library(gplots)
library(beanplot)

# Load a file; apply an inital filter to drop anything with e.g. negative values so logarithms can be used
load.fromfilename = function(filename, range.log = c(-1,9)) { 
	ff = read.FCS(filename, emptyValue = FALSE, alter.names = TRUE)
	gate = make.gate.nonneg(ff, range.log)
	ff.subset = Subset(ff, gate)
	return(ff.subset)
}

# Given a FlowFrame, test a potential DAPI cutoff ('thresh') by histogramming and scatter-plotting
# with the cutoff
# Added 2015-08-04
test.DAPI = function(ff, thresh = 3.0,newPlot = TRUE, dapi.expr = 6, ssc.expr = 2) {
	if(typeof(ff) == 'character') { ff = read.FCS(ff, emptyValue = FALSE, alter.names = TRUE) }
	gate = make.gate.nonneg(ff, range.log=c(1,9), gateSizes = TRUE)
	ff = Subset(ff, gate); ff = update.param.ff(ff)
	l.trans = logTransform(logbase=10)

	ff.names = get.name.rules(ff)
	ff.trans.str = paste('transform(ff,',ff.names['DAPI'],'=l.trans(',ff.names['DAPI'],'))')
	ff = eval(parse(text=ff.trans.str))
	ff = transform(ff,SSC.A = l.trans(SSC.A))
	ff = update.param.ff(ff)
	if(newPlot) { par(mfrow=c(1,2)) }
	hist(ff@exprs[,dapi.expr],breaks=200); abline(v=thresh,col='red')
	plot(ff@exprs[,c(ssc.expr,dapi.expr)],pch='.'); abline(h=thresh,col='red')
}

# DAPI-screen a FlowFrame; return a subsetted FlowFrame with the lower of two KMeans-assigned DAPI clusters
# Modified 2015-08-03 to use one-size-fits-all threshold (10^cutoff); previously used kmeans-based approach
# which turns out to not always be totally accurate
filter.DAPI = function(ff, range.log = c(-1,9), gate, cutoff) {
	if(missing(gate)) { gate = make.gate.nonneg(ff, range.log) }
	ff.subset = Subset(ff, gate)
	ff.DAPI.name = as.character(get.name.rules(ff)['DAPI'])
	if(missing(cutoff) || cutoff <= 0) {
		ff.DAPI.val = list(c('Low','High'))
		names(ff.DAPI.val) = ff.DAPI.name
		ff.subset = split(ff.subset, kmeansFilter(ff.DAPI.val))$Low
	}
	else {
		DAPIgate.list = list(10^(c(range.log[1],cutoff))); names(DAPIgate.list) = ff.DAPI.name
		DAPIgate = rectangleGate(filterId='Fluorescence Region',DAPIgate.list)
		ff.subset = Subset(ff.subset, DAPIgate)
	}
	ff.subset = update.param.ff(ff.subset)		
	return (ff.subset)
}

# Return a matrix with each row a cell, each channel a column
extract.data = function(ff, vals = NULL) { 
	if(length(vals) == 0) { return(ff@exprs) }
	ff.rules = get.name.rules(ff)

	ans = list()
	for(i in 1:length(vals)) { ans[[i]] = ff@exprs[,ff.rules[ vals[i] ]] }
	ans.mat = matrix( nrow = length(ans), ncol = max(unlist(lapply(ans,length))) )
	for(i in 1:length(ans)) { ans.mat[i,] = ans[[i]] }
	rownames(ans.mat) = vals
	return (t(ans.mat))
}

# Given a single 'rats' table or a list of them, plot 'size' of each
# Issue: points from later tables overlay previous - would be better to plot them randomly. Current hack:
# choose size so that things aren't totally obscured but distribution still clear (doesn't always work out)
# Support for setting colors ('col') added 2018-03-08
basicPlot = function(datTable, lines = FALSE, size = 0, main = '', xlim = c(3,5), ylim = c(2,5),col = NA, log = TRUE,...) {
	if(typeof(datTable) != 'list') { datTable = list(datTable) }
  if(is.na(col[1])) { col = 1:length(datTable) }
	
	if(size > 0) { for(i in 1:length(datTable)) { if(nrow(datTable[[i]]) > size) {
			datTable[[i]] = datTable[[i]][sample(x=1:nrow(datTable[[i]]),size=size,replace = FALSE),]
	}}}
	if(log) { datTable = lapply(datTable,log10) }
	datTable.1 = datTable[[1]]
	plot(datTable.1,pch = '.', main = main, xlim = xlim, ylim = ylim,col=col[1],...)
	if(lines) {for (i in seq(-3,3,by = 0.5)) {abline(a=i,b=1,col = 'blue',lty = 2)}}
	if( length(datTable) > 1 ) { for(i in 2:length(datTable)) { points(datTable[[i]],pch = '.', col = col[i]) }}
}

# Nice beanplots of inducible promoter experiment
clean.beanplot = function(x,out.frac = 0.02, ...) {
  beanplot(sapply(x, drop.outliers, out.frac = out.frac), 
           what = c(FALSE, TRUE, FALSE, FALSE), bw="nrd0", col = cols.bean,
           ylab = 'TEF-normalized strength (log10)',
           ...)
}
get.rats = function(dat.in) {
  dat.in = log10(dat.in)
  dat.in = dat.in[,2] - dat.in[,1]
  return(dat.in)
}
trim.name = function(n, prefix) {
  n = strsplit(n,'.',fixed = TRUE)[[1]][1]
  n = strsplit(n,prefix)[[1]][2]
  return(n)
}
drop.outliers = function(dat.in, out.frac) {
  lb = quantile(dat.in, out.frac)
  ub = quantile(dat.in, 1 - out.frac)
  dat.in = dat.in[intersect(which(dat.in > lb), which(dat.in < ub))]
  return(dat.in)
}




###########################################################################################################
# Helper functions for the above

###########################################################################################################

l.trans = logTransform(logbase = 10)

# apply over a char.list of flowsets (names), reassigns flowset objects.
updateParameters = function(flowsets) {
  for (p in flowsets) {
    tmp = get(p)
    assign(p, update.param.fs(tmp), envir=.GlobalEnv)
  }
}

# apply over a single flowset, returns flowset (unnamed)
update.param.fs = function(flowset) {
 fsApply(flowset, update.param.ff)
}

# apply over a single flowframe, returns flowframe (unnamed)
update.param.ff = function(ff) {
  channel.names = colnames(ff)
  minValues = c()
  for (p in channel.names) { 
    minValues[p] = min(ff@exprs[,p])
  }
  maxValues = c()
  for (p in channel.names) { 
    maxValues[p] = max(ff@exprs[,p])
  }
  params = parameters(ff)
  params$minRange = minValues
  params$maxRange = maxValues
  parameters(ff) = params
  return(ff)
}

# Given a flowframe, find its naming convention for channels
# compare to "targets", given in the order mCherry, GFP, DAPI
get.name.rules = function(ff, targets = rbind(c('FL4.A','FL7.A','FL1.A'),c('Y2.A','B1.A','V1.A'),
					c('mCherry.A','FITC.A','Pacific.Blue.A'),c('mCherry-A','FITC-A','Pacific Blue-A'))) {
	colnames(targets) = c('mCherry','GFP','DAPI')
	names.tmp = colnames(ff@exprs)
	ans = names.tmp
	for (i in 1:length(ans)) {
		names(ans)[i] = ans[i]
		for(j in 1:ncol(targets)) { if(ans[i] %in% targets[,j]) {
			names(ans)[i] = colnames(targets)[j]
		}}
	}
	return(ans) #vector with names = human-readable name, value = machine-readable name of channel
}

# make a gate with names appropriate to a flowframe using get.name.rules
make.gate.nonneg = function(ff, range.log = c(-1,9), gateSizes = FALSE) {
	gate.range = 10^range.log
	ff.names = get.name.rules(ff)
	ff.mC = as.character(ff.names['mCherry'])
	ff.GFP = as.character(ff.names['GFP'])
	ff.DAPI = as.character(ff.names['DAPI'])
	ff.list = list(gate.range,gate.range,gate.range)
	names(ff.list) = c(ff.mC, ff.GFP, ff.DAPI)
	if(gateSizes) { 	
		ff.list = list(gate.range,gate.range,gate.range,gate.range,gate.range)
		names(ff.list) = c(ff.mC, ff.GFP, ff.DAPI,'FSC.A','SSC.A')
	}
	ans = rectangleGate(filterId="Fluorescence Region", ff.list)
	return(ans)
}


