# General functions



# Takes a vector and returns a vector of equal length containing all trues (used for selecting all of a given vector)
trues = function(vec) {
	rep(TRUE,length(vec))
}

# Returns the fractional part of a vector
fpart = function(vec) {
	ret = vec - as.integer(vec)
	ret
}

# Returns odds/evens from a vector
odds=function(vec) {
	stopifnot(class(vec)=="integer")
	ret = vec[fpart(vec/2)!=0]
	ret
}
evens=function(vec) {
	stopifnot(class(vec)=="integer")
	ret = vec[fpart(vec/2)==0]
	ret
}

# Round a numeric vector down to the nearest roundvec
roundnear <- function(vec,roundvec) {
  .Deprecated("round_any","plyr","Deprecated. Gave inaccurate results and duplicated functionality already available.")
}

# Takes a dataframe and replicates the chosen observations n times
expandDF = function(df,obs,numtimes=1) {
	dfnew=df
	for(i in 1:numtimes) {
		dfnew=rbind(dfnew,df[obs,])
	}
	return(dfnew)
}

# Takes a dataframe and splits it into a bunch of data.frames held in a list, according to one variable
splitDF = function(df,splitvar) {
	lvls = levels(as.factor(df[[splitvar]]))
	splitdfs = list()
	for(lvl in lvls) {
		splitdfs[[length(splitdfs)+1]] = subset(df,get(splitvar)==lvl)
		names(splitdfs)[length(splitdfs)] <- lvl
	}
	return(splitdfs)
}

# Takes a list of data.frames produced by splitDF and returns them as one appended data.frame
unsplitDF = function(splitdfs) {
	returndf <- splitdfs[[1]]
	if(length(splitdfs) > 1) {
		for(i in 2:length(splitdfs)) {
			returndf <- rbind(returndf,splitdfs[[i]])
		}
	}
	return(returndf)
}


# Returns whether a vector is homogenous or not
homogenous <- function(vec) {
	return(sd(as.integer(as.factor(vec)))==0)
}


# return a vector containing the locations (or T/F for) of the middle of a group
middle.group=function(vec,type="tf") {
	# type is either "tf" for a T/F vector or "loc" for a locations numeric vector

	group.lengths <- rle(vec)$lengths
	group.lengths.shift <- c(0,group.lengths[1:(length(group.lengths)-1)])
	locations=cumsum(group.lengths.shift)+floor(group.lengths/2+1)
	if(type=="loc") 	return(locations)
	
	tf=rep(FALSE,length(vec))
	tf[locations] <- TRUE
	if(type=="tf") 	return(tf)
	stop("Improper type specified")
}

# Make a table by group
# Usage:
#   print(latex.table.by(test.df), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)
#   then add \usepackage{multirow} to the preamble of your LaTeX document
#   for longtable support, add ,tabular.environment='longtable' to the print command (plus add in ,floating=FALSE), then \usepackage{longtable} to the LaTeX preamble
latex.table.by = function(df,num.by.vars=1,...) {
	require(xtable)
	# first num.by.vars groups must be sorted and in descending order of priority
	if(!is.numeric(num.by.vars) | length(num.by.vars)!=1) {
		stop("num.by.vars must be a number")
	}
	# Create a by.vars vector
	by.vars=1:num.by.vars
	
	numcols=length(colnames(df))
	df.original=df

	# Initialize our clines variable (gives the start column of the cline for each row)
	clines = rep(num.by.vars+1,length(df[[1]]))
	# - Make grouping columns multirow - #
	for(b in rev(by.vars)) {
		
		# Create a groups variable for all by.vars up to the current one
		groups=rep("",length(df[[b]]))
		for(by.vars.index in 1:b) {
			groups = paste(groups,df.original[[by.vars.index]],sep="")
		}
		# Add multirow code to current column according to the groups pattern
		df[[b]] <- as.character(df[[b]])
		rle.lengths <- rle(groups)$lengths
		first <- !duplicated(groups)
		df[[b]][!first] <- ""
		df[[b]][first] <- paste("\\multirow{", rle.lengths, "}{*}{", df[[b]][first], "}")
		
		# Store this by.var's information in the clines variable
		clines[first]=b
	}
	
	# Specify horizontal lines wherever all combinations of grouping variables change
	df[[1]]<-paste("\\cline{",clines,"-",numcols,"}",df[[1]],sep="")
	
	
	align.by.vars = sapply(list(rep("|c", (length(by.vars)+1) )),paste,collapse="")
	align.other.vars = sapply(list(rep("r|", (length(colnames(df))-length(by.vars)) )),paste,collapse="")
	align.df = paste("|", align.by.vars , "|" , align.other.vars ,sep="")

	xt=xtable(df, align = align.df,...)
	
	
	return(xt)
	
}


bytable = function(datavec,indices,ops=c(quote(mean)),ops.desc=list(mean="Mean"),na.rm=TRUE) {
	groups=as.character()
	combinations.others=c()

	# indices should be a list of grouping vectors, just like you would pass to -by-, but with sensible names for each vector
	if(!is.list(indices)) { 
		stop("indices needs to be a list")
	}
	# Create a selector variable from the indices given as a list
	if(length(indices) > 1) {
		for(indexnum in length(indices):1) {
			groups=paste(groups,indices[[indexnum]],sep="")
		}
	}	
	if(length(indices)==1) {
		groups=indices[[1]]
	}
	first=!duplicated(groups)
	
	# Initialize data frame with grouping variables (indices)
	bynames=dimnames(by(datavec,indices,function(x) x=1)) # run a dummy by statement to get the name order out...highly inefficient...could use indices.levels=lapply(indices,function(x) x[!duplicated(x)]) instead, as long as we're sure the ordering is the same
	for(indexnum in length(indices):1) {
		# get the number of combinations of other index levels after this one (e.g. the number of replicates we need to make of each one in this index)
		others.selector=rep(TRUE,length(indices))
		others.selector[length(indices):indexnum]=FALSE
		numcombinations.others = prod(unlist(subset(lapply(bynames,length),others.selector)))
		# Replicate each level of this index the number of existing combinations of other indices
		newcolumn=rep(bynames[[indexnum]],each=numcombinations.others)
		
		if(indexnum==length(indices)) { # first run
			by.df=data.frame(newcolumn)
		}
		if(indexnum!=length(indices)) {
			# newcolumn is too short by some multiple so we have to fix that
			newcolumn=rep(newcolumn, length(rownames(by.df))/length(newcolumn) )
			# now attach our new column
			by.df=cbind(by.df,newcolumn)
		}
	}
	
	colnames(by.df)<-rev(names(indices))

	

	# Run -by- for each operation
	for(op in ops) {
		by.df[[deparse(op)]]=as.numeric(by(datavec,indices,eval(op)))
		colnames(by.df)[ colnames(by.df)==deparse(op) ] = ops.desc[[deparse(op)]]
	}
	
	if(na.rm) {
		#this assumes that the NA's in the last one will be the same as the NA's in all ops
		by.df=subset(by.df,!is.na(by.df[[length(by.df)]]))
	}
	
	return(by.df)
}

# Convert a `by` object to a data.frame (reducing dimensionality and adding repetition as necessary)
as.data.frame.by <- function( x, row.names=NULL, optional=FALSE, colnames=paste("IDX",seq(length(dim(x))),sep="" ), na.rm=TRUE, ... ) {
  num.by.vars <- length(dim(x))
	res <- melt(unclass(x))
  if(na.rm) { res <- na.omit(res) }
	colnames(res)[seq(num.by.vars)] <- colnames
  if(!is.null(row.names)) { row.names(res) <- row.names }
	res <- res[ do.call(order,res[ , seq(num.by.vars)] ) , ] # Sort the results by the by vars in the heirarchy given
	res
}


# Convert a data.frame's variables to character if they are factor or ordered
remove.factors = function(df) {
	for(varnum in 1:length(df)) {
		if("factor" %in% class(df[,varnum])) {
			df[varnum]=as.character(df[,varnum])
		}
	}
	return(df)
}

# Shift a vector over by n spots
# wrap adds the entry at the beginning to the end
# pad does nothing unless wrap is false, in which case it specifies whether to pad with NAs
shift <- function(x,...) {
	require(plyr)
	UseMethod("shift",x)
}
shift.default <- function(x,n=1,wrap=TRUE,pad=FALSE,...) {
	if(length(x)<abs(n)) { 
		#stop("Length of vector must be greater than the magnitude of n \n") 
	}
	if(n==0) { 
		return(x) 
	} else if(length(x)==n) { 
		# return empty
		length(x) <- 0
		return(x)
	} else if(n>0) {
		returnvec <- x[seq(n+1,length(x) )]
		if(wrap) {
			returnvec <- c(returnvec,x[seq(n)])
		} else if(pad) {
			returnvec <- c(returnvec,rep(NA,n))
		}
	} else if(n<0) {
		returnvec <- x[seq(1,length(x)-abs(n))]
		if(wrap) {
			returnvec <- c( x[seq(length(x)-abs(n)+1,length(x))], returnvec )
		} else if(pad) {
			returnvec <- c( rep(NA,abs(n)), returnvec )
		}
		
	}
	return(returnvec)
}
shift.data.frame <- function(x,...) {
  colwiseShift <- colwise(shift.default)
  colwiseShift(x,...)
}

# Classify values into groups based on which numbers they're between
between = function(vec,cutpoints) {
	n <- length(cutpoints) - 1
	cutpoints_hi <- cutpoints[seq(2,length(cutpoints))]
	cutpoints_lo <- cutpoints[seq(length(cutpoints)-1)]
	tweened <- rep(NA,length(vec))
	for(i in seq(n)) {
		tweened[vec >= cutpoints_lo[i] & vec < cutpoints_hi[i]] <- i
	}
	tweened[vec >= cutpoints_hi[n]] <- n
	
	return(tweened)
}

# Bin values into n equally wide groups
bin = function(vec, n=10) {
	cutpoints <- quantile(vec,probs=seq(0,1,1/n))
	cutpoints_hi <- cutpoints[seq(2,length(cutpoints))]
	cutpoints_lo <- cutpoints[seq(length(cutpoints)-1)]
	
	binned <- rep(NA,length(vec))
	for(i in seq(n)) {
		binned[vec >= cutpoints_lo[i] & vec < cutpoints_hi[i]] <- cutpoints_lo[i]
	}
	binned[vec >= cutpoints_hi[n]] <- cutpoints_lo[n]
	
	return(binned)
}


# Kludgy horizontal histogram function (really should just fix the lattice equivalent)
hist_horiz = function(formula, data,n=20) {
	# Import data from formula
	parsed <- latticeParseFormula(formula,data=data)
	dt <- parsed$right
	gp <- as.integer(as.factor(parsed$condition[[1]]))
	num_gps <- length(table(parsed$condition[[1]]))
	cutpoints <- quantile(dt,probs=seq(0,1,1/n))
	dt.tweened <- between(dt,cutpoints)
	
	# Check inputs for validity
	if(!is.null(parsed$left)) {
		stop("Left side of formula is not null.\n")
	}
	if(length(parsed$condition)!=1) {
		stop("NO higher-level ordering supported yet.\n")
	}
	# - Plot - #
	# Overall range of plot
	#plot.window(xlim = c(1,num_gps*150), ylim=range(dt))
	par(mfcol=c(1,num_gps))
	
	
	# Loop through groups and do barplots
	for(g in seq(num_gps)) {
		dt.gp <- subset(dt.tweened,gp==g)
		barplot(table(dt.gp),horiz=TRUE,axes=FALSE,ylim=c(1,n))
	}
}



# panel function for xyplot to create lattice plots of the empirical CDF
panel.ecdf <- function(x,y,lines=TRUE,...) {
	require(grid)
	require(lattice)
	if(length(x)!=0 & length(y) != 0) {
		if(lines==FALSE) {
			panel.xyplot(x,ecdf(y)(y),...)  #ecdf() returns a function which we then have to feed a vector back into to get the ecdf
		} else {
			# Sort them by x so the lines come out in the correct order
			sorted <- rbind(x,y)[,order(x,y)]
			panel.lines(sorted[1,],ecdf(sorted[2,])(sorted[2,]),...)
		}
	}
}

# Panel function for densityplot to add in descriptives as text
panel.densityplot.enhanced <- function(x,...) {
	require(grid)
	require(lattice)
	if(length(x)!=0) {
		panel.densityplot(x,...)
		# - Add in mean and SD
		# Compute best locations
		text.x <- unit(0.5,"npc")
		text.y <- unit(.95,"npc")
		text.y.sd <- text.y - unit(9,"points")
		# Now draw them
		grid.text(x=text.x,y=text.y,label=paste("Mean:",round(mean(x),digits=1),sep=""),gp=gpar(fontsize=7))
		grid.text(x=text.x,y=text.y.sd,label=paste("SD:",round(sd(x),digits=1),sep=""),gp=gpar(fontsize=7))
	}
}

# Bar plot divided by three groupings
compareplot <- function(formula, data.frame, show.outlines=FALSE,main="",x.label="",div.axis.major = 10,div.axis.minor = 20,log.x=FALSE,colors.plot=c("salmon","mediumblue","olivedrab","cyan","brown","darkgreen","purple"),panel="panel.tuftebox",box.width.large.scale = .4,box.width.small.scale = .25,box.show.mean=TRUE,box.show.box=FALSE,box.show.whiskers=FALSE,...) {
	require(grid)
	require(lattice)
	grid.newpage()
	# -- Initialize variables and configure -- #
	gp1.titlesize = 9 # Size in points for the titles of group 1 variables
	gp2.titlesize = 7 # Size in points for the titles of group 1 variables
	titlesize.padding = 3 # Number of points to add (half on top half on bottom) to the vertical padding of the title cells
	main.titlesize = 12
	
	x.padding = 25 # Spacing on the top/bottom of our plotting areas, in points
	y.scale = .7 # Scale the (rotated) y axis...this affects how much space there is between gp2 windows
	
	# initialize variables for later
	densities.x = densities.y = as.numeric()
	
	# -- Verify proper conditions -- #
	if(!exists("data") | !exists("formula")) {
		stop("Must include all required parameters")
	}
	if(panel!="panel.tuftebox") {
		stop("Only panel.tuftebox is currently supported")
	}
	
	# -- Interpret formula, etc. etc. -- #
	# - Parse formula - #
	lpf <- latticeParseFormula(formula,data=data.frame)
	if(length(lpf$condition)!= 3) { stop("Must provide three conditions") }
	x <- lpf$right
	if(log.x==TRUE) { x <- log10(x) }
	if(any(is.na(x))) {
		stop("NA's not allowed in your data vector")
	}
	gp <- lpf$condition
	# Confirm they're all factors
	for(i in seq(length(gp))) { 	if(!is.factor(gp[[i]])) { stop("All grouping variables must coerce to factors") }	}
	# Store the descriptives for later use
	levels.gp1 <- levels(gp[[1]])
	levels.gp2 <- levels(gp[[2]])
	levels.gp3 <- levels(gp[[3]])
	num.gp1 <- length(levels.gp1)
	num.gp2 <- length(levels.gp2)
	num.gp3 <- length(levels.gp3)
	
	# -- Draw -- #
	# - Divide into main sections - #
	# Draw main layout (3 main sections: Title, Graphs, Legend)
	pushViewport(viewport(layout=grid.layout(4,2,heights=unit(.99*c(.05,.8,.05,.1),"npc"),widths=unit(.95*c(.1,.9),"npc") ), name="Main"))
	# Draw Title viewport
	seekViewport("Main")
	pushViewport(viewport(layout.pos.col=2,layout.pos.row=1,name="Title"))
	if(show.outlines) {grid.rect() }
	# Draw Graphs viewport - layout for group 1 titles, plus the upper x.padding (lower x.padding is in the group 1 viewport)
	seekViewport("Main")
	pushViewport(viewport(layout.pos.col=2,layout.pos.row=2,name="Graphs",layout=grid.layout(3,num.gp1,heights=unit(c(gp1.titlesize+titlesize.padding,x.padding,1),c("points","points","null")),widths=unit(1/num.gp1,"npc") )))
	if(show.outlines) {grid.rect() }
	# Draw Legend viewport
	seekViewport("Main")
	legend.n.col=ifelse(num.gp3>4,ceiling(num.gp3/2),num.gp3)
	legend.n.row=ifelse(num.gp3>4,2,1)
	pushViewport(viewport(layout.pos.col=2,layout.pos.row=4,name="Legend",layout=grid.layout(legend.n.row,legend.n.col),width=unit(.95,"npc"),height=unit(.95,"npc")))
	if(show.outlines) {grid.rect() }
	# Draw Axis viewport - layout group 1&2 titles, plus x.padding on top/bottom
	seekViewport("Main")
	pushViewport(viewport(layout.pos.col=1,layout.pos.row=2,name="Axis",layout=grid.layout(5,1,heights=unit(c(gp1.titlesize+titlesize.padding,x.padding,1,x.padding,gp2.titlesize+titlesize.padding),c("points","points","null","points","points") ))   )) # layout adjusts for the gp1 and gp2 titles
	if(show.outlines) {grid.rect() }
	
	# - Figure out the max density for the scale - #
	x.range <- range(x)
	y.range <- c(0,1)
	y.range.scaled <- y.range + c(diff(y.range)*(1-y.scale)/2,-diff(y.range)*(1-y.scale)/2) # Scale our range in ways that avoid distortion around 0
	
	# - Divide into num.gp1 graph sections - #
	for(gp1.i in seq(num.gp1)) {
		# - Title areas for gp1 - #
		seekViewport("Graphs")
		pushViewport(viewport(layout.pos.col=gp1.i,layout.pos.row=1, name=paste("gp1.",gp1.i,".title",sep="")))
		grid.text(label=as.character(levels.gp1[gp1.i]), gp=gpar(fontsize=gp1.titlesize))
		grid.rect() # We want this one to show since we're not including padding
		# - Divide into num.gp2 graph sections - #
		seekViewport("Graphs")
		# Graph areas to be divided according to gp2, plus x.padding
		pushViewport(viewport(layout.pos.col=gp1.i,layout.pos.row=3, name=paste("gp1.",gp1.i,sep=""),layout=grid.layout(3,num.gp2,heights=unit(c(1,x.padding,gp2.titlesize+titlesize.padding),c("null","points","points")),widths=unit(1/num.gp2,"npc")) ))
		for(gp2.i in seq(num.gp2)) {
			# - Title areas for GP2 - #
			seekViewport(paste("gp1.",gp1.i,sep=""))
			pushViewport(viewport(layout.pos.col=gp2.i,layout.pos.row=3, name=paste("gp1.",gp1.i,"_gp2.",gp2.i,".title",sep="")))
			# Label in rotated viewport
			pushViewport(viewport(angle=-90))
			grid.text(label=as.character(levels.gp2[gp2.i]), gp=gpar(fontsize=gp2.titlesize))
			if(show.outlines) {grid.rect() }
			# - graph areas for GP2 - #
			seekViewport(paste("gp1.",gp1.i,sep=""))
			pushViewport(viewport(layout.pos.col=gp2.i,layout.pos.row=1, name=paste("gp1.",gp1.i,"_gp2.",gp2.i,sep=""), yscale=unit(x.range,"native") ))
			if(show.outlines) {grid.rect() }
			#pushViewport(viewport(angle=-90,width=convertUnit(unit(1,"npc"),"npc","y","dimension","x","dimension"),height=convertUnit(unit(1,"npc"),"npc","x","dimension","y","dimension"),xscale=unit(x.range,"native"),yscale=unit(y.range.scaled,"native") ))
			for(gp3.i in seq(num.gp3)) {
				# - Create our subsetted (GP3) data gp3.n times on the same gp2 plot viewport - #
				x.gp1gp2gp3 <- subset(x,as.numeric(gp[[1]])==gp1.i & as.numeric(gp[[2]])==gp2.i & as.numeric(gp[[3]])==gp3.i)	
				if(length(x.gp1gp2gp3>0)) { # Handle missing panels
					if(panel=="panel.tuftebox") { # Draw Tufte boxplots if specified
						# Calculate things
						loc.y = y.range.scaled[1]+(gp3.i-.5)*(1/num.gp3)*diff(y.range.scaled) # y coordinate just shifts based on how many things we're plotting in this viewport
						quantiles = quantile(x.gp1gp2gp3)
						iqr = diff(quantiles[c("25%","75%")])
						box.width.tiny = unit(1,"points")
						box.width.small = box.width.small.scale*(1/num.gp3)*diff(y.range.scaled)
						box.width.large = box.width.large.scale*(1/num.gp3)*diff(y.range.scaled)
						# Min/max line (actually goes to 1.5IQR past Q1 or Q3)
						min.reduced = max(quantiles["25%"]-1.5*iqr,min(x.gp1gp2gp3)) # Use true min if 1.5*iqr exceeds it
						max.reduced = min(quantiles["75%"]+1.5*iqr,max(x.gp1gp2gp3)) # Use true max if 1.5*iqr exceeds it
						grid.lines(y=unit(c(min.reduced,quantiles["25%"]),"native"),x=loc.y,default.units="native",gp=gpar(col=colors.plot[gp3.i])) # Min line
						grid.lines(y=unit(c(max.reduced,quantiles["75%"]),"native"),x=loc.y,default.units="native",gp=gpar(col=colors.plot[gp3.i])) # Max line						
						if(box.show.whiskers==TRUE) { # Draw "whiskers" on the min/max
							grid.lines(y=unit(min.reduced,"native"),x=unit(loc.y,"native")+(c(1,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i])) # Min whisker
							grid.lines(y=unit(max.reduced,"native"),x=unit(loc.y,"native")+(c(1,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i])) # Max whisker
						}
						# Q1-Q3 line, shifted just slightly
						if(box.show.mean==FALSE) { # Only show if we're not cluttering it up with the mean/SD diamond already
							# Vertical line
							grid.lines(y=unit(quantiles[c("25%","75%")],"native"),x=unit(loc.y,"native")-box.width.tiny,default.units="native",gp=gpar(col=colors.plot[gp3.i]))
							if(box.show.box==TRUE) { # Show the right side of the box also, if specified
								grid.lines(y=unit(quantiles[c("25%","75%")],"native"),x=unit(loc.y,"native")+box.width.tiny,default.units="native",gp=gpar(col=colors.plot[gp3.i]))
								grid.lines(y=unit(quantiles["25%"],"native"),x=unit(loc.y,"native")+(c(1,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
								grid.lines(y=unit(quantiles["75%"],"native"),x=unit(loc.y,"native")+(c(1,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
							} else {
								# Small horizontal lines to connect it in -- draw only the one to the left half if we're not drawing the full box
								grid.lines(y=unit(quantiles["25%"],"native"),x=unit(loc.y,"native")+(c(0,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
								grid.lines(y=unit(quantiles["75%"],"native"),x=unit(loc.y,"native")+(c(0,-1)*box.width.tiny),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
							}
						}
						# Outliers as points
						outliers=subset(x.gp1gp2gp3,x.gp1gp2gp3<min.reduced | x.gp1gp2gp3>max.reduced )
						if(length(outliers)>0) {
							grid.points(y=unit(outliers,"native"),x=rep(loc.y,length(outliers)),default.units="native",gp=gpar(col=colors.plot[gp3.i],cex=.2),pch=4 )
						}
						# Quartiles 1 and 3
						if(box.show.mean==TRUE) {
							grid.lines(y=unit(rep(quantiles[c("25%")],2),"native"),x=c(loc.y-box.width.small,loc.y+box.width.small),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
							grid.lines(y=unit(rep(quantiles[c("75%")],2),"native"),x=c(loc.y-box.width.small,loc.y+box.width.small),default.units="native",gp=gpar(col=colors.plot[gp3.i]))
						}
						# Median
						grid.points(y=unit(median(x.gp1gp2gp3),"native"),x=loc.y,default.units="native",gp=gpar(col=colors.plot[gp3.i],cex=.3),pch=15)
						# Mean (+/- SD)
						if(box.show.mean==TRUE) {
							meanlines.x = c(mean(x.gp1gp2gp3),mean(x.gp1gp2gp3)-sd(x.gp1gp2gp3),mean(x.gp1gp2gp3),mean(x.gp1gp2gp3)+sd(x.gp1gp2gp3),mean(x.gp1gp2gp3) ) # start at the mean on the left and loop around
							meanlines.y = c(loc.y-box.width.large,loc.y,loc.y+box.width.large,loc.y,loc.y-box.width.large)
							grid.lines(y=meanlines.x,x=meanlines.y,default.units="native",gp=gpar(col=colors.plot[gp3.i]) )
						}
					}
				}
			}
		}
	}
	
	# - Draw in title - #
	seekViewport("Title")
	grid.text(label=main, gp=gpar(fontsize=main.titlesize))

	# - Draw in axis - #
	seekViewport("Axis")
	pushViewport(viewport(layout.pos.col=1,layout.pos.row=3,name="Axis.actual",yscale=unit(x.range,"native") ))
	if(show.outlines) {grid.rect()}
	if(show.outlines) {grid.rect() }
	# Major axis tick marks
	x.range.magnitude <- diff(x.range)
	x.seq=x.range[1]+(1/div.axis.major)*seq(0,div.axis.major)*x.range.magnitude
	mat.x <- rbind(rep(.85,div.axis.major+1),rep(1,div.axis.major+1))
	mat.y <- matrix(rep(x.seq,each=2),nrow=2)
	grid.polyline(x=mat.x,y=mat.y,id.lengths=rep(2,div.axis.major+1),default.units="native")
	# Major axis value labels (depends on value of log.x)
	if(log.x==FALSE) {
		round.digits=-floor(log10(x.range.magnitude))+1
		x.labels <- round(x.seq,round.digits) # Amount of rounding adapts to our range
	} else {
		x.seq.label=10^x.seq
		x.labels=round_sigfig(10^x.seq,1)
	}
	grid.text(label=as.character(x.labels),x=.8,y=x.seq,gp=gpar(fontsize=9),just=c("right","center"),default.units="native")
	# Minor axis tick marks
	mat.x <- rbind(rep(.925,div.axis.minor+1),rep(1,div.axis.minor+1))
	mat.y <- matrix(rep(x.range[1]+x.range.magnitude*(1/div.axis.minor)*seq(0,div.axis.minor),each=2),nrow=2)
	grid.polyline(x=mat.x,y=mat.y,id.lengths=rep(2,div.axis.minor+1),default.units="native")
	# Line on right (the axis itself)
	grid.lines(x=c(1,1),y=c(0,1))
	# Axis title
	upViewport(0)
	pushViewport(viewport(angle=90,x=unit(0.01,"npc"),width=convertUnit(unit(1,"npc"),"npc","y","dimension","x","dimension"),height=convertUnit(unit(.3,"npc"),"npc","x","dimension","y","dimension"))) # viewport rotated 90 degrees
	grid.text(label=x.label,just=c("center","top"),gp=gpar(fontsize=10))
	
	# - Draw in legend - #
	seekViewport("Legend")
	for(gp3.i in seq(num.gp3)) {
		# Handle placement for multirow scenarios
		gp3.col=ifelse(gp3.i>legend.n.col,gp3.i-legend.n.col,gp3.i)
		gp3.row=ifelse(gp3.i>legend.n.col,2,1)
		pushViewport(viewport(layout.pos.col=gp3.col,layout.pos.row=gp3.row,width=unit(.9,"npc"),height=unit(.9,"npc") ))
		legend.colorbox.width=convertUnit(unit(.3,"npc"),"npc","y","dimension","x","dimension")
		# Text
		grid.text(x=unit(.15,"npc")+legend.colorbox.width,hjust=0,label=levels.gp3[gp3.i],gp=gpar(col=colors.plot[gp3.i],fontface="bold"))
		# Color box
		grid.rect(x=unit(.1,"npc"),y=unit(.5,"npc"),hjust=0,height=unit(.3,"npc"),width=legend.colorbox.width,gp=gpar(fill=colors.plot[gp3.i]))
		# Color box outline
		grid.rect(x=unit(.1,"npc"),y=unit(.5,"npc"),hjust=0,height=unit(.3,"npc"),width=legend.colorbox.width)
		popViewport()
	}
}


# Round a vector to n significant digits
round_sigfig <- function(vec,digits=2) {
	#Check inputs
	if(min(digits)<1) {
		stop("Minimum significant figure digits is 1")
	}
	# Make our vector and digits the same length
	if(length(vec)<length(digits)) {
		stop("vec should be longer than or of equal length to digits")
	}
	if(length(vec)>length(digits)) {
		digits <- rep(digits,ceiling(length(vec)/length(digits)))
		if(length(vec)==(length(digits)-1) ) { # Handle odd ratios of length(vec)/length(digits)
			digits <- digits[seq(length(vec))]
		}
	}
	vec.rounded <- round(vec,-floor(log10(vec))+digits-1)
	return(vec.rounded)
}

# Return a vector of the days of the week, in order
daysofweek <- function(start.day="Monday") {
	wkdays <- c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday')
	wkdays <- rep(wkdays,2)
	selector <- unlist(lapply(wkdays,function(x) x==start.day)) #selects start day and one past end day
	selector.numeric <- seq(length(wkdays))[selector]
	selector.numeric[2] <- selector.numeric[2]-1 #Move to last day
	return(wkdays[selector.numeric[1]:selector.numeric[2] ] )
}

# Insert a title page containing the given text.  Good for breaking up sections of plot PDFs.
title.page.new <- function(title.text="") {
	plot.new()
	text(.5,.7,title.text)
}


# Convenience function to return the last/first element of a vector
last <- function(vec) {
	return(vec[length(vec)])
}
first <- function(vec) {
	return(vec[1])
}

# Create a data.frame of quantiles for feeding into e.g. categorize()
quantile_cutpoints <- function(vec,probs) {
	qtiles <- quantile(vec,probs=probs)
	hi <- shift(qtiles,n=1,wrap=FALSE)
	lo <- qtiles[seq(length(hi))]
	deciles <- data.frame(low=lo,high=hi)
	rownames(deciles) <- paste(names(lo),names(hi),sep="-")
	return(deciles)
}

# Categorize a vector based on a data.frame with two columns, the low and high end points of each category
categorize <- function(vec,cutpoints.df,match.min=TRUE,names=TRUE) {
	# Categorize a single point; used with apply below
	cat.one <- function(x,cutpoints.df,names) {
		if (length(x)!=1 & class(x) != "numeric" ) {	stop("x must be a single number.") }
		# Subtract a little from our minimum so it matches
		if(match.min) { cutpoints.df[1,1] <- cutpoints.df[1,1]-.00001 }
		selector <- (x > cutpoints.df[,1]) & (x <= cutpoints.df[,2])
		if("TRUE" %in% names(table(selector)) ) {
			if ( table(selector)[["TRUE"]] != 1 ) { warning(x,"matched more than one category.");return(NA) }
		} else { # If there were no TRUEs then we had 0 match.
			warning(x,"matched zero categories")
			return(NA)
		}
		# Return names or row numbers
		if(names) { 
			return( rownames(cutpoints.df)[selector] )
		} else {
			return(seq(nrow(cutpoints.df))[selector])
		}
		
	}
	sapply(vec,cat.one,cutpoints.df=cutpoints.df,names=names)
}

# Add in methods to handle LME objects in xtable
xtable.lme <- function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, display = NULL, beta.names = NULL, ...) {
	require(xtable)
	return(xtable.summary.lme(summary(x), caption = caption, label = label, align = align, digits = digits, display = display, beta.names = beta.names))
}
xtable.summary.lme <- function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, display = NULL, beta.names=NULL, ...) {
	require(xtable)
	# Grab our data
	x <- data.frame(x$tTable[,-3], check.names = FALSE)
	# Update beta names if specified
	if(!is.null(beta.names)) {
		if(length(beta.names) != nrow(x))	stop(paste("beta.names must have",nrow(x),"elements."))
		rownames(x) <- beta.names
	}
	# Set attributes and return for xtable to deal with
	class(x) <- c("xtable", "data.frame")
	caption(x) <- caption
	label(x) <- label
	align(x) <- switch(1 + is.null(align), align, c("r", "r", "r", "r", "r"))
	digits(x) <- switch(1 + is.null(digits), digits, c(0, 4, 4, 2, 4))
	display(x) <- switch(1 + is.null(display), display, c("s", "f", "f", "f", "f"))
	return(x)
}

# xyplot panel function with rug plots on x and y axes
panel.xyplot_rug <- function(x,y,rug.color="grey",...) {
	panel.xyplot(x,y,...)
	grid.segments(x, unit(0, "npc"), x, unit(3, "mm"),default.units="native",gp=gpar(col=rug.color))
	grid.segments(unit(0, "npc"),y, unit(3, "mm"),y, default.units="native",gp=gpar(col=rug.color))
}

# Returns number of distinct observations in each column of a data frame or in a vector
distinct <- function(input,na.rm=TRUE) {
	if(na.rm!=TRUE) {	
		exclude=c() 
	} else {	
		exclude=c(NA,NaN) 
	}
	return(switch(class(input),
		data.frame=unlist(lapply(input,function(x) length(table(x,exclude=exclude)) )),
		numeric=length(table(input,exclude=exclude)),
		integer=length(table(input,exclude=exclude)),
		character=length(table(input,exclude=exclude))
	))
}

# Convert a by object into a matrix (usage: as.matrix(by(...)) )
	#! only tested on a 2d object
as.matrix.by <- function(x, ...) {
	if(class(x)!= "by") { stop("Must input a by object") }
	ul.x <- unlist(x)
	by.mat <- matrix(data=ul.x,ncol=length(x),nrow=length(ul.x)/length(x) )
	colnames(by.mat) <- names(x)
	rownames(by.mat) <- names(x[[1]])
	
	return(by.mat)
}

# Create a vector that starts with a given number and widens out
searchPattern <- function(center=0,length=5,interval=1) {
	require(gdata)
	vec.up <- seq(center+interval,center+interval*length,interval)
	vec.down <- seq(center-interval,center-interval*length,-interval)
	return(c(center,as.numeric(interleave(vec.up,vec.down))))
}

# Replicate elements of vectors and lists to match another vector
rep_along <- function( x, along.with ) {
  rep( x, times=length(along.with) )
}

# Convert a character string to numeric, dropping any irrelevant characters
destring <- function(x,keep="0-9.") {
  return( as.numeric(gsub(paste("[^",keep,"]+",sep=""),"",x)) )
}

# reshapeasy: Version of reshape with way, way better syntax
# Written with the help of the StackOverflow R community
# x is a data.frame to be reshaped
# direction is "wide" or "long"
# vars are the names of the (stubs of) the variables to be reshaped (if omitted, defaults to everything not in id or vary)
# id are the names of the variables that identify unique observations
# vary is the variable that varies.  Going to wide this variable will cease to exist.  Going to long it will be created.
# omit is a vector of characters which are to be omitted if found at the end of variable names (e.g. price_1 becomes price in long)
# ... are options to be passed to stats::reshape
reshapeasy <- function( data, direction, id=(sapply(data,is.factor) | sapply(data,is.character)), vary=sapply(data,is.numeric), omit=c("_","."), vars=NULL, ... ) {
  if(direction=="wide") data <- stats::reshape( data=data, direction=direction, idvar=id, timevar=vary, ... )
  if(direction=="long") {
    varying <- which(!(colnames(data) %in% id))
    data <- stats::reshape( data=data, direction=direction, idvar=id, varying=varying, timevar=vary, ... )
  }
  colnames(data) <- gsub( paste("[",paste(omit,collapse="",sep=""),"]$",sep=""), "", colnames(data) )
  return(data)
}

# Judicious apply: Apply function to only the specified columns
# Takes a data.frame and returns a data.frame with only the specified columns transformed
japply <- function(df, sel, FUN=function(x) x, ...) {
  df[,sel] <- sapply( df[,sel], FUN, ... )
  df
}

# Stacks lists of data.frames (e.g. from replicate() )
stack.list <- function( x, label=FALSE, ... ) {
  ret <- x[[1]]
  if(label) { ret$from <- 1 }
  if(length(x)==1) return(ret)
  for( i in seq(2,length(x)) ) {
    new <- x[[i]]
    if(label) { new$from <- i }
    ret <- rbind(ret,new)
  }
  return(ret)
}

# Cleans a Mechanical Turk file in a pretty standard way
   # post.process is a user-defined function run on the data.frame after it is cleaned
   # drop.duplicates has three modes.  
      # If false, no duplicates are dropped
      # If true or character use built-in function (drop only if all (TRUE) or the listed (character) "Answer." columns identical).  
      # Or can be a user-defined function that takes in a data.frame and returns a logical vector corresponding to whether to drop a row (TRUE) or not (FALSE).
cleanTurked <- function( dat, drop.duplicates=TRUE , post.process=identity ) {
  inputs <- colnames(dat)[ grep( "Input\\.", colnames(dat) ) ]
  answers <- colnames(dat)[ grep( "Answer\\.", colnames(dat) ) ]
  # Eliminate rejected
  dat <- subset( dat, ! grepl("[xX]",Reject) )
  # Deduplicate
  if( is.logical(drop.duplicates) | is.character(drop.duplicates) ) {
    if(is.character(drop.duplicates) || drop.duplicates ) {
      if(is.logical(drop.duplicates)) { drop.duplicates <- answers } # if it was a logical, drop only rows where all columns were duplicates
      all.identical.byCol <- ddply( dat[,c("HITId",answers)], .(HITId), function(x) sapply(subset(x,select=-c(HITId)), function(v) length(unique(v))==1 ) )
      identical <- apply( all.identical.byCol[,drop.duplicates,drop=FALSE], 1, all ) # Only drop rows where all the desired columns are identical within a group
      ident.df <- data.frame( HITId=all.identical.byCol$HITId, identical=identical )
      dat <- ddply( merge(dat, ident.df), .(HITId), function(x) {
        if(x$identical[1]) return(x[1,])
        else return(x)
      } )
      dat <- subset(dat,select=c(-identical))
    }
  } else if(is.function(drop.duplicates)) { # Use the user's function
    stop("Not yet implemented.  Please use post.process for now.\n")
  } else { stop("drop.duplicates must be a logical or function.\n")}
  # Eliminate all the MT-added columns
  ret <- subset(dat, select=c(inputs,answers) )
  colnames(ret) <- sub( "(Input|Answer)\\.", "", colnames(ret), perl=TRUE )
  ret <- post.process(ret)
  ret
}

# Autoplot method for microbenchmark objects: Prettier graphs for microbenchmark using ggplot2
autoplot.microbenchmark <- function(object, ..., y_max=max(by(object$time,object$expr,uq)) * 1.05 ) {
  uq <- function(x) { quantile(x,.75) }  
  lq <- function(x) { quantile(x,.25) }
  y_min <- 0
  p <- ggplot(object,aes(x=expr,y=time)) + coord_cartesian(ylim = c( y_min , y_max )) 
  p <- p + stat_summary(fun.y=median,fun.ymin = lq, fun.ymax = uq, aes(fill=expr))
  return(p)
}