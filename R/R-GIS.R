# Miscellaneous functions to do GIS-like manipulation of spatial objects
# source("E:/xdrive/projects/Ghana_air/dofiles/R-GIS.R")

# Ari Friedman
# Version 0.9.0
# 8/21/2007

#library(fields)
#library(maps)
#library(geoR)
#library(sp)
#library(maptools)
#library(shapefiles)
#library(spatial)
#library(rgdal)
#library(mapproj)
#library(akima)
#library(mgcv)
#library(lattice)
#library(caTools)
#library(pspline)
#library(grid)
#library(RODBC)

# -*- Define functions for use by other programs -*-

# Simple distance function
simpledist = function(points) {
	# points is a 2x2 matrix, where columns are x,y and rows are the two points
	distance = sqrt( ( points[1,1] - points[2,1] )^2 + ( points[1,2] - points[2,2] )^2 )
	return(distance)
}

# Convert a SpatialLinesDataFrame to a single line matrix with associated segment information
SLDFtoLine = function(lineDF,orderXY=FALSE,segments=TRUE) {
	# Initialize variables
	X = c()
	Y = c()
	Segment = c()
	# Assemble the complete line segment
	numsegments = length(lineDF@lines)
	
	for(seg in 1:numsegments) {
		xnew=lineDF@lines[[seg]]@Lines[[1]]@coords[,1]
		ynew=lineDF@lines[[seg]]@Lines[[1]]@coords[,2]
	
		X = append(X,xnew)
		Y = append(Y,ynew)
		if(segments) {
			Segment = append(Segment,rep.int(seg,length(xnew)))
		}
		else {
			Segment = append(Segment,rep.int(1,length(xnew)))
		}
	}

	# Order the points in increasing X,Y order
	if(orderXY) {
		Y = Y[order(X)]
		Segment = Segment[order(X)]
		X = X[order(X)]
	}
	

	line = data.frame(X,Y,Segment)
	line

}

# Smooth the line segments in a SpatialLinesDataFrame
smoothLines=function(lineDF) {
	# Initialize variables
	#lineX = c()
	#lineY = c()
	#lineSeg = c()
	
	# Assemble the complete line segment
	numsegments = length(lineDF@lines)
	#-Not using this currently
	#for(seg in 1:numsegments) {
	#	lineX = append(lineX,lineDF@lines[[seg]]@Lines[[1]]@coords[,1])
	#	lineY = append(lineY,lineDF@lines[[seg]]@Lines[[1]]@coords[,2])
	#	lineSeg = append(lineSeg,rep.int(seg,length(lineDF@lines[[seg]]@Lines[[1]]@coords[,2])))
	#}
	#spline=sreg(lineX[lineSeg==seg],lineY[lineSeg==seg])
	#lineX[lineSeg==seg]=spline$predicted$x
	#lineY[lineSeg==seg]=spline$predicted$y
	
	# Break the line segment at key points (borrow this from the SNAKES techniques)
	
	# Take our new line segments and smooth them
	for(seg in 1:numsegments) {
		#Order them for the splines
		plot(x,y,main=paste(seg))
		x = lineDF@lines[[seg]]@Lines[[1]]@coords[,1][order(lineDF@lines[[seg]]@Lines[[1]]@coords[,1])]
		y = lineDF@lines[[seg]]@Lines[[1]]@coords[,2][order(lineDF@lines[[seg]]@Lines[[1]]@coords[,1])]
		#Fit and return
		if(length(x)>4) {
			fit <- smooth.Pspline(x, y, method=3)
		}
		fit
		#spline=sreg(lineDF@lines[[seg]]@Lines[[1]]@coords[,1],lineDF@lines[[seg]]@Lines[[1]]@coords[,2])
		#lineDF@lines[[seg]]@Lines[[1]]@coords[,1] = spline$predicted$x
		#lineDF@lines[[seg]]@Lines[[1]]@coords[,2] = spline$predicted$y
	}
}

# Take a grid of regularly spaced points (such as those output by the centroids of Arc's fishnet function) and convert it to various grid data types
pointgrid2SpatialPolygons=function(df,type) {
	# Check that df is a SpatialPointsDataFrame
	if (class(df)[[1]] != "SpatialPointsDataFrame" & class(df)[[1]] !="SpatialPixelsDataFrame") {
		return(-99)
	}
	
	#Output SpatialGrid
	if (type=="SpatialGrid") {
		return(SpatialGrid(points2grid(df)));
	}
	
	#Output SpatialPolygons(DataFrame)
	df_sp=as.SpatialPolygons.GridTopology(points2grid(df))
	if(type=="SpatialPolygons") {
		return (df_sp)
	}
	if(type=="SpatialPolygonsDataFrame") {
		data=as.data.frame(rep(0,length(df_sp@plotOrder)))
		rownames(data)<-getSpPPolygonsIDSlots(df_sp)
		return(SpatialPolygonsDataFrame(df_sp,data))
	}
	
	#Otherwise
	return(-98)
}


# Subset SpatialPolygonsDataFrame or SpatialPointsDataFrame
subsetSPDF = function(SPDF,tf,...) {
	selected_data <- subset(SPDF@data, tf)
	if(class(SPDF)=="SpatialPolygonsDataFrame") {
		SPDF_selected <- subset(SPDF@polygons, tf)
		centroids <- getSpPPolygonsLabptSlots(as.SpatialPolygons.PolygonsList(SPDF_selected))
		x <- centroids[,1]
		y <- centroids[,2]
		export <- SpatialPolygonsDataFrame(as.SpatialPolygons.PolygonsList(SPDF_selected),data=data.frame(x=x, y=y,row.names=getSpPPolygonsIDSlots(as.SpatialPolygons.PolygonsList(SPDF_selected)),selected_data),...)
	}
	if(class(SPDF)=="SpatialPointsDataFrame") {
		export <- SpatialPointsDataFrame(coordinates(SPDF)[tf,],data=selected_data,coords.nrs=SPDF@coords.nrs,proj4string=SPDF@proj4string,...)
	}
	return(export)
}

# Get the areas stored in the polygons and return them in the dataframe slot
SPDFareas = function(SPDF,colname="AREA") {
	numPolys = length(SPDF@polygons)
	areas = data.frame(rep(NA,numPolys))
	colnames(areas) <- c(colname)
	for(polyNum in 1:numPolys) {
		areas[[colname]][polyNum] <- SPDF@polygons[[polyNum]]@area
	}
	returnSPDF = SpatialPolygonsDataFrame(polygons(SPDF),data=cbind(SPDF@data,areas),match.ID=TRUE)
  return(returnSPDF)
}

# Count points in polygon
# Overlays points on polygons and create a new polygon dataset with the count of the points in that polygon
countPointsInPolys = function(points,polys,density=FALSE) {
  grid_count = overlay(x=points,y=polys) # this returns a vector of the same length as points saying which polygon ID they're in
  #-Count points in each polygon (rownames are polygon ID)
  pointcount=tapply(grid_count,grid_count,length)
  #Shift the rownames so they start at 0 instead of 1 like the other functions return
  rownames(pointcount) <- as.character(as.numeric(rownames(pointcount))-1)
  #- Merge with polygons data.frame
  pointcount.df = data.frame(pointcount,rownames=rownames(pointcount),stringsAsFactors=FALSE)
  polys.df = data.frame(polys@data,rownames=rownames(polys@data),stringsAsFactors=FALSE)
  DF = merge(polys.df,pointcount.df,all=TRUE,by.x="rownames")
  DF = DF[order(as.numeric(DF$rownames)),] # sort it
  rownames(DF) <- DF$rownames
  DF$pointcount <- as.numeric(DF$pointcount)
  returnSPDF = SpatialPolygonsDataFrame(polygons(polys),data=DF,match.ID=TRUE)
  return(returnSPDF)
}

# Find closest point to a given point's coordinates
closestPoint = function(point,points) {
	distance = sqrt( ( points[,1] - point[[1]] )^2 + ( points[,2] - point[[2]] )^2 )
	mindist_selector = (distance==min(distance))
	if( length(mindist_selector[mindist_selector==TRUE]) != 1) { # Return NA if there are 2+ points at exactly the same distance
		return(NA)
	} else {
		return(mindist_selector)
	}
}


# Reshape a spatialLinesDataFrame into a series of points with associated information (less efficient because all the segment data gets replicated over each point)
reshapeSLDF = function(SLDF,shape="long") {
	if(shape=="long"){
		# Loop over each segment
		for(seg in 1:length(SLDF@lines)) {
			numRecords = length(SLDF@lines[[seg]]@Lines[[1]]@coords[,1])
			# Longit
			longit = SLDF@lines[[seg]]@Lines[[1]]@coords[,1]
			# Latit
			latit = SLDF@lines[[seg]]@Lines[[1]]@coords[,2]
			# Data
			# Initialize variables
			SLDFdata = as.data.frame(replicate(length(SLDF@data),rep(NA,numRecords)))
			colnames(SLDFdata) <- colnames(SLDF@data)
			for(datacolnum in 1:length(SLDF@data)) {	
				SLDFdata[[datacolnum]]=rep(SLDF@data[[datacolnum]][seg],numRecords)
			}
			
			# Store our results
			if(seg==1) {
				SLDFpoints=cbind(longit,latit,SLDFdata)
			}	else {
				SLDFpoints=rbind(SLDFpoints,cbind(longit,latit,SLDFdata))
			}
		}
		return(SLDFpoints)
	}
	
}


# Interpolate points along a path
interpolatePathpoints = function(pathpoints,dens,tolerance.min=1.2,tolerance.max=50) {
	# dens actually inverse density and is in the units of the x and y in pathpoints (e.g. 1 point per density meters)
	# tolerance.min is in the proportion of the density (e.g. 1.2 means we'll fill in gaps 20% greater than the density size)

	# - Compute distance to next point - #
	#(shift the entire matrix down one row and calculate the distance all at once)
	pathpoints$xNext=c(pathpoints$x[2:length(pathpoints$x)],NA)
	pathpoints$yNext=c(pathpoints$y[2:length(pathpoints$y)],NA)
	pathpoints$distToNext = sqrt((pathpoints$x-pathpoints$xNext)^2 + (pathpoints$y-pathpoints$yNext)^2)
	
	# -- If distance to next point exceeds the dens by the tolerance, interpolate -- #
	# - First group them into chunks to add on our segments to the end of - #
	interpFlag=rep(0,length(pathpoints$x))
	interpFlag[pathpoints$distToNext>dens*tolerance.min & pathpoints$distToNext<dens*tolerance.max]=1
	# Shift the flag one back so our points are at the end of the group
	interpolate=rep(0,length(interpFlag))
	interpolate[2:length(interpFlag)]=interpFlag[1:length(interpFlag)-1]
	# Code the groups
	pathpoints$interpolateGroup=rep(1,length(pathpoints$x))
	for(rownum in 2:length(interpolate)) {
		pathpoints$interpolateGroup[rownum]=pathpoints$interpolateGroup[rownum-1]+interpolate[rownum]
	}
	# Flag the last group so it can be skipped for interpolation #
	numIGroups=max(pathpoints$interpolateGroup)
	# - Add our interpolated points on to the end - #
	bypp=by(pathpoints,pathpoints$interpolateGroup,function(pp) {
		if(pp$interpolateGroup[1]!=numIGroups) { # Skip the last group since it's last record is NA
			numRows = length(pp[[1]]) # get the index of the last row
			pointsToAdd = as.integer(floor(pp$distToNext[numRows] / dens)) # number of points to add
			xIncrement = ((pp$xNext - pp$x) / (pointsToAdd+1))[numRows]
			yIncrement = ((pp$yNext - pp$y) / (pointsToAdd+1))[numRows]
			pp <- expandDF(pp,numRows,pointsToAdd)

			for(pointNum in 1:pointsToAdd) {
				pp[numRows+pointNum,"x"] = pp[numRows+pointNum-1,"x"] + xIncrement
				pp[numRows+pointNum,"y"] = pp[numRows+pointNum-1,"y"] + yIncrement
			}
		}
		return(pp)
	})
	# Now assemble all and return
	returnpp=bypp[[1]]
	for(igroup in 2:numIGroups) {
		returnpp=rbind(returnpp,bypp[[igroup]])
	}
	
	#subset(,select)
	return(returnpp)
}


# Standardize latitude/longitude coordinates
cleanLatLon = function(vec) {
	if(!is.character(vec))		stop("Input vector must be of type character")
	vec = sub("[NnEe]","",vec)
	vec = sub("[WwSs]","-",vec)
	vec = sub("^ +","",vec)
	vec = sub(" +$","",vec)
	vec = sub("^0+","",vec)
	if(length(grep(" ",vec))!=0)		stop("Your vector appears to be in decimal degrees, which are not yet implemented")
	return(as.numeric(vec))
}


# Calculate cumulative distance along a matrix of x,y coordinates
cumDist = function(coords) {
	# coords should be an [n,2] matrix of coordinates
	
	if(any(is.na(coords))) {
		return(NA)
	}
	
	if(dim(coords)[1] == 2) {
		return(simpledist(rbind(coords[1,],coords[2,])))
	}
	
	prevCoords = coords[1:dim(coords)[1]-1,]
	currCoords = coords[2:dim(coords)[1],]
	distToPrev = rep(NA,dim(coords)[1])
	for(rowNum in 2:dim(coords)[1]) {
		distToPrev[rowNum] = simpledist(rbind(prevCoords[(rowNum-1),],currCoords[(rowNum-1),]))
	}
	return(sum(distToPrev,na.rm=TRUE))
}

# Line distance in SpatialLinesDataFrame
#   varname is the name to store the distances in in the SLDF's dataframe
lineDist = function(SLDF, varname="distances") {
	numSegments = length(SLDF@lines)
	Dists = rep(NA,numSegments)
	for(segnum in 1:numSegments) {
		Dists[segnum] <- cumDist(SLDF@lines[[segnum]]@Lines[[1]]@coords)
	}
	SLDF@data[[varname]] = Dists
	return(SLDF)
}

# Create all pairwise distances of points from a SpatialPointsDataFrame
#   names is the variable name in the SPDF's dataframe used to label each point in the resulting matrix
pointDistPairwise = function(SPDF, names = "name") {
	crds <- SPDF@coords #Coordinate data
	# Set up our matrix
	numpoints <- dim(crds)[1]
	pw.mat <- matrix(rep(NaN,numpoints^2),ncol=numpoints)
	for(i in seq(numpoints)) {
		for(j in seq(numpoints)) {
			if(i>=j) {
				pw.mat[i,j] <- simpledist(rbind(crds[i,],crds[j,]))
			}
		}
	}
	colnames(pw.mat) <- SPDF@data[[names]]
	rownames(pw.mat) <- SPDF@data[[names]]
	return(pw.mat)
}

# Convert SpatialPointsDataFrame to a regular data.frame with the coordinates as "x" and "y"
SPDFtoPointsDF <- function(SPDF) {
	coords <- SPDF@coords
	colnames(coords) <- c("x","y")
	return(cbind(SPDF@data,coords))
}