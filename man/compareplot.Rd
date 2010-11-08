\name{compareplot}
\alias{compareplot}

\title{
 Bar plot divided by three groupings
}
\description{
 Bar plot divided by three groupings
}
\usage{
compareplot(formula, data.frame, show.outlines = FALSE, main = "", x.label = "", div.axis.major = 10, div.axis.minor = 20, log.x = FALSE, colors.plot = c("salmon", "mediumblue", "olivedrab", "cyan", "brown", "darkgreen", "purple"), panel = "panel.tuftebox", box.width.large.scale = 0.4, box.width.small.scale = 0.25, box.show.mean = TRUE, box.show.box = FALSE, box.show.whiskers = FALSE, ...)
}
\arguments{
  \item{formula}{
  	  Plot formula
}
  \item{data.frame}{
	Data.frame containing data
}
  \item{show.outlines}{
	Whether to include boxes around plots or leave it open
}
  \item{main}{
	Plot text
}
  \item{x.label}{
	X axis label
}
  \item{div.axis.major}{
	How many major axis ticks to use
}
  \item{div.axis.minor}{
	How many minor axis ticks to use
}
  \item{log.x}{
	Log transform the x data?
}
  \item{colors.plot}{
	Plot colors
}
  \item{panel}{
	Panel function to use
}
  \item{box.width.large.scale}{
%%     ~~Describe \code{box.width.large.scale} here~~
}
  \item{box.width.small.scale}{
%%     ~~Describe \code{box.width.small.scale} here~~
}
  \item{box.show.mean}{
%%     ~~Describe \code{box.show.mean} here~~
}
  \item{box.show.box}{
%%     ~~Describe \code{box.show.box} here~~
}
  \item{box.show.whiskers}{
%%     ~~Describe \code{box.show.whiskers} here~~
}
  \item{\dots}{
	Other arguments to pass to lattice function
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
function (formula, data.frame, show.outlines = FALSE, main = "", 
    x.label = "", div.axis.major = 10, div.axis.minor = 20, log.x = FALSE, 
    colors.plot = c("salmon", "mediumblue", "olivedrab", "cyan", 
        "brown", "darkgreen", "purple"), panel = "panel.tuftebox", 
    box.width.large.scale = 0.4, box.width.small.scale = 0.25, 
    box.show.mean = TRUE, box.show.box = FALSE, box.show.whiskers = FALSE, 
    ...) 
{
    grid.newpage()
    gp1.titlesize = 9
    gp2.titlesize = 7
    titlesize.padding = 3
    main.titlesize = 12
    x.padding = 25
    y.scale = 0.7
    densities.x = densities.y = as.numeric()
    if (!exists("data") | !exists("formula")) {
        stop("Must include all required parameters")
    }
    if (panel != "panel.tuftebox") {
        stop("Only panel.tuftebox is currently supported")
    }
    lpf <- latticeParseFormula(formula, data = data.frame)
    if (length(lpf$condition) != 3) {
        stop("Must provide three conditions")
    }
    x <- lpf$right
    if (log.x == TRUE) {
        x <- log10(x)
    }
    if (any(is.na(x))) {
        stop("NA's not allowed in your data vector")
    }
    gp <- lpf$condition
    for (i in seq(length(gp))) {
        if (!is.factor(gp[[i]])) {
            stop("All grouping variables must coerce to factors")
        }
    }
    levels.gp1 <- levels(gp[[1]])
    levels.gp2 <- levels(gp[[2]])
    levels.gp3 <- levels(gp[[3]])
    num.gp1 <- length(levels.gp1)
    num.gp2 <- length(levels.gp2)
    num.gp3 <- length(levels.gp3)
    pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(0.99 * 
        c(0.05, 0.8, 0.05, 0.1), "npc"), widths = unit(0.95 * 
        c(0.1, 0.9), "npc")), name = "Main"))
    seekViewport("Main")
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1, 
        name = "Title"))
    if (show.outlines) {
        grid.rect()
    }
    seekViewport("Main")
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2, 
        name = "Graphs", layout = grid.layout(3, num.gp1, heights = unit(c(gp1.titlesize + 
            titlesize.padding, x.padding, 1), c("points", "points", 
            "null")), widths = unit(1/num.gp1, "npc"))))
    if (show.outlines) {
        grid.rect()
    }
    seekViewport("Main")
    legend.n.col = ifelse(num.gp3 > 4, ceiling(num.gp3/2), num.gp3)
    legend.n.row = ifelse(num.gp3 > 4, 2, 1)
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 4, 
        name = "Legend", layout = grid.layout(legend.n.row, legend.n.col), 
        width = unit(0.95, "npc"), height = unit(0.95, "npc")))
    if (show.outlines) {
        grid.rect()
    }
    seekViewport("Main")
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2, 
        name = "Axis", layout = grid.layout(5, 1, heights = unit(c(gp1.titlesize + 
            titlesize.padding, x.padding, 1, x.padding, gp2.titlesize + 
            titlesize.padding), c("points", "points", "null", 
            "points", "points")))))
    if (show.outlines) {
        grid.rect()
    }
    x.range <- range(x)
    y.range <- c(0, 1)
    y.range.scaled <- y.range + c(diff(y.range) * (1 - y.scale)/2, 
        -diff(y.range) * (1 - y.scale)/2)
    for (gp1.i in seq(num.gp1)) {
        seekViewport("Graphs")
        pushViewport(viewport(layout.pos.col = gp1.i, layout.pos.row = 1, 
            name = paste("gp1.", gp1.i, ".title", sep = "")))
        grid.text(label = as.character(levels.gp1[gp1.i]), gp = gpar(fontsize = gp1.titlesize))
        grid.rect()
        seekViewport("Graphs")
        pushViewport(viewport(layout.pos.col = gp1.i, layout.pos.row = 3, 
            name = paste("gp1.", gp1.i, sep = ""), layout = grid.layout(3, 
                num.gp2, heights = unit(c(1, x.padding, gp2.titlesize + 
                  titlesize.padding), c("null", "points", "points")), 
                widths = unit(1/num.gp2, "npc"))))
        for (gp2.i in seq(num.gp2)) {
            seekViewport(paste("gp1.", gp1.i, sep = ""))
            pushViewport(viewport(layout.pos.col = gp2.i, layout.pos.row = 3, 
                name = paste("gp1.", gp1.i, "_gp2.", gp2.i, ".title", 
                  sep = "")))
            pushViewport(viewport(angle = -90))
            grid.text(label = as.character(levels.gp2[gp2.i]), 
                gp = gpar(fontsize = gp2.titlesize))
            if (show.outlines) {
                grid.rect()
            }
            seekViewport(paste("gp1.", gp1.i, sep = ""))
            pushViewport(viewport(layout.pos.col = gp2.i, layout.pos.row = 1, 
                name = paste("gp1.", gp1.i, "_gp2.", gp2.i, sep = ""), 
                yscale = unit(x.range, "native")))
            if (show.outlines) {
                grid.rect()
            }
            for (gp3.i in seq(num.gp3)) {
                x.gp1gp2gp3 <- subset(x, as.numeric(gp[[1]]) == 
                  gp1.i & as.numeric(gp[[2]]) == gp2.i & as.numeric(gp[[3]]) == 
                  gp3.i)
                if (length(x.gp1gp2gp3 > 0)) {
                  if (panel == "panel.tuftebox") {
                    loc.y = y.range.scaled[1] + (gp3.i - 0.5) * 
                      (1/num.gp3) * diff(y.range.scaled)
                    quantiles = quantile(x.gp1gp2gp3)
                    iqr = diff(quantiles[c("25\%", "75\%")])
                    box.width.tiny = unit(1, "points")
                    box.width.small = box.width.small.scale * 
                      (1/num.gp3) * diff(y.range.scaled)
                    box.width.large = box.width.large.scale * 
                      (1/num.gp3) * diff(y.range.scaled)
                    min.reduced = max(quantiles["25\%"] - 1.5 * 
                      iqr, min(x.gp1gp2gp3))
                    max.reduced = min(quantiles["75\%"] + 1.5 * 
                      iqr, max(x.gp1gp2gp3))
                    grid.lines(y = unit(c(min.reduced, quantiles["25\%"]), 
                      "native"), x = loc.y, default.unit = "native", 
                      gp = gpar(col = colors.plot[gp3.i]))
                    grid.lines(y = unit(c(max.reduced, quantiles["75\%"]), 
                      "native"), x = loc.y, default.unit = "native", 
                      gp = gpar(col = colors.plot[gp3.i]))
                    if (box.show.whiskers == TRUE) {
                      grid.lines(y = unit(min.reduced, "native"), 
                        x = unit(loc.y, "native") + (c(1, -1) * 
                          box.width.tiny), default.unit = "native", 
                        gp = gpar(col = colors.plot[gp3.i]))
                      grid.lines(y = unit(max.reduced, "native"), 
                        x = unit(loc.y, "native") + (c(1, -1) * 
                          box.width.tiny), default.unit = "native", 
                        gp = gpar(col = colors.plot[gp3.i]))
                    }
                    if (box.show.mean == FALSE) {
                      grid.lines(y = unit(quantiles[c("25\%", 
                        "75\%")], "native"), x = unit(loc.y, "native") - 
                        box.width.tiny, default.unit = "native", 
                        gp = gpar(col = colors.plot[gp3.i]))
                      if (box.show.box == TRUE) {
                        grid.lines(y = unit(quantiles[c("25\%", 
                          "75\%")], "native"), x = unit(loc.y, 
                          "native") + box.width.tiny, default.unit = "native", 
                          gp = gpar(col = colors.plot[gp3.i]))
                        grid.lines(y = unit(quantiles["25\%"], 
                          "native"), x = unit(loc.y, "native") + 
                          (c(1, -1) * box.width.tiny), default.unit = "native", 
                          gp = gpar(col = colors.plot[gp3.i]))
                        grid.lines(y = unit(quantiles["75\%"], 
                          "native"), x = unit(loc.y, "native") + 
                          (c(1, -1) * box.width.tiny), default.unit = "native", 
                          gp = gpar(col = colors.plot[gp3.i]))
                      }
                      else {
                        grid.lines(y = unit(quantiles["25\%"], 
                          "native"), x = unit(loc.y, "native") + 
                          (c(0, -1) * box.width.tiny), default.unit = "native", 
                          gp = gpar(col = colors.plot[gp3.i]))
                        grid.lines(y = unit(quantiles["75\%"], 
                          "native"), x = unit(loc.y, "native") + 
                          (c(0, -1) * box.width.tiny), default.unit = "native", 
                          gp = gpar(col = colors.plot[gp3.i]))
                      }
                    }
                    outliers = subset(x.gp1gp2gp3, x.gp1gp2gp3 < 
                      min.reduced | x.gp1gp2gp3 > max.reduced)
                    if (length(outliers) > 0) {
                      grid.points(y = unit(outliers, "native"), 
                        x = rep(loc.y, length(outliers)), default.unit = "native", 
                        gp = gpar(col = colors.plot[gp3.i], cex = 0.2), 
                        pch = 4)
                    }
                    if (box.show.mean == TRUE) {
                      grid.lines(y = unit(rep(quantiles[c("25\%")], 
                        2), "native"), x = c(loc.y - box.width.small, 
                        loc.y + box.width.small), default.unit = "native", 
                        gp = gpar(col = colors.plot[gp3.i]))
                      grid.lines(y = unit(rep(quantiles[c("75\%")], 
                        2), "native"), x = c(loc.y - box.width.small, 
                        loc.y + box.width.small), default.unit = "native", 
                        gp = gpar(col = colors.plot[gp3.i]))
                    }
                    grid.points(y = unit(median(x.gp1gp2gp3), 
                      "native"), x = loc.y, default.unit = "native", 
                      gp = gpar(col = colors.plot[gp3.i], cex = 0.3), 
                      pch = 15)
                    if (box.show.mean == TRUE) {
                      meanlines.x = c(mean(x.gp1gp2gp3), mean(x.gp1gp2gp3) - 
                        sd(x.gp1gp2gp3), mean(x.gp1gp2gp3), mean(x.gp1gp2gp3) + 
                        sd(x.gp1gp2gp3), mean(x.gp1gp2gp3))
                      meanlines.y = c(loc.y - box.width.large, 
                        loc.y, loc.y + box.width.large, loc.y, 
                        loc.y - box.width.large)
                      grid.lines(y = meanlines.x, x = meanlines.y, 
                        default.unit = "native", gp = gpar(col = colors.plot[gp3.i]))
                    }
                  }
                }
            }
        }
    }
    seekViewport("Title")
    grid.text(label = main, gp = gpar(fontsize = main.titlesize))
    seekViewport("Axis")
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3, 
        name = "Axis.actual", yscale = unit(x.range, "native")))
    if (show.outlines) {
        grid.rect()
    }
    if (show.outlines) {
        grid.rect()
    }
    x.range.magnitude <- diff(x.range)
    x.seq = x.range[1] + (1/div.axis.major) * seq(0, div.axis.major) * 
        x.range.magnitude
    mat.x <- rbind(rep(0.85, div.axis.major + 1), rep(1, div.axis.major + 
        1))
    mat.y <- matrix(rep(x.seq, each = 2), nrow = 2)
    grid.polyline(x = mat.x, y = mat.y, id.lengths = rep(2, div.axis.major + 
        1), default.unit = "native")
    if (log.x == FALSE) {
        round.digits = -floor(log10(x.range.magnitude)) + 1
        x.labels <- round(x.seq, round.digits)
    }
    else {
        x.seq.label = 10^x.seq
        x.labels = round.sigfig(10^x.seq, 1)
    }
    grid.text(label = as.character(x.labels), x = 0.8, y = x.seq, 
        gp = gpar(fontsize = 9), just = c("right", "center"), 
        default.unit = "native")
    mat.x <- rbind(rep(0.925, div.axis.minor + 1), rep(1, div.axis.minor + 
        1))
    mat.y <- matrix(rep(x.range[1] + x.range.magnitude * (1/div.axis.minor) * 
        seq(0, div.axis.minor), each = 2), nrow = 2)
    grid.polyline(x = mat.x, y = mat.y, id.lengths = rep(2, div.axis.minor + 
        1), default.unit = "native")
    grid.lines(x = c(1, 1), y = c(0, 1))
    upViewport(0)
    pushViewport(viewport(angle = 90, x = unit(0.01, "npc"), 
        width = convertUnit(unit(1, "npc"), "npc", "y", "dimension", 
            "x", "dimension"), height = convertUnit(unit(0.3, 
            "npc"), "npc", "x", "dimension", "y", "dimension")))
    grid.text(label = x.label, just = c("center", "top"), gp = gpar(fontsize = 10))
    seekViewport("Legend")
    for (gp3.i in seq(num.gp3)) {
        gp3.col = ifelse(gp3.i > legend.n.col, gp3.i - legend.n.col, 
            gp3.i)
        gp3.row = ifelse(gp3.i > legend.n.col, 2, 1)
        pushViewport(viewport(layout.pos.col = gp3.col, layout.pos.row = gp3.row, 
            width = unit(0.9, "npc"), height = unit(0.9, "npc")))
        legend.colorbox.width = convertUnit(unit(0.3, "npc"), 
            "npc", "y", "dimension", "x", "dimension")
        grid.text(x = unit(0.15, "npc") + legend.colorbox.width, 
            hjust = 0, label = levels.gp3[gp3.i], gp = gpar(col = colors.plot[gp3.i], 
                fontface = "bold"))
        grid.rect(x = unit(0.1, "npc"), y = unit(0.5, "npc"), 
            hjust = 0, height = unit(0.3, "npc"), width = legend.colorbox.width, 
            gp = gpar(fill = colors.plot[gp3.i]))
        grid.rect(x = unit(0.1, "npc"), y = unit(0.5, "npc"), 
            hjust = 0, height = unit(0.3, "npc"), width = legend.colorbox.width)
        popViewport()
    }
  }
}
