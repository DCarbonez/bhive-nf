RCircos.Reset.Plot.Parameters = function(new.params) {
    point.chr <- new.params$point.type
    text.color <- new.params$text.color
    heatmap.color <- new.params$heatmap.color
    hist.color <- new.params$hist.color
    line.color <- new.params$line.color
    scatter.color <- new.params$scatter.color
    tile.color <- new.params$tile.color
    bg.color <- new.params$track.background
    grid.color <- new.params$grid.line.color
    params <- unlist(new.params)
    params <- params[-which(params == point.chr)]
    params <- params[-which(params == text.color)]
    params <- params[-which(params == heatmap.color)]
    params <- params[-which(params == hist.color)]
    params <- params[-which(params == line.color)]
    params <- params[-which(params == scatter.color)]
    params <- params[-which(params == tile.color)]
    params <- params[-which(params == bg.color)]
    params <- params[-which(params == grid.color)]
    params <- as.numeric(params)
    if (sum(is.na(params)) > 0) {
        stop("Plot parameters except of point.type must be numeric.")
    }
    if (sum(params < 0) > 0) {
        stop("Plot parameters cannot have negative values")
    }
    old.params <- RCircos.Get.Plot.Parameters()
    cyto.band.data <- RCircos.Get.Plot.Ideogram()
    cyto.band.data$Unit <- NULL
    cyto.band.data$Location <- NULL
    # Overruling conditional overrule...
#    if (new.params$chrom.paddings != 0) {
#        padding.const <- 3000/1e+06
#        genome.lenth <- sum(cyto.band.data$Length)
#        total.units <- genome.lenth/new.params$base.per.unit
#        the.padding <- round(padding.const * total.units, digits = 0)
#        if (new.params$chrom.paddings > the.padding) {
#            cat(paste("\nNote: chrom.padding", new.params$chrom.paddings, 
#                "is too big,", "and was reset to", the.padding, 
#                "\n"))
#            new.params$chrom.paddings <- the.padding
#        }
#    }
    if (new.params$radius.len < 1) {
        cat("\nNote: radius.len needs be at least 1.0\n\n")
        new.params$radius.len <- 1
    }
    new.params$chr.ideog.pos <- new.params$radius.len + 0.1
    new.params$highlight.pos <- new.params$radius.len + 0.25
    new.params$chr.name.pos <- new.params$radius.len + 0.4
    new.params$highlight.width = round(new.params$radius.len, 
        digits = 0)
    # Why erase user's parameters and replace them by default?
    # There are things I simply cannot understand...
#    new.params$track.in.start <- new.params$radius.len + 0.05
    new.params$track.out.start <- new.params$radius.len + 0.5
    differ <- old.params$plot.radius - old.params$radius.len
    new.params$plot.radius <- new.params$radius.len + differ
    RCircosEnvironment <- NULL
    RCircosEnvironment <- get("RCircos.Env", envir = globalenv())
    RCircosEnvironment[["RCircos.PlotPar"]] <- NULL
    RCircosEnvironment[["RCircos.Cytoband"]] <- NULL
    RCircosEnvironment[["RCircos.Base.Position"]] <- NULL
    RCircosEnvironment[["RCircos.PlotPar"]] <- new.params
    RCircos.Set.Cytoband.data(cyto.band.data)
    RCircos.Set.Base.Plot.Positions()
}



RCircos.Chromosome.Ideogram.Plot = function () {
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    right.side <- nrow(RCircos.Pos)/2
    outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width
    inner.location <- RCircos.Par$chr.ideog.pos
    chroms <- unique(RCircos.Cyto$Chromosome)
    for (a.chr in 1:length(chroms)) {
        the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr],]
        start <- the.chr$Location[1] - the.chr$Unit[1] + 1
        end <- the.chr$Location[nrow(the.chr)]
        mid <- round((end - start + 1)/2, digits = 0) + start
        # This overflow of color is not so cool.
#        chr.color <- the.chr$ChrColor[nrow(the.chr)]
        pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, 
            RCircos.Pos[end:start, 1] * inner.location)
        pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, 
            RCircos.Pos[end:start, 2] * inner.location)
#        polygon(pos.x, pos.y)
        polygon(pos.x, pos.y, col="#EEC90060")
        chr.name <- sub(pattern = "chr", replacement = "", chroms[a.chr])
        text(RCircos.Pos[mid, 1] * RCircos.Par$chr.name.pos, 
            RCircos.Pos[mid, 2] * RCircos.Par$chr.name.pos, label = chr.name)
        # Remove the rotation.
#            srt = RCircos.Pos$degree[mid])
#        lines(RCircos.Pos[start:end, ] * RCircos.Par$highlight.pos, 
#            col = chr.color, lwd = RCircos.Par$highlight.width)
#            col = "grey50", lwd = 2)
    }
    for (a.band in 1:nrow(RCircos.Cyto)) {
        a.color <- RCircos.Cyto$BandColor[a.band]
#        if (a.color == "white") {
#            next
#        }
        if (a.color != "red") next
        a.color = "black"
        start <- RCircos.Cyto$Location[a.band] - RCircos.Cyto$Unit[a.band] + 
            1
        end <- RCircos.Cyto$Location[a.band]
        pos.x <- c(RCircos.Pos[start:end, 1] * outer.location, 
            RCircos.Pos[end:start, 1] * inner.location)
        pos.y <- c(RCircos.Pos[start:end, 2] * outer.location, 
            RCircos.Pos[end:start, 2] * inner.location)
        polygon(pos.x, pos.y, col = a.color, border = NA)
    }
}




Tile.Plot = function(tile.data, track.num, side, col="black") {
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    tile.data <- RCircos.Get.Plot.Data(tile.data, "plot")
    # Pre-process data and arrange tiles in layers.
    the.layer <- 1
    the.chr <- tile.data[1, 1]
    start <- tile.data[1, 2]
    end <- tile.data[1, 3]
    tile.layers <- rep(1, nrow(tile.data))
    # The piece of code commented out below was used to put tiles on
    # different layers in case they overlap. In this version, the
    # tiles are put on the same layer.
#    for (a.row in 2:nrow(tile.data)) {
#        if (tile.data[a.row, 2] >= end) {
#            the.layer <- 1
#            start <- tile.data[a.row, 2]
#            end <- tile.data[a.row, 3]
#        }
#        else if (tile.data[a.row, 1] != the.chr) {
#            the.layer <- 1
#            the.chr <- tile.data[a.row, 1]
#            start <- tile.data[a.row, 2]
#            end <- tile.data[a.row, 3]
#        }
#        else {
#            the.layer <- the.layer + 1
#            if (tile.data[a.row, 3] > end) {
#                end <- tile.data[a.row, 3]
#            }
#        }
#        tile.layers[a.row] <- the.layer
#    }
    locations <- RCircos.Track.Positions(side, track.num)
    out.pos <- locations[1]
    in.pos <- locations[2]
    layer.height <- RCircos.Par$track.height/RCircos.Par$max.layers
    num.layers <- max(tile.layers)
    if (num.layers > RCircos.Par$max.layers) {
        if (side == "in") {
            in.pos <- out.pos - layer.height * num.layers
        }
        else {
            out.pos <- in.pos + layer.height * num.layers
        }
        cat(paste("Tiles plot will use more than one track.", 
            "Please select correct area for next track.\n"))
    }
    if (num.layers < RCircos.Par$max.layers) {
        layer.height <- RCircos.Par$track.height/num.layers
    }
    # tile.colors <- RCircos.Get.Plot.Colors(tile.data, RCircos.Par$tile.color)
    # Do not plot outline, just the tiles.
    # RCircos.Track.Outline(out.pos, in.pos, num.layers)
    the.loc <- ncol(tile.data)
    for (a.row in 1:nrow(tile.data)) {
        tile.len <- tile.data[a.row, 3] - tile.data[a.row, 2]
        tile.range <- round(tile.len/RCircos.Par$base.per.unit/2, 
            digits = 0)
        start <- tile.data[a.row, the.loc] - tile.range
        end <- tile.data[a.row, the.loc] + tile.range
        layer.bot <- in.pos + layer.height * (tile.layers[a.row] - 
            1)
        layer.top <- layer.bot + layer.height * 0.8
        polygon.x <- c(RCircos.Pos[start:end, 1] * layer.top, 
            RCircos.Pos[end:start, 1] * layer.bot)
        polygon.y <- c(RCircos.Pos[start:end, 2] * layer.top, 
            RCircos.Pos[end:start, 2] * layer.bot)
        # polygon(polygon.x, polygon.y, col = tile.colors[a.row])
        polygon(polygon.x, polygon.y, col=col, border=col)
    }
}

