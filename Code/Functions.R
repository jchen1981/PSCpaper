require(matrixStats)
require(qvalue)
require(CpGassoc)
require(vegan)
require(ggplot2)
require(ade4)
require(reshape2)
require(RColorBrewer)
require(gplots)
#################################################################################
B2M <- function(x) {
	x[x==0] <- min(x[x!=0])
	x[x==1] <- max(x[x!=1])
	log(x / (1-x))
}

#################################################################################
fastDist <- function(X) {
	temp <- colSums(X^2)  
	D <- outer(temp, temp, "+") - 2 * t(X) %*% X
	diag(D) <- 0
	sqrt(D)
}

#################################################################################
fastLM <- function(Y, M) {
	Y <- as.matrix(Y)
	XXI <- solve(t(M) %*% M)
	dof <- ncol(Y) - ncol(M)
	est <- XXI %*% t(M) %*% t(Y)
	resid <- t(Y) - M %*% est
	sigma <- sqrt(colSums(resid^2)/dof)
	Pvec <- 2*pt(-abs(t(est/(sqrt(diag(XXI))))/sigma), dof)
	return(Pvec)
}

#################################################################################

vkey2 <- function (map, title = NA, side = 2, stretch = 1.4, x, y, 
		wh)
{
	opar <- par(xpd = NA)
	on.exit(par(opar))
	n <- length(map$breaks) + 1
	dy <- strheight("A")
	aspect <- diff(grconvertX(1:2, from = "inches"))/diff(grconvertY(1:2,
					from = "inches"))
	dx <- dy * aspect
	if (missing(wh)) {
		
		wh <- 1:(n-1)
		
	}
	labs <- format(map$breaks[wh])
	maxlabwidth <- max(strwidth(labs))
	if (missing(x)) {
		x <- grconvertX(1, from = "nfc") - (2 * dx)
		if (side == 4)
			x <- x - maxlabwidth - dx
	}
	else {
		if (is.list(x)) {
			y <- x$y
			x <- x$x
		}
	}
	if (missing(y))
		y <- par("usr")[3] + dy
	ybord <- y + ((0:(n - 1)) * dy * stretch)
	rect(x, ybord[-n], x + dx, ybord[-1], col = map$colors, border = NA)
	if (side == 4) {
		xtext <- x + dx
		text(x = x, y = ybord[n] + (1.5 * dy), title, adj = c(0,
						0))
	}
	if (side == 2) {
		xtext <- x
		text(x = x + dx, y = ybord[n] + (1.5 * dy), title, adj = c(1,
						0))
	}
	text(x = xtext, y = ybord[wh] + 0.5 * dy, labels = labs, pos = side)
}


heatmap.3 <- function(x,
		Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
		distfun = dist,
		hclustfun = hclust,
		dendrogram = c("both","row", "column", "none"),
		symm = FALSE,
		scale = c("none","row", "column"),
		na.rm = TRUE,
		revC = identical(Colv,"Rowv"),
		add.expr,
		breaks,
		symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
		col = "heat.colors",
		colsep,
		rowsep,
		sepcolor = "white",
		sepwidth = c(0.05, 0.05),
		cellnote,
		notecex = 1,
		notecol = "cyan",
		na.color = par("bg"),
		trace = c("none", "column","row", "both"),
		tracecol = "cyan",
		hline = median(breaks),
		vline = median(breaks),
		linecol = tracecol,
		margins = c(5,5),
		ColSideColors = NULL,
		RowSideColors = NULL,
		side.height.fraction=0.3,
		cexRow = 0.2 + 1/log10(nr),
		cexCol = 0.2 + 1/log10(nc),
		labRow = NULL,
		labCol = NULL,
		key = TRUE,
		keysize = 1.5,
		density.info = c("none", "histogram", "density"),
		denscol = tracecol,
		symkey = max(x < 0, na.rm = TRUE) || symbreaks,
		densadj = 0.25,
		main = NULL,
		xlab = NULL,
		ylab = NULL,
		lmat = NULL,
		lhei = NULL,
		lwid = NULL,
		NumColSideColors = 1,
		NumRowSideColors = 1,
		KeyValueName="Value",...){
	
	invalid <- function (x) {
		if (missing(x) || is.null(x) || length(x) == 0)
			return(TRUE)
		if (is.list(x))
			return(all(sapply(x, invalid)))
		else if (is.vector(x))
			return(all(is.na(x)))
		else return(FALSE)
	}
	
	x <- as.matrix(x)
	scale01 <- function(x, low = min(x), high = max(x)) {
		x <- (x - low)/(high - low)
		x
	}
	retval <- list()
	scale <- if (symm && missing(scale))
				"none"
			else match.arg(scale)
	dendrogram <- match.arg(dendrogram)
	trace <- match.arg(trace)
	density.info <- match.arg(density.info)
	if (length(col) == 1 && is.character(col))
		col <- get(col, mode = "function")
	if (!missing(breaks) && (scale != "none"))
		warning("Using scale=\"row\" or scale=\"column\" when breaks are",
				"specified can produce unpredictable results.", "Please consider using only one or the other.")
	if (is.null(Rowv) || is.na(Rowv))
		Rowv <- FALSE
	if (is.null(Colv) || is.na(Colv))
		Colv <- FALSE
	else if (Colv == "Rowv" && !isTRUE(Rowv))
		Colv <- FALSE
	if (length(di <- dim(x)) != 2 || !is.numeric(x))
		stop("`x' must be a numeric matrix")
	nr <- di[1]
	nc <- di[2]
	if (nr <= 1 || nc <= 1)
		stop("`x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2)
		stop("`margins' must be a numeric vector of length 2")
	if (missing(cellnote))
		cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
	if (!inherits(Rowv, "dendrogram")) {
		if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
					c("both", "row"))) {
			if (is.logical(Colv) && (Colv))
				dendrogram <- "column"
			else dedrogram <- "none"
			warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
					dendrogram, "'. Omitting row dendogram.")
		}
	}
	if (!inherits(Colv, "dendrogram")) {
		if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
					c("both", "column"))) {
			if (is.logical(Rowv) && (Rowv))
				dendrogram <- "row"
			else dendrogram <- "none"
			warning("Discrepancy: Colv is FALSE, while dendrogram is `",
					dendrogram, "'. Omitting column dendogram.")
		}
	}
	if (inherits(Rowv, "dendrogram")) {
		ddr <- Rowv
		rowInd <- order.dendrogram(ddr)
	}
	else if (is.integer(Rowv)) {
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd))
			stop("row dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Rowv)) {
		Rowv <- rowMeans(x, na.rm = na.rm)
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd))
			stop("row dendrogram ordering gave index of wrong length")
	}
	else {
		rowInd <- nr:1
	}
	if (inherits(Colv, "dendrogram")) {
		ddc <- Colv
		colInd <- order.dendrogram(ddc)
	}
	else if (identical(Colv, "Rowv")) {
		if (nr != nc)
			stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
		if (exists("ddr")) {
			ddc <- ddr
			colInd <- order.dendrogram(ddc)
		}
		else colInd <- rowInd
	}
	else if (is.integer(Colv)) {
		hcc <- hclustfun(distfun(if (symm)
									x
								else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd))
			stop("column dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Colv)) {
		Colv <- colMeans(x, na.rm = na.rm)
		hcc <- hclustfun(distfun(if (symm)
									x
								else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd))
			stop("column dendrogram ordering gave index of wrong length")
	}
	else {
		colInd <- 1:nc
	}
	retval$rowInd <- rowInd
	retval$colInd <- colInd
	retval$call <- match.call()
	x <- x[rowInd, colInd]
	x.unscaled <- x
	cellnote <- cellnote[rowInd, colInd]
	if (is.null(labRow))
		labRow <- if (is.null(rownames(x)))
					(1:nr)[rowInd]
				else rownames(x)
	else labRow <- labRow[rowInd]
	if (is.null(labCol))
		labCol <- if (is.null(colnames(x)))
					(1:nc)[colInd]
				else colnames(x)
	else labCol <- labCol[colInd]
	if (scale == "row") {
		retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
		x <- sweep(x, 1, rm)
		retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
		x <- sweep(x, 1, sx, "/")
	}
	else if (scale == "column") {
		retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
		x <- sweep(x, 2, rm)
		retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
		x <- sweep(x, 2, sx, "/")
	}
	if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
		if (missing(col) || is.function(col))
			breaks <- 16
		else breaks <- length(col) + 1
	}
	if (length(breaks) == 1) {
		if (!symbreaks)
			breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
					length = breaks)
		else {
			extreme <- max(abs(x), na.rm = TRUE)
			breaks <- seq(-extreme, extreme, length = breaks)
		}
	}
	nbr <- length(breaks)
	ncol <- length(breaks) - 1
	if (class(col) == "function")
		col <- col(ncol)
	min.breaks <- min(breaks)
	max.breaks <- max(breaks)
	x[x < min.breaks] <- min.breaks
	x[x > max.breaks] <- max.breaks
	if (missing(lhei) || is.null(lhei))
		lhei <- c(keysize, 4)
	if (missing(lwid) || is.null(lwid))
		lwid <- c(keysize, 4)
	if (missing(lmat) || is.null(lmat)) {
		lmat <- rbind(4:3, 2:1)
		
		if (!is.null(ColSideColors)) {
			#if (!is.matrix(ColSideColors))
			#stop("'ColSideColors' must be a matrix")
			if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
				stop("'ColSideColors' must be a matrix of nrow(x) rows")
			lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
			#lhei <- c(lhei[1], 0.2, lhei[2])
			lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
		}
		
		if (!is.null(RowSideColors)) {
			#if (!is.matrix(RowSideColors))
			#stop("'RowSideColors' must be a matrix")
			if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
				stop("'RowSideColors' must be a matrix of ncol(x) columns")
			lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
			#lwid <- c(lwid[1], 0.2, lwid[2])
			lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
		}
		lmat[is.na(lmat)] <- 0
	}
	
	if (length(lhei) != nrow(lmat))
		stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
	if (length(lwid) != ncol(lmat))
		stop("lwid must have length = ncol(lmat) =", ncol(lmat))
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	
	if (!is.null(RowSideColors)) {
		if (!is.matrix(RowSideColors)){
			par(mar = c(margins[1], 0, 0, 0.5))
			image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
			box()
		} else {
			par(mar = c(margins[1], 0, 0, 0.5))
			rsc = t(RowSideColors[,rowInd, drop=F])
			rsc.colors = matrix()
			rsc.names = names(table(rsc))
			rsc.i = 1
			for (rsc.name in rsc.names) {
				rsc.colors[rsc.i] = rsc.name
				rsc[rsc == rsc.name] = rsc.i
				rsc.i = rsc.i + 1
			}
			rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
			image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
			
			if (length(rownames(RowSideColors)) > 0) {
				axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
			}
			box()
		}
	}
	
	if (!is.null(ColSideColors)) {
		
		if (!is.matrix(ColSideColors)){
			par(mar = c(0.5, 0, 0, margins[2]))
			image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
			box()
		} else {
			par(mar = c(0.5, 0, 0, margins[2]))
			csc = ColSideColors[colInd, , drop=F]
			csc.colors = matrix()
			csc.names = names(table(csc))
			csc.i = 1
			for (csc.name in csc.names) {
				csc.colors[csc.i] = csc.name
				csc[csc == csc.name] = csc.i
				csc.i = csc.i + 1
			}
			csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
			image(csc, col = as.vector(csc.colors), axes = FALSE)
			if (length(colnames(ColSideColors)) > 0) {
				axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
			}
			box()
		}
	}
	
	par(mar = c(margins[1], 0, 0, margins[2]))
	x <- t(x)
	cellnote <- t(cellnote)
	if (revC) {
		iy <- nr:1
		if (exists("ddr"))
			ddr <- rev(ddr)
		x <- x[, iy]
		cellnote <- cellnote[, iy]
	}
	else iy <- 1:nr
	image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
	retval$carpet <- x
	if (exists("ddr"))
		retval$rowDendrogram <- ddr
	if (exists("ddc"))
		retval$colDendrogram <- ddc
	retval$breaks <- breaks
	retval$col <- col
	if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
		mmat <- ifelse(is.na(x), 1, NA)
		image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
				col = na.color, add = TRUE)
	}
	axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
			cex.axis = cexCol)
	if (!is.null(xlab))
		mtext(xlab, side = 1, line = margins[1] - 1.25)
	axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
			cex.axis = cexRow)
	if (!is.null(ylab))
		mtext(ylab, side = 4, line = margins[2] - 1.25)
	if (!missing(add.expr))
		eval(substitute(add.expr))
	if (!missing(colsep))
		for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
	if (!missing(rowsep))
		for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
	min.scale <- min(breaks)
	max.scale <- max(breaks)
	x.scaled <- scale01(t(x), min.scale, max.scale)
	if (trace %in% c("both", "column")) {
		retval$vline <- vline
		vline.vals <- scale01(vline, min.scale, max.scale)
		for (i in colInd) {
			if (!is.null(vline)) {
				abline(v = i - 0.5 + vline.vals, col = linecol,
						lty = 2)
			}
			xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
			xv <- c(xv[1], xv)
			yv <- 1:length(xv) - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (trace %in% c("both", "row")) {
		retval$hline <- hline
		hline.vals <- scale01(hline, min.scale, max.scale)
		for (i in rowInd) {
			if (!is.null(hline)) {
				abline(h = i + hline, col = linecol, lty = 2)
			}
			yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
			yv <- rev(c(yv[1], yv))
			xv <- length(yv):1 - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (!missing(cellnote))
		text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
				col = notecol, cex = notecex)
	par(mar = c(margins[1], 0, 0, 0))
	if (dendrogram %in% c("both", "row")) {
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	}
	else plot.new()
	par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
	if (dendrogram %in% c("both", "column")) {
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	}
	else plot.new()
	if (!is.null(main))
		title(main, cex.main = 1.5 * op[["cex.main"]])
	if (key) {
		par(mar = c(5, 4, 2, 1), cex = 0.75)
		tmpbreaks <- breaks
		if (symkey) {
			max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
			min.raw <- -max.raw
			tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
			tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
		}
		else {
			min.raw <- min(x, na.rm = TRUE)
			max.raw <- max(x, na.rm = TRUE)
		}
		
		z <- seq(min.raw, max.raw, length = length(col))
		image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
				xaxt = "n", yaxt = "n")
		par(usr = c(0, 1, 0, 1))
		lv <- pretty(breaks)
		xv <- scale01(as.numeric(lv), min.raw, max.raw)
		axis(1, at = xv, labels = lv)
		if (scale == "row")
			mtext(side = 1, "Row Z-Score", line = 2)
		else if (scale == "column")
			mtext(side = 1, "Column Z-Score", line = 2)
		else mtext(side = 1, KeyValueName, line = 2)
		if (density.info == "density") {
			dens <- density(x, adjust = densadj, na.rm = TRUE)
			omit <- dens$x < min(breaks) | dens$x > max(breaks)
			dens$x <- dens$x[-omit]
			dens$y <- dens$y[-omit]
			dens$x <- scale01(dens$x, min.raw, max.raw)
			lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
					lwd = 1)
			axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
			title("Color Key\nand Density Plot")
			par(cex = 0.5)
			mtext(side = 2, "Density", line = 2)
		}
		else if (density.info == "histogram") {
			h <- hist(x, plot = FALSE, breaks = breaks)
			hx <- scale01(breaks, min.raw, max.raw)
			hy <- c(h$counts, h$counts[length(h$counts)])
			lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
					col = denscol)
			axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
			title("Color Key\nand Histogram")
			par(cex = 0.5)
			mtext(side = 2, "Count", line = 2)
		}
		else title("Color Key")
	}
	else plot.new()
	retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
			high = retval$breaks[-1], color = retval$col)
	invisible(retval)
}

generate_heatmap <- function (data.obj,  meta.info, sam.ord=NULL, data.type='P',
		row.col.dat=NULL,  sepwidth=0.01, colsep=NULL, rowsep=NULL, taxa.as.row = TRUE,
		colFnsC=NULL, colFnsF=NULL, Rowv=T, Colv=T, dendrogram='both', margins=c(5, 15),
		in.grid=F, sepcolor='black', is.labCol=T, cexCol=1, cexRow=NULL,
		omas=c(1, 1, 1, 8), width=12, height=6, file.name='Heatmap', return.obj=FALSE, pdf=TRUE, ...) {
	
	colsep0 <- colsep
	rowsep0 <- rowsep
	df <- data.obj$meta.dat
	
	# Determine the col/rowside color
	# Determine the col/rowside color
	if (is.null(colFnsC)) {
		# colFnsC <- c(colorRampPalette(c('black', 'yellow', 'red'), colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue'))))
		# https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
		colFunsC <- c(colorRampPalette(colors = brewer.pal(9,"RdBu")), 
				colorRampPalette(colors = brewer.pal(9,"PRGn")), colorRampPalette(colors = brewer.pal(9,"RdYlBu")),  
				colorRampPalette(colors = brewer.pal(9,"RdGy")), colorRampPalette(colors = brewer.pal(9,"PuOr")))
	} 
	if (is.null(colFnsF)) {
		
		colFnsF <- function (x) {
			if (x <= 6) {
				return(brewer.pal(7, "Set2")[1:x])
			} else {
				return(colorRampPalette(colors = jet(8))(x))
			}
		}
		colFnsF <- c(colFnsF)
	}
	prop <- data.obj$data
	
	# Sort row and column
	if (!is.null(sam.ord)) {
		prop <- prop[, sam.ord, drop=FALSE]
	}
	
	if (is.labCol) {
		labCol <- colnames(prop)
	}  else {
		labCol <- ''
	}
	if (data.type == 'LogC') {
		col.scheme <- c('white', colorpanel(9, 'aliceblue', 'lightblue', 'blue'))
		breaks <- c(-0.1, 0.1, 2, 4, 6, 8, 10, 12, 14, 16, 18)
	}
	# Deal with zeros
	if (data.type == 'B') {
		col.scheme <- c("lightyellow", "red")
		breaks <- c(-0.01, 0.01, 1.1)
	} 
	if (data.type == 'P'){
		#
		#col.scheme = brewer.pal(11, "Spectral")
		col.scheme <- bluered(11)
		breaks <- c(seq(min(prop), max(prop), len=12))
	}
	if (data.type == 'R'){
		col.scheme <- c(white, colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(ncol(prop)-1))
		breaks <- seq(0, ncol(prop), len=ncol(prop)+1)
	}
	
	if (nrow(prop) > 1500) {
		labRow <- ''
		cexRow <- 1
	} else {
		labRow <- rownames(prop)
		if (is.null(cexRow)) {
			cexRow <- ifelse(0.5 * 60 / nrow(prop) > 1, 1, 0.5 * 60 / nrow(prop))
		}
		
	}
	
	if (is.null(colsep0)) {
		if (in.grid == T) {
			colsep <- c(0:ncol(prop))  
		} else {
			colsep <- c(0, ncol(prop))  
		}
	} else {
		colsep <- colsep0
	}
	
	if (is.null(rowsep0)) {
		if (in.grid == T) {
			rowsep <- c(0:nrow(prop))
		} else {
			rowsep <- c(0, nrow(prop))
		}
	} else {
		rowsep <- rowsep0
	}
	
	key.list <- list()
	colsidecol <- NULL
	
	i <- 0
	j <- 0
	
	for (keyID in meta.info) {
		if (is.null(sam.ord)) {
			x <- df[, keyID]
		} else {
			x <- df[sam.ord, keyID]
		}		
		if (is.factor(x)) {
			key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)	
			i <- (i + 1) %% length(colFnsF)
			colsidecol <- cbind(colsidecol, key.list[[keyID]]$colors[x])
		} else {
			key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
			j <- (j + 1) %% length(colFnsC)
			colsidecol <- cbind(colsidecol, cmap(x, key.list[[keyID]]))
		}
	}
	colnames(colsidecol) <- meta.info
	
	# Rev: 2016_09_13
	colsidecol[is.na(colsidecol)] <- 'white'
	
	# add aunbdance key 
	if (data.type == 'B') {
		prop.cmap <- list(breaks=c("0", "1"), colors=c("lightyellow", "red"), base=NA, col.na=NA, right=F, include.lowest=F) 
		KeyName <- 'Value'
	} 
	
	if (data.type == 'P'){
		KeyName <- 'Value'	
	}
	if (data.type == 'R'){
		KeyName <- 'Rank'
	}
	
	if (data.type == 'LogC'){
		KeyName <- 'Log2(count+1)'
	}
	
	# add row col key
	if (!is.null(row.col.dat)) {
		
		phy <- df[, row.col.dat]
		phy <- factor(phy)
		
		rowsidecol<- rainbow(nlevels(phy))[phy]
		rowsidecol <- rbind(rowsidecol, rowsidecol)
		rownames(rowsidecol) <- c('', '')
		phy.cmap <- list(breaks=levels(phy), colors=rainbow(nlevels(phy)), base=NA, col.na=NA, right=F, include.lowest=F) 
	} else {
		rowsidecol <- NULL
	}
	if(pdf == TRUE) {
		pdf(paste0(file.name, '.pdf'), width=width, height=height)
	}
	
	par(oma = omas)
	
	# Pearson correlation distance
	if (taxa.as.row) {
		prop <- prop
	} else {
		prop <- t(prop)
	}
	obj <- heatmap.3(prop, 
			Rowv=Rowv, 
			Colv=Colv, 
#				distfun = dist2,
			dendrogram=dendrogram,
			scale='none',
			col=col.scheme, 
			breaks=breaks, 
			symbreaks=F,
			trace='none',
			margins= margins, 
			colsep = colsep,  
			rowsep = rowsep,  
			sepcolor= sepcolor, 
			sepwidth=c(sepwidth, sepwidth),
			ColSideColors=colsidecol,
			RowSideColors=rowsidecol,
			cexRow=cexRow,
			labRow=labRow,
			labCol=labCol,
			cexCol=cexCol,
			key=(data.type != 'B'), density.info='none', symkey=F, KeyValueName=KeyName,	
			NumColSideColors= 0.5 *length(meta.info),
			NumRowSideColors= 0.5,
			...
	)
	
	par(cex=0.75)
	par(oma=c(0, 0, 1, 0)) 
	
	if (!is.null(rowsidecol)) {
		y.cord <- (1/(length(meta.info)+2)) * (0:((length(meta.info)+2) - 1))
		vkey2(phy.cmap, rowsidecol, x=0, y=-0.2, stretch=1.2)
		
	}
	
	y.cord <- (1/(length(meta.info))) * (0:(length(meta.info) - 1))
	k <- 1
	xs <- c(0.95, 1)
	for (keyID in meta.info) {
		x <- df[, keyID]
		if (is.factor(x)) {
			vkey2(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
		} else {
			vkey(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
		}
		k <- k + 1
	}	
	
	if (data.type == 'B') {
		vkey2(prop.cmap, KeyName, x=xs[(k-1) %% 2 + 1], y=1, stretch=1.2)		
	} 
	
	if (pdf == TRUE) {
		dev.off()	
	}
	
	if (return.obj == TRUE) {
		return(obj)
	}
}
