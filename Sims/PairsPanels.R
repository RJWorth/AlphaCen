# Pretty panels plots with histograms, correlation lines and coefficients
# from http://musicroamer.com/blog/2011/01/16/r-tips-and-tricks-modified-pairs-plot/

panel.cor.scale <- function(x, y, digits=2, prefix="", cex.cor)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r = (cor(x, y,use="pairwise"))
txt <- format(c(r, 0.123456789), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
text(0.5, 0.5, txt, cex = cex * abs(r))
}


panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r = (cor(x, y,use="pairwise"))
txt <- format(c(r, 0.123456789), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
text(0.5, 0.5, txt, cex = cex )
}


panel.hist <- function(x, ...)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(usr[1:2], 0, 1.5) )
h <- hist(x, plot = FALSE)
breaks <- h$breaks; nB <- length(breaks)
y <- h$counts; y <- y/max(y)
rect(breaks[-nB], 0, breaks[-1], y, ...)
}


pairs.panels <- function (x,y,pchs=20,cols='black',smooth=FALSE,scale=FALSE)
{if (smooth ){
if (scale) {
pairs(x,diag.panel=panel.hist,lower.panel=panel.cor.scale,upper.panel=panel.smooth)
}
else {pairs(x,diag.panel=panel.hist,lower.panel=panel.cor,upper.panel=panel.smooth)
} #else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
}
else #smooth is not true
{ if (scale) {pairs(x,diag.panel=panel.hist,lower.panel=panel.cor.scale)
} else {pairs(x, pch=pchs, col=cols,
	diag.panel=panel.hist,lower.panel=panel.cor) }
} #end of else (smooth)
} #end of function
