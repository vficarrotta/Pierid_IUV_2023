### Plotting Phylogeny and ASR of IUV in Pieridae
### 'C:\Program Files\R\R-4.1.0\bin'

# Libraries and functions
library('ape');library(geiger);library('phytools');library('corHMM')

# Required libraries
library('ggplot2');library('caper')

# gplots is needed for the following functions: rich.colors
library('gplots')

# plotrix is needed for the following functions: draw.circle
library('plotrix')

# diversitree is needed for the internal plotting functions extracted below
library('diversitree') # last run using version 0.9-3
plot2.phylo <- diversitree:::plot2.phylo
radial.text <- diversitree:::radial.text

# required functionss

group.label.tip.rad3 <- function (obj, lab, col.bar, col.lab, lwd = 1, offset.bar = 0,
                                  offset.lab = 0, cex = 1, font = 1, check = FALSE, quiet = FALSE, lend="butt", arc.bar.width=5,
                                  ...)
{
  n.taxa <- obj$n.taxa
  np <- obj$Ntip
  n <- obj$n.spp
  if (is.null(n.taxa))
    dt <- 1/6/n * 2 * pi
  else dt <- (n.taxa/2 - 0.5 + 1/6)/n * 2 * pi
  theta <- obj$xy$theta[seq_len(obj$Ntip)] + 0.0225#0.01125 for large tree
  t0 <- tapply(theta - dt, lab, min)
  t1 <- tapply(theta + dt, lab, max)
  str <- names(t0)
  if (check) {
    i <- order(t0)
    t0 <- t0[i]
    t1 <- t1[i]
    str <- str[i]
    g <- integer(length(t0))
    g[1] <- j <- 1
    end <- t1[1]
    for (i in seq_along(g)[-1]) {
      if (t0[i] > end) {
        j <- j + 1
        end <- t1[i]
      }
      else {
        end <- max(end, t1[i])
      }
      g[i] <- j
    }
    tg <- table(g)
    if (any(tg > 1)) {
      if (!quiet) {
        err <- sapply(which(tg != 1), function(x) paste(str[g ==
                                                              x], collapse = ", "))
        warn <- c("Collapsing non-monophyletic groups:",
                  sprintf("\t%s", err))
        warning(paste(warn, collapse = "\n"))
      }
      t0 <- tapply(t0, g, min)
      t1 <- tapply(t1, g, max)
      str <- as.character(tapply(str, g, collapse))
    }
  }
  tm <- (t0 + t1)/2
  r.bar <- rep(max(obj$xx) + offset.bar, length(t0))
  r.lab <- rep(max(obj$xx) + offset.lab, length(t0))
  arcs(t0, t1, r.bar, col = col.bar, lwd = lwd, lend=lend, np=np, arc.bar.width=arc.bar.width)
  if (any(!is.na(col.lab)))
    radial.text(r.lab, tm, str, col = col.lab, font = font,
                cex = cex, ...)
}


arcs <- function (theta0, theta1, r, col, lty = 1, lwd = par("lwd"),
                  np=np, lend="butt", arc.bar.width)
{
  if (length(r) != length(theta0) || length(r) != length(theta1))
    stop("theta0, theta1 and r must be of same length")
  if (any(lty[!is.na(lty)] != 1))
    warning("lwd != 1 will probably be ugly for arcs")
  theta0 <- sort(theta0)
  theta1 <- sort(theta1)
  dx <- max(r) * 2 * pi/np
  nn <- pmax(2, ceiling((theta1 - theta0) * r/2/dx))
  tmp <- lapply(seq_along(nn), function(i) cbind(seq(theta0[i],
                                                     theta1[i], length = nn[i]), rep(r[i], nn[i])))
  tmp0 <- do.call(rbind, lapply(seq_along(nn), function(i) rbind(tmp[[i]][-nn[i],
                                                                          ])))
  tmp1 <- do.call(rbind, lapply(seq_along(nn), function(i) rbind(tmp[[i]][-1,
                                                                          ])))
  if (length(lwd) > 1)
    lwd <- rep(rep(lwd, length = length(theta0)), nn - 1)
  if (length(col) > 1)
    col <- rep(rep(col, length = length(theta0)), nn - 1)
  {
    segments(x0=tmp0[, 2] * cos(tmp0[, 1]), y0=tmp0[, 2] * sin(tmp0[,
                                                                    1]), x1=tmp0[, 2] * cos(tmp0[, 1])*arc.bar.width, y1=tmp0[, 2] * sin(tmp0[, 1])*arc.bar.width, lwd = lwd, col = col, lend=lend)
    }
}

color.legend <- function (xl, yb, xr, yt, legend, rect.col, cex = 1, align = "lt",
                          gradient = "x", border="grey", lwd=0.1,...) {

  oldcex <- par("cex")
  par(xpd = TRUE, cex = cex)
  gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col), gradient = gradient, border=border, lwd=lwd)
  if (gradient == "x") {
    xsqueeze <- (xr - xl)/(2 * length(rect.col))
    textx <- seq(xl + xsqueeze, xr - xsqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      texty <- yb - 0.2 * strheight("O")
      textadj <- c(0.5, 1)
    }
    else {
      texty <- yt + 0.2 * strheight("O")
      textadj <- c(0.5, 0)
    }
  }
  else {
    ysqueeze <- (yt - yb)/(2 * length(rect.col))
    texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      textx <- xr + 0.2 * strwidth("O")
      textadj <- c(0, 0.5)
    }
    else {
      textx <- xl - 0.2 * strwidth("O")
      textadj <- c(1, 0.5)
    }
  }
  text(textx, texty, labels = legend, adj = textadj, ...)
  par(xpd = FALSE, cex = oldcex)
}

gradient.rect <- function (xleft, ybottom, xright, ytop, reds, greens, blues, col = NULL, nslices = 50, gradient = "x", border = par("fg"), lwd) {
  if (is.null(col))
    col <- color.gradient(reds, greens, blues, nslices)
  else nslices <- length(col)
  nrect <- max(unlist(lapply(list(xleft, ybottom, xright, ytop),
                             length)))
  if (nrect > 1) {
    if (length(xleft) < nrect)
      xleft <- rep(xleft, length.out = nrect)
    if (length(ybottom) < nrect)
      ybottom <- rep(ybottom, length.out = nrect)
    if (length(xright) < nrect)
      xright <- rep(xright, length.out = nrect)
    if (length(ytop) < nrect)
      ytop <- rep(ytop, length.out = nrect)
    for (i in 1:nrect) gradient.rect(xleft[i], ybottom[i],
                                     xright[i], ytop[i], reds, greens, blues, col, nslices,
                                     gradient, border = border)
  }
  else {
    if (gradient == "x") {
      xinc <- (xright - xleft)/nslices
      xlefts <- seq(xleft, xright - xinc, length = nslices)
      xrights <- xlefts + xinc
      rect(xlefts, ybottom, xrights, ytop, col = col, lty = 0, lwd=lwd)
      rect(xlefts[1], ybottom, xrights[nslices], ytop,
           border = border, lwd=lwd)
    }
    else {
      yinc <- (ytop - ybottom)/nslices
      ybottoms <- seq(ybottom, ytop - yinc, length = nslices)
      ytops <- ybottoms + yinc
      rect(xleft, ybottoms, xright, ytops, col = col, lty = 0, lwd=lwd)
      rect(xleft, ybottoms[1], xright, ytops[nslices], border = border, lwd=lwd)
    }
  }
  invisible(col)
}

# Data
scales <- read.csv('D:/Pierid_phylogenetics/uv_mat_syn-edit-V2_2-22-22.csv')  # _2 does not have NA species
tree <- ladderize(read.tree('D:/Pierid_Phylogenetics/Pieridae_7114_final_names_trim_dates_2-22-2022.tre'))
physpc <- read.csv('D:/Pierid_Phylogenetics/full-tree_species.csv')

uv <- scales[,4:6];rownames(uv) <- scales[,3]           #male/female presence/absence
uv.data <- treedata(tree, uv, sort = T, warnings = T)   #match tree and data
phy <- uv.data$phy                                      #trimmed tree
data <- uv.data$data                                    #trimmed data

### Read in corHMM output 2-22-2022

    # *^& reinsert *^&

####
# Quick Plot
####
# male.ER.2cat
plotRECON(phy,male.ER.2cat$states,no.margin=TRUE,show.tip.label = FALSE)
# binary.ER.fix3, both presence
plotRECON(phy, binary.ER.fix3$states, no.margin=TRUE, show.tip.label = FALSE)

##########################
# Full phylogeny - Black #
##########################
tree <- ladderize(read.tree('D:/Pierid_Phylogenetics/Pieridae_7114_final_names_trim_dates_2-22-2022.tre'))
treespc <- read.csv('D:/Pierid_Phylogenetics/full-tree_species.csv')

fam1 <- treespc$Subfamily
tribe1 <- treespc$Tribe

data.colors <- rep("#000000",dim(tree$edge)[1]) # black, 0

th <- max(branching.times(tree))
gr.col <- gray.colors(2, start=0.90, end=0.95)

plot(tree, show.tip.label=FALSE, type="f", edge.width=1, no.margin=T, root.edge=TRUE, edge.color="white", plot=FALSE)
draw.circle(0,0, radius=th, col=gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-20, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-40, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-60, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-80, col =gr.col[1], border=gr.col[1])

par(new=T)
obj <- plot2.phylo(tree, show.tip.label=F, type="f", edge.width=2, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(data.colors))
add.scale.bar(x=-5, y=25, length=10)

par(new=T)
group.label.tip.rad3(obj, fam1, c( "black","light grey"), col.lab="black", offset.bar=.75, offset.lab=7, cex=1.0, lwd=15, arc.bar.width=1.025)

par(new=T)
group.label.tip.rad3(obj, tribe1, c("black","light grey"), col.lab="black", offset.bar=3.75, offset.lab=7, cex=1.0, lwd=15, arc.bar.width=1.025)

####################################################
# Full PHylogeny - IUV A/P, dropped taxa colored   #
####################################################
# tree is full tree
# phy is reduced for analysis

# Time tree, black branches with data, grey without
tree <- ladderize(read.tree('D:/Pierid_Phylogenetics/Pieridae_7114_final_names_trim_dates_2-22-2022.tre'))

# IUV data
scales <- read.csv('D:/Pierid_phylogenetics/uv_mat_syn-edit-V2_2-22-22.csv')

# Full tree species labeling (subfam, tribe, genus, species)
treespc <- read.csv('D:/Pierid_Phylogenetics/full-tree_species.csv')

# Mark tips dropped from tree in light grey, remove unwanted tips in illustrator
com <- phy$tip.label                        # Common Taxa
fulltree <- tree$tip.label                  # Full Tree
treeonly <- setdiff(fulltree, com)          # Tree only
ptreeonly <- scales$Species[which(scales$male_ultrastructure_dorsal==1)]

fam1 <- treespc$Subfamily
tribe1 <- treespc$Tribe

data.colors <- rep("#ffc04d",dim(tree$edge)[1]) # orange, 0
    # #ffdb99 # lighter orange
    # #ffc04d # darker orange
    # data.colors <- rep("#FFFFFF",dim(tree$edge)[1]) # white, 0
    # data.colors <- rep("#000000",dim(tree$edge)[1]) # black, 0

data.colors[match(match(scales$Species[which(scales$male_ultrastructure_dorsal==1)],tree$tip.label),tree$edge[,2])] <- "hotpink" #hotpink,1
data.colors[match(match(treeonly, tree$tip.label), tree$edge[,2])] <- '#6c8993' # darker grey blue, NA

th <- max(branching.times(tree))
gr.col <- gray.colors(2, start=0.90, end=0.95)

# Tree
    # pdf(file="Pieridae_7114_SHL_names_trim_dates.pdf", width=10, height=10, onefile=TRUE)

plot(tree, show.tip.label=FALSE, type="f", edge.width=1, no.margin=T, root.edge=TRUE, edge.color="white", plot=FALSE)
draw.circle(0,0, radius=th, col=gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-20, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-40, col =gr.col[1], border=gr.col[1])
draw.circle(0,0, radius=th-60, col=gr.col[2], border=gr.col[2])
draw.circle(0,0, radius=th-80, col =gr.col[1], border=gr.col[1])

par(new=T)
obj <- plot2.phylo(tree, show.tip.label=F, type="f", edge.width=2, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(data.colors))
add.scale.bar(x=-5, y=25, length=10)

par(new=T)
group.label.tip.rad3(obj, fam1, c( "black","light grey"), col.lab="black", offset.bar=.75, offset.lab=7, cex=1.0, lwd=15, arc.bar.width=1.025)

par(new=T)
group.label.tip.rad3(obj, tribe1, c("black","light grey"), col.lab="black", offset.bar=3.75, offset.lab=7, cex=1.0, lwd=15, arc.bar.width=1.025)

    # dev.off()

######################
# A/P Color settings #
######################
# male.ER.2cat
# tree used for analysis is phy
# scalesc is the overlapped data

# keep data common between full dataset and full tree using scales
scalesc <- data.frame()
for(i in 1:length(scales$Species)){
    if(scales$Species[i] %in% phy$tip.label == T){
        scalesc <- rbind(scales[i,], scalesc)
    }
}

# keep data common between full dataset and full tree using treespc
scalesc1 <- data.frame()
for(i in 1:length(treespc$Species)){
    if(treespc$Species[i] %in% phy$tip.label == T){
        scalesc1 <- rbind(treespc[i,], scalesc1)
    }
}

fam2 <- scalesc1$Subfamily[match(phy$tip.label, scalesc1$Species)]  # EX: data$clade[match(tree$tip.label,data$species)]
tribe2 <- scalesc1$Tribe[match(phy$tip.label, scalesc1$Species)]

# plotting & tabling
th <- max(branching.times(phy))
gr.col <- gray.colors(2, start=0.90, end=0.95)

nodes <- max.col(male.ER.2cat$states)
node.colors <- c('#ffc04d','cyan','black','magenta')[nodes]
tip.colors <- c("#ffc04d","cyan")[a[,2][match(phy$tip.label,a[,1])]+1]
model.colors <- c(tip.colors,node.colors)[phy$edge[,2]]

# alt tips
gg <- node.colors[phy$edge[,1][match(1:length(phy$Nnode),phy$edge[,2])]-length(phy$Nnode)]
xx <- cbind(tip.colors,tip.colors)
xx[,2][which(xx[,1] %in% "#ffc04d")] <- 'black'
xx[,2][which(xx[,1] %in% "cyan")] <- 'magenta'
df <- as.matrix(cbind(xx,gg))
alt.tip.colors <- apply(df,1,function(x) if(max(table(x))==2){names(which.max(table(x)))}else x[1])
model.colors.AP <- c(alt.tip.colors,node.colors)[phy$edge[,2]]
    # cyan = presence slow , 1,R1
    # orange #ffc04d = absence slow, 2,R1
    # magenta = presence fast, 1,R2
    # black = absence fast, 2,R2

    # pdf(file="Pieridae_male_UV_ER2.pdf", width=10, height=10, onefile=TRUE)

plot(phy, show.tip.label=FALSE, type="f", edge.width=1, no.margin=T, root.edge=TRUE, edge.color="white", plot=FALSE)
draw.circle(0,0, radius=th, col=gr.col[1], border=gr.col[1]); draw.circle(0,0, radius=th-20, col=gr.col[2], border=gr.col[2]); draw.circle(0,0, radius=th-40, col =gr.col[1], border=gr.col[1]); draw.circle(0,0, radius=th-60, col=gr.col[2], border=gr.col[2]); draw.circle(0,0, radius=th-80, col =gr.col[1], border=gr.col[1])

par(new=T)
obj <- plot2.phylo(phy, show.tip.label=FALSE, type="f", edge.width=2.4, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(model.colors.AP))
add.scale.bar(x=-5, y=15, length=10)

par(new=T)
group.label.tip.rad3(obj, fam2, c( "black","light grey"), col.lab="black", offset.bar=1, offset.lab=4, cex=1.2, lwd=25, arc.bar.width=1.025)

par(new=T)
group.label.tip.rad3(obj, tribe2, c("black","light grey"), col.lab="black", offset.bar=3.75, offset.lab=7, cex=1.0, lwd=15, arc.bar.width=1.025)

    # dev.off()

########################
# Mutual ornamentation #
########################
# binary.ER.fix3
# tree used for analysis is phy
# scalesc is the overlapped data

# plotting & tabling
th <- max(branching.times(tree))
gr.col <- gray.colors(2, start=0.90, end=0.95)

nodes <- max.col(binary.ER.fix3$states)
node.colors <- c('#ffc04d','cyan','magenta')[nodes]
tip.colors <- c("#ffc04d","cyan",'magenta')[c[,2][match(phy$tip.label,c[,1])]+1]
model.colors.mut <- c(tip.colors,node.colors)[phy$edge[,2]]

    # Orange #ffc04d = absence 1
    # Cyan = presence m 2
    # Hot pink #fe019a = presence both 3

### Tree plot
    # pdf(file="Pieridae_binary.ARD2.pdf", width=10, height=10, onefile=TRUE)

plot(phy, show.tip.label=FALSE, type="f", edge.width=1, no.margin=T, root.edge=TRUE, edge.color="white", plot=FALSE)
draw.circle(0,0, radius=th, col=gr.col[1], border=gr.col[1]); draw.circle(0,0, radius=th-20, col=gr.col[2], border=gr.col[2]); draw.circle(0,0, radius=th-40, col =gr.col[1], border=gr.col[1]); draw.circle(0,0, radius=th-60, col=gr.col[2], border=gr.col[2]); draw.circle(0,0, radius=th-80, col =gr.col[1], border=gr.col[1])

par(new=T)
obj <- plot2.phylo(phy, show.tip.label=FALSE, type="f", edge.width=2.4, no.margin=TRUE, root.edge=TRUE, edge.color=as.matrix(model.colors.mut))
add.scale.bar(x=-5, y=15, length=10)

par(new=T)
group.label.tip.rad3(obj, fam2, c( "black","light grey"), col.lab="black", offset.bar=1, offset.lab=4, cex=1.2, lwd=25, arc.bar.width=1.025)

par(new=T)
group.label.tip.rad3(obj, tribe2, c( "black","light grey"), col.lab="black", offset.bar=3.75, offset.lab=4, cex=1.2, lwd=25, arc.bar.width=1.025)

    # dev.off()
