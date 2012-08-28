library(latticeExtra)
library(reshape)

n.positive <- 33

read.mrp <- function(filename) {
  data <- read.table(filename, header=TRUE, check.names=FALSE)
  colnames(data) <- c('FLD', 'Forest', 'GERP', 'NNet', 'SVMmap')
  data <- data[,c(1,4,5,2,3)]

  # Remove last row
  data <- data[-nrow(data),]

  melt(data)
}

panel.mrp <- function(horizontal, ...) {
  panel.xyplot(..., horizontal=horizontal)
  panel.bwplot(..., horizontal=FALSE, 
    notch=TRUE, box.ratio=1, pch="|", do.out=FALSE)
}

plot.mrp <- function(data) {
  col <- c('black')
  pch <- 1
  xyplot(value ~ variable, data, 
    horizontal=TRUE, jitter.x=TRUE, factor=2.5,
    scales=list(alternating=FALSE, y=list(limits=c(-0.01, 0.8), at=c(0, 0.2, 0.4, 0.6))), group=variable,
    strip=strip.custom(par.strip.text=list(cex=1),bg='transparent'),
    panel=panel.performance,
    xlab=NULL, ylab="R-Precision (17 positive instances)",
    par.settings=list(superpose.symbol=list(col=col, pch=pch), 
      box.rectangle=list(col='black'),
      box.umbrella=list(col='black')))
}


read.loo <- function(filename) {
  data <- read.table(filename, header=TRUE, check.names=FALSE)
  colnames(data) <- c('FLD', 'Forest', 'GERP', 'NNet', 'SVMmap')
  data <- data[,c(1,4,5,2,3)]

  threshs=1:30
  cdf <- matrix(nrow=length(threshs), ncol=ncol(data))

  colnames(cdf) <- colnames(data)

  for (i in threshs) {
    cdf[i,] <- colSums(data <= i)
  }
  cdf <- cdf / n.positive
  cdf <- as.data.frame(cdf)
  cdf$thresh <- threshs
  cdf <- melt(cdf, id.vars=c('thresh'))

  cdf
}

plot.loo <- function(data, x=value ~ thresh) {
  col <- c('darkgray', 'darkgray', 'darkgray', 'black', 'red') 
  lty <- c(2,3,4,1,5)
  xyplot(x, data,
    scales=list(y=list(limits=c(0,max(data$value)*1.05), axs='i', alternating=1, at=c(0, 0.2, 0.4, 0.6)),
                x=list(limits=c(0,max(threshs)+0.5), alternating=1, axs='i', at=c(5, 10, 15, 20))),
    groups=variable, type='l',
    xlab="Rank cutoff", ylab="Fraction of 33 diseases",
    strip=strip.custom(par.strip.text=list(cex=1),bg='transparent'),
    auto.key=list(lines=TRUE,points=FALSE,columns=3,cex=0.9,padding.text=1.2),
    par.settings=list(superpose.line=list(col=col, lty=lty, lwd=2), strip.border=list(col=NA))
  )
}





model.curve.areas <- function(clip=TRUE) {
  suffix <- if (clip) { 'average.clip' } else { 'average' }
  models <- matrix(c('forest', 'RF', 'gerp', 'GERP', 'svm', 'SVM', 'nn8', 'NN', 'fld', 'FLD'), ncol=2, byrow=TRUE)
  datasets <- matrix(c('rare', 'Rare (AF < 5%)', 'common', 'Common (5% < AF < 50%)'), ncol=2, byrow=TRUE)
  max.x <- 25

  for (i in 1:nrow(models)) {
    for (j in 1:nrow(datasets)) {
      fn <- paste(models[i,1],datasets[j,1],suffix,sep='.')
      fp <- paste('data', fn, sep="/")
      data <- read.table(fp, header=FALSE, check.names=FALSE)
      # Remove last row
      data <- data[-nrow(data),]
      colnames(data) <- c('x', 'y')

      data$y <- data$y * n.positive
      data$x <- data$x * switch(j, n.negative.rare, n.negative.common)

      # add interpolated point
      interp <- approx(data$x, data$y, max.x)
      data <- rbind(data[data$x < max.x,], interp)

      # calculate area
      trap.area <- sum(diff(data$x) * (head(data$y, n=-1) + tail(data$y, n=-1)) / 2)
      print(paste(models[i,2], datasets[j,2], trap.area))
    }
  }
  for (j in 1:nrow(datasets)) {
    rand.slope <- n.positive / switch(j, n.negative.rare, n.negative.common)
    rand.y <- rand.slope * max.x
    trap.area <- rand.y * max.x / 2
    print(paste("random", datasets[j,2], trap.area))
  }
}



get.model.curves <- function(clip=TRUE) {
  suffix <- if (clip) { 'average.clip' } else { 'average' }
  models <- matrix(c('forest', 'RF', 'gerp', 'GERP', 'svm', 'SVM', 'nn8', 'NN', 'fld', 'FLD'), ncol=2, byrow=TRUE)
  datasets <- matrix(c('rare', 'Rare (AF < 5%)', 'common', 'Common (5% < AF < 50%)'), ncol=2, byrow=TRUE)
  data <- c()
  for (i in 1:nrow(models)) {
    for (j in 1:nrow(datasets)) {
      fn <- paste(models[i,1],datasets[j,1],suffix,sep='.')
      fp <- paste('data', fn, sep="/")
      data <- rbind(data, read.curve(models[i,2], datasets[j,2], fp))
    }
  }

  data$group <- factor(data$group, levels=models[,2])
  data$dataset <- factor(data$dataset, levels=datasets[,2])

  data
}

get.feature.curves <- function(clip=TRUE, model='forest') {
  suffix <- if (clip) { 'average.clip' } else { 'average' }

  #forest.tt.d6.rare.average.clip
  features <- matrix(c(
      'all', 'all', 
      'd1', 'no GERP',
      'd5', 'no splice',
      'd4', 'no ESX',
      'd2', 'no RSCU',
      'd3', 'no CpG',
      'd6', 'no ∆G'
    ), ncol=2, byrow=TRUE)
  datasets <- matrix(c('rare', 'Rare (AF < 5%)', 'common', 'Common (5% < AF < 50%)'), ncol=2, byrow=TRUE)
  data <- c()
  for (i in 1:nrow(features)) {
    for (j in 1:nrow(datasets)) {
      fn <- paste(model,'tt',features[i,1],datasets[j,1],suffix,sep='.')
      fp <- paste('data', fn, sep="/")
      data <- rbind(data, read.curve(features[i,2], datasets[j,2], fp))
    }
  }
  data$group <- factor(data$group, levels=features[,2])
  data$dataset <- factor(data$dataset, levels=datasets[,2])

  data
}

plot.model.curves <- function(data) {
  col <- c('black', 'red', 'darkgray', 'darkgray', 'darkgray') 
  lty <- c(1,5,2,3,4)
  xyplot(y ~ x | dataset, data,
    groups=group, type='l', aspect='iso', layout=c(2,1),
    xlab="False Positive Rate (FPR)", ylab="True Positive Rate (TPR)",
    strip=strip.custom(par.strip.text=list(cex=1),bg='transparent'),
    auto.key=list(lines=TRUE,points=FALSE,columns=3,cex=0.9,padding.text=1.2),
    par.settings=list(superpose.line=list(col=col, lty=lty))
  )
}



plot.model.curves.unclipped <- function(data) {
  data$y <- data$y * n.positive
  rows.rare <- data$dataset == "Rare (AF < 5%)"
  data$x[rows.rare] <- data$x[rows.rare] * n.negative.rare
  data$x[!rows.rare] <- data$x[!rows.rare] * n.negative.common

  col <- c('black', 'red', 'darkgray', 'darkgray', 'darkgray') 
  lty <- c(1,5,2,3,4)
  xyplot(y ~ x | dataset, data,
    scales=list(y=list(limits=c(0,5.9), axs='i', alternating=1),
                x=list(limits=c(0,26), alternating=1, axs='i')),
    groups=group, type='l', layout=c(2,1),
    panel=panel.ref,
    xlab="Number of false positives", ylab="Number of true positives",
    strip=strip.custom(par.strip.text=list(cex=1), bg='transparent'),
    auto.key=list(lines=TRUE,points=FALSE,columns=2,cex=0.9,padding.text=1.2),
    par.settings=list(superpose.line=list(col=col, lty=lty), strip.border=list(col=NA))
  )
}

plot.feature.curves <- function(data) {
  col <- c('black', 'red', gray.colors(nlevels(data$group)))
  lty <- c(1,5,2,3,4,6,7)
  xyplot(y ~ x | dataset, data,
    groups=group, type='l', aspect='iso', layout=c(2,1),
    xlab="False Positive Rate (FPR)", ylab="True Positive Rate (TPR)",
    strip=strip.custom(par.strip.text=list(cex=1),bg='transparent'),
    auto.key=list(lines=TRUE,points=FALSE,columns=4,cex=0.9,padding.text=1.2),
    par.settings=list(superpose.line=list(col=col, lty=lty))
  )
}

plot.feature.curves.unclipped <- function(data) {
  data$y <- data$y * n.positive
  rows.rare <- data$dataset == "Rare (AF < 5%)"
  data$x[rows.rare] <- data$x[rows.rare] * n.negative.rare
  data$x[!rows.rare] <- data$x[!rows.rare] * n.negative.common

  col <- c('black', 'red', gray.colors(nlevels(data$group)))
  lty <- c(1,5,2,3,4,6,7)
  xyplot(y ~ x | dataset, data,
    scales=list(y=list(limits=c(0, 5.9), axs='i', alternating=1),
                x=list(alternating=1, limits=c(0,26), axs='i')),
    groups=group, type='l', layout=c(2,1),
    panel=panel.ref,
    xlab="Number of false positives", ylab="Number of true positives",
    strip=strip.custom(par.strip.text=list(cex=1),bg='transparent'),
    auto.key=list(lines=TRUE,points=FALSE,columns=2,cex=0.9,padding.text=1.2),
    par.settings=list(superpose.line=list(col=col, lty=lty), strip.border=list(col=NA))
  )
}




read.cv.perf <- function(filename) {
  data <- read.delim(filename, check.names=FALSE, header=FALSE, comment.char="#")
  data <- data[,c(1,3,6,7,8)]
  colnames(data) <- c('model', 'dataset', '25', '10', '5')
  data <- melt.data.frame(data, id.vars=c('model', 'dataset'))
  ##data$variable <- as.numeric(levels(data$variable))[data$variable]
  old <- rev(c('25', '10', '5'))
  new <- rev(c('In top 25', 'In top 10', 'In top 5'))
  data$variable <- factor(data$variable, levels=old, labels=new)

  old <- c('forest', 'svm', 'nn8', 'fld', 'gerp')
  new <- c('RF', 'SVM', 'NN', 'FLD', 'GERP')
  data$model <- factor(data$model, levels=old, labels=new)

  data <- subset(data, dataset=="rare")
  data$dataset <- factor(data$dataset)

  data
}

get.performance <- function() {
  read.cv.perf('data/cv_model_performance.txt')
}

panel.performance <- function(horizontal, ...) {
  panel.xyplot(..., horizontal=horizontal)
  panel.bwplot(..., horizontal=FALSE, 
    notch=TRUE, box.ratio=1, pch="|", do.out=FALSE)
}

plot.performance <- function(data) {
  col <- c('black')
  pch <- 1
  xyplot(value ~ model | variable, data, 
    groups=dataset, horizontal=TRUE, jitter.x=TRUE, factor=1.0,
    layout=c(3, 1),
     scales=list(alternating=FALSE),
    strip=strip.custom(par.strip.text=list(cex=1),bg='transparent'),
    panel=panel.performance,
    xlab="Method used to rank SNVs", ylab="# of held-out SNVs ranked in top 5, 10, 25",
    par.settings=list(superpose.symbol=list(col=col, pch=pch), 
      box.rectangle=list(col='black'),
      box.umbrella=list(col='black')))
}

#data.melt <- melt.data.frame(data, id.vars=c('model', 'features'))
#data.melt$variable <- as.numeric(levels(data.melt$variable))[data.melt$variable]

#data.all <- subset(data.melt, features=="all")
#data.all$model <- factor(data.all$model)

#data.rf <- subset(data.melt, model=="random forest")

#data.plot1 <- xyplot(value ~ variable, 
#  data.all, groups=model, type='l', 
#  auto.key=list(points=FALSE, lines=TRUE, columns=2), 
#  xlab="Top N variants", ylab="Fraction of CV iterations",
#  xlim=rev(extendrange(data.all$variable)),
#  ylim=c(0,20),
#  scales=list(x=list(at=c(100,50,25,10,5))))
