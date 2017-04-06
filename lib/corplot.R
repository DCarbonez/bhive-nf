corplot <- function(dataset,cor.use="na.or.complete") {
  if (!is.data.frame(dataset)) {
     stop('data must be of type data.frame')
     return()
  }
  data = dataset
  cols = ncol(data)
  rows = cols
  D = 0.8/cols
  maxvals = unlist(lapply(data,max))
  par(mfrow=c(rows,cols))
  # Remove infinite values -> NaN
  for (i in 1:length(data)) {
      data[is.infinite(data[,i]),i] = NaN
  }
  corvals = cor(data, use=cor.use)
  for (j in 1:cols) {
    for (i in 1:rows) {
      figlims = 0.1+c((j-1)*D,j*D,(i-1)*D,i*D)
      margins = c(0.5,0.5,0.5,0.5)
      if (i==1 && j==1) { par(fig=figlims,mar=margins,new=F) }
      else {par(fig=figlims,mar=margins,new=T)}
      if (j <= i) {
        # prepare plot for text
        plot(c(0, 1), c(0, 1), ann = F, type = 'n', xaxt='n',yaxt='n')
        if (j < i) {
          # show correlation
            text(x = 0.5, y = 0.5, gettextf("%.4f",corvals[i,j]), cex = 1.6, col = "black")
        } else {
          # plot names(data)[i]
          text(x = 0.5, y = 0.5, names(data)[i], cex = 1.6, col = "black")
        }
      } else {
        plot(data[,j], data[,i], xlim=c(0,maxvals[j]), ylim=c(0,maxvals[i]), pch=1, ann=F, xaxt='n', yaxt='n',cex=1)
        axis(side=1, tick=F, labels=F)
        grid()
        if (i == 1) { axis(side=1) }
        if (j == 1) { axis(side=2) }
        if (i == rows) { axis(side=3) }
        if (j == cols) { axis(side=4) }
      }
    }
  }
}

