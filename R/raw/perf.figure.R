
# timingOd <- read.csv("./data/timingOd.csv")
# timingO2 <- read.csv("./data/timingO2.csv")
# timingOx <- read.csv("./data/timingOx.csv")
# 
# timing <- data.frame(its=timingOd$its, sseOd=timingOd$sseOd, Od=timingOd$Od,
#                      sseO2=timingO2$sseO2, O2=timingO2$O2,
#                      sseOx=timingOx$sseOx, Ox=timingOx$Ox)

timing <- read.csv("R/data/performance.tests/test.test.csv")



plot.timing <- function(plot.title="Сравнение производительности",
                        x.label="Размерность входных данных", y.label="Время (с)"){
  steps = seq(1, length(timing$its), length(timing$its)/5)
  
  xlimit = c(0, max(timing$its))
  maxy = 0
  for (i in 2:length(timing))
    maxy = max(maxy, max(timing[,i]))
  ylimit = c(0, maxy)
  
  mode.col = c(2:length(timing)-1)
  
  plot(0, 0, col="white", xlim = xlimit, ylim=ylimit, xlab=x.label, ylab=y.label)
  for (i in 2:length(timing)){
    points(timing[,1][steps], timing[,i][steps], col=mode.col[i-1], pch=17)
    lines(timing[,1], timing[,i], col=mode.col[i-1], type="l")
  }
  m = length(timing)-1
  #first value - upper left corner, second - lower rigth
  y.legend.corner.coords = c(10, 13)
  legend(c(0, 5), y.legend.corner.coords, legend=c(names(timing)[-1]), col=mode.col, 
         pch=rep(19, times=m), cex=rep(.8, times=m), y.intersp=.5, x.intersp=.5)
  
  title("Сравнение производительности")
  
}