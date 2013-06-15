
# timingOd <- read.csv("./data/timingOd.csv")
# timingO2 <- read.csv("./data/timingO2.csv")
# timingOx <- read.csv("./data/timingOx.csv")
# 
# timing <- data.frame(its=timingOd$its, sseOd=timingOd$sseOd, Od=timingOd$Od,
#                      sseO2=timingO2$sseO2, O2=timingO2$O2,
#                      sseOx=timingOx$sseOx, Ox=timingOx$Ox)

# timing <- read.csv("R/data/report2013/perf.csv")
# 
# SSOR <- timing$SSOR
# SSOR.par <- timing$par_SSOR
# 
# timing <- data.frame(its=c(0, 15, timing$its[-1]), core2=c(0, 0.67, (timing$core2/timing$core2.par)[-1]), i5=c(0, 0.4, (timing$i5/timing$i5.par)[-1]),
#                      i7=c(0, 0.5, (timing$i7/timing$i7.par)[-1]))
timing <- read.csv("R/data/ololo.csv")

plot.timing <- function(plot.title="Сравнение производительности",
                        x.label="Размерность матрицы", y.label="Время выполнения (мс.)"){
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
  y.legend.corner.coords = c(10, 7)
  legend(c(0, 800), y.legend.corner.coords, legend=c(names(timing)[-1]), col=mode.col, 
         pch=rep(19, times=m), cex=rep(.8, times=m), y.intersp=.5, x.intersp=.5)
  
  title(plot.title)
  
}


# plot.div <- function(){
#   plot(timing$its[-1], (timing$SSOR/timing$par_SSOR)[-1], type="l", 
#        xlab="Количество разбиений", ylab="Отношение количества итераций")
#   title("График соотношения количества итераций при использовании
#         последовательного и параллельного предобуславливатей \n с изменением количества разбиений")
#   print(timing$SSOR/timing$par_SSOR)
# }
# 
# s <- read.csv("R/data/report2013/solve4.csv")
# 
# temperaturegram <- function(I, J, a, b, h1=a/(I-1), h2=b/(J-1)){
#   x <- seq(from=0, to=a, by=h1)
#   y <- seq(from=0, to=b, by=h2)
#   sl <- matrix(s$solve, (I-2), (J-2), byrow=FALSE)
#   image(x[2:(I-1)], y[2:(J-1)], sl, xlab="", ylab="")
#   title("Решение 4")
# }
