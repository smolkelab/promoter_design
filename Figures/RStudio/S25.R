library(data.table)
setwd('D:/Promoter Design Data/Sequence Alignment Scores/')
fns = dir()

fn.ord = c(1,3,2,4)
plot.names = c('Forward alignment scores',# - pGPD',
               'Reverse alignment scores',#' - pGPD',
               'Forward alignment scores',#' - pZEV',
               'Reverse alignment scores')#' - pZEV')

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S25.png', 
    units = 'in', res = 144, width = 6, height = 6)

par(mfrow = c(2,2))
for(i in 1:length(fns)) {
  x = fread(fns[fn.ord[i]])$V1
  hist(x, main = '', xlab = plot.names[i])
  abline(v=192, lty = 2, lwd = 2, col = 'red')
}

dev.off()

