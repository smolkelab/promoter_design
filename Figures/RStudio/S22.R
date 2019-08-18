# Figure 6A: demonstrate analysis by single and double mutagenesis for sequence 0_0

library(data.table)
library(ggplot2)
library(reshape2)
library(viridis)

setwd('D:/Promoter Design Data/Characterize Mutagenesis/')
png.head = 'D:/Promoter Design Data/Figures/PNGs/'

diagonal.zoom.melted = function(x,w) {
  ans = matrix(nrow = nrow(x)*(2*w+1), ncol = 3)
  colnames(ans) = c('Position','Offset','Score')
  offset = -1*w
  position = 0
  for(i in 1:nrow(ans)) {
    if(position == nrow(x)) {position = 0; offset = offset + 1}
    position = position + 1
    comp = position + offset
    val = 0
    if(comp > 0 && comp < nrow(x)) {val = x[comp, position]}
    ans[i,] = c(position, offset, val)
  }
  return(data.frame(ans))
}


# double mutagenesis analysis of whole 0_0
fn_l2_mat = '0_0/0_GPD_strong_nofilter_mean_screen_0.45_selected_0_double_l2.csv'
fn_l2 = paste0(png.head, 'S22_full.png', collapse = '')
png(filename = fn_l2, units = 'cm', width = 13.5, height = 10, res = 600)
x = as.matrix(read.csv(fn_l2_mat, header = FALSE))
colnames(x) = 1:ncol(x)
xd = data.frame(melt(x))
names(xd) = c('Position.1','Position.2','Score')
ggplot(data = xd, aes(x=Position.1,y=Position.2,fill=Score)) + geom_tile() + theme_bw() + scale_fill_viridis() +
  theme(legend.key.size = unit(0.5,'cm')) + xlab('Position') + ylab('Position')
dev.off()


fn_l2_diag = paste0(png.head, 'S22_diag.png', collapse = '')
png(filename = fn_l2_diag, units = 'cm', width = 11.5, height = 4, res = 600)
yd = diagonal.zoom.melted(x,20)
ggplot(data = yd, aes(x=Position,y=Offset,fill=Score)) + geom_tile() + theme_bw() + scale_fill_viridis() +
  theme(legend.position='none')
dev.off()










