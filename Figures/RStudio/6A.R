# Figure 6A: demonstrate analysis by single and double mutagenesis for sequence 0_0

library(data.table)
library(ggplot2)

setwd('D:/Promoter Design Data/Characterize Mutagenesis/')
png.head = 'D:/Promoter Design Data/Figures/PNGs/'

get.seq = function(x.m, n.bases) {
  if(length(n.bases) == 1) { n.bases = strsplit(n.bases, '')[[1]]}
  n.pos = 1
  seq.min = min(x.m$pos); seq.max = max(x.m$pos)
  ans = rep('',seq.max - seq.min + 1)
  for(i in 1:length(ans)) {
    posx = (i-1) + seq.min
    x.tmp = x.m[x.m$pos == posx,]
    x.tmp = x.tmp[abs(x.tmp$value) == min(abs(x.tmp$value)),]
    if(nrow(x.tmp) == 4) { ans[i] = n.bases[n.pos]; n.pos = n.pos + 1 } #{ ans[i] = 'N'}
    else { ans[i] = as.character(x.tmp$variable) }
  }
  return(ans)
}

# double mutagenesis analysis of whole 0_0
#fn_l2_mat = '0_0/0_GPD_strong_nofilter_mean_screen_0.45_selected_0_double_l2.csv'
#fn_l2 = paste0(png.head, '6A_double.png', collapse = '')
#png(filename = fn_l2, units = 'cm', width = 6.6, height = 6.6, res = 600)
#x = as.matrix(read.csv(fn_l2_mat, header = FALSE))
#heatmap(x, Rowv = NA, Colv = NA, labRow = NA, labCol = NA,
#        scale = 'none', col = gray.colors(256, start = 0, end = 1))
#dev.off()

#x = fread('char_mut_0_0.csv', header = FALSE)
x = fread('0_0/0_GPD_strong_nofilter_mean_screen_0.45_selected_0_single.csv')
names(x) = c('A','C','G','T')
x$pos = 1:nrow(x)
# single mutagenesis analysis of whole 0_0
x.m = melt(x, measure.vars = c('A','C','G','T'))
png(filename = paste0(png.head, '6A_single.png', collapse = ''),
    units = 'cm', width = 8, height = 3.5, res = 600)
p = ggplot(data = x.m, aes(x = pos, y = value, color = variable)) + geom_line(lwd = 0.2, alpha = 0.7) + 
  # fix alpha washing out legend, cf. https://stackoverflow.com/questions/5290003/how-to-set-legend-alpha-with-ggplot2
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlim(c(0,363)) +
  labs(x = 'Position', y = 'Score Diff.', color = 'Base')
print(p)
dev.off()

# single mutagenesis analysis of Gcn4p site

x.m = melt(x[26:50,], measure.vars = c('A','C','G','T'))

png(filename = paste0(png.head, '6A_single_Gcn4p.png', collapse = ''),
    units = 'cm', width = 8, height = 3.5, res = 600)
p = ggplot(data = x.m, aes(x = pos, y = value, color = variable)) + geom_line(lwd = 0.5, alpha = 0.7) + 
  # fix alpha washing out legend, cf. https://stackoverflow.com/questions/5290003/how-to-set-legend-alpha-with-ggplot2
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  theme_bw() +
  annotate('text', 38:44, rep(0.005, 7), label = strsplit('TGASTCA','')[[1]], size = 1.6) + 
  geom_hline(yintercept = 0, lty = 2) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm')) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(x = 'Position', y = 'Score Diff.', color = 'Base')
print(p)
dev.off()

