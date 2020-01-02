library(data.table)
library(ggplot2)
# use 3, 14, 25
setwd('D:/Promoter Design Data/Characterize Mutagenesis/summarized_use')
png.head = 'D:/Promoter Design Data/Figures/PNGs/'

line.width = 1

# gray out the constant regions
const.gpd = data.frame(x=c(41,50,50,41), y=c(0,0,1,1), t='a')
const.gpd = rbind(const.gpd, data.frame(x=c(62,69,69,62), y=c(0,0,1,1), t='b'))
const.gpd = rbind(const.gpd, data.frame(x=c(95,102,102,95), y=c(0,0,1,1), t='c'))
const.gpd = rbind(const.gpd, data.frame(x=c(111,120,120,111), y=c(0,0,1,1), t='d'))
const.gpd = rbind(const.gpd, data.frame(x=c(195,204,204,195), y=c(0,0,1,1), t='e'))
const.gpd = rbind(const.gpd, data.frame(x=c(270,276,276,270), y=c(0,0,1,1), t='f'))

const.zev = data.frame(x=c(41,49,49,41), y=c(0,0,1,1), t='a')
const.zev = rbind(const.zev, data.frame(x = c(57, 65, 65, 57), y = c(0,0,1,1), t= 'b'))
const.zev = rbind(const.zev, data.frame(x = c(75, 83, 83, 75), y = c(0,0,1,1), t= 'c'))
const.zev = rbind(const.zev, data.frame(x = c(121, 138, 138, 121), y = c(0,0,1,1), t= 'd'))
const.zev = rbind(const.zev, data.frame(x = c(204, 210, 210, 204), y = c(0,0,1,1), t= 'e'))

prep.one.summarized = function(fn.in, setname, start, end) {
  x = fread(fn.in)
  x = x[start:end,]
  setnames(x, 'V1', 'Score')
  x$Score = x$Score/max(x$Score)
  x$Position = 1:nrow(x)
  x$Set = setname
  return(x)  
}

fns = c('3.csv', '14.csv', '25.csv')
ns = c('pGPD', 'pZEV-I', 'pZEV-AR')
starts = c('26', '59', '59')
ends = c('337', '304', '304')
preps = list()
for(i in 1:length(fns)) {
  preps[[i]] = prep.one.summarized(fns[i], ns[i], starts[i], ends[i])
}

preps = rbindlist(preps)

ylab.txt = 'Position\nimportance score'

png(filename = paste0(png.head, '5A_pGPD.png', collapse = ''),
    units = 'cm', width = 8, height = 3.5, res = 600)
ggplot(data = preps[preps$Set == 'pGPD',], aes(x = Position, y = Score)) + #geom_point() + 
  geom_line(size=line.width) + theme_bw() + xlab('Position (pGPD)') + ylab(ylab.txt) +
  geom_polygon(data = const.gpd, mapping=aes(x=x, y=y, group=t), fill = '#888888') +
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm'))
dev.off()

png(filename = paste0(png.head, '5A_pZEV-I.png', collapse = ''),
    units = 'cm', width = 8, height = 3.5, res = 600)
ggplot(data = preps[preps$Set == 'pZEV-I',], aes(x = Position, y = Score)) + #geom_point() + 
  geom_line(size=line.width) + theme_bw() + xlab('Position (pZEV - induced)') + ylab(ylab.txt) + 
  geom_polygon(data = const.zev, mapping=aes(x=x, y=y, group=t), fill = '#888888') +
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm'))
dev.off()

png(filename = paste0(png.head, '5A_pZEV-AR.png', collapse = ''),
    units = 'cm', width = 8, height = 3.5, res = 600)
ggplot(data = preps[preps$Set == 'pZEV-AR',], aes(x = Position, y = Score)) + #geom_point() + 
  geom_line(size=line.width) + theme_bw() + xlab('Position (pZEV - activation ratio)') + ylab(ylab.txt) +
  geom_polygon(data = const.zev, mapping=aes(x=x, y=y, group=t), fill = '#888888') + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.2,'cm'))
dev.off()

get_seq = function(fn, xmin = 26, xmax = 65) {
  mat = fread(fn)
  mat = mat[xmin:xmax,]
  names(mat) = c('A','C','G','T')
  ans = vector('character', nrow(mat))
  for (i in 1:nrow(mat)) {
    x = abs(mat[i,])
    y = which(x == 0)
    if(length(y) != 1) {ans[i] = 'X'}
    else {ans[i] = names(mat)[y]}
  }
  return(paste0(ans, collapse = ''))
}
get_seq('14_ZEV_induced_nofilter_1sd_evolve-thresh_1.6_selected_0_single.csv',302,302)

# what positions are represented at -3?
get.seq.onepos = function(dir_in, header, pos) {
  currdir = getwd()
  setwd(dir_in)
  fns = dir()
  use = startsWith(fns, header)
  fns = fns[use]
  ans = vector(length = length(fns), mode = 'character')
  for(i in 1:length(fns)) {
    ans[i] = get_seq(fns[i], pos, pos)
  }
  setwd(currdir)
  return(ans)
}

seqs.gpd = table(get.seq.onepos('D:/Promoter Design Data/Characterize Mutagenesis/collected',
                               '3_', 335))
# A   C   T 
# 103   9   1
seqs.zev.i = table(get.seq.onepos('D:/Promoter Design Data/Characterize Mutagenesis/collected',
                                '14_', 302))
# A   C   T 
# 110   5   1
seqs.zev.ar = table(get.seq.onepos('D:/Promoter Design Data/Characterize Mutagenesis/collected',
                                  '25_', 302))
# A  C  T 
# 97 10  8 
sum(seqs.gpd) + sum(seqs.zev.i) + sum(seqs.zev.ar) # 344 - vs 10 total w/ T (2.91%, < 3%)

