library(data.table)
library(ggplot2)
png.head = 'D:/Promoter Design Data/Figures/PNGs/'
setwd('D:/Promoter Design Data/Characterize Mutagenesis/collected')
fns = dir()
fns = fns[startsWith(fns, '25_ZEV')]

DNA = c('A', 'C', 'G', 'T')
pos = list(c(108,111), c(124,127), c(142,145))
x = list(); mins = list(); calls = list()
for(i in 1:length(pos)) {
x[[i]] = lapply(fns, function(x) as.matrix(fread(x))[pos[[i]][1]:pos[[i]][2],] )
mins[[i]] = lapply(x[[i]], function(x) apply(x,1,min) )
calls[[i]] = lapply(x[[i]], function(x) DNA[apply(x, 1, function(y) which(y == 0))])
}

df = list()
for(i in 1:length(pos)) {
  df[[i]] = data.frame(Score=sapply(mins[[i]], median),
                       Seq=sapply(calls[[i]], function(x) paste0(x,collapse='')),
                       Pos = i)
}
df = rbindlist(df)
df$is.GCTA = df$Seq == 'GCTA'

for(i in 1:3) {
print(median(df$Score[df$Pos == i & df$is.GCTA]))
      print(median(df$Score[df$Pos == i & !df$is.GCTA]))
            print(wilcox.test(df$Score[df$Pos == i & df$is.GCTA], 
            df$Score[df$Pos == i & !df$is.GCTA], 'two.sided')$p.value)
}

stat_box_data = function(y) { return(data.frame(y = min(y) - 0.02, label = length(y))) }

df$Pos = paste0('Site ', df$Pos)

png(filename = paste0(png.head, '5E.png', collapse = ''),
    units = 'cm', width = 10, height = 6, res = 600)
ggplot(data = df, aes(x = is.GCTA, y = Score, color = is.GCTA)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9, size=3) +
  facet_grid(. ~ Pos) +
  geom_jitter(width = 0.2, size = 0.5) + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.position = 'none') + 
  labs(x = 'GCTA after ZEV ATF binding site\n(pZEV-activation ratio)',
       y='Median score differentials') + 
  ylim(-0.45, 0)
dev.off()


