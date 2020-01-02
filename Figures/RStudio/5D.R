library(data.table)
library(ggplot2)
png.head = 'D:/Promoter Design Data/Figures/PNGs/'
setwd('D:/Promoter Design Data/Characterize Mutagenesis/collected')
fns = dir()
fns = fns[startsWith(fns, '14_ZEV')]

DNA = c('A', 'C', 'G', 'T')

get.ta.reps = function(x) {
  y = rep(0, length(x))
  for(i in 1:(length(x)-1)) {
    if(i > 1) {nextval = y[i-1]+1}
    else {nextval = 1}
    if(x[i] == 'T' & x[i+1] == 'A') {y[i] = nextval; y[i+1] = nextval + 1}
  }
  return(y)
}

set.max = function(x) {
  for(i in 2:length(x)) {
    if(x[i] > x[i-1]) {x[(i - x[i] + 1):i] = x[i]   }
  }
  return(x)
}


x = lapply(fns, function(x) as.matrix(fread(x))[59:98,] )
mins = lapply(x, function(x) apply(x,1,min) )
calls = lapply(x, function(x) DNA[apply(x, 1, function(y) which(y == 0))])
tas = lapply(calls, function(x) set.max(get.ta.reps(x)))


for(i in 1:length(x)) { x[[i]] = data.table(score=mins[[i]], fn = fns[i], pos = 59:98,
                                            call = calls[[i]], ta.ct = tas[[i]]) }

x = rbindlist(x)

# overall effect of TA repeat
wilcox.test(x$score[x$ta.ct == 0], x$score[x$ta.ct != 0],'two.sided')$p.value # 9.120946e-221
sum(x$ta.ct == 0)
sum(x$ta.ct != 0)
mean(x$score[x$ta.ct == 0])
sd(x$score[x$ta.ct == 0])
mean(x$score[x$ta.ct != 0])
sd(x$score[x$ta.ct != 0])


# boxplot ALL bases
#boxplot(x$score ~ as.factor(x$ta.ct))

# get just the median score for the longest TA repeat in each seq

get.median = function(one.mins, one.tas) {
  mv = max(one.tas)
  x = one.mins[one.tas == mv]
  return(median(x))
}

get.rejects = function(one.mins, one.tas) {
  mv = max(one.tas)
  x = one.mins[one.tas != mv]
  return(x)
}


meds = vector(mode = 'numeric', length = length(fns))
for(i in 1:length(meds)) {  meds[i] = get.median(mins[[i]], tas[[i]]) }
lens = sapply(tas, max)

lens.cropped = lens
lens.cropped[lens.cropped > 12] = 12

#boxplot(meds ~ as.factor(lens.cropped))

# get all scores outside these regions
rejects = list()
for(i in 1:length(meds)) {  rejects[[i]] = get.rejects(mins[[i]], tas[[i]]) }
rejects = do.call(c, rejects)

box.df = data.frame(meds = meds, lens = lens.cropped)
rej.df = data.frame(meds = rejects, lens = 0)
box.df = rbind(box.df, rej.df)
box.df$lens[box.df$lens == 12] = '12+'
box.df$lens = as.factor(box.df$lens)
box.df$lens = factor(box.df$lens, levels = levels(box.df$lens)[c(1,4,5,2,3)])

box.df.nozero = box.df[box.df$lens != 0,]
box.zero = box.df[box.df$lens == 0,]
box.zero$lens = 'Non-TA'
box.zero$dummy = 'a'
box.df.nozero$dummy = ''
box.df.nozero$dummy[box.df.nozero$lens == '6'] = 'b'
box.df.nozero$dummy[box.df.nozero$lens == '8'] = 'c'
box.df.nozero$dummy[box.df.nozero$lens == '10'] = 'd'
box.df.nozero$dummy[box.df.nozero$lens == '12+'] = 'e'

stat_box_data = function(y) { return(data.frame(y = min(y) - 0.005, label = length(y))) }

png(filename = paste0(png.head, '5D.png', collapse = ''),
    units = 'cm', width = 9, height = 6, res = 600)
g = ggplot(data = box.df.nozero, aes(x = lens, y = meds, color = dummy)) + 
  geom_boxplot(outlier.shape = NA) + 
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 0.9, size=3) +
  geom_jitter(width = 0.2, size = 0.5) +
  geom_violin(data = box.zero, aes(fill = lens), alpha = 0.5) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
        legend.title = element_blank(), legend.position = 'none') +
  labs(x = 'TA repeat lengths\n(pZEV-induced)', y='Median score differentials')  +
  scale_x_discrete(limits=c('Non-TA', '6', '8', '10', '12+')) + 
  ylim(-0.18, 0)
print(g)
dev.off()


l6 = box.df.nozero$meds[box.df.nozero$lens == '6']
l8 = box.df.nozero$meds[box.df.nozero$lens == '8']
l10 = box.df.nozero$meds[box.df.nozero$lens == '10']
l12 = box.df.nozero$meds[box.df.nozero$lens == '12+']

wilcox.test(l6, l8, 'two.sided')$p.value
wilcox.test(l8, l10, 'two.sided')$p.value
wilcox.test(l10, l12, 'two.sided')$p.value



