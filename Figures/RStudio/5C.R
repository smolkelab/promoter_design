# identify significant motifs across pGPD
library(data.table)
library(ggplot2)
setwd('D:/Promoter Design Data/Code/')
source('yeastract_motifs.R')
setwd('D:/Promoter Design Data/Characterize Mutagenesis/motifs')
png.head = 'D:/Promoter Design Data/Figures/PNGs/'

get.motifs = function(dir_in, header) {
  currdir = getwd()
  setwd(dir_in)
  fns = dir()
  use = startsWith(fns, header)
  fns = fns[use]
  ans = list()
  for(i in 1:length(fns)) {
    ans[[i]] = fread(fns[i], header = TRUE)
    ans[[i]]$fn = fns[i]
  }
  setwd(currdir)
  return(rbindlist(ans))
}

motifs.pgpd = get.motifs('D:/Promoter Design Data/Characterize Mutagenesis/motifs',
                         '3_')

get.positions = function(table.in) {
  ans = rep(0, max(table.in$End))
  for(i in 1:nrow(table.in)) {
    for(j in table.in$Start[i]:table.in$End[i]) {ans[j] = ans[j] +1}
  }
  return(ans)
}
#pos.pgpd = get.positions(motifs.pgpd)
#barplot(pos.pgpd)

motifs.pgpd = find.yeastract(motifs.pgpd)
motifs.pgpd = merge.fr(motifs.pgpd)

motifs.tmp = copy(motifs.pgpd)
motifs.tmp$V1 = NULL
motifs.tmp$End = NULL
motifs.tmp$Start = NULL
motifs.tmp$Strengths = NULL
motifs.tmp$Seq = NULL
motifs.tmp$fn = NULL
motifs.tmp$Seq.formatted = NULL
motifs.tmp = as.matrix(motifs.tmp)
represented = apply(motifs.tmp,1,sum) > 0
motifs.pgpd$has.tf = represented

#pos.pgpd.rep = get.positions(motifs.pgpd[motifs.pgpd$represented])
#pos.pgpd.norep = get.positions(motifs.pgpd[!motifs.pgpd$represented])
#barplot(pos.pgpd.rep)
#barplot(pos.pgpd.norep)

# non-Yeastract
#motifs.pgpd$AT4 = grepl('[AT]{4,}', motifs.pgpd$Seq)
#motifs.pgpd$GC4 = grepl('[GC]{4,}', motifs.pgpd$Seq)
#motifs.pgpd$has.bias = motifs.pgpd$AT4 | motifs.pgpd$GC4
#motifs.pgpd$Category = ''
#motifs.pgpd$Category[motifs.pgpd$has.tf & motifs.pgpd$has.bias] = 'Both'
#motifs.pgpd$Category[motifs.pgpd$has.tf & !motifs.pgpd$has.bias] = 'TF'
#motifs.pgpd$Category[!motifs.pgpd$has.tf & motifs.pgpd$has.bias] = 'Bias'
#motifs.pgpd$Category[!motifs.pgpd$has.tf & !motifs.pgpd$has.bias] = 'Neither'

#png(filename = paste0(png.head, 'new_5B_pGPD-motifs.png', collapse = ''),
#    units = 'cm', width = 8, height = 6, res = 600)
#ggplot(data = motifs.pgpd, aes(x = '', y = Category, color = Category, fill = Category)) + 
#  #geom_histogram(bins=20) +
#  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + 
#  theme_bw() + 
#  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
#  theme(axis.text = element_text(size=6), axis.title = element_text(size=8, face='bold'),
#        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
#        legend.key.size = unit(0.2,'cm'))
#dev.off()

tf.cts = apply(as.matrix(motifs.tmp),2,sum)
other.cts = sum(tf.cts[tf.cts < 10])
tf.cts = tf.cts[tf.cts >= 10]
tf.cts[11] = other.cts
names(tf.cts)[11] = 'Other'

motifs.pgpd$Length = nchar(motifs.pgpd$Seq)
motifs.pgpd$Total = motifs.pgpd$Strengths*motifs.pgpd$Length

sum(motifs.pgpd$Length)/(113*259)


# pie chart

motifs.tmp = data.table(motifs.tmp)
x = apply(motifs.tmp, 2, sum)
x = x[x > 0]

names.use = c('Hac1p', 'Ino2p', 'Hap1p', 'Gcn4p', 'Rtg1p', 'Stp1p', 'Stb5p')
y = x[names(x) %in% names.use]
#Stp1p Stb5p Hac1p Hap1p Ino2p Gcn4p Rtg1p 
#22    36     4    11     2    16    17 
names(y) = c('Stp1p', 'Stb5p/Oaf1p/Haa1p', 'Hac1p/Swi5p/Ace2p', 'Hap1p', 'Ino2/4p',
             'Gcn4p/Bas1p, Yap1/3p', 'Rtg1/3p')

# Source Data
write.csv(data.frame(Name=names(y), Counts=y), 'D:/Promoter Design Data/Source Data/5C.csv', quote = FALSE, row.names = FALSE)

png(filename = paste0(png.head, '5C.png', collapse = ''),
    units = 'cm', width = 11, height = 5, res = 600)
ggplot(data.frame(Name=names(y), Counts=y),
       aes(x='', y=Counts, fill=Name))+
  #aes(x=Name, y = Counts, fill=Name))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.4, size = 8, face='bold')) +
  theme(axis.text = element_text(size=6), axis.title = element_blank(),
        legend.text = element_text(size=6), legend.title = element_text(size=8, face='bold'),
        legend.key.size = unit(0.3,'cm')) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(), panel.border = element_blank()) + 
  labs(fill = 'Transcription factor') + 
  scale_fill_viridis_d(option='cividis')
dev.off()

# are TFs enriched? Compare against a random population.
# all TFs appear to be in fully randomized sequences.
# conclusion: TFs are indeed enriched (at least for 6-11 bp; not many longer)

#table(motifs.pgpd$Length)
#6  7  8  9 10 11 12 14 16 
#32 55 26 30 16 10  6  1  1 
true.totals = c(32,55,26,30,16,10,6,0,1,0,1)

#table(motifs.pgpd$Length[motifs.pgpd$has.tf])
#6  7  8  9 10 11 12 14 16 
#14 20 11 20 13  9  5  1  1 

true.hits = c(14,20,11,20,13,9,5,0,1,0,1)

# how many sequences of a given length will have a hit?
test.random.seqs = function(seq.len, num.seqs = 1e4) {
  DNA = c('A','C','G','T')
  seqs = sapply(1:num.seqs, function(x) paste0(sample(DNA, seq.len, replace = TRUE), collapse = ''))
  yeastract = find.yeastract(data.table(Seq=seqs))
  yeastract$Seq = NULL
  return(sum(apply(yeastract, 1, sum) > 0))
}
set.seed(42)

# 1075 1609 2312 2849 3365 3913 4297 4793 5210 5711 5908
cts = sapply(6:16, test.random.seqs)


p.vals = vector(mode='numeric', length = 16 - 6 + 1)
for(i in 1:length(p.vals)) {
  p.vals[i] = fisher.test(rbind(c(true.totals[i] - true.hits[i], true.hits[i]),
                                c(1e4 - cts[i], cts[i])), alternative = "two.sided")$p.value
}
names(p.vals) = 6:16
p.vals
#6            7            8            9           10           11           12 
#2.080517e-06 2.827857e-04 3.276758e-02 1.654516e-05 1.313685e-04 1.401907e-03 9.078892e-02 
#13           14           15           16 
#1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00



