# Get the read tables for MiSeq and Nextseq, for both GPD and ZEV
# Get the number of reads for each sequence; join tables as in original analysis;
# plot # reads in both against each other.
setwd('D:/Promoter Design Data/Read Table Comparison/')
library(data.table)

miseq.GPD = 'filtered_read_table_miseq_GPD.txt'
nextseq.GPD = 'final_merged_reads_nextseq_GPD.csv'
miseq.ZEV = 'filtered_read_table_miseq_ZEV.txt'
nextseq.ZEV = 'final_merged_reads_nextseq_ZEV.csv'
key.len = 35

process.one.pair = function(miseq, nextseq, key.len) {
  miseq.table = fread(miseq); miseq.table$V2 = NULL # drop groups
  nextseq.table = fread(nextseq, header = TRUE); nextseq.table$V1 = NULL
  # get keys in miseq (first 35 bp); drop repeated rows
  miseq.table$keys = sapply(miseq.table$V1, function(x) substr(x, 1, key.len))
  dups = table(miseq.table$keys); dups = dups[dups > 1]; dups = names(dups)
  miseq.table = miseq.table[!miseq.table$keys %in% dups,]
  # drop seqs
  miseq.table$V1 = NULL; seqs = miseq.table$keys; miseq.table$keys = NULL
  read.cts = apply(miseq.table, 1, sum)
  miseq.table = data.table(seqs = seqs, miseq.cts = read.cts)  
  
  nextseq.table$keys = sapply(nextseq.table$Seq, function(x) substr(x, 1, key.len))
  dups = table(nextseq.table$keys); dups = dups[dups > 1]; dups = names(dups)
  nextseq.table = nextseq.table[!nextseq.table$keys %in% dups,]
  nextseq.table$Seq = NULL; seqs = nextseq.table$keys; nextseq.table$keys = NULL
  read.cts = apply(nextseq.table, 1, sum)
  nextseq.table = data.table(seqs = seqs, nextseq.cts = read.cts)
  setkey(miseq.table, seqs); setkey(nextseq.table, seqs)
  ans = miseq.table[nextseq.table, on = 'seqs', nomatch = 0]
  return(ans)
}

dat.GPD = process.one.pair(miseq.GPD, nextseq.GPD, key.len)
dat.ZEV = process.one.pair(miseq.ZEV, nextseq.ZEV, key.len)

png(filename = 'D:/Promoter Design Data/Figures/PNGs/S27.png',
    units = 'in', res = 144, width = 7, height = 5.25)
par(mfrow = c(2,3))
hist(log10(dat.GPD$miseq.cts), breaks = 100, main = '', xlab = 'Read Counts (log10) - GPD MiSeq', xlim = c(0,4))
hist(log10(dat.GPD$nextseq.cts), breaks = 100, main = '', xlab = 'Read Counts (log10) - GPD NextSeq', xlim = c(0,4))
plot(log10(dat.GPD$miseq.cts), log10(dat.GPD$nextseq.cts), pch = '.', main = '', xlab = 'Read Counts (log10) - GPD Miseq', ylab = 'Read Counts (log10) - GPD Nextseq')
hist(log10(dat.ZEV$miseq.cts), breaks = 100, main = '', xlab = 'Read Counts (log10) - ZEV MiSeq', xlim = c(0,4))
hist(log10(dat.ZEV$nextseq.cts), breaks = 100, main = '', xlab = 'Read Counts (log10) - ZEV NextSeq', xlim = c(0,4))
plot(log10(dat.ZEV$miseq.cts), log10(dat.ZEV$nextseq.cts), pch = '.', main = '', xlab = 'Read Counts (log10) - ZEV Miseq', ylab = 'Read Counts (log10) - ZEV Nextseq')
dev.off()



