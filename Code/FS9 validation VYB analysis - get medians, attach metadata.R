setwd('D:/Promoter Design Data/Final Flow Validation FCS/')
source('D:/Promoter Design Data/Code/FCS file analysis.R')
library(data.table)

# fn.key is based on the plate maps provided by Twist
fn.key = 'plate_key.csv'
fn.design = 'D:/Promoter Design Data/Validation/FS9_validation_design_summary.csv'

subdirs = c('R1','R2','R3')
table.list = list()
for(i in 1:length(subdirs)) {
  x = subdirs[i]
  fns = dir(paste0(getwd(),'/',x, collapse=''))
  fns = paste0(x,'/',fns)
  well.id = tstrsplit(fns, 'P')[[2]] # P for plate ID
  well.id = paste0('P', well.id)
  well.id = tstrsplit(well.id, '.', fixed = TRUE)[[1]]
  table.list[[i]] = data.table(Filename=fns, Replicate=x, Condition = well.id)
}
dat = rbindlist(table.list)

# Extract the medians
get.one.median = function(fn) {
  x = extract.data(load.fromfilename(fn))
  rats = log10(x[,10]/x[,8])
  return(median(rats))
}
dat$Median = sapply(dat$Filename, get.one.median)

# Separate out the induction condition
dat$Induced = !is.na(tstrsplit(dat$Condition, '_')[[2]])
dat$Condition = tstrsplit(dat$Condition, '_')[[1]]
# Correct some inconsistent filenames from an error in naming files
# in the Induced condition
for(i in 1:nrow(dat)) {
  x = dat$Condition[i]
  if(nchar(x) == 4) {
    x = paste0(substr(x,1,3),'0',substr(x,4,4), collapse = '')
  dat$Condition[i] = x
}}

# Recast the table: one row per design, not per data file
dat.cast = dcast(dat, Condition + Induced ~ Replicate, value.var = 'Median')
dat.cast = dcast(dat.cast, Condition ~ Induced, value.var = c('R1','R2','R3'))

# Map the 'Condition' well to the insert sequences
table.key = fread(fn.key)
setkey(dat.cast, 'Condition'); setkey(table.key, 'Condition')
dat.joined = dat.cast[table.key, nomatch = 0]

# Use the insert sequences to join with the table 'FS9 validation design summary.csv'
table.designed = fread(fn.design)
# need to extract original sequences from dat.joined
get.seq = function(x) {
  if(is.na(x)) {return(NA)}
  stopifnot(nchar(x) %in% c(300, 366))
  if(nchar(x) == 366) { x = substr(x,29,341); x = paste0('TACGTAAATAATTAATAGTAGTGAC',x,'TGTCTAAAGGTGAAGAATTATTCAC',collapse = '')  }
  else { x = substr(x,29,275); x = paste0('TTTATCATTATCAATACTCGCCATTTCAAAGAATACGTAAATAATTAATAGTAGTGAC',x,'TGTCTAAAGGTGAAGAATTATTCACTGGTGTTGTCCCAATTTTGGTTGAATTAGATGG',collapse = '') }
  return(x)
}

dat.joined$Seqs = sapply(dat.joined$Insert.sequence, get.seq)
setkey(dat.joined, Seqs); setkey(table.designed, Seqs)

dat.all = merge(dat.joined, table.designed, all.x = TRUE)
write.csv(dat.all, 'D:/Promoter Design Data/Validation/validation_output.csv')




