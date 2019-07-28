# Select sequences for validation experiments

# Selection parameters
# Total # points to measure in each range
range.pts = 8; num.pts = 3; tolerance = 0.1

# Threshold # sequences to measure for each subpool
subpool.ct.thresh = 5

# Threshold # reads required for a sequence in each replicate to trust the measurement
read.thresh = 10

# Number of off-scale sequences to measure for each endpoint-subpool intersection
offscale.ct.a.min = 3
offscale.ct.a.max = 0
offscale.ct.b.min = 0
offscale.ct.b.max = 3

# Number of promising sequences to pick for each target
num.promising = 10

# cloning homology to affix to all sequences: use 28/25 bp here, and bring up to 40 total by PCR for gap-repair
# 28 on the front, because Twist minimum sequence length is 300 bp
homology = c('GAATACGTAAATAATTAATAGTAGTGAC','TGTCTAAAGGTGAAGAATTATTCAC')

# load analysis
setwd('D:/Promoter Design Data/Validation/')
require(data.table)
tables = fread('merged_FS9_results.csv'); tables$V1 = NULL
fn.summary = 'FS9_validation_design_summary.csv'
fn.fasta = 'FS9_validation_designs_from_R.fa'
header.fasta = '>FS9_validation_'

# First, though, are these means reliable? How many reads per sequence?
# Download the file 'merged_read_table.csv' from Instance 3, ~/facs-seq_test/design_testing
merged_reads = fread('merged_read_table.csv', header = TRUE)
merged_reads$`Unnamed: 0` = NULL
seqs = merged_reads$Seq
merged_reads$Seq = NULL
merged_reads_mat = as.matrix(merged_reads)
rownames(merged_reads_mat) = seqs
merged_reads_a = merged_reads_mat[,as.numeric(colnames(merged_reads_mat)) < 12]
merged_reads_b = merged_reads_mat[,as.numeric(colnames(merged_reads_mat)) >= 12]
merged_reads$Reads_A = apply(merged_reads_a, 1, sum)
merged_reads$Reads_B = apply(merged_reads_b, 1, sum)
merged_reads$Reads_Min = apply(cbind(merged_reads$Reads_A, merged_reads$Reads_B), 1, min)
merged_reads$Seqs = seqs
setnames(merged_reads, 'Seqs', 'Seqs.orig')
setkey(merged_reads, Seqs.orig)
setkey(tables, Seqs.orig)
tables = merge(tables, merged_reads, all.x = TRUE)
tables$Reads_Min[is.na(tables$Reads_Min)] = 0 # sequences that aren't represented at all in the read table

# Range validation
# identify points which can be used for validation
range_a = seq(from = min(tables$Means_A[tables$Reads_Min >= read.thresh], na.rm = TRUE), to = max(tables$Means_A[tables$Reads_Min >= read.thresh], na.rm = TRUE), length.out = range.pts)
range_b = seq(from = min(tables$Means_B[tables$Reads_Min >= read.thresh], na.rm = TRUE), to = max(tables$Means_B[tables$Reads_Min >= read.thresh], na.rm = TRUE), length.out = range.pts)

tables$range_a = NA; tables$range_b = NA
for(i in 1:range.pts) {
  tables$range_a[which(abs(tables$Means_A - range_a[i]) < tolerance)] = i
  tables$range_b[which(abs(tables$Means_B - range_b[i]) < tolerance)] = i
}

# Missing subpool correction
# identify subpools that need more sequences measured
u = unique(tables$Experiment.designed)
tables$Subpool.undersampled = NA
for(i in 1:length(u)) {
  if(sum(tables$Experiment.designed == u[i] & !is.na(tables$Scores) & tables$Reads_Min >= read.thresh) < subpool.ct.thresh) {
    tables$Subpool.undersampled[tables$Experiment.designed == u[i]] = u[i]
  }
}

# Off-scale correction
tables$Scale_A = 0
tables$Scale_B = 0
tables$Scale_A[tables$Means_A == min(tables$Means_A, na.rm = TRUE)] = -1
tables$Scale_A[tables$Means_A == max(tables$Means_A, na.rm = TRUE)] = 1
tables$Scale_B[tables$Means_B == min(tables$Means_B, na.rm = TRUE)] = -1
tables$Scale_B[tables$Means_B == max(tables$Means_B, na.rm = TRUE)] = 1

# Promising sequence identification: include strong controls if these are dominant
# Select best subpools w/r/t score, GC content, and diversity (considered together):
# focus on evolution strategies with GC filter.


# For GPD, use pool 23 ('GPD gcfilter evolve 0.6, 1sd')
tables$GPD_cand = tables$Experiment.designed == 23 & tables$Reads_Min >= read.thresh
# For ZEV Induced, use pool 34 ('ZEV induced gcfilter evolve 1.6, 1sd')
tables$ZEV_I_cand = tables$Experiment.designed == 34 & tables$Reads_Min >= read.thresh
  #   pick some with high AR!
tables$ZEV_I_AR_cand = tables$Target == 'Induced' & tables$AR > 2.5 & !is.na(tables$AR)
# For ZEV AR use pool 44 ('ZEV AR gcfilter evolve 1.85')
tables$ZEV_AR_cand = tables$Experiment.designed == 44 & tables$Reads_Min >= read.thresh
tables$ZEV_AR_cand[is.na(tables$ZEV_AR_cand)] = FALSE
# Also, for ZEV AR, synthesize promising off-scale seqs
tables$ZEV_AR_offscale_cand = tables$Means_A == min(tables$Means_A, na.rm = TRUE) & 
  tables$Means_B > 1.0 & tables$Reads_Min >= read.thresh
tables$ZEV_AR_offscale_cand[is.na(tables$ZEV_AR_offscale_cand)] = FALSE

# Lastly, pick 5-10 interesting outliers? Or not - just leaving this here in case I change my mind
tables$hand_selected = FALSE
hand.sel.rows = numeric(0)
tables$hand_selected[hand.sel.rows] = TRUE


# Identify sequences to potentially synthesize
tables$synthesis_score = as.numeric(!is.na(tables$range_a)) + 
  as.numeric(!is.na(tables$range_b)) +
  as.numeric(!is.na(tables$Subpool.undersampled)) + 
  as.numeric((tables$Scale_A != 0)) + as.numeric((tables$Scale_B != 0)) + 
  as.numeric(tables$GPD_cand) + as.numeric(tables$ZEV_I_cand) + 
  as.numeric(tables$ZEV_I_AR_cand) + as.numeric(tables$ZEV_AR_cand) + 
  as.numeric(tables$ZEV_AR_offscale_cand) + as.numeric(tables$hand_selected)


# Choose sequences to synthesize
tables$synthesize = FALSE

# Reject sequences that will be rejected by Twist 
# seqs that contain repeat regions of 20 bp or more
check.internal.repeats = function(seq, w=20) {
  substrs = sapply(1:(nchar(seq) - w+1), function(start, seq) substr(seq, start, start+(w-1)), seq = seq)
  return(max(table(substrs)) == 1)
}

# seqs that have a difference in GC content of > 52% between the lowest and highest 50-bp windows
check.gc.windows = function(seq,w=50, difference = 0.52) {
  substrs = sapply(1:(nchar(seq) - w+1), function(start, seq) substr(seq, start, start+(w-1)), seq = seq)  
  gcs = sapply(substrs, function(x) sum(strsplit(x, '')[[1]] %in% c('G','C')))
  gcs = gcs/w
  return(max(gcs) - min(gcs) < difference)
}

# seqs with a homopolymer of 10 (or more) bp
check.homopolymer = function(seq, w=10) {
  substrs = sapply(1:(nchar(seq) - w+1), function(start, seq) substr(seq, start, start+(w-1)), seq = seq)  
  longests = sapply(substrs, function(x) max(table(strsplit(x,'')[[1]]))  )
  longest = max(longests)
  return(longest < w)
}

tables$synthesizable = sapply(tables$Seqs, function(x) check.internal.repeats(x) & check.gc.windows(x) & check.homopolymer(x))
stopifnot(sum(is.na(tables$synthesizable)) == 0)

is.eligible = function(tables, idx) {
  r = tables[idx,]
  print(idx)
  # reject this sequence if it can't be synthesized, or if an identical sequence is already to be ordered (there are e.g. some duplicate controls)
  if(!(tables$synthesizable[i])) { return(FALSE) }
  if(tables$Seqs[i] %in% tables$Seqs[tables$synthesize]) { return(FALSE) }
  # range validation - don't use sequences that went off-scale in the other pass
  useful.range.a = !is.na(r$range_a) & r$Scale_B == 0 & r$Reads_Min >= read.thresh & sum(tables$synthesize[!is.na(tables$range_a) & tables$range_a == r$range_a & tables$Scale_B == 0]) < num.pts
  stopifnot(!is.na(useful.range.a))
  useful.range.b = !is.na(r$range_b) & r$Scale_A == 0 & r$Reads_Min >= read.thresh & sum(tables$synthesize[!is.na(tables$range_b) & tables$range_b == r$range_b & tables$Scale_A == 0]) < num.pts
  stopifnot(!is.na(useful.range.b))
  # collect missing subpools
  useful.missing.subpool = !is.na(r$Subpool.undersampled) & r$Reads_Min < read.thresh & sum(tables$synthesize[!is.na(tables$Subpool.undersampled) & tables$Subpool.undersampled == r$Subpool.undersampled]) < subpool.ct.thresh
  stopifnot(!is.na(useful.missing.subpool))
  # collect off-scale sequences
  useful.offscale.a = (r$Scale_A == -1 & sum(tables$synthesize[tables$Experiment.designed == r$Experiment.designed & tables$Scale_A == r$Scale_A]) < offscale.ct.a.min) |
    (r$Scale_A == 1 & sum(tables$synthesize[tables$Experiment.designed == r$Experiment.designed & tables$Scale_A == r$Scale_A]) < offscale.ct.a.max)
  
  useful.offscale.b = (r$Scale_B == -1 & sum(tables$synthesize[tables$Experiment.designed == r$Experiment.designed & tables$Scale_B == r$Scale_B]) < offscale.ct.b.min) |
      (r$Scale_B == 1 & sum(tables$synthesize[tables$Experiment.designed == r$Experiment.designed & tables$Scale_B == r$Scale_B]) < offscale.ct.b.max)

  stopifnot(!is.na(useful.offscale.a))
  stopifnot(!is.na(useful.offscale.a))
  # collect promising sequences
  useful.gpd = r$GPD_cand & sum(tables$synthesize[tables$GPD_cand]) < num.promising
  useful.zev.i = r$ZEV_I_cand & sum(tables$synthesize[tables$ZEV_I_cand]) < num.promising
  useful.zev.i.ar = r$ZEV_I_AR_cand & sum(tables$synthesize[tables$ZEV_I_AR_cand]) < num.promising
  useful.zev.ar = r$ZEV_AR_cand & sum(tables$synthesize[tables$ZEV_AR_cand]) < num.promising
  useful.zev.ar.offscale = r$ZEV_AR_offscale_cand & sum(tables$synthesize[tables$ZEV_AR_offscale_cand]) < num.promising
  useful.handpick = r$hand_selected & sum(tables$synthesize[tables$hand_selected]) < length(hand.sel.rows)
  stopifnot(!is.na(useful.gpd))
  stopifnot(!is.na(useful.zev.i))
  stopifnot(!is.na(useful.zev.i.ar))
  stopifnot(!is.na(useful.zev.ar))
  stopifnot(!is.na(useful.zev.ar.offscale))
  stopifnot(!is.na(useful.handpick))

  return(sum(c(useful.range.a, useful.range.b, useful.missing.subpool, useful.offscale.a, useful.offscale.b,
             useful.gpd, useful.zev.i, useful.zev.i.ar, useful.zev.ar, useful.zev.ar.offscale, useful.handpick)) > 0)
}

# Approach: reorder 'tables' to put 'most useful' sequences (those which can contribute to the most of my goals)
# on top, after randomizing the row order.
# Check whether each sequence is 'eligible' - that is, it can contribute to a goal we don't have enough sequences for yet - 
# and if yes, set it to be synthesized.

# reorder 'tables' for sequence selection
set.seed(2017)
tables = tables[sample(nrow(tables)),]
tables = tables[order(-tables$synthesis_score)]
for(i in 1:nrow(tables)) { tables$synthesize[i] = is.eligible(tables, i) }

sum(tables$synthesize)
plot(tables$Means_A[tables$Reads_Min >= read.thresh], tables$Means_B[tables$Reads_Min >= read.thresh], pch = 20)
points(tables$Means_A[tables$synthesize], tables$Means_B[tables$synthesize], pch = 20, col = 'red')

# Save a summary CSV of the sequence design process
write.csv(tables[tables$synthesize,],file = fn.summary, quote = FALSE, row.names = TRUE)


# Generate a FASTA file with these sequences; names are a header plus the row number
# Need to use tables$Seqs, but it's padded to 363 bp - remove padding so real sequence
# has 40 bp of homology on both sides
tables.synthesize = tables[tables$synthesize]
tables.synthesize$final.seq = ''
for(i in 1:nrow(tables.synthesize)) {
  if(tables.synthesize$Promoter[i] == 'GPD') {tables.synthesize$final.seq[i] = paste0(homology[1],substr(tables.synthesize[i]$Seqs,26,338),homology[2], collapse = '')}
  else {tables.synthesize$final.seq[i] = paste0(homology[1],substr(tables.synthesize[i]$Seqs,59,305),homology[2], collapse = '') }
}

# Test whether design goals actually met!
# Test range.pts
for(i in 1:range.pts) {
  print(c('Range position', i))
  print(c('A', sum(tables$range_a[tables$synthesize] == i, na.rm = TRUE)))
  print(c('B', sum(tables$range_b[tables$synthesize] == i, na.rm = TRUE)))
}

# Test subpool.ct
print('Undersampled subpools - total')
print(table(tables$Subpool.undersampled))
print('Undersampled subpools - to synthesize')
print(table(tables$Subpool.undersampled[tables$synthesize]))
print('Undersampled subpools - total after synthesis')
print(table(tables$Subpool.undersampled[tables$synthesize | (!is.na(tables$Means_A) & tables$Reads_Min >= read.thresh)]))

# Test offscales
print('Testing offscales')
print('Offscales by pool (total) - A min.')
offscales = vector(mode = 'numeric', length = max(tables$Experiment.designed) + 1)
print(table(tables$Experiment.designed[tables$Scale_A == -1]))
print('Offscales by pool (to synthesize) - A min.')
print(table(tables$Experiment.designed[tables$Scale_A == -1 & tables$synthesize]))
#print('Offscales by pool (total) - A max.')
#print(table(tables$Experiment.designed[tables$Scale_A == 1]))
#print('Offscales by pool (to synthesize) - A max.')
#print(table(tables$Experiment.designed[tables$Scale_A == 1 & tables$synthesize]))
#print('Offscales by pool (total) - B min.')
#print(table(tables$Experiment.designed[tables$Scale_B == -1]))
#print('Offscales by pool (to synthesize) - B min.')
#print(table(tables$Experiment.designed[tables$Scale_B == -1 & tables$synthesize]))
print('Offscales by pool (total) - B max.')
print(table(tables$Experiment.designed[tables$Scale_B == 1]))
print('Offscales by pool (to synthesize) - B max.')
print(table(tables$Experiment.designed[tables$Scale_B == 1 & tables$synthesize]))

# Test promising sets
print('Testing promising')
print('GPD')
print(sum(tables$GPD_cand[tables$synthesize]))
print('ZEV induced')
print(sum(tables$ZEV_I_cand[tables$synthesize]))
print('ZEV induced - AR')
print(sum(tables$ZEV_I_AR_cand[tables$synthesize]))
print('ZEV - AR')
print(sum(tables$ZEV_AR_cand[tables$synthesize]))
print('ZEV - AR, offscale')
print(sum(tables$ZEV_AR_offscale_cand[tables$synthesize]))

# Lastly, reorder the table to group GPD and ZEV sequences together - 
# this will make experiments downstream a little easier, as we don't need to test inducible conditions
# for GPD.
tables.synthesize = tables.synthesize[order(tables.synthesize$Promoter, decreasing = TRUE)]

lines.out = vector(mode='character', length = 2*nrow(tables.synthesize))
for(i in 1:nrow(tables.synthesize)) {
  lines.out[2*i - 1] = paste0(header.fasta, i, collapse = '')
  lines.out[2*i] = tables.synthesize$final.seq[i]
  
}

# control promoters: traditional control promoters and selected ZEV promoters
# Leaving these here for reference - but can save a good amount of money by just cloning them out of existing plasmids

control.names = c('>FS9_Validation_GPD',
             '>FS9_Validation_TEF',
             '>FS9_Validation_ADH',
             '>FS9_Validation_PGK',
             '>FS9_Validation_TPI',
             '>FS9_Validation_CYC',
             '>FS9_Validation_ZEV_Pr3',
             '>FS9_Validation_ZEV_Pr4',
             '>FS9_Validation_ZEV_Pr8')

# GPD from pBK18 (this version was modified to match WT; others may have cloning sites, etc., still)
# TEF from pBK18 
# ADH1 from pCS2660
# PGK1 from pCS2663
# TPI1 from pCS2661
# CYC1 from pCS2659
# ZEV_Pr3 from pTH3
# ZEV_Pr4 from pTH4
# ZEV_Pr8 from pTH5

control.seqs = c('GAGTTTATCATTATCAATACTCGCCATTTCAAAGAATACGTAAATAATTAATAGTAGTGATTTTCCTAACTTTATTTAGTCAAAGAATTAGCCTTTTAATTCTGCTGTAACCCGTACATGCCCAAAATAGGGGGCGGGTTACACAGAATATATAACATCGTAGGTGTCTGGGTGAACAGTTTATTCCTGGCATCCACTAAATATAATGGAGCCCGCTTTTTAAGCTGGCATCCAGAAAAAAAAAGAATCCCAGCACCAAAATATTGTTTTCTTCACCAACCATCAGTTCATAGGTCCATTCTCTTAGCGCAACTACAGAGAACAGGGGCACAAACAGGCAAAAAACGGGCACAACCTCAATGGAGTGATGCAACCTGCCTGGAGTAAATGATGACACAAGGCAATTGACCCACGCATGTATCTATCTCATTTTCTTACACCTTCTATTACCTTCTGCTCTCTCTGATTTGGAAAAAGCTGAAAAAAAAGGTTGAAACCAGTTCCCTGAAATTATTCCCCTACTTGACTAATAAGTATATAAAGACGGTAGGTATTGATTGTAATTCTGTAAATCTATTTCTTAAACTTCTTAAATTCTACTTTTATAGTTAGTCTTTTTTTTAGTTTTAAAACACCAAGAACTTAGTTTCGAataaacacacataaacaaacaaa',
                 'ATAGCTTCAAAATGTTTCTACTCCTTTTTTACTCTTCCAGATTTTCTCGGACTCCGCGCATCGCCGTACCACTTCAAAACACCCAAGCACAGCATACTAAATTTCCCCTCTTTCTTCCTCTAGGGTGTCGTTAATTACCCGTACTAAAGGTTTGGAAAAGAAAAAAGAGACCGCCTCGTTTCTTTTTCTTCGTCGAAAAAGGCAATAAAAATTTTTATCACGTTTCTTTTTCTTGAAAATTTTTTTTTTTGATTTTTTTCTCTTTCGATGACCTCCCATTGATATTTAAGTTAATAAACGGTCTTCAATTTCTCAAGTTTCAGTTTCATTTTTCTTGTTCTATTACAACTTTTTTTACTTCTTGCTCATTAGAAAGAAAGCATAGCAATCTAATCTAAGTTTTGCCGCGGGAAATA',
                 'gcatgcaacttcttttctttttttttcttttctctctcccccgttgttgtctcaccatatccgcaatgacaaaaaaatgatggaagacactaaaggaaaaaattaacgacaaagacagcaccaacagatgtcgttgttccagagctgatgaggggtatctcgaagcacacgaaactttttccttccttcattcacgcacactactctctaatgagcaacggtatacggccttccttccagttacttgaatttgaaataaaaaaaagtttgctgtcttgctatcaagtataaatagacctgcaattattaatcttttgtttcctcgtcattgttctcgttccctttcttccttgtttctttttctgcacaatatttcaagctataccaagcatacaattataca',
                 'aggcatttgcaagaattactcgtgagtaaggaaagagtgaggaactatcgcatacctgcatttaaagatgccgatttgggcgcgaatcctttattttggcttcaccctcatactattatcagggccagaaaaaggaagtgtttccctccttcttgaattgatgttaccctcataaagcacgtggcctcttatcgagaaagaaattaccgtcgctcgtgatttgtttgcaaaaagaacaaaactgaaaaaacccagacacgctcgacttcctgtcttcctattgattgcagcttccaatttcgtcacacaacaaggtcctagcgacggctcacaggttttgtaacaagcaatcgaaggttctggaatggcgggaaagggtttagtaccacatgctatgatgcccactgtgatctccagagcaaagttcgttcgatcgtactgttactctctctctttcaaacagaattgtccgaatcgtgtgacaacaacagcctgttctcacacactcttttcttctaaccaagggggtggtttagtttagtagaacctcgtgaaacttacatttacatatatataaacttgcataaattggtcaatgcaagaaatacatatttggtcttttctaattcgtagtttttcaagttcttagatgctttctttttctcttttttacagatcatcaaggaagtaattatctactttttacaacaaatatataca',
                 'cggaccttaatacattcagacacttctgcggtatcaccctacttattcccttcgagattatatctaggaacccatcaggttggtggaagattacccgttctaagacttttcagcttcctctattgatgttacacctggacaccccttttctggcatccagtttttaatcttcagtggcatgtgagattctccgaaattaattaaagcaatcacacaattctctcggataccacctcggttgaaactgacaggtggtttgttacgcatgctaatgcaaaggagcctatatacctttggctcggctgctgtaacagggaatataaagggcagcataatttaggagtttagtgaacttgcaacatttactattttcccttcttacgtaaatatttttctttttaattctaaatcaatctttttcaattttttgtttgtattcttttcttgcttaaatctataactacaaaaaacacatacataaactaaaatataca',
                 'gagcgttggttggtggatcaagcccacgcgtaggcaatcctcgagcagatccgccaggcgtgtatatatagcgtggatggccaggcaactttagtgctgacacatacaggcatatatatatgtgtgcgacgacacatgatcatatggcatgcatgtgctctgtatgtatataaaactcttgttttcttcttttctctaaatattctttccttatacattaggacctttgcagcataaattactatacttctatagacacacaaacacaaatacacacactaaattaatatataca',
                 'ttatattgaattttcaaaaattcttactttttttttggatggacgcaaagaagtttaataatcatattacatggcattaccaccatatacatatccatatctaatcttacttatatgttgtggaaatgtaaagagccccattatcttagcctaaaaaaacctgcgtgggcgtctctttggaactttcagtaatacgcttgcgtgggcgaactgctcattgctatattgaagtccgtgcgtcctcgtcttcaccggtcgcgttcctgaaacgcagatgtgcctaacaataaagattctagcgtgggcgcaatactagcttttatggttatgagcgtgggcgagaggaaaaattggcagtaaccgcgtgggcgtggccccacaaaccttcaaattaacgaatcaaattaacgcgtgggcgaaccataggatgataatgcgattagttttttagccttatttctggggtaattaatcagcgaagcgatgatttttgatctattaacagatatataaatggaaaagctgcataaccactttaactaatactttcaacattttcagtttgtattacttcttattcaaatgtcataaaagtatcaacaaaaaattgttaatatacctctatactttaacgtcaaggagaaaaaactataTTTGCCGCCCAAAGAGCCGAACGAACTACCTTAACA',
                 'ttatattgaattttcaaaaattcttactttttttttggatggacgcaaagaagtttaataatcatattacatggcattaccaccatatacatatccatatctaatcttacttatatgttgtggaaatgtaaagagccccattatcttagcctaaaaaaaccttctctttggaactttcagtaatacgcttaactgctcattgctatattgaagtccgtgcgtcctcgtcttcaccggtcgcgttcctgaaacgcagatgtgcctaacaataaagattctacaatactagcttttatggttatgaagaggaaaaattggcagtaacctggccccacaaaccttcaaattaacgaatcaaattaagcggccgcgtgggcgTTACTCAAGgcgtgggcgtgcgtgggcgggcgtgggcgtgcgtgggcgtctagacaaccataggatgataatgcgattagttttttagccttatttctggggtaattaatcagcgaagcgatgatttttgatctattaacagatatataaatggaaaagctgcataaccactCtaactaCtactGtcaacattCtcagtGtgtattGcttcttattcaaatgtcataCaagtatcaacaaCaaattgttaatatacctctatactGtaacgtcaaggagaaaaaactataTTTGCCGCCCAAAGAGCCGAACGAACTACCTTAACA',
                 'gcgtgggcgaattggtgcgtgggcgccaattggtgcgtgggcgtcgagcagatccgccaggcgtgtatatatagcgtggatggccaggcaactttagtgctgacacatacaggcatatatatatgtgtgcgacgacacatgatcatatggcatgcatgtgctctgtatgtatataaaactcttgttttcttcttttctctaaatattctttccttatacattaggacctttgcagcataaattactatacttctatagacacacaaacacaaatacacacactaaattaataTTTGCCGCCCAAAGAGCCGAACGAACTACCTTAACA')

control.seqs = toupper(control.seqs)
controls = vector(mode='character',length=length(control.names)*2)
for(i in 1:length(control.names)) {
  controls[2*i - 1] = control.names[i]
  # add adapters to control seqs - my designed promoters always have the 'A' in ATG already,
  # these control seqs need it added
  controls[2*i] = paste0(homology[1], control.seqs[i], 'A', homology[2], collapse = '')
}

# don't actually synthesize the controls - cheaper this way
#lines.out = c(lines.out, controls)

writeLines(lines.out, fn.fasta)

# Twist considers two of the Pool 2 (GPD range) 'problematic' and won't synthesize them (rows 142 and 147 in the final table)
# Not a priority to recapitulate this.
# Hand-edit the output file to drop these sequences, altering the numbering (for only 4 sequences: 143 -> 142 through 146 -> 145).
# Save as 'FS9_validation_designs_to_order.fa'
# 1-86 are ZEV, 87-145 are GPD
