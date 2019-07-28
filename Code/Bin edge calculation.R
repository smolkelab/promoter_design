source('D:/Promoter Design Data/Code/FCS file analysis.R')

# directory for raw flow data
setwd('D:/Promoter Design Data/Bin Edges/')
### Needed functions ###

# get equivalents of 'vals', mapping from 'rats.from' to 'rats.to'
map.vals = function(vals, rats.from, rats.to) {
  map.frame = data.frame(from=rats.from, to=rats.to)
  map.model = lm(to ~ from, data = map.frame)
  pred.frame = data.frame(from=vals)
  return(predict(map.model, pred.frame))
}

get.ratios = function(fn, num.pts, mc.min = 0, gfp.min = 0) {
  x = extract.data(load.fromfilename(fn), vals = c('mCherry', 'GFP'))
  x = x[x[,1] > log10(mc.min),]
  x = x[x[,2] > log10(gfp.min),]
  x = log10(x[,2]) - log10(x[,1])
  x = sample(x, num.pts, replace = FALSE)
  x = sort(x, decreasing = TRUE)
  print(length(x))
  return(x)
}

### pZEV experiment ###

ctrl_zev = '1A_008.fcs' # Induced pZEV library, measured with (reference) conditions used for 'Uninduced' sort in pZEV experiment
dat_zev = '1A_001_009.fcs' # Induced pZEV library, measured with conditions used for 'Induced' sort in pZEV experiment

#to_rats = get.ratios(ctrl_zev, 10000, mc.min = 2.8, gfp.min = 3.8, FALSE, FALSE)
#from_rats = get.ratios(dat_zev, 10000, 2.7, 2.6, FALSE, FALSE)
to_rats = get.ratios(ctrl_zev, 10000, mc.min = 2.8, gfp.min = 3.8)
from_rats = get.ratios(dat_zev, 10000, 2.7, 2.6)

pzev.orig = c(-0.39794001,-0.30980392,-0.22184875,-0.13667714,-0.04575749,0.04139269,0.12710480,0.21484385,0.30103000,0.38916608,0.47712125)
pzev.cal = map.vals(pzev.orig, from_rats, to_rats)


### Validation experiment ###
# here we have to calibrate both the uninduced and induced condition, 
# and we have multiple samples collected under the experimental condition for each

### Calibrate Uninduced ###
ctrl_0 = 'Specimen_001_0R_004.fcs'
rats_0 = get.ratios(ctrl_0, 5000)
dat_0 = list('Specimen_001_0A_002.fcs', 'Specimen_001_0B_010.fcs', 'Specimen_001_0C_003.fcs')
output_0 = lapply(dat_0, get.ratios, num.pts=5000)

edges.orig = seq(from = log10(0.04), to = log10(3), length.out = 11)
final.frame = data.frame(orig=edges.orig)
final.frame$Uninduced_A = map.vals(edges.orig, output_0[[1]], rats_0)
final.frame$Uninduced_B = map.vals(edges.orig, output_0[[2]], rats_0)
final.frame$Uninduced_C = map.vals(edges.orig, output_0[[3]], rats_0)
final.frame$Uninduced_Mean = apply(cbind(final.frame$Uninduced_A, final.frame$Uninduced_B, final.frame$Uninduced_C), 1, mean)

### Calibrate induced ###
ctrl_1 = 'Specimen_001_1R_005.fcs'
rats_1 = get.ratios(ctrl_1, 5000)
dat_1 = list('Specimen_001_1A_006.fcs', 'Specimen_001_1B_007.fcs', 'Specimen_001_1C_008.fcs')
output_1 = lapply(dat_1, get.ratios, num.pts=5000)

final.frame$Induced_A = map.vals(edges.orig, output_1[[1]], rats_1)
final.frame$Induced_B = map.vals(edges.orig, output_1[[2]], rats_1)
final.frame$Induced_C = map.vals(edges.orig, output_1[[3]], rats_1)
final.frame$Induced_Mean = apply(cbind(final.frame$Induced_A, final.frame$Induced_B, final.frame$Induced_C), 1, mean)

