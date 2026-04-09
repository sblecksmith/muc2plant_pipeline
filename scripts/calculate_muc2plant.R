# title: calculate_muc2plant.R
# author: Sarah Blecksmith
# purpose: calculate the mucin to plant cazyme ratio using the unrounded aggregated CAZyme family table
# and the families unique to one substrate output from CAZy_family_substrates.R uses substrate annotation scheme 
# from DOI: 10.1126/science.aan4834


library(tidyverse)
library(argparse)

parser$add_argument('--input', '-i', help= 'This is the file of aggregated cazyme abundances from dbcan')
parser$add_argument('--families', '-i', help= 'This is the file of cazyme substrate annotations from Smits, et al.')
parser$add_argument('--output', '-o', help= 'This is the output file containing muc2plant for each sample')


xargs<- parser$parse_args()

cazymes <- read.table(xargs$input)
muc2plant_file <- read.table(xargs$output)
all_substrates_sep <- read.table(xargs$families)

# Unique families
# Collect the plant families
plant <- filter(all_substrates_sep, substrate == 'plant')
plant_fams <- data.frame()
for (fam in plant$cazyme_family){
  fams <- data.frame()
  # for family GH11, we want to catch all the subfamilies but not
  # GH114 or GH117.  So this grep matches on GH11 alone and all the GH_X
  fams <- filter(cazymes, grepl(paste0("^", fam, "$|^", fam, "_"), Family))
  plant_fams <- rbind(plant_fams, fams)
}

# Collect the mucin families
mucin <- filter(all_substrates_sep, substrate == 'mucin')
mucin_fams <- data.frame()
for (fam in mucin$cazyme_family){
  fams <- data.frame()
  # for family GH11, we want to catch all the subfamilies but not
  # GH114 or GH117.  So this grep matches on GH11 alone and all the GH_X
  fams <- filter(cazymes, grepl(paste0("^", fam, "$|^", fam, "_"), Family))
  mucin_fams <- rbind(mucin_fams, fams)
}

# Sum the mucin and plant families and make new dataframe
plant_sum <- data.frame(metagenome_id = colnames(plant_fams[,-1]),plant_family_total=colSums(plant_fams[, -1]))
mucin_sum <- data.frame(metagenome_id = colnames(mucin_fams[,-1]),mucin_family_total=colSums(mucin_fams[, -1]))


merged = plant_sum %>%
  left_join(mucin_sum, by = "metagenome_id") 

# calculate unique muc2plant
merged$muc2plant_unique <- merged$mucin_family_total/merged$plant_family_total

write.csv(merged, muc2plant_file, row.names = FALSE)

