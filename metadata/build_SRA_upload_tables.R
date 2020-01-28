#build_SRA_upload_tables.R
library(tidyverse)
rm(list=ls())
bdat = read_tsv('metadata/sample_names.txt', col_names='sample_name')

#make sample titles
snames = bdat$sample_name
stitles = sub('NB', 'Nesticus.barri_', snames)
stitles = sub('PH', 'Ptomaphagus.hatchi_', stitles)
stitles = sub('BT', 'Buggytop_', stitles)
stitles = sub('SB', 'SewaneeBlowhole_', stitles)
stitles = sub('GV', 'Grapevine_', stitles)
stitles = sub('ST', 'SolomonsTemple_', stitles)
bdat$sample_title=stitles

#set up empty accession
bdat$bioproject_accession = ''

#set up organism
bdat$organism = if_else(grepl('^NB', snames),
                   'Nesticus barri',
                   'Ptomaphagus hatchi')


#set up breed as species and cave
bdat$isolate = substr(snames, 5,10)
cave = sapply(stitles, function(x) strsplit(x, '_')[[1]][2])
species = sapply(stitles, function(x) strsplit(x, '_')[[1]][1])
bdat$breed = paste(bdat$organism, cave, sep=' from ')
bdat$host=''
bdat$isolation_source=paste(cave, 'cave')
bdat$collection_date=if_else(cave=="Grapevine",
                             '2018-09',
                             '2018-10')
bdat$geo_loc_name = 'USA: Franklin County TN'
bdat$tissue = if_else(grepl('^NB', snames),
                      'cephalothorax',
                      'entire beetle')




#write out and fill in rest maually
bdat %>% 
  write_csv(path='metadata/Biosample_attributes.csv')
