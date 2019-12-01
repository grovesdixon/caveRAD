#rename_caveRAD.R

library('tidyverse')
bcFile = 'metadata/barcoding_table.csv'
bcdat = read_csv(bcFile) %>% 
  mutate(barcodePair=paste(ligBarcodeSeq,PCRbarcode,sep='_'))
bcdat
    


infileName='metadata/trimmedFiles.txt'
fdat = read_tsv(infileName, col_names=c('fileName')) %>% 
  mutate(fn=sub('Pool-', '', fileName),
         fn=sub('.trim', '', fn)) %>% 
  separate(fn, into=c('PCRbarcode', 'sNum', 'ligBarcodeSeq')) %>% 
  mutate(barcodePair = paste(ligBarcodeSeq,PCRbarcode,sep='_'))
fdat

mdat = bcdat %>% 
  left_join(fdat, by = 'barcodePair') %>% 
  select(sample, fileName)
nrow(bcdat)
nrow(fdat)
nrow(mdat)
sum(!is.na(mdat$fileName))
length(unique(mdat$fileName))
  
nrow(na.omit(mdat))

write_tsv(na.omit(mdat), path='metadata/mergedBCs.tsv', col_names=F)
