#!/usr/bin/env Rscript
#angsd_fst_summary.R
#gather summary data from Angsd command: realSFS fst print

library(tidyverse)
library(plotrix)



# READ IN AND CALCULATE PER SITE FST VALUES -------------------------------
#the per site are the unweighted Fst values (http://www.popgen.dk/angsd/index.php/Fst)
#see section ANGSD ALLELE FREQUENCIES AND FST in cave_RAD_processing_walkthrough.txt

args = commandArgs(trailingOnly=TRUE)
abFile = args[1] #probably allAB.tsv
print(paste('Reading in infile ', abFile, '...', sep=''))
abdat = read_tsv(abFile)

#pairs that share a watershed
same = c('BT.GV',
         'SB.ST')

#get AB ratio
print('Calculating A/B ratio and calling watershed sharing...')
dat = abdat %>% 
  mutate(ABratio = A/B,
         shared = if_else(pair %in% same,
                          'shared',
                          'unshared'))


# GATHER BOXPLOT STATS FOR THE PERSITE FST VALUES -------------------------
#the mean of all these values together is the unweighted Fst
#so these boxplots summarize the distribution of unweighted Fst accross the dataset

boxList = list()
for (ip in unique(dat$pair)){
  print(paste(ip,'...', sep=''))
  sub = dat %>% 
    filter(pair==ip)
  boxstats = boxplot.stats(sub$ABratio, coef = 1.5, do.conf = TRUE, do.out = TRUE)$stat
  print(boxstats)
  boxList[[ip]]=boxstats
}
b=data.frame(t(data.frame(boxList)))
boxCols = c('ymin', 'lower', 'middle', 'upper', 'ymax')
colnames(b)=boxCols
b$pair = rownames(b)
b

#write out boxplot data
b %>% 
  write_tsv(path='./unweighted_fst_boxStats.tsv')



# RE-CALCULATE THE OVERALL RESULTS AS SANITY CHECK ------------------------
#these should match the results given by realSFS fst stats

rdat = dat %>% 
  group_by(pair) %>% 
  filter(!is.na(ABratio)) %>% 
  summarize(wieghted = sum(A)/sum(B),
            unwieghted = mean(ABratio),
            uw.med = median(ABratio),
            uw.sd = sd(ABratio),
            uw.stderr = std.error(ABratio))

#write out
rdat %>% 
  write_tsv(path='./overall_fsts_from_AB.tsv')



# GET STATS BY SHARED/UNSHARED WATERSHED ----------------------------------

#FIRST GET OVERALL
wres = dat %>% 
  group_by(shared) %>% 
  summarize(weighted = sum(A)/sum(B),
            unweighted = mean(ABratio, na.rm=TRUE))

wres %>% 
  write_tsv(path='watershed_overall_fst_fromAB.tsv')

#GET BOXPLOT STATS FOR UNWEIGHTED 
sdat = dat %>% 
  filter(shared=='shared')
udat = dat %>% 
  filter(shared=='unshared')

uwShared = boxplot.stats(sdat$ABratio, coef = 1.5, do.conf = TRUE, do.out = TRUE)$stat
uwUnshared = boxplot.stats(udat$ABratio, coef = 1.5, do.conf = TRUE, do.out = TRUE)$stat

wsBs = rbind(uwShared,
             uwUnshared) %>% 
  data.frame()
colnames(wsBs) = boxCols

wsBs %>% 
  write_tsv(path='watershed_unweighted_fst_boxStats.tsv')




# GET WEIGHTED FST ACCROSS PSEUDOCHROMOSOMES ------------------------------

#get overall
print('Gettting weighted Fst accross pseudochromosomes')
wChr = dat %>% 
  group_by(chr, pair, shared) %>% 
  summarize(weighted = sum(A)/sum(B))

wChr %>% 
  write_tsv('weighted_fst_by_pseudochromosome.tsv')

# 
# print('Gettting weighted Fst accross pseudochromosomes split by watershed')
# wChrWs = dat %>% 
#   group_by(chr, shared) %>% 
#   summarize(weighted = sum(A)/sum(B))
# 
# wChrWs %>% 
#   write_tsv('watershed_weighted_fst_by_pseudochromosome.tsv')
# 
# 
# #then get by cave
# print('Getting weighted Fst by cave and shared watershed')
# caves = c('BT','GV', 'SB', 'ST')
# caveRes = list()
# for (icave in caves){
#   caveSub = dat %>% 
#     filter(grepl(icave, pair))
#   wsChr = caveSub %>% 
#     group_by(chr, shared) %>% 
#     summarize(weighted = sum(A)/sum(B)) %>% 
#     mutate(cave=icave)
#   caveRes = rbind(caveRes, wsChr)
# }
# caveRes
# 
# caveRes %>% 
#   write_tsv(path='cave_watershed_weighted_fst_by_pseudochromosome.tsv')

# #GET WEIGHTED FST ACCROSS LOCI (this takes too long keeping commented for refernce)
# #this is following equation 4 given Fumagalli et al. 2013
# #use the cdh_alltags_cc.tab file output during locus generation
# print('Reading in locus coordinates...')
# cdh = read_tsv('cdh_alltags_cc.tab',
#                col_names = c('tag', 'chr', 'start', 'end'))
# cdh
# utags = cdh$tag
# ntags = length(utags)
# uchrs = unique(cdh$chr)
# 
# 
# resList = list()
# for (ichr in uchrs){
#   csub = dat %>% 
#     filter(chr==ichr)
#   tsub = cdh %>% 
#     filter(chr==ichr) %>% 
#     data.frame()
#   rownames(tsub)=tsub$tag
#   ctags = tsub$tag
#   print('----')
#   print(paste(ichr,'...', sep=''))
#   print(paste(length(ctags), 'tags'))
#   for (itag in ctags){
#     start = tsub[itag,'start']
#     end = tsub[itag,'end']
#     tagSub = csub %>% 
#       filter(pos >= start & pos <= end)
#     tagRes = tagSub %>% 
#       group_by(shared) %>% 
#       summarize(weighted = sum(A)/sum(B)) %>% 
#       pull(weighted)
#     if (length(tagRes) < 2){
#       next
#     }
#     names(tagRes) = c('shared', 'unshared')
#     resList[[itag]]=tagRes
#   }
# }
# resDf = data.frame(t(data.frame(resList)))
# 
# weighted.shared.bs = boxplot.stats(resDf$shared, coef = 1.5, do.conf = TRUE, do.out = TRUE)$stat
# weighted.unshared.bs = boxplot.stats(resDf$unshared, coef = 1.5, do.conf = TRUE, do.out = TRUE)$stat
# 
# save(unweighted.shared.bs,
#      unweighted.unshared.bs,
#      weighted.shared.bs,
#      weighted.unshared.bs,
#      file='shared_and_unshared_boxplotStats.Rdata')

# 
# 
# 
# # GET RATIOS FOR SHARED AND UNSHARED WATERSHEDS ---------------------------
# 
# caves = c('BT',
#           'GV',
#           'SB',
#           'ST')
# same = c('BT.GV',
#          'SB.ST')
# 
# 
# allRatios = c()
# rboxList = list()
# for (c in caves){
#   print('------------------')
#   print(paste(c, '...', sep=''))
#   csub = dat %>% 
#     filter(grepl(c, pair))
#   wsub = csub %>% 
#     select(-A, -B) %>% 
#     pivot_wider(names_from = pair, values_from = ABratio) %>% 
#     data.frame()
#   shared = wsub[,colnames(wsub) %in% same]
#   unshared = wsub[,!colnames(wsub) %in% same]
#   print('subset:')
#   print(head(wsub))
#   print(head(shared))
#   print('unshared:')
#   print(head(unshared))
#   r1 = shared/unshared[,3]
#   r2 = shared/unshared[,4]
#   bs = boxplot.stats(append(r1,r2), coef = 1.5, do.conf = TRUE, do.out = TRUE)$stat
#   allRatios = append(allRatios, append(r1,r2))
#   print(bs)
#   rboxList[[c]] = bs
# }
# res=data.frame(t(data.frame(rboxList)))
# colnames(res)=c('ymin', 'lower', 'middle', 'upper', 'ymax')
# res$cave = rownames(res)
# all = boxplot.stats(allRatios, coef = 1.5, do.conf = TRUE, do.out = TRUE)$stat
# allRow = append(all, 'all')
# res2 = rbind(res, allRow)
# 
# #write out
# res2 %>% 
#   write_tsv(path='fst_waterShedRatio.tsv')
#   
# 
# 
# 
# 
# 
# # REPEAT WITHOUT GRAPEVILLE  ----------------------------------------------
# 
# caves=c('SB', 'ST')
# dat = dat %>% 
#   filter(!grepl('GV',pair))
# 
# rboxList = list()
# allRatios = c()
# for (c in caves){
#   print('------------------')
#   print(paste(c, '...', sep=''))
#   csub = dat %>% 
#     filter(grepl(c, pair))
#   wsub = csub %>% 
#     select(-A, -B) %>% 
#     pivot_wider(names_from = pair, values_from = ABratio) %>% 
#     data.frame()
#   shared = wsub[,colnames(wsub) %in% same]
#   unshared = wsub[,!colnames(wsub) %in% same]
#   print('subset:')
#   print(head(wsub))
#   print(head(shared))
#   print('unshared:')
#   print(head(unshared))
#   r1 = shared/unshared[,3]
#   # r2 = shared/unshared[,4]
#   bs = boxplot.stats(r1, coef = 1.5, do.conf = TRUE, do.out = TRUE)$stat
#   print(bs)
#   rboxList[[c]] = bs
#   # allUnshared = append(unshared[,1], unshared[,2])
#   # mnShared = mean(shared, na.rm=TRUE)
#   # mnUnshared = mean(allUnshared, na.rm=TRUE)
#   # sharedArr = append(sharedArr, mnShared)
#   # unsharedArr = append(unsharedArr, mnUnshared)
#   # ratio = mnShared/mnUnshared
#   # print(c(mnShared, mnUnshared,  ratio))
# }
# 
# res=data.frame(t(data.frame(rboxList)))
# colnames(res)=c('ymin', 'lower', 'middle', 'upper', 'ymax')
# res$cave = rownames(res)
# res
# 
# #write out
# res %>% 
#   write_tsv(path='fst_waterShedRatio_NOGRAPE.tsv')
# 
