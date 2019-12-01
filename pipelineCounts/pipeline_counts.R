#pipeline_counts.R
library(tidyverse)
library(plotrix)


# RAW ---------------------------------------------------------------------

raw = read_tsv('pipelineCounts/raw_read_counts.tsv',
               col_names = c('pool', 'value', 'stat'))
tot.raw = sum(raw$value)
tot.raw
paste(tot.raw/1e6, 'million reads')



# PER SAMPLE --------------------------------------------------------------

trim = read_tsv('pipelineCounts/trimmed_read_counts.txt',
                col_names = c('sample', 'value', 'stat')) %>% 
  mutate(spp=substr(sample, start=1, stop=2),
         cave = substr(sample, start=3, stop=4))

#total trim
tot.trim = sum(trim$value)
tot.trim / 1e6

tot.trim/tot.raw

#summarize
trim %>% 
  group_by(spp, cave) %>% 
  summarize(mnCount = mean(value)/1e6,
            seCount = std.error(value)/1e6)

#summarize
trim %>% 
  group_by(spp) %>% 
  summarize(mnCount = mean(value)/1e6,
            seCount = std.error(value)/1e6)

#stats
nb = trim %>% 
  filter(spp=='NB') 
lm.n = lm(nb$value ~ nb$cave)
summary(lm.n)

fit <- aov(value ~ cave, data=nb)
summary(fit)




# SAMPLE SIZES ------------------------------------------------------------

ns = read.table('ibs_pca_admix/nesticus/bams') %>% 
  mutate(cave=substr(V1, start=3, stop=4),
         spp=substr(V1, start=1, stop=2))

ps = read.table('ibs_pca_admix/ptomaphagus/bams') %>% 
  mutate(cave=substr(V1, start=3, stop=4),
         spp=substr(V1, start=1, stop=2))

res=rbind(ns, ps) %>% 
  group_by(spp, cave) %>% 
  summarize(N=n())
res

res %>% 
  write_tsv(path='metadata/sampleSizes.tsv')










