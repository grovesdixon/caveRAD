#plot_pi.R


library(tidyverse)
library(cowplot)
library(gtools)
theme_set(theme_cowplot())
rm(list=ls())



# LOAD CAVE DATA ----------------------------------------------------------

ll=load('metadata/cave_data.Rdata')
ll


# READ IN ALLELE FREQUENCIES ----------------------------------------------

#here we get expected heterozygosity h as:
# h = (n / (n-1))(1-sum(pi*pi))
#where n is the number of sequences in the sample, and pi is the frequency of the ith allele at the site
read_in_get_h = function(a.file, species){
  dat = read_tsv(a.file) %>% 
    mutate(homo.sum = p1*p1 + p2*p2,
           h = (N_CHR/(N_CHR-1)) * (1-homo.sum),
           spp=species)
}


n.adat = read_in_get_h('nucleotide_diversity/nesticus/alleleFrequencies.tsv', 'Nesticus')
p.adat = read_in_get_h('nucleotide_diversity/ptomaphagus/alleleFrequencies.tsv', 'Ptomaphagus')
adat = rbind(n.adat, p.adat) %>% 
  rename(cave=pop) %>% 
  left_join(cdat, by='cave')

# CALCULATE NE ------------------------------------------------------------

#basing this on Hahn 2019 page 53
#Under a neutral model, we expect:
#E(pi) = theta = 4Neu
#where E(pi) = the expected average number of nucleotide differences between sequences per site
#pi = the sum of site heterozygosities
#so sum pi from vcftools accross all variants, then divide by the number of sites sequenced


#set up variables for Ne
n.Ntags = 130061 #get this from cdh_alltags_cc.tab
p.Ntags = 264486 #get this from cdh_alltags_cc.tab
n.Nsites = n.Ntags * 26 #bcg1 cut site NOT including 2-mer sticky ends
p.Nsites = p.Ntags * 26
u =  2.8e-9#using 4.65E-9 mutations / bp from (Keightley et al. 2014)


#GET PI BASED ON h (page 44): pi is sum of h accross sites divided by total sites
#AND Ne BASED ON PI (page 53)
ndat = adat %>% 
  group_by(spp, cave, length_m) %>% 
  summarize(sumPi = sum(h, na.rm=TRUE)) %>% 
  mutate(perSitePi = if_else(spp=='Nesticus',
                             sumPi/n.Nsites,
                             sumPi/p.Nsites),
         Ne = perSitePi / (4*u))
  
#min and max
ndat[ndat$perSitePi==min(ndat$perSitePi),]
ndat %>% 
  filter(cave!='all') %>% 
  pull(perSitePi) %>% 
  max()

#plot Ne
ndat %>% 
  ggplot() +
  geom_bar(aes(x=cave, y=Ne, fill=spp), stat='identity', position='dodge')

#plot pi
ndat %>% 
  ggplot() +
  geom_bar(aes(x=cave, y=perSitePi, fill=spp), stat='identity', position='dodge')


#get medians
ndat %>% 
  group_by(spp, cave) %>% 
  summarize(median = median(Ne))
  


#correlate pi with cave length
unlog = function(x){
  return(round(10^x, digits=0))
}

#plotting df
noAll = ndat %>% 
  filter(cave !='all') %>% 
  mutate(logl = log(length_m, 10),
         spp2=if_else(spp=='Nesticus',
                     'N. barri',
                     'P. hatchi'),
         spp2 = factor(spp2, levels=c('N. barri', 'P. hatchi')))

#plotting variables
buffDenom = 10
minl = min(noAll$logl)
maxl = max(noAll$logl)
xbuff= (maxl - minl)/buffDenom
minpi = min(noAll$perSitePi)
maxpi = max(noAll$perSitePi)
ybuff = (maxpi - minpi)/buffDenom
lmpi = lm(noAll$perSitePi~noAll$logl)
summary(lmpi)
YLIMs = c(1.1e-4, 3.55e-4)
r2=round(summary(lmpi)$r.squared, 2)
r2str = paste(r2, '*', sep='')

#build plot
noAll %>% 
  ggplot(aes(x=logl, y=perSitePi)) +
  geom_smooth(se=FALSE, color='black') +
  geom_point(aes(shape=spp2, color=cave), size=6) +
  labs(x='cave length (m)',
       y=expression(pi),
       shape='species',
       subtitle=bquote('R'^2~"="~.(r2str))) +
  scale_x_continuous(breaks = seq(minl, maxl, length.out=4),
                     labels = unlog,
                     limits = c(minl-xbuff, maxl+xbuff)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),
                     limits=YLIMs,
                     breaks=c(seq(1.5e-4, 3.5e-4, 1e-4))) +
  theme(axis.title.y = element_text(angle=0, vjust=0.5, size=20)) 


#PLOT AGAIN WITH NE ESTIMATES INSTEAD
ndat %>% 
  filter(cave !='all') %>% 
  data.frame() %>% 
  mutate(spp=factor(spp, levels=c('Ptomaphagus', 'Nesticus'))) %>% 
  ggplot(aes(x=log(length_m, 10), y=Ne)) +
  geom_smooth(method='lm', se=TRUE, color='black') +
  # geom_smooth(aes(shape=spp), se=FALSE, lwd=0.5, color='black') +
  geom_point(aes(shape=spp, color=cave), size=6) +
  labs(x='cave length (m)', y=bquote(N[e]), shape='species') +
  scale_x_continuous(breaks = c(2.5, 3, 3.5), limits=c(2.5, 3.5), labels = unlog)




# ADD IN THE ANGSD THETA ESTIMATES ----------------------------------------

read_pestPG = function(pestPG_file, cave, spp){
  pestPG = pestPG_file
  read_tsv(pestPG) %>% 
    mutate(tPps = tP/nSites,
           cave = cave,
           spp = spp)
}

npgList = list(read_pestPG('nucleotide_diversity/nesticus/BT_outFold.thetas.idx.pestPG', 'BT', 'Nesticus'),
               read_pestPG('nucleotide_diversity/nesticus/GV_outFold.thetas.idx.pestPG', 'GV', 'Nesticus'),
               read_pestPG('nucleotide_diversity/nesticus/SB_outFold.thetas.idx.pestPG', 'SB', 'Nesticus'),
               read_pestPG('nucleotide_diversity/nesticus/ST_outFold.thetas.idx.pestPG', 'ST', 'Nesticus'),
               read_pestPG('nucleotide_diversity/ptomaphagus/BT_outFold.thetas.idx.pestPG', 'BT', 'Ptomaphagus'),
               read_pestPG('nucleotide_diversity/ptomaphagus/GV_outFold.thetas.idx.pestPG', 'GV', 'Ptomaphagus'),
               read_pestPG('nucleotide_diversity/ptomaphagus/SB_outFold.thetas.idx.pestPG', 'SB', 'Ptomaphagus'),
               read_pestPG('nucleotide_diversity/ptomaphagus/ST_outFold.thetas.idx.pestPG', 'ST', 'Ptomaphagus'))

noAll = purrr::reduce(npgList, rbind) %>% 
  group_by(cave, spp) %>% 
  summarize(angsd.pi = mean(tPps)) %>% 
  left_join(ndat, by = c('spp', 'cave')) %>% 
  mutate(logl = log(length_m, 10),
         spp2=if_else(spp=='Nesticus',
                      'N. barri',
                      'P. hatchi'),
         spp2 = factor(spp2, levels=c('N. barri', 'P. hatchi')),
         angsd.Ne = angsd.pi / (4*u)) %>% 
  arrange(angsd.pi)

#plotting variables
buffDenom = 10
minl = min(noAll$logl)
maxl = max(noAll$logl)
xbuff= (maxl - minl)/buffDenom
minpi = min(noAll$angsd.pi)
maxpi = max(noAll$angsd.pi)
ybuff = (maxpi - minpi)/buffDenom
lmpi = lm(noAll$angsd.pi~noAll$logl)
summary(lmpi)
YLIMs = c(1.1e-4, 3.55e-4)
r2=round(summary(lmpi)$r.squared, 3)
pval=summary(lmpi)$coefficients[2,4]
stars.pval(pval)
starString = stars.pval(pval)[1]
r2str = paste(r2, starString, sep='')

#build plot
noAll %>% 
  ggplot(aes(x=logl, y=angsd.pi)) +
  geom_smooth(se=FALSE, color='black', method='lm') +
  geom_point(aes(shape=spp2, color=cave), size=6) +
  labs(x='cave length (m)',
       y=expression(pi),
       shape='species',
       subtitle=bquote('R'^2~"="~.(r2str))) +
  scale_x_continuous(breaks = seq(minl, maxl, length.out=4),
                     labels = unlog,
                     limits = c(minl-xbuff, maxl+xbuff)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.title.y = element_text(angle=0, vjust=0.5, size=20)) 


#LOOK AT STATS
#for nesticus
nest = noAll %>% 
  filter(spp=='Nesticus') 
nest

lmn = lm(nest$angsd.pi~nest$logl)
summary(lmn)

#for ptomaphagus
ptom = noAll %>% 
  filter(spp=='Ptomaphagus') 
ptom


lmp = lm(ptom$angsd.pi~ptom$logl)
summary(lmp)


#Write out as table for simplicity
noAll %>% 
  select(spp2, cave, length_m, angsd.pi, angsd.Ne, perSitePi, Ne) %>% 
  write_csv(path='nucleotide_diversity/resultsTable.csv')







