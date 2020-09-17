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
lmpi = lm(angsd.pi~logl, data=noAll)
summary(lmpi)
YLIMs = c(1.1e-4, 3.55e-4)
r2=round(summary(lmpi)$r.squared, 3)
pval=summary(lmpi)$coefficients[2,4]
stars.pval(pval)
starString = stars.pval(pval)[1]
# starString = '' #get rid of stars to show we're not doing a statistical test
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


#run multiple regression
lmm = lm(noAll$angsd.pi ~ noAll$spp + noAll$length_m + noAll$cave)
summary(lmm)

#LOOK AT STATS AND ASSUMPTIONS FOR REGRESSION

#build boxplot
boxplot(noAll$angsd.pi ~ noAll$spp,
        xlab='species', ylab='Pi')
boxplot(noAll$logl,ylab='log cave length')
#no points flagged as outliers

#check normality
#individually
noAll %>% 
  ggplot(aes(x=angsd.pi)) +
  geom_density() +
  scale_x_continuous(breaks = c(0, 0.0015, 0.003), limits = c(0,0.004)) +
  facet_wrap(~spp)

#together
noAll %>% 
  ggplot(aes(x=angsd.pi)) +
  geom_density() +
  lims(x=c(0,0.004))

#length
noAll %>% 
  ggplot(aes(x=logl)) +
  geom_density() +
  lims(x=c(1,5))

#normality of residuals for length alone
noAll$pred = predict(lmpi)
noAll$resids = noAll$angsd.pi - noAll$pred
noAll %>% 
  ggplot(aes(x=resids)) +
  geom_density() +
  lims(x=c(-8e-4, 8e-4))

#normality of residuals for length alone
noAll$pred = predict(lmm)
noAll$resids = noAll$angsd.pi - noAll$pred
noAll %>% 
  ggplot(aes(x=resids)) +
  geom_density() +
  lims(x=c(-8e-4, 8e-4))


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

#look at plot without log transform
#build plot
noAll %>% 
  ggplot(aes(x=length_m, y=angsd.pi)) +
  geom_smooth(se=FALSE, color='black') +
  geom_point(aes(shape=spp2, color=cave), size=6) +
  labs(x='cave length (m)',
       y=expression(pi),
       shape='species') +
  # scale_x_continuous(breaks = seq(minl, maxl, length.out=4),
  #                    labels = unlog,
  #                    limits = c(minl-xbuff, maxl+xbuff)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.title.y = element_text(angle=0, vjust=0.5, size=20)) 

#get observed correlations
n = noAll %>% 
  filter(spp=='Nesticus')
p = noAll %>% 
  filter(spp == 'Ptomaphagus')
n_actual = cor(n$angsd.pi, n$length_m)
p_actual = cor(p$angsd.pi, p$length_m)
all_actual = cor(noAll$angsd.pi, noAll$length_m)


# do some additional stats ------------------------------------------------
#we need extra confidence because there are so few data points.
require(gtools)

#PERMUTATIONS OF LENGTHS
#nesticus
target = 'Nesticus'
pi = noAll %>% 
  filter(spp == target)
n_actual = cor(pi$angsd.pi, pi$length_m)
perms = data.frame(t(permutations(n = 4, r = 4, v = pi$length_m)))
cors = data.frame(x=t(cor(pi$angsd.pi, perms))) 
cors %>% 
  ggplot(aes(x=x)) +
  geom_density() +
  geom_vline(xintercept = n_actual, color='red')
sum(cors$x > n_actual) / nrow(cors)

#ptomaphagus
target = 'Ptomaphagus'
pi = noAll %>% 
  filter(spp == target)
p_actual = cor(pi$angsd.pi, pi$length_m)
perms = data.frame(t(permutations(n = 4, r = 4, v = pi$length_m)))
cors = data.frame(x=t(cor(pi$angsd.pi, perms))) 
cors %>% 
  ggplot(aes(x=x)) +
  geom_density() +
  geom_vline(xintercept = p_actual, color='red')
sum(cors$x > p_actual) / nrow(cors)


# simulate resampling -----------------------------------------------------
#simulate cave lenths and Pi estimates as normal curves based on observed means and sds
#randomly sample from these to see how improbable our correlations were

set.seed(123)

#set up simulated curve of cave lengths
cave_lengths = noAll %>% 
  filter(spp == target) %>% 
  pull(length_m)
norm_length = rnorm(1000, mean = mean(cave_lengths), sd = sd(cave_lengths))
norm_length_dens = density(norm_length, n=10000)

#SET UP SIMULATED CURVES OF PI ESTIMATES FOR EACH SPECIES

simulate_pi = function(target){
  #nesticus
  obs_pi = noAll %>% 
    filter(spp == target) %>% 
    pull(angsd.pi)
  norm_pi <- rnorm(1000, mean = mean(obs_pi), sd = sd(obs_pi))
  obs_dens = density(obs_pi, n=10000)
  norm_dens = density(norm_pi, n=10000)
  
  #trace observed and simulated densities
  d=data.frame(obs_x = obs_dens$x,
               obs_y = obs_dens$y,
               norm_x = norm_dens$x,
               norm_y = norm_dens$y)
  plt=d %>% 
    ggplot() +
    geom_line(aes(x=obs_x, y=obs_y)) +
    geom_line(aes(x=norm_x, y=norm_y), color='red', lty=2) +
    labs(x='Pi', y='density', title=target)
  plot(plt)
  return(norm_pi)
}

#nesticus
n_norm_pi = simulate_pi('Nesticus')
p_norm_pi = simulate_pi('Ptomaphagus')


#sample from hypothetical and get potential correlations with set cave lengths
n_iter = 10000
n_test_cors = c()
p_test_cors = c()
for (i in 1:n_iter){
  n_test_pi = sample(n_norm_pi, 4)
  p_test_pi = sample(p_norm_pi, 4)
  n_test_cors = append(n_test_cors, cor(n_test_pi, cave_lengths))
  p_test_cors = append(p_test_cors, cor(p_test_pi, cave_lengths))
}

#plot
plot_sim = function(test_cors, actual_pi){
  res = tibble(i = 1:n_iter,
               cor = test_cors)
  plt=res %>% 
    ggplot() +
    geom_density(aes(x=cor)) +
    geom_vline(xintercept = actual_pi)
  plot(plt)
  return(res)
}
n_res = plot_sim(n_test_cors, n_actual)
p_res = plot_sim(p_test_cors, p_actual)

#stat
sum(n_res$cor >= n_actual)
sum(n_res$cor >= n_actual) / nrow(n_res)
sum(p_res$cor >= p_actual)
sum(p_res$cor >= p_actual) / nrow(p_res)

#repeat for hypothetical cave lengths
n_iter = 10000
n_test_cors = c()
p_test_cors = c()
for (i in 1:n_iter){
  test_length = sample(norm_length, 4)
  n_test_pi = sample(n_norm_pi, 4)
  p_test_pi = sample(p_norm_pi, 4)
  n_test_cors = append(n_test_cors, cor(n_test_pi, test_length))
  p_test_cors = append(p_test_cors, cor(p_test_pi, test_length))
}


#plot
n_res = plot_sim(n_test_cors, n_actual)
p_res = plot_sim(p_test_cors, p_actual)

#stat
sum(n_res$cor >= n_actual)
sum(n_res$cor >= n_actual) / nrow(n_res)
sum(p_res$cor >= p_actual)
sum(p_res$cor >= p_actual) / nrow(p_res)

#getting positive correlation
sum(n_res$cor >= 0.5)
sum(n_res$cor >= 0.5) / nrow(n_res)
sum(p_res$cor >= 0.5)
sum(p_res$cor >= 0.5) / nrow(p_res)


#PERFORM THE DOUBLE SIMULATION, AS FOR THE TWO SPECIES
n_iter = 10000
test_cors = c()
for (i in 1:n_iter){
  test_length = sample(norm_length, 4)
  n_test_pi = sample(n_norm_pi, 4)
  p_test_pi = sample(p_norm_pi, 4)
  c_length = append(test_length, test_length)
  c_pi = append(n_test_pi, p_test_pi)
  test_cors = append(test_cors, cor(c_pi, c_length))
}

#plot combined
res = plot_sim(test_cors, actual)

#stat
sum(res$cor >= actual)
sum(res$cor >= actual) / nrow(res)
