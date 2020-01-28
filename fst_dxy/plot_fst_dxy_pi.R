#plot_fst_dxy.R


library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
rm(list=ls())



# PLOT WEIGHTED FST BY PSEUDOCHROMOSOME -----------------------------------

read_fst_by_pseudochrom = function(fileName){
  adat = read_tsv(fileName) %>% 
    arrange(shared,
            weighted) %>% 
    mutate(shared = if_else(shared=='shared',
                            'same',
                            'different'),
           shared = factor(shared, levels=c('same', 'different')))
  adat$pair = sub('.', '-', adat$pair, fixed=TRUE)
  adat$pair = factor(adat$pair, levels = unique(adat$pair))
  return(adat)
}
ndat = read_fst_by_pseudochrom('fst_dxy/nesticus/weighted_fst_by_pseudochromosome.tsv')
pdat = read_fst_by_pseudochrom('fst_dxy/ptomaphagus/weighted_fst_by_pseudochromosome.tsv')

plot_overall_watershed = function(fstdat, TITLE=''){
  fstdat %>% 
    ggplot(aes(x=shared,y=weighted, fill=shared)) + 
    geom_boxplot() +
    labs(title=TITLE)
}

nOver = plot_overall_watershed(ndat)
pOver = plot_overall_watershed(pdat)

plot_watershed_by_pair = function(fstdat, TITLE){
  fstdat %>% 
    ggplot(aes(x=pair,y=weighted,fill=shared)) + 
    geom_boxplot() +
    labs(title=TITLE,
         fill='watershed',
         y=bquote(weighted~F[ST])) +
    theme(plot.title = element_text(face='italic',
                                    hjust=0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=10,
                                     angle=10),
          axis.title.y = element_text(angle=90,
                                      vjust=0.5))
}

YLIM = c(0.25, 0.58)
MARG = 3
lpairs = plot_watershed_by_pair(ndat, 'LEGEND')
l = cowplot::get_legend(lpairs)
nPairs = plot_watershed_by_pair(ndat, 'N. barri') + 
  theme(legend.position = 'none') + 
  lims(y=YLIM)
pPairs = plot_watershed_by_pair(pdat, 'P. hatchi') + 
  theme(axis.title.y = element_blank(),
        legend.position = c(0.65, 0.8),
        legend.box.background = element_rect(color="black", size=0.2),
        legend.box.margin = margin(MARG, MARG, MARG, MARG)) +
  lims(y=YLIM)
plts = plot_grid(nPairs, pPairs, rel_widths = c(1,0.96))
xlab = ggdraw() + draw_label('cave pair')
top = plot_grid(plts, xlab, nrow=2, rel_heights = c(14,1))
top



#PLOT PAIRWISE BY CAVE
fstdat = ndat
fstdat = pdat

plot_by_cave = function(fstdat){
  splt1 = fstdat %>% 
    separate(pair, into=c('cave1', 'cave2'))
  splt2 =splt1
  splt2$cave1 = splt1$cave2
  splt2$cave2 = splt1$cave1
  sdat = rbind(splt1,splt2)
  
  sdat %>% 
    ggplot(aes(x=cave1,y=weighted,fill=cave2, color=shared)) +
    geom_boxplot() +
    scale_color_manual(values=c('blue', 'black'))
}

nbyCave = plot_by_cave(ndat)
pbyCave = plot_by_cave(pdat)

plot_grid(nbyCave, pbyCave)


# UPLOAD HARD GENOTYPE CALL DATA -------------------------------------------------------------

#for caves
ll=load('metadata/cave_data.Rdata') #see pair_cave_dat.R in metadata/
pcdat=pcdat %>% 
  mutate(group=paste(cave1,cave2,sep='-')) %>% 
  select(-cave1,-cave2)


#read in datasets (see BASIC STATS WITH VCFTOOLS and DXY in cave_RAD_processing_walkthrough.txt)
LEVELS = c('BT-GV',
           'SB-ST',
           'BT-SB',
           'BT-ST',
           'GV-SB',
           'GV-ST')

read_in = function(fileName, stat, spp){
  d=read_tsv(fileName) %>% 
    filter(group!='ALL_CAVES') %>% 
    mutate(group=sub('_', '-', group),
           stat=stat, 
           spp=spp) %>% 
    left_join(pcdat, by='group') %>% 
    mutate(group=factor(group, levels=LEVELS))
  colnames(d)[3]='value'
  return(d)
}
nf.dat = read_in('fst_dxy/nesticus/all_fst_results.tsv', 'fst', 'n')
pf.dat = read_in('fst_dxy/ptomaphagus/all_fst_results.tsv', 'fst', 'p')
nd.dat = read_in('fst_dxy/nesticus/all_dxy_results.tsv', 'dxy', 'n')
pd.dat = read_in('fst_dxy/ptomaphagus/all_dxy_results.tsv', 'dxy', 'p')
dat = list(nf.dat, pf.dat, nd.dat, pd.dat) %>% 
  purrr::reduce(rbind) 
unique(dat$group)


#SUMMARIZE
#by cave fst
dat %>% 
  group_by(group, spp, stat) %>% 
  summarize(mn = mean(value, na.rm=TRUE)) %>% 
  filter(stat=='fst') %>% 
  mutate(M=1/mn-1) %>% #where M = 4Ne*m; see page 94 in Molecular Population Genetics (Hahn 2018)
  arrange(mn)

#estimate 4NeM
dat %>% 
  group_by(group, spp, stat) %>% 
  summarize(mn = mean(value, na.rm=TRUE)) %>% 
  filter(stat=='fst') %>% 
  arrange(mn) %>% 
  mutate()


#by cave dxy
dat %>% 
  group_by(group, spp, stat) %>% 
  summarize(mn = mean(value, na.rm=TRUE))

#pairwise only
dat %>% 
  filter(group != "ALL-CAVES") %>% 
  group_by(spp, stat) %>% 
  summarize(mn = mean(value, na.rm=TRUE))


#plot using facet grid
dat %>% 
  filter(group != "ALL-CAVES") %>% 
  mutate(stat=factor(stat, levels=c('fst', 'dxy'))) %>% 
  ggplot(aes(x = group, y = value)) + 
  facet_grid(stat~spp) + 
  geom_boxplot(aes(fill=group)) +
  labs(y='', x='cave')


#BUILD CUSTOM FACET GRID
plot_box = function(dat, ycol, ylab, ylim){
  dat %>% 
    filter(group!='ALL-CAVES') %>% 
    ggplot(aes_string(x='group', y='value', fill='group')) +
    geom_boxplot() +
    labs(y=ylab) +
    lims(y=ylim) +
    theme(legend.position = 'none',
          axis.title.x = element_blank())
}


build_spp_panels = function(f.plt, d.plt, main, labels){
  plts = plot_grid(f.plt, d.plt, nrow=2)
  main = ggdraw() + draw_label(main, fontface='italic')
  r = plot_grid(main, plts, nrow=2, rel_heights = c(1,12))
  return(r)
}


#build individual plots
YLIM = c(0,1)
nf = plot_box(nf.dat, 'value', bquote(F[ST]), YLIM) +
  theme(axis.text.x = element_blank())
pf = plot_box(pf.dat, 'value', bquote(F[ST]), YLIM) + 
  theme(axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank(),
  axis.text.x = element_blank())
nd = plot_box(nd.dat, 'value', bquote(d[XY]), YLIM)
pd = plot_box(pd.dat, 'value', bquote(d[XY]), YLIM) + 
  theme(axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank())

#assemble
l=cowplot::get_legend(nf + theme(legend.position='right'))
n.pan = build_spp_panels(nf, nd, 'Nesticus', c('A','C'))
p.pan = build_spp_panels(pf, pd, 'Ptomaphagus', c('B','D'))
pans = plot_grid(n.pan, p.pan, nrow=1, rel_widths=c(1,0.8))
xlab = ggdraw() + draw_label('cave pair')
final = plot_grid(pans,xlab,nrow=2,rel_heights = c(12,1))
final


# ASSEMBLE TABLE OF DIVERGENCE ESTIMATES ----------------------------------

#GET THE VCFTOOLS RESULTS
#for nesticus
nhard = dat %>% 
  filter(spp=='n') %>% 
  group_by(group, stat) %>% 
  summarize(mn = mean(value, na.rm=TRUE)) %>% 
  arrange(mn) %>% 
  pivot_wider(names_from = stat,
              values_from = mn)

#ptomaphagus
phard = dat %>% 
  filter(spp=='p') %>% 
  group_by(group, stat) %>% 
  summarize(mn = mean(value, na.rm=TRUE)) %>% 
  arrange(mn) %>% 
  pivot_wider(names_from = stat,
              values_from = mn)
  
#GET THE ANGSD RESULTS
#for nesticus
nangsd = read.table('fst_dxy/nesticus/overall_fst.tsv',
                    header = TRUE,
                    stringsAsFactors = FALSE) %>% 
  mutate(group=sub('.', '-', pair, fixed=TRUE)) %>% 
  select(group, weighted, unweighted) %>% 
  as_tibble()

#for ptomaphagus
pangsd = read.table('fst_dxy/ptomaphagus/overall_fst.tsv',
                    header = TRUE,
                    stringsAsFactors = FALSE) %>% 
  mutate(group=sub('.', '-', pair, fixed=TRUE)) %>% 
  select(group, weighted, unweighted) %>% 
  as_tibble()

#MERGE THEM

#nesticus
ndist = nangsd %>% 
  full_join(nhard, by = 'group')
ndist
ndist %>% 
  write_csv('fst_dxy/nesticus/resultsTable.csv')

#ptom
pdist = pangsd %>% 
  full_join(phard, by = 'group')
pdist
pdist %>% 
  write_csv('fst_dxy/ptomaphagus/resultsTable.csv')



#CHECK AGREEMENT
ndist$species = 'N. barri'
pdist$species = 'P. hatchi'
ddat = rbind(ndist,
             pdist)


ddat %>% 
  pivot_longer(unweighted:fst) %>% 
  ggplot(aes(x=weighted,
             y=value,
             color=species)) +
  geom_smooth(method='lm',
              se=FALSE) +
  geom_point() +
  facet_wrap(~ name, scales = "free") +
  labs(x=bquote(weighted~F[ST]),
       y='')


#barplots of weighted
barplot_fst = function(fstdat, diffStat){
  share = c('BT-GV',
            'SB-ST')
  splt1 = fstdat %>% 
    mutate(watershed = if_else(group %in% share,
                            'same',
                            'different')) %>% 
    separate(group, into=c('cave1', 'cave2'))
  splt2 =splt1
  splt2$cave1 = splt1$cave2
  splt2$cave2 = splt1$cave1
  sdat = rbind(splt1,splt2)
  
  sdat %>% 
    ggplot(aes_string(x='cave1',y=diffStat,fill='cave2', color='watershed')) +
    geom_bar(stat='identity', position='dodge') +
    scale_color_manual(values=c('grey', 'black'))
}

nweight = barplot_fst(ndist, 'weighted') + 
  labs(title='N. barri', x=bquote(weighted~F[ST])) +
  theme(plot.title = element_text(face='italic',hjust=0.5),
        legend.position = 'right')
pweight = barplot_fst(pdist, 'weighted') + 
  labs(title='P. hatchi', x=bquote(weighted~F[ST])) +
  theme(plot.title = element_text(face='italic',hjust=0.5))
nfst = barplot_fst(ndist, 'fst') +
  labs(y='Weir')
pfst = barplot_fst(pdist, 'fst') +
  labs(y='Weir')
ndxy = barplot_fst(ndist, 'dxy') +
  labs(y=bquote(d[XY]))
pdxy = barplot_fst(pdist, 'dxy') +
  labs(y=bquote(d[XY]))

l=cowplot::get_legend(nweight)
pltList = list(nweight,pweight,nfst,pfst,ndxy,pdxy)
modList = lapply(pltList, function(x) x+theme(legend.position = 'none',
                                              axis.title.x = element_blank()))
plts=plot_grid(plotlist = modList, nrow=3)
xlab = ggdraw() + draw_label('cave1')
top=plot_grid(plts, xlab, nrow=2, rel_heights=c(20,1))
plot_grid(top,l, nrow=1, rel_widths = c(6,1))

# PLOT AGAINST CAVE DATA --------------------------------------------------
#BY WATERSHED
#fst
dat %>% 
  filter(stat=='fst') %>% 
  mutate(watershed = if_else(ws1==ws2,
                             'Same',
                             'Different')) %>% 
  ggplot() +
  geom_boxplot(aes(x=spp, y=value, fill=watershed))

#dxy
dat %>% 
  filter(stat=='dxy') %>% 
  mutate(watershed = if_else(ws1==ws2,
                             'Same',
                             'Different')) %>% 
  ggplot() +
  geom_boxplot(aes(x=spp, y=value, fill=watershed))

#CAVE DISTANCE
sdat = dat %>% 
  group_by(group,spp, stat) %>% 
  summarize(dist=mean(dist),
            mnStat = mean(value, na.rm=TRUE))

#fst
sdat %>% 
  filter(stat=='fst') %>% 
  ggplot(aes(x=dist,y=mnStat,color=group, shape=spp)) +
  geom_point(size=5) +
  geom_smooth(se=FALSE)

#dxy
sdat %>% 
  filter(stat=='dxy') %>% 
  ggplot(aes(x=dist,y=mnStat,color=group, shape=spp)) +
  geom_point(size=5) +
  geom_smooth(se=FALSE)



# GET DXY FROM ANSD ALLELE FREQUENCIES ------------------------------------
#didn't end up using this, but keeping for reference

pops = c('BT','GV','SB','ST')
spp = 'nesticus'

#Funciton to read in and merge minor allele frequency data from each population
#these maf files are generated under the ANGSD DEMOGRAPHIC ANALYSIS section in cave_RAD_processing_walkthrough.txt
read_in_mafs = function(pops, spp){
  mafList = list()
  for (ipop in pops){
    print(paste('Reading in ', ipop, '...', sep=''))
    inFile = paste('large_ignored/', spp, '/', ipop, '.mafs', sep='')
    pmaf = read.table(inFile, header=TRUE) %>% 
      dplyr::select(chromo,position,minor,ref,anc,knownEM) %>% 
      set_names(c('chromo',
                  'position',
                  paste(ipop, 'minor', sep='.'),
                  'ref',
                  'anc',
                  paste(ipop, 'freq', sep='.')))
    mafList[[ipop]] = pmaf
  }
  mafs = purrr::reduce(mafList, left_join, by=c('chromo','position','ref','anc'))
  return(mafs)
}

nmaf = read_in_mafs(pops, 'nesticus') %>% 
  na.omit()
pmaf = read_in_mafs(pops, 'ptomaphagus') %>% 
  na.omit()



get_ansd_dxy = function(mafDat, pops){
  dxyRes = data.frame()
  for (c1 in pops){
    otherPops = pops[pops != c1]
    for (c2 in otherPops){
      pair=paste(c1,c2,sep='-')
      print(paste('Running dxy for pair ',pair,'...', sep=''))
      ma.x = mafDat[,paste(c1,'minor',sep='.')] #minor allele for population x
      maf.x = mafDat[,paste(c1,'freq',sep='.')] 
      maj.x = 1-maf.x
      ma.y = mafDat[,paste(c2,'minor',sep='.')]
      maf.y = mafDat[,paste(c2,'freq',sep='.')]
      maj.y = 1-maf.y
      k = ma.x==ma.y
      
      dxyDat = data.frame(ma.x,
                          maf.x,
                          maj.x,
                          ma.y,
                          maf.y,
                          maj.y,
                          k)
      dxyDat$dxy = 
        maf.x*maf.y*as.numeric(!k) +  #multiply by 0 if minor allele is the same or 1 if different
        maf.x*maj.y*as.numeric(k) +  #multiply by 1 if minor allele is the same or 0 if different
        maj.x*maj.y*as.numeric(!k) +
        maj.x*maf.y*as.numeric(k)
      mnDxy = mean(dxyDat$dxy,
                   na.rm=TRUE)
      subRes = data.frame(group=pair,
                          dxy = mnDxy)
      dxyRes = rbind(dxyRes, subRes)
    }
  }
  return(dxyRes)
}

ndxy = get_ansd_dxy(nmaf, pops)
pdxy = get_ansd_dxy(pmaf, pops)

