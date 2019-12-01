#plot_fst_dxy.R


library(tidyverse)
library(cowplot)
rm(list=ls())


#read in datasets (see BASIC STATS WITH VCFTOOLS and DXY in cave_RAD_processing_walkthrough.txt)
read_in = function(fileName, stat, spp){
  d=read_tsv(fileName) %>% 
    mutate(group=sub('_', '-', group),
           stat=stat, 
           spp=spp)
  colnames(d)[3]='value'
  return(d)
}
nf.dat = read_in('mpileup_analyses/nesticus/all_fst_results.tsv', 'fst', 'n')
pf.dat = read_in('mpileup_analyses/ptomaphagus/all_fst_results.tsv', 'fst', 'p')
nd.dat = read_in('mpileup_analyses/nesticus/all_dxy_results.tsv', 'dxy', 'n')
pd.dat = read_in('mpileup_analyses/ptomaphagus/all_dxy_results.tsv', 'dxy', 'p')
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


#GET STATS

