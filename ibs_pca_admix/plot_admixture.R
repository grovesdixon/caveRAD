#plot_admixture.R
library(tidyverse)
library(rlang)
library(cowplot)
rm(list=ls())
source('ibs_pca_admix/admixture_functions.R')



# GLOBAL VARS -------------------------------------------------------------

KPOPS=5
PREFIX='mydata_'
XLAB = 'cave'
YLAB = 'ancestry'


# BUILD FOR NESTICUS ------------------------------------------------------

#read in the data
dfList = read_admixture(samplesFile = 'ibs_pca_admix/nesticus/samples.txt',
                          qoptDir = 'ibs_pca_admix/nesticus/',
                          kpops = KPOPS,
                          prefix = PREFIX)

#apply the plotting function to the list of dataframes
n.splots=lapply(names(dfList), function(x) plot_admix_samples(x))
n.splots[[2]] = n.splots[[2]] + theme(axis.text.y = element_blank())
n.splots[[4]] = n.splots[[4]] + theme(axis.text.y = element_blank())
names(n.splots)=c('k2', 'k3', 'k4', 'k5')

# #plot the list of plots
# plts=plot_grid(plotlist=n.splots, nrow=2)  
# ylab = ggdraw() + draw_label(YLAB, angle=90)
# xlab = ggdraw() + draw_label(XLAB)
# top0=plot_grid(ylab, plts, nrow=1, rel_widths=c(1, 25))
# title=ggdraw() + draw_label('Nesticus', fontface='italic')
# top=plot_grid(title, top0, nrow=2, rel_heights=c(1,15))
# n.plt = plot_grid(top, xlab, nrow=2, rel_heights=c(16,1), rel_widths = c(1,0.9))
# # n.plt


#plot the list of plots
plts=plot_grid(plotlist=n.splots, nrow=2, rel_widths=c(1, 0.75))
ylab = ggdraw() + draw_label(YLAB, angle=90)
top0=plot_grid(ylab, plts, nrow=1, rel_widths=c(1, 25))
title=ggdraw() + draw_label('N. barri', fontface='italic')
n.plt=plot_grid(title, top0, nrow=2, rel_heights=c(1,15))
# n.plt




# BUILD FOR PTOMAPHAGUS ---------------------------------------------------

#read in the data
dfList = read_admixture(samplesFile = 'ibs_pca_admix/ptomaphagus/samples.txt',
                        qoptDir = 'ibs_pca_admix/ptomaphagus/',
                        kpops = KPOPS,
                        prefix = PREFIX)

#apply the plotting function to the list of dataframes
p.splots=lapply(names(dfList), function(x) plot_admix_samples(x))
p.splots = lapply(p.splots, function(x) return(x+theme(axis.text.y = element_blank())))
names(n.splots)=c('k2', 'k3', 'k4', 'k5')

# #plot the list of plots
# plts=plot_grid(plotlist=p.splots, nrow=2)  
# ylab= ggdraw() + draw_label(YLAB, angle=90)
# xlab = ggdraw() + draw_label(XLAB)
# top0=plot_grid(ylab, plts, nrow=1, rel_widths=c(1,16))
# title=ggdraw() + draw_label('Ptomaphagus', fontface='italic')
# top=plot_grid(title, top0, nrow=2, rel_heights=c(1,15))
# p.plt = plot_grid(top, xlab, nrow=2, rel_heights=c(16,1))
# # p.plt

#plot the list of plots
plts=plot_grid(plotlist=p.splots, nrow=2)  
title=ggdraw() + draw_label('P. hatchi', fontface='italic')
p.plt=plot_grid(title, plts, nrow=2, rel_heights=c(1,15))
# p.plt


f.top = plot_grid(n.plt, p.plt, nrow=1, rel_widths = c(1, 0.9))
xlab = ggdraw() + draw_label(XLAB)
final = plot_grid(f.top, xlab, nrow=2, rel_heights = c(15,1))
final


