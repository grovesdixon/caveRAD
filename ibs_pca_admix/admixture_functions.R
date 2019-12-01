#funciton to read in admixture results from a directory
read_admixture = function(samplesFile, qoptDir, kpops, prefix, main){
  samples = read.table(samplesFile)$V1
  caves=substr(samples, start=3, stop=4)
  dfList = list()
  for (i in 2:kpops){
    infile=paste(c(prefix,'k',i,'.qopt'), collapse='')
    inPath = paste(qoptDir, infile, sep='')
    indat = read.table(inPath)
    colnames(indat) = paste('pop', 1:i, sep='')
    indat$sample = samples
    indat$cave = substr(samples, start=3, stop=4)
    head(indat)
    dfName = paste('k=',i,sep='')
    dfList[[dfName]] = as_tibble(indat) 
  }
  return(dfList)
}

#funciton to get coordinates for divider lines in barplots
get_segment_coords = function(sites){
  ns = table(sites)
  xs = c(ns[1],
         ns[1]+ns[2],
         ns[1]+ns[2]+ns[3],
         NA) + 0.5
  lab.xs = c(ns[1]/2,
             ns[1]+ns[2]/2,
             ns[1]+ns[2]+ns[3]/2,
             ns[1]+ns[2]+ns[3]+ns[4]/2) + 0.5
  seg.df = data.frame('x1'=xs,
                      'x2'=xs,
                      'y1'=0,
                      'y2' = 1,
                      'lab.xs'=lab.xs,
                      'site'=unique(sites))
  return(seg.df)
}

#funciton to sort
sort_by_pop_props = function(dat){
  sdat = data.frame()
  sites = unique(dat$cave)
  for (s in sites){
    sub=dat %>% 
      filter(cave==s)
    popSums = sub %>% 
      gather(key='pop', value='proportion', grep('pop', colnames(sub))) %>%
      group_by(pop) %>% 
      summarize(popSum=sum(proportion)) %>% 
      arrange(popSum)
    # maxPop = popSums %>% 
    #   filter(popSum==max(popSum)) %>% 
    #   pull(pop)
    sigPops = popSums %>% 
      filter(popSum > 0.1) %>% 
      arrange(pop) %>% 
      pull(pop) %>% 
      rev()
    s.sub = sub %>% 
      arrange(!!sym(sigPops[1]))
    sdat=rbind(sdat,s.sub)
  }
  return(sdat)
}


#function to plot with sample names
plot_admix_samples = function(dfName){
  dat=dfList[[dfName]]
  seg.df = get_segment_coords(dat$cave)
  dat %>%
    sort_by_pop_props() %>% 
    mutate(num=1:nrow(dat)) %>% 
    gather(key='pop', value='proportion', grep('pop', colnames(dat))) %>%
    ggplot() +
    geom_bar(aes(x=num, y=proportion, fill=pop), stat="identity", width=1.0) +
    labs(x='Sample', y='Admixture', subtitle=dfName) +
    scale_x_continuous(breaks = seg.df$lab.xs, labels=seg.df$site) +
    geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2), data=seg.df, lwd=0.5) +
    theme(legend.position='none',
          axis.title = element_blank())
}
