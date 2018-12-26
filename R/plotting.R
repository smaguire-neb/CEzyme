# Make peak plots, probably send to PDF 6 or so / page
colors<-function(channel){
  case_when(channel == "ch1" ~ "royalblue3",
            channel == "ch2" ~ "springgreen3",
            channel == "ch3" ~ "goldenrod2",
            channel == "ch4" ~ "firebrick2",
            channel == "ch5" ~ "orange")
}

basic.plot<-function(data,x,y,ymin,ymax,fill){
  ggplot(data,aes_string(x=x,y=y)) + geom_line() +
    geom_ribbon(aes_string(ymin=ymin,ymax=y),fill=fill)
}


add.segments<-function(basicplot,peaks,seg.col,ar = arrow(angle = 90,length = unit(0.1, "inches"),
                                                          ends = "both", type = "open")){
  if(any(peaks$height > 0)){
    ymax<-max(peaks$height) + 2000
    basicplot+
      geom_segment(data=peaks,aes_string(x="start.bps",xend="end.bps"),
                   y=0,yend=0,col = seg.col,arrow = ar) +
      ylab("RFU") + xlab("Size (bp)") +
      geom_label(data=peaks,aes(x= size.bps,y= height+1000,
                                label=paste0("size: ",round(size.bps,digits = 2),"\n",
                                             "area: ",round(area,digits = 2))),
                 size=2) +
      ylim(c(-50,ymax))
  } else{
    basicplot +
      ylab("RFU") + xlab("Size (bp)") +
      geom_label(x=50,y=0,label="No Peaks Found",
                 size=2)

  }
}
explore.plots<-function(data,peaks,channel=c("ch1","ch2","ch3","ch4"),xlims=c(0,120),
                        color = NULL, file.name = "rplots.pdf"){
  #require(friendlyeval)
  require(ggplot2)
  require(cowplot)
  if(is.null(color)){
    color <- colors(channel)
  }

  data <- map(data, function(x) filter(x,size.bps >= 0 & size.bps <= xlims[2])) %>%
    map(function(x) group_by(x,position) %>%
          mutate(ymin = `if`(!!rlang::ensym(channel) > 0, return(0), return(!!rlang::ensym(channel)))) %>%
          group_by())
  basicplots<-flatten(list(map(data,basic.plot,x="size.bps",y=channel,ymin="ymin",ymax=channel,fill=color)))
  fullplots<-flatten(list(map2(basicplots,peaks,function(z,w) add.segments(basicplot = z,peaks = w,seg.col = "red"))))
  plot.names<-as.list(names(fullplots))
  fullplots<-map2(fullplots,plot.names, function(x,y) x+ggtitle(y))

  pdf(file = file.name)
  i<-1
  while(i <= length(fullplots)){
    i.max<-ifelse(i+3<=length(fullplots),i+3,length(fullplots))
    print(plot_grid(plotlist = fullplots[i:i.max],ncol=2,nrow=2))
    i<-i+3
  }
  dev.off()
}

plot.specific<-function(data,peaks,channel=c("ch1","ch2","ch3","ch4"),xlims=c(0,120),
                        color = NULL,labels=F,selected){
  #require(friendlyeval)
  require(ggplot2)
  if(is.null(color)){
    color <- colors(channel)
  }

  data <- map(data, function(x) filter(x,size.bps >= 0 & size.bps <= xlims[2])) %>%
    map(function(x) group_by(x,position) %>%
          mutate(ymin = `if`(!!rlang::ensym(channel) > 0, return(0), return(!!rlang::ensym(channel)))) %>%
          group_by())
  basicplots<-flatten(list(map(data,basic.plot,x="size.bps",y=channel,ymin="ymin",ymax=channel,fill=color)))
  if(labels ==T){
    basicplots<-flatten(list(map2(basicplots,peaks,function(z,w) add.segments(basicplot = z,peaks = w,seg.col = "red"))))
  }
  plot.names<-as.list(names(basicplots))
  fullplots<-map2(basicplots,plot.names, function(x,y) x+ggtitle(y))
  fullplots[[selected]]+
    xlab("Size (bp)") + ylab("RFU") +
    theme(panel.background = element_rect(fill="gray95")) %>%
    return()
}

