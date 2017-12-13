library(ggplot2)
library(scales)
library(grid)


#filename="/home/aylingm/testgenome/patronol/contigs_46_100k.1_clean.fa.degrees"
#filename="/home/aylingm/testgenome/mc_consensus_3k/contigs.fa.degrees"

cmd_args = commandArgs(TRUE);
filename=cmd_args[1];
#print(filename)

format="png"
# these need to be fixed! use dirname(), and check for directory existence, mkdir if not there
graphsdir=paste(dirname(filename), "/graphs", sep="");
if(!dir.exists(graphsdir)){
  dir.create(graphsdir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

#print(graphsdir)
  
types = c("I", "Y", "X","E");
min_graph=100 # need to have this as a variable that can be passed to the script
scaling_factor=100
bin_size=0.0001
title_font=7
decimal_points=4

degrees <- read.delim(filename)
if (min(degrees$total)<min_graph){
  min_graph<-min(degrees$total)-1
}
degrees <- degrees[degrees$total>min_graph,]


for (t in 1:4) {
  x_high<-c(0.95,1)
  x_low<-c(0,0.05)
  type = types[t];
  if (type=="I"){
    Is<-data.frame(degrees$for.1.rev.1.)
    degree_dist<-data.frame(Is)
    x_min<-x_high[1]
    x_max<-x_high[2]
  }
  else if (type=="Y"){
    Ys<-data.frame(degrees$for.1.rev.2., degrees$for.1.rev.3., degrees$for.1.rev.4., degrees$for.2.rev.1., degrees$for.3.rev.1., degrees$for.4.rev.1.)
    degree_dist<-data.frame(rowMeans(Ys))
    x_min<-x_low[1]
    x_max<-x_low[2]
  }
  else if (type=="X"){
    Xs<-data.frame(degrees$for.1.rev.2., degrees$for.1.rev.3., degrees$for.1.rev.4., degrees$for.2.rev.1., degrees$for.3.rev.1., degrees$for.4.rev.1.,degrees$for.2.rev.2., degrees$for.2.rev.3., degrees$for.2.rev.4., degrees$for.3.rev.3., degrees$for.3.rev.2., degrees$for.3.rev.4., degrees$for.4.rev.4., degrees$for.4.rev.2., degrees$for.4.rev.3.)
    degree_dist<-data.frame(rowMeans(Xs))
    x_min<-x_low[1]
    x_max<-x_low[2]
  }
  else if (type=="E"){
    Es<-data.frame(degrees$for.0.rev.0., degrees$for.0.rev.1., degrees$for.0.rev.2., degrees$for.0.rev.3., degrees$for.0.rev.4., degrees$for.1.rev.0., degrees$for.2.rev.0., degrees$for.3.rev.0., degrees$for.4.rev.0.)
    degree_dist<-data.frame(rowMeans(Es))
    x_min<-x_low[1]
    x_max<-x_low[2]
  }
    
  degrees_png <- paste(graphsdir, "/", type,"_degrees.png", sep="");
  png(degrees_png, width=1200, height=800)
  degree_dist<-degree_dist*scaling_factor
  sd_dd<-sd(degree_dist[,1])
  mean_dd<-mean(degree_dist[,1])
  print(ggplot(degree_dist, aes(degree_dist)) + ggtitle(paste(type ,"nodes probability,", "mean=",round(mean_dd,decimal_points),"sd=",round(sd_dd,decimal_points))) + 
          xlab("probability (%)") + ylab("Subgraph Count") + geom_histogram(binwidth=bin_size*scaling_factor) + 
          scale_x_continuous(breaks=seq(0, 1*scaling_factor, bin_size*10*scaling_factor), limits=c(x_min*scaling_factor, x_max*scaling_factor))) +
          theme(plot.title = element_text(size=title_font))
  garbage <- dev.off()
}

type="E"
degrees.not_meaned<-data.frame(degrees[degrees$total>min_graph,-26])
degrees.m<-data.frame("FOR"=rep(c(0:4),each = 5),"REV"=rep(0:4,5), data.frame(colMeans(degrees[degrees$total>min_graph,-26]*scaling_factor))) # manually 'melt' the average for each point
names(degrees.m)[3]<-"prob"
degrees_png <- paste(graphsdir, "/", type,"_heatmap.png", sep="")
png(degrees_png, width=1200, height=800)
print(ggplot(degrees.m, aes(x=FOR, y=REV)) + geom_tile(aes(fill = prob), colour = "white") + ggtitle(paste(filename, "_", type,"_heatmap")) + theme(plot.title = element_text(size=title_font)) + scale_fill_gradient(low = "white", high = "steelblue", trans='log'))
garbage <- dev.off()

degrees.not_meaned<-data.frame(degrees[degrees$total>min_graph,-26])
degrees.not_meaned$for.1.rev.1.<-0
degrees.m<-data.frame("FOR"=rep(c(0:4),each = 5),"REV"=rep(0:4,5), data.frame(colMeans(degrees.not_meaned*scaling_factor))) # manually 'melt' the average for each point
names(degrees.m)[3]<-"prob"
degrees_png <- paste(graphsdir, "/", type,"_heatmap_noI.png", sep="")
png(degrees_png, width=1200, height=800)
print(ggplot(degrees.m, aes(x=FOR, y=REV)) + geom_tile(aes(fill = prob), colour = "white")+ ggtitle(paste(filename, "_", type,"_heatmap (ex. I nodes)")) + theme(plot.title = element_text(size=title_font))+ scale_fill_gradient(low = "white", high = "steelblue", trans='log'))
garbage <- dev.off()
