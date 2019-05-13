rm(list=ls())
plot.new()
dev.off()

library("ggplot2")
library("reshape2")
library("gridExtra")
library("data.table")

library(grid)
library(lattice)
library(ggpmisc)

### general set of useful functions###

# calculate slope of regression
slope <- function(x,y){
  return(cor(x,y)*sqrt(var(y))/sqrt(var(x)))
}

#calculate Rsquared for linear regression
rsquared <-function(x,y){
  return(cor(x, y) ^ 2)
}

# extract legend from ggplot for multiplotting
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


### read_data from HPC runs
path = "../Results/SIA_Results/"
extend <- function(v,max_length,len){
  len = length(v)
  if(len < max_length){finalfreq = v[len]
  addfreq = rep(finalfreq,max_length - len)
  v = as.vector(as.numeric(c(v,addfreq)))
  }
  return(v)
}
read_data_frequency <- function(seed,task){
  Invasion_path = paste(path,"Invasion_Frequency_",sep="")
  Inv_freq <- strsplit(readLines(paste(Invasion_path,seed,"_", task,".csv",sep="")),split=",") 
  max_length <- max(unlist(lapply(Inv_freq,length)))
  if(max_length< 4000){max_length = 4000}
  Inv_freq <- lapply(Inv_freq,extend,max_length = max_length)
  matrix <- t(matrix(unlist(Inv_freq), nrow=length(Inv_freq), byrow=T))
  suppressWarnings(storage.mode(matrix) <- "numeric")
  matrix <- as.data.frame(cbind(matrix,seq(1:nrow(matrix))))
  colnames(matrix) <- c(seq(1:(ncol(matrix)-1)),"Generation")
  if(max_length > 4000){matrix <- matrix[1:4000,]}
  dataframe <- melt(matrix,id.vars = c("Generation"),variable.name = "Run", value.name = "Invasion_Frequency")
  return(matrix)
}

read_data_frequency2 <- function(seed,task,maxlength){
  Invasion_path = paste(path,"Invasion_Frequency_",sep="")
  Inv_freq <- strsplit(readLines(paste(Invasion_path,seed,"_", task,".csv",sep="")),split=",") 
  max_length <- max(unlist(lapply(Inv_freq,length)))
  if(max_length< maxlength){max_length = maxlength}
  Inv_freq <- lapply(Inv_freq,extend,max_length = max_length)
  matrix <- t(matrix(unlist(Inv_freq), nrow=length(Inv_freq), byrow=T))
  suppressWarnings(storage.mode(matrix) <- "numeric")
  matrix <- as.data.frame(cbind(matrix,seq(1:nrow(matrix))))
  colnames(matrix) <- c(seq(1:(ncol(matrix)-1)),"Generation")
  if(max_length > maxlength){matrix <- matrix[1:maxlength,]}
  dataframe <- melt(matrix,id.vars = c("Generation"),variable.name = "Run", value.name = "Invasion_Frequency")
  return(matrix)
}


# stop()
plotdata= fread('Summary_Data.csv')
s = unique(plotdata$selection_strength)
u = unique(plotdata$mutation_rate)
p = unique(plotdata$propagule_pressure)
i <- unique(plotdata$resident_popsize)
r <- unique(plotdata$invader_popsize)
g <- unique(plotdata$invader_generations)

params = expand.grid(s,u,p,i,r,g)
param_inv_freq <- function(a,b,c,d,e,f){
  row = which(plotdata$selection_strength==s[a] &
                plotdata$mutation_rate == u[b] &
                plotdata$propagule_pressure == p[c] &
                plotdata$resident_popsize == i[d] &
                plotdata$invader_popsize == r[e]&
                plotdata$invader_generations== g[f])
  seed = plotdata$seed[row][1]
  task = plotdata$task_number[row][1]
  dataframe <- read_data_frequency(seed,task)
  return(dataframe)
}


Fig2 <- function(){
# 
#   df_freq1 <- param_inv_freq(2,3,10,3,2,5)
#   fwrite(df_freq1,'Timeseries/Fig2a.csv')
#   df_freq2 <- param_inv_freq(2,3,10,3,6,5)
#   fwrite(df_freq2,'Timeseries/Fig2b.csv')
#   df_freq3 <- param_inv_freq(2,3,10,3,10,5)
#   fwrite(df_freq3,'Timeseries/Fig2c.csv')
#   
  df_freq1 = fread('Timeseries/Fig2a.csv')
  df_freq2 = fread('Timeseries/Fig2b.csv')
  df_freq3 = fread('Timeseries/Fig2c.csv')
  
  df_freq1 <- melt(df_freq1,id.vars='Generation',variable.name = "Run",value.name = "Invasion_Frequency")
  p1 <- ggplot(df_freq1, aes(x = Generation, y = Invasion_Frequency,group=as.factor(Run))) + geom_line(size=1,col='#F8766D') + ylab("Invasion_Frequency") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none") +
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    xlim(0,4000) + ylim(0,1) + ylab("Invader Frequency")
  df_freq1_4 = df_freq1[df_freq1$Generation == 4,]
  df_freq1_40 = df_freq1[df_freq1$Generation == 40,]
  df_freq1_400 = df_freq1[df_freq1$Generation == 400,]
  df_freq1_4000 = df_freq1[df_freq1$Generation == 4000,]

  df_freq2 <- melt(df_freq2,id.vars='Generation',variable.name = "Run",value.name = "Invasion_Frequency")
  p2 <- ggplot(df_freq2, aes(x = Generation, y = Invasion_Frequency,group=as.factor(Run))) + geom_line(size=1,col='#00BA38') + ylab("Invasion_Frequency") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none") +
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    xlim(0,4000) + ylim(0,1) + ylab("Invader Frequency")
  df_freq2_4 = df_freq2[df_freq2$Generation == 4,]
  df_freq2_40 = df_freq2[df_freq2$Generation == 40,]
  df_freq2_400 = df_freq2[df_freq2$Generation == 400,]
  df_freq2_4000 = df_freq2[df_freq2$Generation == 4000,]
  
  df_freq3 <- melt(df_freq3,id.vars='Generation',variable.name = "Run",value.name = "Invasion_Frequency")
  p3 <- ggplot(df_freq3, aes(x = Generation, y = Invasion_Frequency,group=as.factor(Run))) + geom_line(size=1,col='#619CFF') + ylab("Invasion_Frequency") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none") +
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    xlim(0,4000) + ylim(0,1) + ylab("Invader Frequency")
  df_freq3_4 = df_freq3[df_freq3$Generation == 4,]
  df_freq3_40 = df_freq3[df_freq3$Generation == 40,]
  df_freq3_400 = df_freq3[df_freq3$Generation == 400,]
  df_freq3_4000 = df_freq3[df_freq3$Generation == 4000,]
  
  df_4 = rbind(cbind(df_freq1_4,Outcome = 'Extinction'),
                cbind(df_freq2_4,Outcome = 'Both'),
                cbind(df_freq3_4,Outcome = 'Fixation'))
  df_40 = rbind(cbind(df_freq1_40,Outcome = 'Extinction'),
                 cbind(df_freq2_40,Outcome = 'Both'),
                 cbind(df_freq3_40,Outcome = 'Fixation'))
  df_400 = rbind(cbind(df_freq1_400,Outcome = 'Extinction'),
                 cbind(df_freq2_400,Outcome = 'Both'),
                 cbind(df_freq3_400,Outcome = 'Fixation'))
  df_4000 = rbind(cbind(df_freq1_4000,Outcome = 'Extinction'),
                 cbind(df_freq2_4000,Outcome = 'Both'),
                 cbind(df_freq3_4000,Outcome = 'Fixation'))
  p4 <- ggplot(df_4) + geom_density(aes(x = Invasion_Frequency,fill=as.factor(Outcome)),alpha=0.5,bw=0.005) + ylab("Count") + theme_classic() + 
    scale_x_continuous(expand=c(0,0),limits=c(0,1)) + scale_y_continuous(expand=c(0,0)) + ggtitle('t = 4') + xlab("Invader Frequency") +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none") +
    theme(panel.border= element_blank()) +     theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))
  p5 <- ggplot(df_40) + geom_density(aes(x = Invasion_Frequency,fill=as.factor(Outcome)),alpha=0.5,bw=0.05) + ylab("Count") + theme_classic() + 
    scale_x_continuous(expand=c(0,0),limits=c(0,1)) + scale_y_continuous(expand=c(0,0)) + ggtitle('t = 40')  + xlab("Invader Frequency") +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none")  +
    theme(panel.border= element_blank())+ theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))
  p6 <- ggplot(df_400) + geom_density(aes(x = Invasion_Frequency,fill=as.factor(Outcome)),alpha=0.5,bw=0.05) + ylab("Count") + theme_classic() + 
    scale_x_continuous(expand=c(0,0),limits=c(0,1)) + scale_y_continuous(expand=c(0,0)) + ggtitle('t = 400')  + xlab("Invader Frequency") +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none")  +
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))
  p7 <- ggplot(df_4000) + geom_density(aes(x = Invasion_Frequency,fill=as.factor(Outcome)),alpha=0.5,bw=0.05) + ylab("Count") + theme_classic() + 
    scale_x_continuous(expand=c(0,0),limits=c(0,1)) + scale_y_continuous(expand=c(0,0)) + ggtitle('t = 4000')  + xlab("Invader Frequency") +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none")  +
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))
  p8 <- ggplot(df_4) + geom_density(aes(x = Invasion_Frequency,fill=as.factor(Outcome)),alpha=0.5,bw=0.002) + ylab("") + theme_classic() + 
    scale_x_continuous(expand=c(0,0),limits=c(0,.1),breaks=c(0,0.1)) + scale_y_continuous(expand=c(0,0),breaks=c(0,60)) + xlab('') +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none")  +
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))
  p1 = arrangeGrob(p1,top=textGrob('A)',x=0.05))
  p2 = arrangeGrob(p2,top=textGrob('B)',x=0.05))
  p3 = arrangeGrob(p3,top=textGrob('C)',x=0.05))
  d4 = grid.arrange(p4,p8,layout_matrix=rbind(c(1,1,1,1,1,1,1),
                                              c(1,1,2,2,2,2,1),
                                              c(1,1,2,2,2,2,1),
                                              c(1,1,2,2,2,2,1),
                                              c(1,1,2,2,2,2,1),
                                              c(1,1,1,1,1,1,1),
                                              c(1,1,1,1,1,1,1)))
  p4 = arrangeGrob(d4,top=textGrob('D)',x=0.05))
  p5 = arrangeGrob(p5,top=textGrob('E)',x=0.05))
  p6 = arrangeGrob(p6,top=textGrob('F)',x=0.05))
  p7 = arrangeGrob(p7,top=textGrob('G)',x=0.05))
  ggsave(file = '../Plots/Fig2.png', grid.arrange(arrangeGrob(p1),arrangeGrob(p2),arrangeGrob(p3),arrangeGrob(p4),arrangeGrob(p5),arrangeGrob(p6),arrangeGrob(p7),
                                                             layout_matrix=rbind(c(1,1,1,1,2,2,2,2,3,3,3,3),c(4,4,4,5,5,5,6,6,6,7,7,7))),width = 10, height = 6)
  # ggsave(filename = "../Plots/InvasionModel/OneRun.pdf",p1)
  return()
}

Fig3 <- function(){
  row1 = which(plotdata$selection_strength==s[3] & 
                 plotdata$mutation_rate == u[3] &
                 plotdata$propagule_pressure == p[10] & 
                 plotdata$resident_popsize ==i[3])
  plotdf1 <- plotdata[row1]
  # plotdf1$Success = plotdf1$Frequency1000
  # we are not actually meaning over anything are we?
  sum_df  =plotdf1[,list(mean=mean(Adaptation_Rate_Invader),sd=sd(Adaptation_Rate_Invader)),by=list(invader_popsize,invader_generations)]
  sum_df2  =plotdf1[,list(mean=mean(Frequency40),sd=sd(Frequency40)),by=list(invader_popsize,invader_generations,propagule_pressure)]

  leg <- g_legend(ggplot(sum_df, aes(x = invader_popsize, y =mean,color = as.factor(invader_generations))) +  geom_line(size=2) + geom_point(size=4) +
    xlab("Community Size") +ylab("Adaptation Rate") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())  + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(color='Generations') +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")))
  p1 <- ggplot(sum_df, aes(x = invader_popsize, y =mean,color = as.factor(invader_generations))) +  geom_line(size=2) + geom_point(size=4) +
    xlab("Invader Community Size") +ylab("Invader Adaptation Rate") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())  + guides(color=FALSE) + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(color='Generations') +  
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  # plotdf2 = plotdf1[plotdf1$invader_popsize == plotdf1$resident_popsize]
  p2 <- ggplot(plotdf1, aes(x = Adaptation_Rate_Invader, y = Relative_Propagule_Fitness,color = as.factor(invader_generations))) + 
    geom_point(size=0.5) + 
    xlab("Invader Adaptation Rate") +ylab("Relative Invader Fitness") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())  + guides(color=FALSE) + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(color='Generations') + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  p3 <- ggplot(plotdf1, aes(x = Relative_Propagule_Fitness, y =Frequency40,color = as.factor(invader_generations))) +  
    geom_point(size=0.5) + 
    xlab("Relative Invader Fitness") +ylab("Invasion Success") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + guides(color=FALSE) + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(color='Generations') + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  p4 <- ggplot(sum_df2, aes(x = invader_popsize, y = mean,color = as.factor(invader_generations))) +  geom_line(size=2) + geom_point(size=4) +
    xlab("Invader Community Size") +ylab("Invasion Success") + theme_classic() + geom_vline(xintercept = i[3],linetype=2) +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())  + theme(legend.key = element_blank()) + scale_x_log10() + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(color='Generations') + guides(fill=FALSE) + guides(color=guide_legend(override.aes=list(fill=NA))) + 
    # geom_errorbar(aes(ymin = mean - 1.96 * sd/sqrt(100), ymax = mean +1.96 * sd/sqrt(100),fill=as.factor(invader_generations)), alpha = 0.5,show.legend=FALSE) +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) + guides(color=FALSE)
  p1 = arrangeGrob(p1,top=textGrob('A)',x=0.05,gp=gpar(fontsize=15)))
  p2 = arrangeGrob(p2,top=textGrob('B)',x=0.05,gp=gpar(fontsize=15)))
  p3 = arrangeGrob(p3,top=textGrob('C)',x=0.05,gp=gpar(fontsize=15)))
  p4 = arrangeGrob(p4,top=textGrob('D)',x=0.05,gp=gpar(fontsize=15)))
  ggsave(file = '../Plots/Fig3.png',grid.arrange(p1,p2,p3,p4,leg
                                                           ,layout_matrix = rbind(c(1,1,1,2,2,2,NA),c(1,1,1,2,2,2,5),c(3,3,3,4,4,4,5),c(3,3,3,4,4,4,NA))),width = 13, height = 10)
  return( list(p1,p2,p3,p4))
}


Fig4 <- function(){
  row1 = which(plotdata$selection_strength==s[3] & 
                 plotdata$mutation_rate == u[3] &
                 plotdata$resident_popsize == i[3])
  plotdf1 <- plotdata[row1]
  options("scipen"=10) 
  tdf1 =  plotdf1[which(plotdf1$invader_generations == g[10] & plotdf1$invader_popsize %in% r[10]),]
  tdf2 =  plotdf1[which(plotdf1$invader_generations == g[2] & plotdf1$invader_popsize %in% r[2]),]
  tdf3 =  plotdf1[which(plotdf1$invader_generations == g[2] & plotdf1$invader_popsize %in% r[10]),]
  tbdf = rbind(tdf1,tdf2,tdf3)
  sum_df  =tbdf[,list(mean=mean(Frequency40),sd=sd(Frequency40)),by=list(invader_popsize,invader_generations,propagule_pressure,resident_popsize)]
  sum_df$treatment = paste('I = ', sum_df$invader_popsize, ' G = ',sum_df$invader_generations,sep = '') 
  p1 <-  ggplot(sum_df,aes(x=as.numeric(propagule_pressure), y=mean, shape=treatment,
                           group =treatment,colour=treatment)) + 
    theme_classic()+
    geom_point(size= 4) + geom_line(size = 2)  +  
    theme(legend.position=c(0.7,0.45)) + scale_x_log10() + 
    theme(panel.border= element_blank()) + 
    labs(x = 'Propagule Pressure (P)') +  
    theme(legend.position = c(0.7,0.6)) +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(y = 'Invasion Success',group='',shape='',col='') + geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd,fill=as.factor(treatment)), alpha = 0.5,show.legend=FALSE)+ 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm")) 

  #Mean is does not actuall change the valueas they are identical across factors
  sum_df2 = plotdf1[,list(pval = mean(PPpval),mean=mean(PP),sd=sd(PP),Relative_Propagule_Fitness = mean(Relative_Propagule_Fitness)),by=list(invader_popsize,invader_generations,resident_popsize)]
  sum_df2$correctedpval = sum_df2$pval*nrow(sum_df2)
  sum_df2[correctedpval<1e-16,]$correctedpval = 1e-16
  #NB if all values are significant you will need to change manual colours 
  p2 <- ggplot(sum_df2, aes(x = Relative_Propagule_Fitness,y = mean,col=correctedpval)) +  geom_point(size=2) + 
    xlab("Relative Invader Fitness") +ylab("Effect of Propagule Pressure") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + geom_vline(xintercept = 1,linetype=2,color='Red') +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    # guides(col=FALSE)+
    # scale_colour_manual(values=c('Red','Black')) +
    theme(text = element_text(size=18))   + labs(col='p-value') +
    scale_color_gradient(low='red',high='black',trans='log10',breaks=c(1e-2,1e-9,1e-16))
  
  sum_df3  =tbdf[,list(mean=mean(Frequency40),sd=sd(Frequency40)) ,by=list(invader_popsize,invader_generations,Propagule_Richness,propagule_pressure,resident_popsize)]
  sum_df3$treatment = paste('I = ', sum_df3$invader_popsize, ' G = ',sum_df3$invader_generations,sep = '') 
  # sum_df3$treatment = paste('R = ', sum_df3$invader_popsize, ' G = ',sum_df3$invader_generations,sep = '')
  # sum_df3 = subset(sum_df3,Propagule_Richness <=8)
  sum_df3 = sum_df3[which(sum_df3$propagule_pressure  == max(sum_df3$propagule_pressure)),]
  p3 <-  ggplot(sum_df3,aes(x=as.numeric(Propagule_Richness), y=mean, shape=treatment,
                           group =treatment,colour=treatment)) + 
    theme_classic()+
    geom_point(size= 4) + geom_line(size = 2)  +  
    theme(legend.position=c(0.7,0.45)) + scale_x_log10(breaks=c(3,10,30),limits=c(3,35)) +
    theme(panel.border= element_blank()) + 
    labs(x = 'Number of Invading Genotypes') +
    theme(legend.position = c(0.7,0.6)) +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(y = 'Invasion Success',group='',shape='',col='') + 
    geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd,fill=as.factor(treatment)), alpha = 0.5,show.legend=FALSE)+ 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm")) 
  
  sum_df4 = plotdf1[,list(pval = mean(PRpval),mean=mean(PR),sd=sd(PR),Relative_Propagule_Fitness = mean(Relative_Propagule_Fitness)),by=list(invader_popsize,invader_generations,resident_popsize,propagule_pressure)]
  #NA have arisen if there is no variation in invasion outcome or invader reichness. By defualt we assign this to have no effect so PR =0 and pval = 1
  sum_df4[is.na(sum_df4$mean)]$pval=1
  sum_df4[is.na(sum_df4$mean)]$mean=0
  #correct for multiple comparisons when assigning significance (bonferioni)
  sum_df4 = sum_df4[propagule_pressure==200,]
  sum_df4$correctedpval = sum_df4$pval*nrow(sum_df4)
  sum_df4[correctedpval>1]$correctedpval =1
  sum_df4[correctedpval<1e-6]$correctedpval = 1e-6
  # p4 <- ggplot(sum_df4, aes(x = Relative_Propagule_Fitness,y = mean,col=Signif)) +  geom_point(size=2) + 
  #   xlab("Relative Invader Fitness") +ylab("Effect of Invading Genotype Number") + theme_classic() + 
  #   theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
  #   theme(panel.border= element_blank())   + geom_vline(xintercept = 1,linetype=2,color='Red') +
  #   theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
  #   guides(col=FALSE)+ scale_colour_manual(values=c('Black','Red')) +
  #   theme(text = element_text(size=18)) 
  p4 <- ggplot(sum_df4, aes(x = Relative_Propagule_Fitness,y = mean,col=correctedpval)) +  geom_point(size=2) + 
    xlab("Relative Invader Fitness") +ylab("Effect of Invading Genotype Number") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + geom_vline(xintercept = 1,linetype=2,color='Red') +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    # guides(col=FALSE)+
    # scale_colour_manual(values=c('Red','Black')) +
    theme(text = element_text(size=18))   + labs(col='p-value') +
    scale_color_gradient(low='red',high='black',trans='log10',
                         labels = c('1','1e-02','1e-04','1e-06'),
                         breaks=c(1,1e-2,1e-4,1e-6))
  p1 = arrangeGrob(p1,top=textGrob('A)',x=0.05,gp=gpar(fontsize=15)))
  p2 = arrangeGrob(p2,top=textGrob('B)',x=0.05,gp=gpar(fontsize=15)))
  p3 = arrangeGrob(p3,top=textGrob('C)',x=0.05,gp=gpar(fontsize=15)))
  p4 = arrangeGrob(p4,top=textGrob('D)',x=0.05,gp=gpar(fontsize=15)))
  ggsave(file = '../Plots/Fig4.png',grid.arrange(arrangeGrob(p1),arrangeGrob(p2),arrangeGrob(p3),arrangeGrob(p4)
                                                          ,layout_matrix = rbind(c(1,2),c(3,4))),width = 12, height = 10)
  return( list(p1,p2,p3,p4))
}



Fig5<- function(){
  row1 = which(plotdata$selection_strength==s[1] &
                 plotdata$mutation_rate == u[3] &
                 plotdata$invader_popsize == r[6] &
                 plotdata$invader_generations == g[6] &
                 plotdata$propagule_pressure==p[10])
  row2 = which(plotdata$selection_strength==s[3] &
                 plotdata$mutation_rate == u[3] &
                 plotdata$invader_popsize == r[6] &
                 plotdata$invader_generations == g[6] &
                 plotdata$propagule_pressure==p[10])
  plotdf1 <- plotdata[row1] 
  plotdf2 <- plotdata[row2]
  plotdf3 = plotdf2[resident_popsize <= i[3],]
  p1 <-ggplot(data=plotdf1,aes(x=as.factor(resident_popsize),y=Frequency40)) + 
    geom_boxplot() + theme_classic() + 
    theme(text = element_text(size=18)) +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + #scale_y_continuous(limits=c(0,1)) +
    labs(x= 'Resident Community Size',y = 'Invasion Success') + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))
  p2 <-ggplot(data=plotdf2,aes(x=as.factor(resident_popsize),y=Frequency40)) + 
    geom_boxplot() + theme_classic() + 
    theme(text = element_text(size=18)) +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + #scale_y_continuous(limits=c(0,1)) +
    labs(x= 'Resident Community Size',y = 'Invasion Success') + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))


  reg_df = data.frame()
  for(i in unique(plotdf3$resident_popsize)){
    mod = lm(Frequency40~Resident_Richness_End,data=plotdf3[resident_popsize==i,])
    reg_df = rbind(reg_df,c(signif(mod$coefficients[2],2),signif(summary(mod)$coefficients[2,4],2),i))
  }
  mod = lm(Frequency40~Resident_Richness_End,data=plotdf3)
  print(summary(mod))
  colnames(reg_df) = c('Slope','P-Value','resident_popsize')
  p3 <-ggplot(data=plotdf3,
              aes(x=Resident_Richness_End,y=Frequency40,col=as.factor(resident_popsize))) + 
    geom_point() + theme_classic() + 
    theme(text = element_text(size=18)) +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())  + 
    geom_smooth(method='lm',se=FALSE) +
    geom_segment(x = 2, xend = 18,y = mod$coefficients[1]  + 2*mod$coefficients[2],
                 yend = mod$coefficients[1]  + 18*mod$coefficients[2],
                 col='Black')+
    annotate(x = 20, y =1, 
             label=paste('slope =', reg_df$Slope[1],', p = ' , reg_df$`P-Value`[1]),
             geom='text',size=5,colour='#F8766D') +
    annotate(x = 20, y =0.96, 
             label=paste('slope =', reg_df$Slope[2],', p = ' , reg_df$`P-Value`[2]),
             geom='text',size=5,colour='#00BA38') +
    annotate(x = 20, y =0.92, 
             label=paste('slope =', reg_df$Slope[3],', p = ' , reg_df$`P-Value`[3]),
             geom='text',size=5,colour='#619CFF') +
    annotate(x = 20, y =0.88, 
               label=paste('slope =', signif(mod$coefficients[2],2),', p <1e-6'),
               geom='text',size=5,colour='Black') +
    labs(x= 'Resident Community Genotypic Richness',y = 'Invasion Success',col='R') + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))


  options("scipen"=10)
  p1 = arrangeGrob(p1,top=textGrob('A)',x=0.05,gp=gpar(fontsize=15)))
  p2 = arrangeGrob(p2,top=textGrob('B)',x=0.05,gp=gpar(fontsize=15)))
  p3 = arrangeGrob(p3,top=textGrob('C)',x=0.05,gp=gpar(fontsize=15)))
  ggsave(file = '../Plots/Fig5.png',grid.arrange(arrangeGrob(p1),arrangeGrob(p2),arrangeGrob(p3),
                                                          layout_matrix = rbind(c(1,3,3),c(2,3,3))),
         width = 14, height = 8)

  return()
}

Fig6 <- function(){
  fig1 =read.csv('Jousset.csv')
  row1 = which(plotdata$selection_strength==s[1] & 
                 plotdata$mutation_rate == u[3] &
                 plotdata$invader_popsize == r[10] &
                 plotdata$invader_generations == g[1] &
                 plotdata$propagule_pressure == p[10])
  plotdf1 <- plotdata[row1]
  plotdf1=plotdf1[plotdf1$resident_popsize<=10000]
  # p1 <- ggplot(plotdf1,aes(x=Resident_Richness_End,y=Frequency100)) + geom_point()+ theme_classic() + scale_y_log10() + 
  #   geom_smooth(method='lm',formula=y~x) + annotate("text",x=10,y=0.3,label=as.character(expression(R^2==0.32)),parse=TRUE) + 
  #   annotate("text",x=13,y=0.3,label=as.character(expression(P<= 0.0001)),parse=TRUE)+ xlab('Genotypic Richness Resident Community') + ylab('Invasion Success (log_10_frequency)') +
  #   theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm"))
  p1a <- ggplot(plotdf1,aes(x=Resident_Richness_End,y=Frequency40*100)) + geom_point()+ theme_classic() +
    scale_y_log10(breaks=c(0.1, 0.316,1,3.16,10,31.6,100),labels = c('-1','-0.5','0','0.5','1','1.5','2'),limits=c(0.1,100)) + 
    scale_x_continuous(breaks=c(1,3,5,7,9,11))  + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    geom_smooth(method='lm',formula=y~x) + #annotate("text",x=11,y=0.1,label=as.character(expression(R^2==0.43)),parse=TRUE) +
    annotate("text",x=9,y=0.1,label=as.character(expression(R^2*'='*0.59*'    '*P<= 1e-6)),parse=TRUE)+
    labs(x='Resident Community Genotypic Richness',y=expression(atop('Invasion Sucess', '(log'[10]*'frequency)'))) +
    theme(text = element_text(size=15),axis.text = element_text(size=15))
  p1b <- ggplot(fig1,aes(x=Genotyic.Richness,y=Invader.Success)) + geom_point()+ theme_classic() +
    scale_y_continuous(breaks=c(-1,-0.5,-0,0.5,1,1.5,2),labels = c('-1','-0.5','0','0.5','1','1.5','2')) + 
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8),limits=c(1,8))   + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    geom_smooth(method='lm',formula=y~x) + #annotate("text",x=11,y=0.1,label=as.character(expression(R^2==0.43)),parse=TRUE) +
    annotate("text",x=7,y=-1,label=as.character(expression(R^2*'='*0.16*'    '*P<= 0.01)),parse=TRUE)+
    labs(x='Resident Community Genotypic Richness',y=expression(atop('Invasion Sucess', '(log'[10]*'frequency)'))) +
    theme(text = element_text(size=15),axis.text = element_text(size=15))
  fig2 = read.csv('Acosta.csv')
  row2 = c(plotdata$selection_strength==s[1] & 
             plotdata$mutation_rate == u[3] &
             plotdata$invader_popsize == r[5] &
             plotdata$invader_generations == g[5] &
             plotdata$resident_popsize == i[3]&
             plotdata$propagule_pressure %in% c(1,10,100))
  plotdf2 <- plotdata[row2]
  # freq_df = rbind(cbind(melt(read_data_frequency2(unique(plotdf2$seed)[1],unique(plotdf2$task_number),500),id.vars='Generation',variable.name='Run',value.name='Invasion_Frequency'),'p' = 1),
  #                 cbind(melt(read_data_frequency2(unique(plotdf2$seed)[2],unique(plotdf2$task_number),500),id.vars='Generation',variable.name='Run',value.name='Invasion_Frequency'),'p' = 10),
  #                 cbind(melt(read_data_frequency2(unique(plotdf2$seed)[3],unique(plotdf2$task_number),500),id.vars='Generation',variable.name='Run',value.name='Invasion_Frequency'),'p' = 100))
  # fwrite(freq_df,'Timeseries/Fig6D')
  freq_df = fread('Timeseries/Fig6d.csv')
  ag <- aggregate(. ~ p*Generation, freq_df, function(x) c(mean = mean(x), sd = sd(x)))
  ag$mean = ag$Invasion_Frequency[,1] *i[3]
  ag$sd = ag$Invasion_Frequency[,2]*i[3]
  ag$Treatment = NA
  ag$Treatment[ag$p==1] = '1'
  ag$Treatment[ag$p==10] = '10'
  ag$Treatment[ag$p==100] = '100'
  ag$Treatment = factor(ag$Treatment,levels=c('1','10','100'))
  fig2$Treatment = factor(fig2$Treatment,levels=c('Low','Medium','High'))
  
  fig2$Time=fig2$Time-min(fig2$Time)
  p2a <-  ggplot(ag[ag$Generation%in%c(1,100,200),],aes(x = Generation,y=mean,colour=Treatment)) + 
    theme_classic()+ geom_line(size=2) +
    xlim(0,200) + geom_point(size=4)  +  
    labs(x = 'Time(Generations)',y =expression(atop('Invasion Sucess', '(log'[10]*'Abundance')),colour='P') +
    scale_y_log10(breaks=c(0.1,1,10,100,1000),labels = c('0  ','1  ','2  ','3  ','4  '),limits=c(0.1,1000))+  
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    #+ #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd,fill=as.factor(p)), alpha = 0.5,show.legend=FALSE)+ 
    theme(legend.position = c(0.8,0.58),legend.text = element_text(size=12),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm")) +
    theme(text = element_text(size=15),axis.text = element_text(size=15))
  p2b <-  ggplot(fig2,aes(x = Time,y=Abundance.Log.,colour=as.factor(Treatment))) + 
    theme_classic()+ geom_line(size=2) + geom_point(size=4)+
    labs(x = 'Time(Days)',y = expression(atop('Invasion Sucess', '(log'[10]*'Abundance)')),colour='P') + 
    scale_y_continuous(breaks=c(0,1,2,3,4),labels = c('0  ','1  ','2  ','3  ','4  '))+  
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
  #+ #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd,fill=as.factor(p)), alpha = 0.5,show.legend=FALSE)+ 
    theme(legend.position = c(0.8,0.58),legend.text = element_text(size=12),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm")) +
    theme(text = element_text(size=15),axis.text = element_text(size=15))
  
  
  ggsave(file = '../Plots/Fig6.png',grid.arrange(arrangeGrob(p1b,top=textGrob('A)',x=0.05,gp=gpar(fontsize=15))),
                                                            arrangeGrob(p2b,top=textGrob('C)',x=0.05,gp=gpar(fontsize=15))),
                                                            arrangeGrob(p1a,top=textGrob('B)',x=0.05,gp=gpar(fontsize=15))),
                                                            arrangeGrob(p2a,top=textGrob('D)',x=0.05,gp=gpar(fontsize=15))),ncol=2),
                                                            width = 12, height = 10)
}           
# 

#NBF fig2 requires the raw data to run which is humungous
Fig2()
Fig3()
Fig4()
Fig5()
#Similary fig 6 panel 1b requires a large data file to be accessed.
Fig6()
# Resident_Plots(2,2,2)
