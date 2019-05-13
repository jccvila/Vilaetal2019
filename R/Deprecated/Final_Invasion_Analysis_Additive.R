rm(list=ls())
plot.new()
dev.off()

library("ggplot2")
library("reshape2")
library("gridExtra")
library("data.table")

library(grid)
library(lattice)


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


### read_data from HPC reuns


path = "../Results/SIA_Results/"

extend <- function(v,max_length,len){
  len = length(v)
  if(len < max_length){finalfreq = v[len]
  addfreq = rep(finalfreq,max_length - len)
  v = as.vector(as.numeric(c(v,addfreq)))
  }
  return(v)
}

sum_length <- function(v,len_limit){
  len = length(v)
  if(len > len_limit){v = v[1:len_limit]}
  if(len < len_limit){v = c(v,(v[rep(len,len_limit-len)]))}
  return(sum(as.numeric(v)))
}

end_length <- function(v,len_limit){
  len = length(v)
  if(len >= len_limit){v = v[len_limit]}
  if(len < len_limit){v = v[len]}
  return(v)
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

Invador_Frequency4000 <- function(seed,task){
  Invasion_path = paste(path,"Invasion_Frequency_",sep="")
  Inv_freq <- strsplit(readLines(paste(Invasion_path,seed,"_", task,".csv",sep="")),split=",") 
  return(as.numeric(unlist(lapply(Inv_freq,end_length,len_limit = 4000))))
}
Invador_Frequency400 <- function(seed,task){
  Invasion_path = paste(path,"Invasion_Frequency_",sep="")
  Inv_freq <- strsplit(readLines(paste(Invasion_path,seed,"_", task,".csv",sep="")),split=",") 
  return(as.numeric(unlist(lapply(Inv_freq,end_length,len_limit = 400))))
}

Invador_Frequency40 <- function(seed,task){
  Invasion_path = paste(path,"Invasion_Frequency_",sep="")
  Inv_freq <- strsplit(readLines(paste(Invasion_path,seed,"_", task,".csv",sep="")),split=",") 
  return(as.numeric(unlist(lapply(Inv_freq,end_length,len_limit = 40))))
}

Invador_Frequency4 <- function(seed,task){
  Invasion_path = paste(path,"Invasion_Frequency_",sep="")
  Inv_freq <- strsplit(readLines(paste(Invasion_path,seed,"_", task,".csv",sep="")),split=",") 
  return(as.numeric(unlist(lapply(Inv_freq,end_length,len_limit = 4))))
}

read_data_summary<- function(seed,task){
  Invasion_path = paste(path,"Summary_Data_",seed,'_',task,'.csv',sep="")
  mydf = read.csv(Invasion_path,header=FALSE)
  tmydf = setNames(data.frame(t(mydf[,-1])), mydf[,1])
  # tmydf$Invasion_Success = Invasion_Success(seed,task)
  # paramaters = read_data_paramaters(seed,task)
  # tmydf$invader_popsize = paramaters$Value[8]
  # tmydf$invader_generations = paramaters$Value[9]
  return(tmydf)
}

read_data_paramaters <- function(seed,task){
  params <- read.csv(paste(path,"Paramaters_",seed,"_", task,".csv",sep=""),check.names = FALSE,header=FALSE,col.names = c("Paramater","Value"))
  return(params)
}

read_full <- function(j,i){
  Frequency4 <- Invador_Frequency4(j,i)
  Frequency40 <- Invador_Frequency40(j,i)
  Frequency400 <- Invador_Frequency400(j,i)
  Frequency4000 <- Invador_Frequency4000(j,i)
  # Frequency1 <- Invador_Frequency1(j,i)
  param_data = read_data_paramaters(j,i)
  Paramaters <- t(replicate(length(Frequency4), param_data[,2]))
  colnames(Paramaters) <- param_data[,1]
  Other <- read_data_summary(j,i)
  return(cbind(Paramaters,Other,Frequency4,Frequency40,Frequency400,Frequency4000))
}

# read_seed <- function(j){return(rbindlist(lapply(1:100,read_full,j=j)))}
# a <- rbindlist(lapply(1:90,read_seed))
# print('a')
# b <- rbindlist(lapply(91:180,read_seed))
# print('b')
# c <- rbindlist(lapply(181:270,read_seed))
# print('c')
# d <- rbindlist(lapply(271:360,read_seed))
# print('d')
# e <- rbindlist(lapply(361:450,read_seed))
# print('e')
# plotdata <- rbind(a,b,c,d,e)
# save(plotdata,file = 'Additive_Data.rds')
# load('Additive_Data.rds')
# # #
# plotdata$Relative_Fitness = plotdata$Invader_Mean_End/plotdata$Invaded_Mean_End
# plotdata$Relative_Propagule_Fitness = plotdata$Propagule_Fitness_Mean/plotdata$Invaded_Mean_End
# plotdata$Fitness_Increase_Invaded = plotdata$Invader_Mean_End - plotdata$Invader_Mean_Start
# plotdata$Fitness_Increase_Invader= plotdata$Invaded_Mean_End - plotdata$Invaded_Mean_Start
# plotdata$Adaptation_Rate_Invader = plotdata$Fitness_Increase_Invaded/plotdata$invader_generations
# plotdata$Adaptation_Rate_Invaded = plotdata$Fitness_Increase_Invader/plotdata$invader_generations
# plotdata$Adaptation_Rate_Difference= plotdata$Adaptation_Rate_Invader - plotdata$Adaptation_Rate_Invaded
# plotdata$Relative_Propagule_Fitness = plotdata$Propagule_Fitness_Mean/plotdata$Invaded_Mean_End
# #
# # # generic ###
# # ### single param combinations ###
# #
# s = unique(plotdata$selection_strength)
# u = unique(plotdata$mutation_rate)
# p = unique(plotdata$propagule_pressure)
# i <- unique(plotdata$invaded_popsize)
# r <- unique(plotdata$invader_popsize)
# g <- unique(plotdata$invader_generations)
# 
# #Propagule effect
# #Propagule effect rsquared
# plotdata$PP = 0
# plotdata$PPRS= 0
# #Propagule Richness
# #Propagule richness rsquared
# plotdata$PR = 0
# plotdata$PRRS= 0
# #Invaded Size
# #Invaded Size Rsquared
# plotdata$IS= 0
# plotdata$ISRS= 0
# #Invaded Richness
# #Invaded Richness Rsquared
# plotdata$IR = 0
# plotdata$IRRS= 0
# #
# params = data.table(expand.grid(s,u,p,i,r,g))
# #
# PP = c()
# PPRS = c()
# PR = c()
# PRRS =c()
# IS = c()
# ISRS = c()
# IR = c()
# IRRS = c()
# plotdata[,PP:=slope(propagule_pressure,Frequency40),by=list(selection_strength,mutation_rate,invaded_popsize,invader_popsize,invader_generations)]
# plotdata[,PPRS:=rsquared(propagule_pressure,Frequency40),by=list(selection_strength,mutation_rate,invaded_popsize,invader_popsize,invader_generations)]
# plotdata[,PR:=slope(Propagule_Richness,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,invaded_popsize,invader_popsize,invader_generations)]
# plotdata[,PRRS:=rsquared(Propagule_Richness,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,invaded_popsize,invader_popsize,invader_generations)]
# plotdata[,IS:=slope(invaded_popsize,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,invader_popsize,invader_generations)]
# plotdata[,ISRS:=rsquared(invaded_popsize,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,invader_popsize,invader_generations)]
# plotdata[,IR:=slope(Invaded_Richness_Start,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,invaded_popsize,invader_popsize,invader_generations)]
# plotdata[,IRRS:=rsquared(Invaded_Richness_Start,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,invaded_popsize,invader_popsize,invader_generations)]

# for(j in 1:nrow(params)){
#   print(j)
#   tr = as.matrix(params[j])
#   t1 = plotdata[selection_strength==tr[1]&mutation_rate==tr[2]&invaded_popsize==tr[4]&invader_popsize==tr[5] &invader_generations==tr[6]]
#   t2 = plotdata[selection_strength==tr[1]&mutation_rate==tr[2]&propagule_pressure==tr[3]&invader_popsize==tr[5] &invader_generations==tr[6]]
#   t3 = plotdata[selection_strength==tr[1]&mutation_rate==tr[2]&propagule_pressure==tr[3]&invaded_popsize == tr[4]&invader_popsize==tr[5]&invader_generations==tr[6]]
#   PP = c(PP,slope(t1$propagule_pressure,t1$Frequency40))
#   PPRS =c(PPRS,rsquared(t1$propagule_pressure,t1$Frequency40))
#   PR=c(PR,slope(t3$Propagule_Richness,t3$Frequency40))
#   PRRS=c(rsquared(t3$Propagule_Richness,t3$Frequency40))
#   IS = c(IS,slope(t2$invaded_popsize,t2$Frequency40))
#   ISRS = c(ISRS,rsquared(t2$invaded_popsize,t2$Frequency40))
#   IR = c(IR,slope(t3$Invaded_Richness_Start,t3$Frequency40))
#   IRRS = c(IRRS,rsquared(t3$Invaded_Richness_Start,t3$Frequency40))
# }
# 
# # 
# # # function takes abundance vector and converts it into popsize.
# save(plotdata,file = 'Additive_Data3.rds')
# stop()
load('Additive_Data3.rds')
s = unique(plotdata$selection_strength)
u = unique(plotdata$mutation_rate)
p = unique(plotdata$propagule_pressure)
i <- unique(plotdata$invaded_popsize)
r <- unique(plotdata$invader_popsize)
g <- unique(plotdata$invader_generations)

params = expand.grid(s,u,p,i,r,g)
param_inv_freq <- function(a,b,c,d,e,f){
  row = which(plotdata$selection_strength==s[a] & 
                plotdata$mutation_rate == u[b] &
                plotdata$propagule_pressure == p[c] & 
                plotdata$invaded_popsize == i[d] & 
                plotdata$invader_popsize == r[e]&
                plotdata$invader_generations== g[f])
  seed = plotdata$seed[row][1]
  task = plotdata$task_number[row][1]
  dataframe <- read_data_frequency(seed,task)
  return(dataframe)
}


Time_series_plot <- function(){
  # row = which(plotdata$selection_strength==s[a] & 
  #               plotdata$mutation_rate == u[b] &
  #               plotdata$propagule_pressure == p[c] & 
  #               plotdata$invaded_popsize == i[d] & 
  #               plotdata$invader_popsize == r[e]&
  #               plotdata$invader_generations== g[f])
  df_freq1 <- param_inv_freq(3,3,10,3,2,5)
  df_freq1 <- melt(df_freq1,id.vars='Generation',variable.name = "Run",value.name = "Invasion_Frequency")
  p1 <- ggplot(df_freq1, aes(x = Generation, y = Invasion_Frequency,group=as.factor(Run))) + geom_line(size=1,col='#F8766D') + ylab("Invasion_Frequency") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none") +
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    xlim(0,4000) + ylim(0,1) + ylab("Invader Frequency")
  df_freq1_4 = df_freq1[df_freq1$Generation == 4,]
  df_freq1_40 = df_freq1[df_freq1$Generation == 40,]
  df_freq1_400 = df_freq1[df_freq1$Generation == 400,]
  df_freq1_4000 = df_freq1[df_freq1$Generation == 4000,]
  df_freq2 <- param_inv_freq(3,3,10,3,6,5)
  df_freq2 <- melt(df_freq2,id.vars='Generation',variable.name = "Run",value.name = "Invasion_Frequency")
  p2 <- ggplot(df_freq2, aes(x = Generation, y = Invasion_Frequency,group=as.factor(Run))) + geom_line(size=1,col='#00BA38') + ylab("Invasion_Frequency") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none") +
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    xlim(0,4000) + ylim(0,1) + ylab("Invader Frequency")
  df_freq2_4 = df_freq2[df_freq2$Generation == 4,]
  df_freq2_40 = df_freq2[df_freq2$Generation == 40,]
  df_freq2_400 = df_freq2[df_freq2$Generation == 400,]
  df_freq2_4000 = df_freq2[df_freq2$Generation == 4000,]
  df_freq3 <- param_inv_freq(3,3,10,3,10,5)
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
  ggsave(file = '../../Plots/Timeseries_A.png', grid.arrange(arrangeGrob(p1),arrangeGrob(p2),arrangeGrob(p3),arrangeGrob(p4),arrangeGrob(p5),arrangeGrob(p6),arrangeGrob(p7),
                                                             layout_matrix=rbind(c(1,1,1,1,2,2,2,2,3,3,3,3),c(4,4,4,5,5,5,6,6,6,7,7,7))),width = 10, height = 6)
  # ggsave(filename = "../../Plots/InvasionModel/OneRun.pdf",p1)
  return()
}

Evolution_Plots <- function(a,b,d){
  d=3
  row1 = which(plotdata$selection_strength==s[3] & 
                 plotdata$mutation_rate == u[3] &
                 plotdata$propagule_pressure == p[10] & 
                 plotdata$invaded_popsize ==i[3])
  plotdf1 <- plotdata[row1]
  # plotdf1$Success = plotdf1$Frequency1000
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
  # plotdf2 = plotdf1[plotdf1$invader_popsize == plotdf1$invaded_popsize]
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
    xlab("Invader Community Size") +ylab("Invasion Success") + theme_classic() + geom_vline(xintercept = i[d],linetype=2) +
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
  ggsave(file = '../../Plots/Evolution_A.png',grid.arrange(p1,p2,p3,p4,leg
                                                           ,layout_matrix = rbind(c(1,1,1,2,2,2,NA),c(1,1,1,2,2,2,5),c(3,3,3,4,4,4,5),c(3,3,3,4,4,4,NA))),width = 13, height = 10)
  return( list(p1,p2,p3,p4))
}


Invader_Plots <- function(a,b,d){
  row1 = which(plotdata$selection_strength==s[3] & 
                 plotdata$mutation_rate == u[3] &
                 plotdata$invaded_popsize == i[3])
  plotdf1 <- plotdata[row1]
  options("scipen"=10) 
  tdf1 =  plotdf1[which(plotdf1$invader_generations == g[10] & plotdf1$invader_popsize %in% r[10]),]
  tdf2 =  plotdf1[which(plotdf1$invader_generations == g[2] & plotdf1$invader_popsize %in% r[2]),]
  tdf3 =  plotdf1[which(plotdf1$invader_generations == g[2] & plotdf1$invader_popsize %in% r[10]),]
  tbdf = rbind(tdf1,tdf2,tdf3)
  sum_df  =tbdf[,list(mean=mean(Frequency40),sd=sd(Frequency40)),by=list(invader_popsize,invader_generations,propagule_pressure,invaded_popsize)]
  sum_df$treatment = paste('R = ', sum_df$invader_popsize, ' G = ',sum_df$invader_generations,sep = '') 
  p1 <-  ggplot(sum_df,aes(x=as.numeric(propagule_pressure), y=mean, group =as.factor(treatment),colour=treatment)) + theme_classic()+
    geom_point(size= 4) + geom_line(size = 2)  +   scale_color_discrete(name="") + 
    theme(legend.position=c(0.7,0.45)) + scale_x_log10() + 
    theme(panel.border= element_blank()) + labs(x = 'Propagule Pressure (P)') +  theme(legend.position = c(0.7,0.6)) +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(y = 'Invasion Success',colour='Paramaters') + geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd,fill=as.factor(treatment)), alpha = 0.5,show.legend=FALSE)+ 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm")) 
  sum_df2 = plotdf1[,list(mean=mean(PP),sd=sd(PP),Relative_Propagule_Fitness = mean(Invader_Mean_End/Invaded_Mean_End)),by=list(invader_popsize,invader_generations,invaded_popsize)]
  p2 <- ggplot(sum_df2, aes(x = Relative_Propagule_Fitness,y = mean)) +  geom_point(size=2) + 
    xlab("Relative Invader Fitness") +ylab("Effect of Propagule Pressure") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + geom_vline(xintercept = 1,linetype=2,color='Red') +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(color='Invaded Community Size') + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  sum_df3  =tbdf[,list(mean=mean(Frequency40),sd=sd(Frequency40)) ,by=list(invader_popsize,invader_generations,Propagule_Richness,propagule_pressure,invaded_popsize)]
  sum_df3$treatment = paste('R = ', sum_df3$invader_popsize, ' G = ',sum_df3$invader_generations,sep = '') 
  # sum_df3$treatment = paste('R = ', sum_df3$invader_popsize, ' G = ',sum_df3$invader_generations,sep = '')
  # sum_df3 = subset(sum_df3,Propagule_Richness <=8)
  sum_df3 = sum_df3[which(sum_df3$propagule_pressure  == max(sum_df3$propagule_pressure)),]
  p3 <-  ggplot(sum_df3,aes(x=as.numeric(Propagule_Richness), y=mean, group =treatment,colour=treatment)) + theme_classic()+
    geom_point(size= 4) + geom_line(size = 2)  +   scale_color_discrete(name="") + 
    theme(legend.position=c(0.7,0.45)) + scale_x_log10(breaks=c(3,10,30),limits=c(3,35)) +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    theme(panel.border= element_blank()) + labs(x = 'Number of Invading Genotypes') +  theme(legend.position = c(0.7,0.6)) +
    labs(y = 'Invasion Success',colour='Paramaters') + geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd,fill=as.factor(treatment)), alpha = 0.5,show.legend=FALSE) + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm")) 
  sum_df4 = plotdf1[,list(mean=mean(PR),sd=sd(PR),Relative_Propagule_Fitness = mean(Invader_Mean_End/Invaded_Mean_End)),by=list(invader_popsize,invader_generations,invaded_popsize,propagule_pressure)]
  p4 <- ggplot(sum_df4, aes(x = Relative_Propagule_Fitness,y = mean)) +  geom_point(size=2) + 
    xlab("Relative Invader Fitness") +ylab("Effect of Invading Genotype Number") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + geom_vline(xintercept = 1,linetype=2,color='Red') +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(color='Invaded Community Size') + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  ggsave(file = '../../Plots/Invaders_A.png',grid.arrange(arrangeGrob(p1),arrangeGrob(p2),arrangeGrob(p3),arrangeGrob(p4)
                                                          ,layout_matrix = rbind(c(1,2),c(3,4))),width = 12, height = 10)
  return( list(p1,p2,p3,p4))
}



# Invaded_Plots <- function(a,b,d){
#   row1 = which(plotdata$selection_strength==s[3] & 
#                  plotdata$mutation_rate == u[3] &
#                  plotdata$invader_popsize == r[6] & 
#                  plotdata$invader_generations == g[6] &
#                  plotdata$propagule_pressure==p[5])
#   plotdf1 <- plotdata[row1]
#   sum_df <-aggregate(plotdf1,by=list(invader_popsize),function(x) mean(x))
#   a= plotdf1[,mean(Frequency40),by=list(invaded_popsize)]
#   
#   options("scipen"=10) 
# 
#   ggsave(file = '../../Plots/Invaded_A.png',grid.arrange(arrangeGrob(p1),arrangeGrob(p2),arrangeGrob(p3),arrangeGrob(p4)
#                                                           ,layout_matrix = rbind(c(1,2),c(3,4))),width = 12, height = 10)
#   
#   return( list(p1,p2))
# }

Recreation_plots <- function(){
  fig1 =read.csv('Jousset.csv')
  row1 = which(plotdata$selection_strength==s[1] & 
                 plotdata$mutation_rate == u[3] &
                 plotdata$invader_popsize == r[10] &
                 plotdata$invader_generations == g[1] &
                 plotdata$propagule_pressure == p[10])
  plotdf1 <- plotdata[row1]
  plotdf1=plotdf1[plotdf1$invaded_popsize<=10000]
  # p1 <- ggplot(plotdf1,aes(x=Invaded_Richness_End,y=Frequency100)) + geom_point()+ theme_classic() + scale_y_log10() + 
  #   geom_smooth(method='lm',formula=y~x) + annotate("text",x=10,y=0.3,label=as.character(expression(R^2==0.32)),parse=TRUE) + 
  #   annotate("text",x=13,y=0.3,label=as.character(expression(P<= 0.0001)),parse=TRUE)+ xlab('Genotypic Richness Resident Community') + ylab('Invasion Success (log_10_frequency)') +
  #   theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm"))
  p1a <- ggplot(plotdf1,aes(x=Invaded_Richness_End,y=Frequency40*100)) + geom_point()+ theme_classic() +
    scale_y_log10(breaks=c(0.1, 0.316,1,3.16,10,31.6,100),labels = c('-1','-0.5','0','0.5','1','1.5','2'),limits=c(0.1,100)) + 
    scale_x_continuous(breaks=c(1,3,5,7,9,11))  + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    geom_smooth(method='lm',formula=y~x) + #annotate("text",x=11,y=0.1,label=as.character(expression(R^2==0.43)),parse=TRUE) +
    annotate("text",x=9,y=0.1,label=as.character(expression(R^2*'='*0.59*'    '*P<= 1e-6)),parse=TRUE)+
    labs(x='Genotypic richness resident community',y=expression(atop('Invader Success', '(log'[10]*'frequency)'))) +
    theme(text = element_text(size=15),axis.text = element_text(size=15))
  p1b <- ggplot(fig1,aes(x=Genotyic.Richness,y=Invader.Success)) + geom_point()+ theme_classic() +
    scale_y_continuous(breaks=c(-1,-0.5,-0,0.5,1,1.5,2),labels = c('-1','-0.5','0','0.5','1','1.5','2')) + 
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8),limits=c(1,8))   + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    geom_smooth(method='lm',formula=y~x) + #annotate("text",x=11,y=0.1,label=as.character(expression(R^2==0.43)),parse=TRUE) +
    annotate("text",x=7,y=-1,label=as.character(expression(R^2*'='*0.16*'    '*P<= 0.01)),parse=TRUE)+
    labs(x='Genotypic richness resident community',y=expression(atop('Invader Success', '(log'[10]*'frequency)'))) +
    theme(text = element_text(size=15),axis.text = element_text(size=15))
  fig2 = read.csv('Acosta.csv')
  row2 = c(plotdata$selection_strength==s[1] & 
             plotdata$mutation_rate == u[3] &
             plotdata$invader_popsize == r[5] &
             plotdata$invader_generations == g[5] &
             plotdata$invaded_popsize == i[3]&
             plotdata$propagule_pressure %in% c(1,10,100))
  plotdf2 <- plotdata[row2]
  freq_df = rbind(cbind(melt(read_data_frequency2(unique(plotdf2$seed)[1],unique(plotdf2$task_number),500),id.vars='Generation',variable.name='Run',value.name='Invasion_Frequency'),'p' = 1),
                  cbind(melt(read_data_frequency2(unique(plotdf2$seed)[2],unique(plotdf2$task_number),500),id.vars='Generation',variable.name='Run',value.name='Invasion_Frequency'),'p' = 10),
                  cbind(melt(read_data_frequency2(unique(plotdf2$seed)[3],unique(plotdf2$task_number),500),id.vars='Generation',variable.name='Run',value.name='Invasion_Frequency'),'p' = 100))
  ag <- aggregate(. ~ p*Generation, freq_df, function(x) c(mean = mean(x), sd = sd(x)))
  ag$mean = ag$Invasion_Frequency[,1] *i[3]
  ag$sd = ag$Invasion_Frequency[,2]*i[3]
  ag$Treatment = NA
  ag$Treatment[ag$p==1] = 'Low'
  ag$Treatment[ag$p==10] = 'Medium'
  ag$Treatment[ag$p==100] = 'High'
  ag$Treatment = factor(ag$Treatment,levels=c('Low','Medium','High'))
  fig2$Treatment = factor(fig2$Treatment,levels=c('Low','Medium','High'))
  
  fig2$Time=fig2$Time-min(fig2$Time)
  p2a <-  ggplot(ag[ag$Generation%in%c(1,100,200),],aes(x = Generation,y=mean,colour=Treatment)) + 
    theme_classic()+ geom_line(size=2) +
    xlim(0,200) + geom_point(size=4)  +  
    labs(x = 'Time(Generations)',y =expression(atop('Invader Success', '(log'[10]*'frequency)')),colour='P') +
    scale_y_log10(breaks=c(0.1,1,10,100,1000),labels = c('0  ','1  ','2  ','3  ','4  '),limits=c(0.1,1000))+  
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    #+ #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd,fill=as.factor(p)), alpha = 0.5,show.legend=FALSE)+ 
    theme(legend.position = c(0.8,0.58),legend.text = element_text(size=12),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm")) +
    theme(text = element_text(size=15),axis.text = element_text(size=15))
  p2b <-  ggplot(fig2,aes(x = Time,y=Abundance.Log.,colour=as.factor(Treatment))) + 
    theme_classic()+ geom_line(size=2) + geom_point(size=4)+
    labs(x = 'Time(Days)',y = expression(atop('Invader Success', '(log'[10]*'frequency)')),colour='P') + 
    scale_y_continuous(breaks=c(0,1,2,3,4),labels = c('0  ','1  ','2  ','3  ','4  '))+  
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
  #+ #geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd,fill=as.factor(p)), alpha = 0.5,show.legend=FALSE)+ 
    theme(legend.position = c(0.8,0.58),legend.text = element_text(size=12),legend.key.size = unit(0.5, "cm"),legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.1, unit="cm")) +
    theme(text = element_text(size=15),axis.text = element_text(size=15))
  
  
  ggsave(file = '../../Plots/Recreation_A.png',grid.arrange(arrangeGrob(p1b,top=textGrob('A)',x=0.05,gp=gpar(fontsize=15))),
                                                            arrangeGrob(p2b,top=textGrob('C)',x=0.05,gp=gpar(fontsize=15))),
                                                            arrangeGrob(p1a,top=textGrob('B)',x=0.05,gp=gpar(fontsize=15))),
                                                            arrangeGrob(p2a,top=textGrob('D)',x=0.05,gp=gpar(fontsize=15))),ncol=2),
                                                            width = 12, height = 10)
}           
# 

Time_series_plot()
Evolution_Plots(2,2,2)
Invader_Plots(2,2,2)
Recreation_plots()
# Invaded_Plots(2,2,2)