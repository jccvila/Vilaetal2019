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


#calculate sloep pvalue for linear regression
pval <- function(x,y){
  pv = try(summary(lm(y~x))$coefficients[2, 4],silent=TRUE)
  if(!is.numeric(pv)){
    pv =1
  }
  return(pv)
}

slope2 <- function(x1,x2,y){
  slope = try(summary(lm(y/x2~x1))$coefficients[2, 1],silent=TRUE)
  if(!is.numeric(slope)){
    slope =0
  }
  return(slope)
}

#calculate sloep pvalue for linear regression
pval2 <- function(x1,x2,y){
  pv = try(summary(lm(y/x2~x1))$coefficients[2, 4],silent=TRUE)
  if(!is.numeric(pv)){
    pv =1
  }
  return(pv)
}

# extract legend from ggplot for multiplotting
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


### read_data from HPC reuns


path = "../Results/SIA_Results/"

# extends a vector to a given length by repeat the final element n times (i.e if fixation has achieved all subsequent points)
extend <- function(v,max_length){
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

#The HPC outputs a very large number of small files which i have to merege into an enormous data.set only need to run this once.
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
load('Additive_Data.rds')
#Fix a deprecated issue with naming (invaded to resident)
#Fix naming of resident generation to max invasion time9 to reflect what it actually means
colnames(plotdata) = gsub('resident_generations','max_invasion_time',
                          gsub('Invaded','Resident',
                               gsub('invaded','resident',colnames(plotdata))))

#resident generations
# #
plotdata$Relative_Fitness = plotdata$Invader_Mean_End/plotdata$Resident_Mean_End
plotdata$Relative_Propagule_Fitness = plotdata$Propagule_Fitness_Mean/plotdata$Resident_Mean_End
plotdata$Fitness_Increase_Resident = plotdata$Invader_Mean_End - plotdata$Invader_Mean_Start
plotdata$Fitness_Increase_Invader= plotdata$Resident_Mean_End - plotdata$Resident_Mean_Start
plotdata$Adaptation_Rate_Invader = plotdata$Fitness_Increase_Resident/plotdata$invader_generations
plotdata$Adaptation_Rate_Resident = plotdata$Fitness_Increase_Invader/plotdata$invader_generations
plotdata$Adaptation_Rate_Difference= plotdata$Adaptation_Rate_Invader - plotdata$Adaptation_Rate_Resident
plotdata$Relative_Propagule_Fitness = plotdata$Propagule_Fitness_Mean/plotdata$Resident_Mean_End
#
# # # generic ###
# # ### single param combinations ###
#
s = unique(plotdata$selection_strength)
u = unique(plotdata$mutation_rate)
p = unique(plotdata$propagule_pressure)
r <- unique(plotdata$resident_popsize)
i <- unique(plotdata$invader_popsize)
g <- unique(plotdata$invader_generations)


plotdata[,PP:=slope(propagule_pressure,Frequency40),by=list(selection_strength,mutation_rate,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PPRS:=rsquared(propagule_pressure,Frequency40),by=list(selection_strength,mutation_rate,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PPpval:=pval(propagule_pressure,Frequency40),by=list(selection_strength,mutation_rate,resident_popsize,invader_popsize,invader_generations)]

plotdata[,PR:=slope(Propagule_Richness,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PRRS:=rsquared(Propagule_Richness,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PRpval:=pval(Propagule_Richness,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]


plotdata[,PV:=slope(Propagule_Fitness_Variance,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PVRS:=rsquared(Propagule_Fitness_Variance,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PVpval:=pval(Propagule_Fitness_Variance,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]


plotdata[,PR_PM:=slope(Propagule_Richness,Propagule_Fitness_Mean),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PR_PMRS:=rsquared(Propagule_Richness,Propagule_Fitness_Mean),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PR_PMpval:=pval(Propagule_Richness,Propagule_Fitness_Mean),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]

plotdata[,PR_PV:=slope(Propagule_Richness,Propagule_Fitness_Variance),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PR_PVRS:=rsquared(Propagule_Richness,Propagule_Fitness_Variance),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PR_PVpval:=pval(Propagule_Richness,Propagule_Fitness_Variance),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]

plotdata[,PR2pval:=pval2(Propagule_Richness,Propagule_Fitness_Variance,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PR2slope:=slope2(Propagule_Richness,Propagule_Fitness_Variance,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]

plotdata[,PR3pval:=pval2(Propagule_Richness,Propagule_Fitness_Mean,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,PR3slope:=slope2(Propagule_Richness,Propagule_Fitness_Mean,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]

plotdata[,IS:=slope(resident_popsize,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,invader_popsize,invader_generations)]
plotdata[,ISRS:=rsquared(resident_popsize,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,invader_popsize,invader_generations)]
plotdata[,ISpval:=pval(resident_popsize,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,invader_popsize,invader_generations)]

plotdata[,IR:=slope(Resident_Richness_End,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,IRRS:=rsquared(Resident_Richness_End,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]
plotdata[,IRpval:=pval(Resident_Richness_End,Frequency40),by=list(selection_strength,mutation_rate,propagule_pressure,resident_popsize,invader_popsize,invader_generations)]


#
save(plotdata,file = 'Additive_Data3.rds')
fwrite(plotdata,'Summary_Data.csv')