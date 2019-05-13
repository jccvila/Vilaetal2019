rm(list=ls())
plot.new()
dev.off()

library("ggplot2")
library("reshape2")
library("gridExtra")

library(grid)
library(lattice)
### core functions###

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

## Multiple plot function ###r
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## sum vects missing values are 0
sum_vect <- function(x,y){
  if(length(x) <= length(y)){
    length(x) <- length(y)
    x[is.na(x)] <- 0
  }
  else {
    length(y) <- length(x)
    y[is.na(y)] <- 0} 
  sum <- x + y
  return(sum)
}

## sum vect missing values are extended
sum_vect2 <- function(x,y){
  oldlengthy <- length(y)
  oldlengthx <- length(x)
  if(length(x) <= length(y)){
    length(x) <- length(y)
    x[is.na(x)] <- x[oldlengthx]
  }
  else {
    length(y) <- length(x)
    y[is.na(y)] <- y[oldlengthy]} 
  sum <- x + y
  return(sum)
}

###read _data
path = "../../Results/Invasion/Summary_Invasion/"

read_data_paramaters <- function(seed,task){
  params <- read.csv(paste(path,"Paramaters_",seed,"_", task,".csv",sep=""),check.names = FALSE,header=FALSE,col.names = c("Paramater","Value"))
  return(params)
}


read_data_frequency <- function(seed,task){
  Invasion_path = paste(path,"Invador_Frequency_",sep="")
  Inv_freq <- strsplit(readLines(paste(Invasion_path,seed,"_", task,".csv",sep="")),split=",") 
  lengths <- vector()
  for(i in 1:length(Inv_freq)){lengths <- c(lengths,length(Inv_freq[[i]]))}
  max_length <- max(lengths)
  for(i in 1:length(Inv_freq)){
    if(lengths[i] == max_length){next}
    else{
      finalfreq = Inv_freq[[i]][lengths[i]]
      addfreq = rep(finalfreq,max_length - lengths[i])
      Inv_freq[[i]] = as.vector(as.numeric(c(Inv_freq[[i]],addfreq)))
    }
  }
  matrix <- t(matrix(unlist(Inv_freq), nrow=length(Inv_freq), byrow=T))
  suppressWarnings(storage.mode(matrix) <- "numeric")
  matrix <- as.data.frame(cbind(matrix,seq(1:nrow(matrix))))
  colnames(matrix) <- c(seq(1:(ncol(matrix)-1)),"Generation")
  dataframe <- melt(matrix,id.vars = c("Generation"),variable.name = "Run", value.name = "Invador_Frequency")
  return(dataframe)
}


read_data_start_fitness <- function(seed,task){
  start_fit <- as.numeric(read.csv(paste(path,"Start_Fitness_",seed,"_", task,".csv",sep=""),check.names = FALSE, header = FALSE))
  return(start_fit)
}


read_data_end_fitness <- function(seed,task){
  start_fit <- as.numeric(read.csv(paste(path,"End_Fitness_",seed,"_", task,".csv",sep=""),check.names = FALSE, header = FALSE))
  return(start_fit)
}

read_data_propagule_fitness <- function(seed,task){
  prop_fit <- as.numeric(read.csv(paste(path,"Propagule_Fitness_",seed,"_", task,".csv",sep=""),check.names = FALSE, header = FALSE))
  return(prop_fit)
}

read_data_propagule_richness <- function(seed,task){
  prop_rich <- as.numeric(read.csv(paste(path,"Propagule_Richness_",seed,"_", task,".csv",sep=""),check.names = FALSE, header = FALSE))
  return(prop_rich)
}


# generic ### 

# function takes abundance vector and converts it into popsize.

# plots invaodr frequency for single run
plot_invador_frequency <- function(seed,task,run){
  InvadorFrequency <- read_data_frequency(seed,task)
  subset <- InvadorFrequency[which(InvadorFrequency$Run == run),]
  plot(subset$Generation,subset$Invador_Frequency,type="l",ylim=c(0,1))
}

### load paramaters

loadparamaters <- function(max_seed,max_task){
  paramatersdf <- data.frame()
    counter = 1
    progressbar <- txtProgressBar(min = 0, max = max_seed * max_task, style = 3)
    for (seed in 1:max_seed){
      possibleError <- tryCatch(read_data_frequency(seed,1),error = function(e) e)
      if(inherits(possibleError, "error")){next}
    for (task in 1:max_task){
      possibleError <- tryCatch(read_data_frequency(seed,task),error = function(e) e)
      if(inherits(possibleError, "error")){break}
      params <- read_data_paramaters(seed,task)
      Seed <- params[1,2]
      Task <-params[2,2]
      Selection_strength <- params[3,2]
      Mutation_rate <- params[4,2]
      Propagule_pressure <- params[5,2]
      Home_popsize <- params[6,2]
      Invador_popsize <- params[8,2]
      Invador_generations <- params[9,2]
      Start_fitness <- mean(read_data_start_fitness(seed,task))
      End_fitness <- mean(read_data_end_fitness(seed,task))
      Propagule_fitness <- mean(read_data_propagule_fitness(seed,task)/Propagule_pressure)
      Propagule_richness <- mean(read_data_propagule_richness(seed,task))
      paramatersdf <- rbind(paramatersdf,data.frame(Seed,Task,Home_popsize,Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Start_fitness,End_fitness,Propagule_fitness,Propagule_richness))
      counter <- counter+1
      setTxtProgressBar(progressbar,(seed-1)*max_task + task)
    }
    }
  paramatersdf$Propagule_proportion <- paramatersdf$Propagule_pressure/paramatersdf$Home_popsize
  paramatersdf$Relative_Fitness <- paramatersdf$End_fitness/paramatersdf$Start_fitness
  
  return(paramatersdf)
}

### single param combinations ### 


# sum across runs with same paramaters and returns summarised version data
sum_stat_invador_frequency <- function(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Home_popsize){
  row = which(finaldf$Invador_popsize == Invador_popsize & finaldf$Invador_generations == Invador_generations & finaldf$Mutation_rate == Mutation_rate & finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure & finaldf$Home_popsize == Home_popsize)
  seed = finaldf$Seed[row]
  task = finaldf$Task[row]
  dataframe <- read_data_frequency(seed,task)
  return(dataframe)
}

# calculates mean invador_frequency for first 1000 generations

relative_propagule_fitness <- function(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Home_popsize){}

proportion_calculator <- function(data){
  Endsim <- max(unique(data$Generation))
  EndsimData <- data[which(data$Generation == Endsim),]
  n <- nrow(EndsimData)
  proportion  <- sum(EndsimData$Invador_Frequency)/n
  SE <- 1.96 * sqrt((proportion * (1 - proportion))/n)
  return(data.frame(Mean = proportion, Upper_CI = proportion +SE, Lower_CI = proportion - SE))
}

success_calculator <- function(data){
  if(max(unique(data$Generation)) >= 1000){big = TRUE}
  if(max(unique(data$Generation)) < 1000){big = FALSE}
  mean <- aggregate(Invador_Frequency ~  Generation, data = data, FUN = function (x) mean(x))
  min <- aggregate(Invador_Frequency ~  Generation, data = data, FUN = function (x) min(x))
  max <- aggregate(Invador_Frequency ~  Generation, data = data, FUN = function (x) max(x))
  if(big == TRUE){
    mean <- mean[1:1000,]
    min <- min[1:1000,]
    max <- max[1:1000,]
  }
  if(big == FALSE){
    finalmean = mean[nrow(mean),]
    finalmin = min[nrow(min),]
    finalmax = max[nrow(max),]
    addmean = data.frame(finalmean[rep(1,1000 - nrow(mean)),])
    addmin = data.frame(finalmin[rep(1,1000 - nrow(min)),])
    addmax = data.frame(finalmax[rep(1,1000 - nrow(max)),])
    mean = rbind(mean,addmean)
    min = rbind(min,addmin)
    max = rbind(max,addmax)
  }
  
  return(data.frame(Mean = sum(mean$Invador_Frequency), Max = sum(max$Invador_Frequency), Min = sum(min$Invador_Frequency)))
}

gradient_calculator <- function(data){
  mean <- aggregate(Invador_Frequency ~  Generation, data = data, FUN = function (x) mean(x))
  min <- aggregate(Invador_Frequency ~  Generation, data = data, FUN = function (x) min(x))
  max <- aggregate(Invador_Frequency ~  Generation, data = data, FUN = function (x) max(x))
  mean <- mean[12,]
  min <- min[1:2,]
  max <- max[1:2,]
  return(data.frame(Mean = mean$Invador_Frequency[2] - mean$Invador_Frequenc[1], Max = max$Invador_Frequency[2] - max$Invador_Frequenc[1],Min = min$Invador_Frequency[2] - min$Invador_Frequenc[1]))
}


proportion_updator <- function(dataframe){
  #   Invasion_success <- vector
  Proportion <- data.frame()
  progressbar <- txtProgressBar(min = 0, max = nrow(paramatersdf), style = 3)
  for(i in 1:nrow(dataframe )){
    data =  read_data_frequency(paramatersdf$Seed[i],paramatersdf$Task[i])
    #     Temp_Invasion_success = Invasion_success(data)$Mean
    Temp_Proportion = proportion_calculator(data)
    Proportion = rbind(Proportion,Temp_Proportion)
    setTxtProgressBar(progressbar,i)
    
  }
  dataframe <- cbind(dataframe,Proportion)
  return(dataframe)
}

success_updator <- function(dataframe){
  #   Invasion_success <- 
  Proportion <- data.frame()
  progressbar <- txtProgressBar(min = 0, max = nrow(paramatersdf), style = 3)
  for(i in 1:nrow(dataframe)){
    data =  read_data_frequency(paramatersdf$Seed[i],paramatersdf$Task[i])
    #     Temp_Invasion_success = Invasion_success(data)$Mean
    Temp_Proportion = success_calculator(data)
    Proportion = rbind(Proportion,Temp_Proportion)
    setTxtProgressBar(progressbar,i)
    
  }
  dataframe <- cbind(dataframe,Proportion)
  return(dataframe)
}

gradient_updator <- function(dataframe){
  #   Invasion_success <- 
  Proportion <- data.frame()
  progressbar <- txtProgressBar(min = 0, max = nrow(paramatersdf), style = 3)
  for(i in 1:nrow(dataframe)){
    data =  read_data_frequency(paramatersdf$Seed[i],paramatersdf$Task[i])
    #     Temp_Invasion_success = Invasion_success(data)$Mean
    Temp_Proportion = gradient_calculator(data)
    Proportion = rbind(Proportion,Temp_Proportion)
    setTxtProgressBar(progressbar,i)
    
  }
  dataframe <- cbind(dataframe,Proportion)
  return(dataframe)
}

# invasion succes - area under curve for 1000 generations.
# invasion success - 

invfreqvstime <- function(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Home_popsize){
  df_freq <- sum_stat_invador_frequency(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Home_popsize)
  p1 <- ggplot(df_freq, aes(x = Generation, y = Invador_Frequency,col=as.factor(Run))) + geom_line(size=1) + ylab("Invador_Frequency") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position="none") +
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1)) +
    xlim(0,2000) + ylim(0,1) +  theme(text = element_text(size=26)) + ylab("Invader Frequency")
  ggsave(filename = "../../Plots/InvasionModel/OneRun.pdf",p1)
  return()
}

invfreqvstime2 <- function(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Home_popsize){
  data <- sum_stat_invador_frequency(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Home_popsize)
  mean <- aggregate(Invador_Frequency ~  Generation, data = data, FUN = function (x) mean(x))
  min <- aggregate(Invador_Frequency ~  Generation, data = data, FUN = function (x) min(x))
  max <- aggregate(Invador_Frequency ~  Generation, data = data, FUN = function (x) max(x))
  df_freq = data.frame(Generation = mean$Generation, Mean = mean$Invador_Frequency, Min = min$Invador_Frequency, Max = max$Invador_Frequency)
  p1 <- ggplot(df_freq, aes(x = Generation, y = Mean)) + geom_line(size=1) + ylab("Mean Fitness Category") + theme_classic()+ 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(panel.border= element_blank()) + theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1)) + 
    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = 0.5)
  return(p1)
}

Invasionvsgenerations <- function(Invador_popsize,Mutation_rate,Selection_strength,Propagule_pressure){
  plotdf = finaldf[which(finaldf$Invador_popsize  == Invador_popsize & finaldf$Mutation_rate == Mutation_rate & finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure),]
  p1 <- ggplot(plotdf, aes(x = Invador_generations, y = Proportion,color = as.factor(Home_popsize))) +  geom_point(size=2) + geom_line(size=1) + 
    xlab(expression( "Invader Generations g"[r]) ) +ylab("Proportion Succesful Invasion") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + scale_x_log10() +
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color=expression("J"[d])) + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"))

  p2 <- ggplot(plotdf, aes(x = Invador_generations, y = Invasion_success,color = as.factor(Home_popsize))) +  geom_point(size=2) + geom_line(size=1) + 
    xlab(expression( "Invader Generations g"[r]) ) +ylab("Success") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color=expression("J"[d])) + scale_x_log10() +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"))  
  return(list(p1,p2))
}

Invasionvspopsize <- function(Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure){
  plotdf = finaldf[which(finaldf$Invador_generations  == Invador_generations & finaldf$Mutation_rate == Mutation_rate & finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure),]
  p1 <- ggplot(plotdf, aes(x = Invador_popsize, y = Proportion,color = as.factor(Home_popsize))) +  geom_point(size=2) + geom_line(size=1) + 
    xlab(expression( "Invader Community Size J"[r]) ) +ylab("Proportion Succesful Invasion") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color=expression("J"[d]))+ scale_x_log10() +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"))  
  p3 <- ggplot(plotdf, aes(x = Invador_popsize, y = Invasion_success,color = as.factor(Home_popsize))) +  geom_point(size=2) + geom_line(size=1) + 
    xlab(expression( "Invader Community Size J"[r]) ) +ylab("Success") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color=expression("J"[d])) + scale_x_log10() +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"))  
  return(list(p1,p3))
}

Invasionvspropagulepressure <- function(Invador_generations,Invador_popsize,Mutation_rate,Selection_strength){
  plotdf = finaldf[which(finaldf$Invador_generations  == Invador_generations & finaldf$Mutation_rate == Mutation_rate & finaldf$Selection_strength == Selection_strength & finaldf$Invador_popsize == Invador_popsize),]
  p1 <- ggplot(plotdf, aes(x = Propagule_pressure, y = Proportion,color = as.factor(Home_popsize))) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Propagule Pressure p") +ylab("Proportion Succesful Invasion") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color=expression("J"[d])) + scale_x_log10() +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"))
  p3 <- ggplot(plotdf, aes(x = Propagule_pressure, y = Invasion_success,color = as.factor(Home_popsize))) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Propagule Pressure p")  +ylab("Success") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color=expression("J"[d])) + scale_x_log10() +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"))
  return(list(p1,p3))
}



Invasionvsrelativefitness <- function(Mutation_rate,Selection_strength,Home_popsize){
  plotdf = finaldf[which(finaldf$Mutation_rate == Mutation_rate & finaldf$Selection_strength == Selection_strength & finaldf$Home_popsize == Home_popsize),]
  p1 <- ggplot(plotdf, aes(x = End_fitness, y = Proportion,color = as.factor(Propagule_pressure))) +  geom_point(size=2)  + 
    xlab("Relative Propagule Fitness") +ylab("Proportion Succesful Invasion") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='p') + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  p3 <- ggplot(plotdf, aes(x = End_fitness, y = Invasion_success,color = as.factor(Invador_popsize))) +  geom_point(size=2) + 
    xlab("Relative Propagule Fitness") +ylab("Success") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='p') + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  return( list(p1,p3))
}

Internalplots <- function(Mutation_rate,Selection_strength,Invador_generations,Invador_popsize,Home_popsize,Propagule_pressure){
  paramsdf = finaldf[which(finaldf$Mutation_rate == Mutation_rate & finaldf$Selection_strength == Selection_strength & finaldf$Home_popsize == Home_popsize & finaldf$Invador_generations == Invador_generations & finaldf$Invador_popsize == Invador_popsize),]
  plotdf <- data.frame()
  listproppres <- unique(paramsdf$Propagule_pressure)[c(1,5,10)]
  for (j in listproppres){
   Propagule_pressure = j
   paramsdf2 <- paramsdf[which(paramsdf$Propagule_pressure  == Propagule_pressure),]
   Frequency <- read_data_frequency(paramsdf2$Seed,paramsdf2$Task)
   Frequency <- subset(Frequency,Frequency$Generation <= 1000 )
  SortedDF <- data.frame()
  for(i in 1:length(unique(Frequency$Run))){
    subsetfreq <- subset(Frequency,Frequency$Run == i)
    if(nrow(subsetfreq) <= 1000){SortedDF <- rbind(SortedDF,subsetfreq)}
    else{
      newdf <- data.frame(subsetfreq[rep(1,1000 - nrow(subsetfreq)),])
      SortedDF = rbind(SortedDF,rbind(subsetfreq,newdf))
    }
  }
  Richness <- read_data_propagule_richness(paramsdf2$Seed,paramsdf2$Task)
  Relative_Fitness <- (read_data_propagule_fitness(paramsdf2$Seed,paramsdf2$Task)/Propagule_pressure)/read_data_start_fitness(paramsdf2$Seed,paramsdf2$Task)
  Adaptation_rate <- (read_data_end_fitness(paramsdf2$Seed,paramsdf2$Task)/read_data_start_fitness(paramsdf2$Seed,paramsdf2$Task))/read_data_paramaters(201,1)[[9,2]]
  
  aggregated <- aggregate(Invador_Frequency ~  Run, data = SortedDF, FUN = function (x) sum(x))
  plotdf <- rbind(plotdf,cbind(aggregated,Richness,Relative_Fitness,Propagule_pressure,Adaptation_rate,Richness))
  }
  p1 <- ggplot(plotdf, aes(x = Relative_Fitness, y = Invador_Frequency,color = as.factor(Propagule_pressure))) +  geom_point(size=3)  + 
    xlab("Relative Propagule Fitness") +ylab("Sucess") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='p') + scale_colour_manual(values = c("Yellow","Green","Purple")) +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  p2 <- ggplot(plotdf, aes(x = Adaptation_rate, y = Relative_Fitness,color = as.factor(Propagule_pressure))) +  geom_point(size=3)  + 
    xlab("Adaptation Rate") +ylab("Relative Propagule Fitnesss") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='p') + scale_colour_manual(values = c("Yellow","Green","Purple")) +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  paramsdf = finaldf[which(finaldf$Mutation_rate == Mutation_rate & finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure & finaldf$Invador_generations == Invador_generations & finaldf$Invador_popsize == Invador_popsize),]
  plotdf <- data.frame()
  listhomepop <- unique(paramsdf$Home_popsize)
  for (j in listhomepop){
    Home_popsize = j
    paramsdf2 <- paramsdf[which(paramsdf$Home_popsize == Home_popsize),]
    Frequency <- read_data_frequency(paramsdf2$Seed,paramsdf2$Task)
    Frequency <- subset(Frequency,Frequency$Generation <= 1000 )
    SortedDF <- data.frame()
    for(i in 1:length(unique(Frequency$Run))){
      subsetfreq <- subset(Frequency,Frequency$Run == i)
      if(nrow(subsetfreq) <= 1000){SortedDF <- rbind(SortedDF,subsetfreq)}
      else{
        newdf <- data.frame(subsetfreq[rep(1,1000 - nrow(subsetfreq)),])
        SortedDF = rbind(SortedDF,rbind(subsetfreq,newdf))
      }
    }
    Richness <- read_data_propagule_richness(paramsdf2$Seed,paramsdf2$Task)
    Relative_Fitness <- (read_data_propagule_fitness(paramsdf2$Seed,paramsdf2$Task)/Propagule_pressure)/read_data_start_fitness(paramsdf2$Seed,paramsdf2$Task)
    Adaptation_rate <- (read_data_end_fitness(paramsdf2$Seed,paramsdf2$Task)/read_data_start_fitness(paramsdf2$Seed,paramsdf2$Task))/read_data_paramaters(201,1)[[9,2]]
    aggregated <- aggregate(Invador_Frequency ~  Run, data = SortedDF, FUN = function (x) sum(x))
    gradient <- summary(lm(aggregated$Invador_Frequency~Richness))$coefficients[[2]]
    plotdf <- rbind(plotdf,cbind(aggregated,Richness,Relative_Fitness,Home_popsize,gradient, Adaptation_rate))
  }
  plotdf$Label = paste("m = ", plotdf$gradient)
  p3 <- ggplot(plotdf, aes(x = Richness, y = Invador_Frequency,color = as.factor(Home_popsize))) +  geom_point(size=3)  + 
    xlab("Propagule Richness") +ylab("Sucess") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color=expression("J"[d])) +
    geom_smooth(method="lm",formula = y ~ x,se = FALSE,fullrange = TRUE) +
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  return( list(p1,p2,p3))
}

richness_coefficient_calculator <- function(seed,task){
  Propagule_Richness <- read_data_propagule_richness(seed , task)
  Invasion_freq_df <- read_data_frequency(seed,task)
  Invasion_Outcome <- vector()
  for(j in unique(Invasion_freq_df$Run)){
    subset <- Invasion_freq_df[which(Invasion_freq_df$Run == j),]
    subset_size =nrow(subset)
    if(subset_size >= 1000){
      sum_freq <- sum(subset$Invador_Frequency[1:1000])
      Invasion_Outcome <- c(Invasion_Outcome,sum_freq)
    }
    else{
      sum_freq <- sum(subset$Invador_Frequency[1:subset_size]) + sum(rep(subset$Invador_Frequency[subset_size],1000 - subset_size))
      Invasion_Outcome <- c(Invasion_Outcome, sum_freq)
    }
  }
  coefficient = summary(lm(Invasion_Outcome ~ Propagule_Richness))$coefficients[[2]]
  return(coefficient)
}


Invasionvspropagulerichness <- function(Mutation_rate,Selection_strength,Home_popsize,Propagule_pressure){
  # plotdf = finaldf[which(finaldf$Mutation_rate == Mutation_rate & finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure),]
  # plotdf = finaldf
  # progressbar <- txtProgressBar(min = 0, nrow(plotdf), style = 3)
  # for(i in 1:nrow(plotdf)){
  #   setTxtProgressBar(progressbar,i)
  #   plotdf$Propagule_richness_correlation[i] = richness_coefficient_calculator(plotdf$Seed[i] , plotdf$Task[i])
  # }
  # save(plotdf,file = "RichnessCoefficients.Rdata")
  load(file = "RichnessCoefficients.Rdata")
  plotdf = subset(plotdf,plotdf$Home_popsize != 750000)
  p1 <- ggplot(plotdf,aes(Propagule_richness_correlation,col = as.factor(Home_popsize))) + geom_freqpoly(binwidth = 5,size=1) +
    xlab("Slope") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color=expression("J"[d])) + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  plotdf2 = subset(plotdf,plotdf$Propagule_pressure == c(21,67,222))
  p2 <- ggplot(plotdf2,aes(Propagule_richness_correlation,col = as.factor(Propagule_pressure))) + geom_freqpoly(binwidth = 5,size=1) +
    xlab("Slope") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(colour = "p") + 
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm")) 
  
  return(list(p1,p2))
}


ProportionPlots <- function(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure){
  Home_popsize = unique(finaldf$Home_popsize)[[1]]
  Popsize_Generations = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure & finaldf$Home_popsize == Home_popsize),]
  Popsize_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_generations == Invador_generations & finaldf$Home_popsize == Home_popsize),]
  Generations_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_popsize == Invador_popsize & finaldf$Home_popsize == Home_popsize),]
  PropaguleFitness_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Home_popsize),]
  p1 <-  ggplot(Popsize_Generations,aes(x=as.factor(Invador_generations), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Proportion)) +
    xlab(expression( "Invader Generations g"[r]) ) +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p2 <-  ggplot(Popsize_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Proportion)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p3 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Proportion)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Generations g"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  
  Home_popsize = unique(finaldf$Home_popsize)[[2]]
  Popsize_Generations = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure & finaldf$Home_popsize == Home_popsize),]
  Popsize_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_generations == Invador_generations & finaldf$Home_popsize == Home_popsize),]
  Generations_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_popsize == Invador_popsize & finaldf$Home_popsize == Home_popsize),]
  PropaguleFitness_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Home_popsize),]
  p4 <-  ggplot(Popsize_Generations,aes(x=as.factor(Invador_generations), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Proportion)) +
    xlab(expression( "Invader Generations g"[r]) ) +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p5 <-  ggplot(Popsize_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Proportion)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p6 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Proportion)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Generations g"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  
  Home_popsize = unique(finaldf$Home_popsize)[[3]]
  Popsize_Generations = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure & finaldf$Home_popsize == Home_popsize),]
  Popsize_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_generations == Invador_generations & finaldf$Home_popsize == Home_popsize),]
  Generations_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_popsize == Invador_popsize & finaldf$Home_popsize == Home_popsize),]
  PropaguleFitness_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Home_popsize),]
  p7 <-  ggplot(Popsize_Generations,aes(x=as.factor(Invador_generations), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Proportion)) +
    xlab(expression( "Invader Generations g"[r]) ) +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p8 <-  ggplot(Popsize_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Proportion)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p9 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Proportion)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Generations g"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  
  p10 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Proportion)) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm")) +
    theme(axis.title=element_text(size=18),legend.title=element_text(size=24),legend.title.align = -2)
  
  
  a <-  arrangeGrob(p1, left = textGrob("a)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)),
                    top = textGrob(bquote("J"[d] == .(unique(finaldf$Home_popsize)[1])),
                                   gp=gpar(col="black", fontsize=24,fontface="bold")))
  b <-  arrangeGrob(p2, left = textGrob("b)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  c <-  arrangeGrob(p3, left = textGrob("c)", x = unit(0, "npc"),
                                        y   = unit(1, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  
  d <-  arrangeGrob(p4, left = textGrob("d)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)),
                    top = textGrob(bquote("J"[d] == .(unique(finaldf$Home_popsize)[2])),
                                   gp=gpar(col="black", fontsize=24,fontface="bold")))
  e <-  arrangeGrob(p5, left = textGrob("e)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  f <-  arrangeGrob(p6, left = textGrob("f)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  
  g <-  arrangeGrob(p7, left = textGrob("g)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)),
                    top = textGrob(bquote("J"[d] == .(unique(finaldf$Home_popsize)[3])),
                                   gp=gpar(col="black", fontsize=24,fontface="bold")))
  h <-  arrangeGrob(p8, left = textGrob("h)", x = unit(0, "npc")
                                        , y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  i <-  arrangeGrob(p9, left = textGrob("i)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  
  legend <- g_legend(p10)
  lay = rbind(c(1,1,1,2,2,2,3,3,3,NA),c(1,1,1,2,2,2,3,3,3,NA),c(1,1,1,2,2,2,3,3,3,NA),
              c(4,4,4,5,5,5,6,6,6,10),c(4,4,4,5,5,5,6,6,6,10),c(4,4,4,5,5,5,6,6,6,10),
              c(7,7,7,8,8,8,9,9,9,NA),c(7,7,7,8,8,8,9,9,9,NA),c(7,7,7,8,8,8,9,9,9,NA))
  pdf(paste("../../Plots/InvasionModel/ProportionPlots.pdf"),height=18.5,width=22)
  grid.arrange(a,d,g,b,e,h,c,f,i,legend,layout_matrix = lay,top = "")
  dev.off()
  return()
}


SuccessPlots <- function(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure){
  Home_popsize = unique(finaldf$Home_popsize)[[1]]
  Popsize_Generations = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure & finaldf$Home_popsize == Home_popsize),]
  Popsize_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_generations == Invador_generations & finaldf$Home_popsize == Home_popsize),]
  Generations_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_popsize == Invador_popsize & finaldf$Home_popsize == Home_popsize),]
  PropaguleFitness_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Home_popsize),]
  p1 <-  ggplot(Popsize_Generations,aes(x=as.factor(Invador_generations), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab(expression( "Invader Generations g"[r]) ) +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Invasion Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p2 <-  ggplot(Popsize_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Invasion Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p3 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Generations g"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Invasion Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  
  Home_popsize = unique(finaldf$Home_popsize)[[2]]
  Popsize_Generations = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure & finaldf$Home_popsize == Home_popsize),]
  Popsize_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_generations == Invador_generations & finaldf$Home_popsize == Home_popsize),]
  Generations_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_popsize == Invador_popsize & finaldf$Home_popsize == Home_popsize),]
  PropaguleFitness_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Home_popsize),]
  p4 <-  ggplot(Popsize_Generations,aes(x=as.factor(Invador_generations), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab(expression( "Invader Generations g"[r]) ) +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Invasion Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p5 <-  ggplot(Popsize_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Invasion Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p6 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Generations g"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Invasion Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  
  Home_popsize = unique(finaldf$Home_popsize)[[3]]
  Popsize_Generations = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure & finaldf$Home_popsize == Home_popsize),]
  Popsize_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_generations == Invador_generations & finaldf$Home_popsize == Home_popsize),]
  Generations_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_popsize == Invador_popsize & finaldf$Home_popsize == Home_popsize),]
  PropaguleFitness_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Home_popsize),]
  p7 <-  ggplot(Popsize_Generations,aes(x=as.factor(Invador_generations), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab(expression( "Invader Generations g"[r]) ) +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Invasion Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p8 <-  ggplot(Popsize_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Invasion Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p9 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Generations g"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Invasion Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  
  p10 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Invasion_success)) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) +
    labs(fill='Success') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm")) +
    theme(axis.title=element_text(size=18),legend.title=element_text(size=24),legend.title.align = -2)

  
  
  a <-  arrangeGrob(p1, left = textGrob("a)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)),
                    top = textGrob(bquote("J"[d] == .(unique(finaldf$Home_popsize)[1])),
                                   gp=gpar(col="black", fontsize=24,fontface="bold")))
  b <-  arrangeGrob(p2, left = textGrob("b)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  c <-  arrangeGrob(p3, left = textGrob("c)", x = unit(0, "npc"),
                                        y   = unit(1, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  
  d <-  arrangeGrob(p4, left = textGrob("d)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)),
                    top = textGrob(bquote("J"[d] == .(unique(finaldf$Home_popsize)[2])),
                                   gp=gpar(col="black", fontsize=24,fontface="bold")))
  e <-  arrangeGrob(p5, left = textGrob("e)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  f <-  arrangeGrob(p6, left = textGrob("f)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  
  g <-  arrangeGrob(p7, left = textGrob("g)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)),
                    top = textGrob(bquote("J"[d] == .(unique(finaldf$Home_popsize)[3])),
                                   gp=gpar(col="black", fontsize=24,fontface="bold")))
  h <-  arrangeGrob(p8, left = textGrob("h)", x = unit(0, "npc")
                                        , y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  i <-  arrangeGrob(p9, left = textGrob("i)", x = unit(0, "npc"),
                                        y   = unit(0.95, "npc"), just=c("left","top"),
                                        gp=gpar(col="black", fontsize=24)))
  
  legend <- g_legend(p10)
  lay = rbind(c(1,1,1,2,2,2,3,3,3,NA),c(1,1,1,2,2,2,3,3,3,NA),c(1,1,1,2,2,2,3,3,3,NA),
              c(4,4,4,5,5,5,6,6,6,10),c(4,4,4,5,5,5,6,6,6,10),c(4,4,4,5,5,5,6,6,6,10),
              c(7,7,7,8,8,8,9,9,9,NA),c(7,7,7,8,8,8,9,9,9,NA),c(7,7,7,8,8,8,9,9,9,NA))
  pdf(paste("../../Plots/InvasionModel/SuccessPlots.pdf"),height=18.5,width=22)
  grid.arrange(a,d,g,b,e,h,c,f,i,legend,layout_matrix = lay,top = "")
  dev.off()
  return()
}

MeasurePlots <- function(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Home_popsize){
  Popsize_Generations = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Propagule_pressure == Propagule_pressure & finaldf$Home_popsize == Home_popsize),]
  Popsize_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_generations == Invador_generations & finaldf$Home_popsize == Home_popsize),]
  Generations_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Invador_popsize == Invador_popsize & finaldf$Home_popsize == Home_popsize),]
  PropaguleFitness_PropagulePressure = finaldf[which(finaldf$Mutation_rate == Mutation_rate &finaldf$Selection_strength == Selection_strength & finaldf$Home_popsize),]
  p1 <-  ggplot(Popsize_Generations,aes(x=as.factor(Invador_generations), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab(expression( "Invader Generations g"[r]) ) +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Sum') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p2 <-  ggplot(Popsize_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Sum') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p3 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Invasion_success)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Generations g"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) + guides(fill=FALSE) +
    labs(fill='Sum') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p4 <-  ggplot(Popsize_Generations,aes(x=as.factor(Invador_generations), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Proportion)) +
    xlab(expression( "Invader Generations g"[r]) ) +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p5 <-  ggplot(Popsize_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_popsize))) +  geom_tile(aes(fill=Proportion)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Community Size J"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  p6 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Proportion)) +
    xlab("Propagule Pressure p") +ylab(expression( "Invader Generations g"[r]) ) + theme_classic() + theme(panel.border= element_blank()) + 
    theme(text = element_text(size=18),plot.background= element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) + guides(fill=FALSE) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm"))+theme(axis.title=element_text(size=24))
  #legends
  p7 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Invasion_success)) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 500,limits=c(0,1000)) +
    labs(fill='Sum') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm")) +
    theme(axis.title=element_text(size=18),legend.title=element_text(size=18))
  p8 <-  ggplot(Generations_PropagulePressure,aes(x=as.factor(Propagule_pressure), y=as.factor(Invador_generations))) +  geom_tile(aes(fill=Proportion)) +
    scale_fill_gradient2(low = "yellow", mid = "red", high = "blue",midpoint = 0.5, limits = c(0,1)) +
    labs(fill='Proportion') + theme(legend.text = element_text(size=18),legend.key.size = unit(1, "cm")) +
    theme(axis.title=element_text(size=18),legend.title=element_text(size=18))
  a <-  arrangeGrob(p1, top = textGrob("a)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  b <-  arrangeGrob(p2, top = textGrob("b)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  c <-  arrangeGrob(p3, top = textGrob("c)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  d <-  arrangeGrob(p4, top = textGrob("d)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  e <-  arrangeGrob(p5, top = textGrob("e)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  f <-  arrangeGrob(p6, top = textGrob("f)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  g <- g_legend(p7)
  h <- g_legend(p8)
  lay = rbind(c(1,1,1,2,2,2,3,3,3,7),c(4,4,4,5,5,5,6,6,6,8))
  pdf(paste("../../Plots/InvasionModel/MeasurePlots.pdf"),height=14,width=24)
  grid.arrange(a,b,c,d,e,f,g,h,layout_matrix = lay,top = "")
  dev.off()
  return()
}
multiplots <- function(Home_popsize,Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure){
  invspop <- Invasionvspopsize(Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure)
  invsgen <- Invasionvsgenerations(Invador_popsize,Mutation_rate,Selection_strength,Propagule_pressure)
  invsprop<- Invasionvspropagulepressure(Invador_generations,Invador_popsize,Mutation_rate,Selection_strength)
#   invsfit <- Invasionvsrelativefitness(Mutation_rate,Selection_strength,Home_popsize)
  invsrich <- Invasionvspropagulerichness(Mutation_rate,Selection_strength,Home_popsize,Propagule_pressure)
  invvsint <- Internalplots(Mutation_rate,Selection_strength,Invador_generations,
                            Invador_popsize,Home_popsize,Propagule_pressure)
  a <- invspop[[2]]
#   b <- invsgen[[2]]
  b <- invsprop[[2]]
  
  c <- invvsint[[2]]
  d <- invvsint[[1]]
  
  f <- invsrich[[1]]
  e <- invvsint[[3]]
  
  
  a <-  arrangeGrob(a, top = textGrob("a)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  b <-  arrangeGrob(b, top = textGrob("b)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  c <-  arrangeGrob(c, top = textGrob("c)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  d <-  arrangeGrob(d, top = textGrob("d)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  e <-  arrangeGrob(e, top = textGrob("e)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  f <-  arrangeGrob(f, top = textGrob("f)", x = unit(0, "npc")
                                      , y   = unit(1, "npc"), just=c("left","top"),
                                      gp=gpar(col="black", fontsize=18)))
  lay = rbind(c(1,2),c(3,4),c(5,6))
  pdf(paste("../../Plots/InvasionModel/Mixplots2.pdf"),height=15,width=15)
  grid.arrange(a,b,c,d,e,f,layout_matrix = lay,top = "")
  dev.off()
#   
#   a <- Invasionvspopsize(Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure)[[1]]
#   b <- Invasionvspropagulepressure(Invador_generations,Invador_popsize,Mutation_rate,Selection_strength)[[1]]
#   c <- Invasionvsrelativefitness(Mutation_rate,Selection_strength,Home_popsize)[[1]]
#   d <- Invasionvspropagulerichness(Mutation_rate,Selection_strength,Home_popsize,Invador_generations)[[2]]
#   pdf(paste("../../Plots/InvasionModel/MixplotssProportion2.pdf"),height=6,width=10)
#   multiplot(a,c,b,d,cols=2)
#   dev.off()
  return()
}
# paramatersdf <- suppressWarnings(loadparamaters(400,100))
# proportiondf <- proportion_updator(paramatersdf)
# successdf <- success_updator(paramatersdf)
# # gradientdf <-gradient_updator(paramatersdf)
# finaldf <- cbind(paramatersdf,data.frame(Proportion = proportiondf$Mean, Invasion_success = successdf$Mean))
load("summaryInvasionData.Rdata")

Selection_strength = 0.01
Mutation_rate = 1e-04
Home_popsize = 25000
Invador_popsize = 1000
Invador_generations = 1050
Propagule_pressure = 21
invfreqvstime(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Home_popsize)
Invador_popsize = 10500
MeasurePlots(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure,Home_popsize)
# SuccessPlots(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure)
# ProportionPlots(Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure)
multiplots(Home_popsize,Invador_popsize,Invador_generations,Mutation_rate,Selection_strength,Propagule_pressure)
