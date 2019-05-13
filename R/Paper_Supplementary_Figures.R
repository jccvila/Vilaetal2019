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


plotdata= fread('Summary_Data.csv')
s = unique(plotdata$selection_strength)
u = unique(plotdata$mutation_rate)
p = unique(plotdata$propagule_pressure)
i <- unique(plotdata$resident_popsize)
r <- unique(plotdata$invader_popsize)
g <- unique(plotdata$invader_generations)

params = expand.grid(s,u,p,i,r,g)

SuppFig1 <- function(){
  
  row1 = which(plotdata$resident_popsize==i[3] & plotdata$mutation_rate !=0)

  plotdf1 <- plotdata[row1]

  sum_df  =plotdf1[,list(initial = mean(Resident_Richness_Start),mean=mean(Invader_Richness_End),sd=sd(Invader_Richness_End)),by=list(invader_popsize,invader_generations,mutation_rate,selection_strength)]
  sum_df$mutation_rate = paste('u = ',sum_df$mutation_rate,sep='')
  sum_df$selection_strength = paste('s = ',sum_df$selection_strength,sep='')

  p1 <- ggplot(sum_df, aes(x = invader_generations, y =mean,color = as.factor(invader_popsize))) +  geom_line(size=2) + geom_point(size=4) +
    xlab("Number of Generations") +ylab("Number of Genotypes") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + scale_y_log10() + scale_x_log10() +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(color='Community Size') +  
    theme(text = element_text(size=18),legend.text = element_text(size=15),legend.key.size = unit(0.5, "cm"))  +
    facet_grid(mutation_rate~selection_strength,scales='free_y')
  ggsave(file = '../Plots/Supp1.png',p1,width = 16, height = 8)
}  

SuppFig2 <- function(){
    
  row1 = which(plotdata$selection_strength==s[3] & 
               plotdata$mutation_rate == u[3] &
               plotdata$resident_popsize == i[3])
  plotdf1 <- plotdata[row1]
  options("scipen"=1000) 
  tdf = plotdf1[which(plotdf1$invader_generations == g[2] & plotdf1$invader_popsize %in% r[10]),]
  mod = lm(tdf$Frequency40~tdf$propagule_pressure)
  slope = signif(mod$coefficients[2],2)
  print(summary(mod))
  pvalue = 'p < 1e-6'
  p1 <-  ggplot(tdf,aes(x=as.numeric(propagule_pressure), y=Frequency40)) +
    geom_jitter(size= 4,width=0.02,height=0,colour='grey') + 
    geom_smooth(method='lm',colour='red') +
    theme_classic() + 
    labs(x = 'Propagule Pressure (P)') +  
    theme(legend.position = c(0.7,0.6)) +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    labs(y = 'Invasion Success',group='',shape='',col='') +
    annotate(geom='text',x = 20,y=0.25,label=paste('\t\t\t\t Slope =', slope,', ' , pvalue),size=5)
  mod2 = lm(tdf$Frequency40~tdf$Propagule_Richness)
  slope2 = signif(mod2$coefficients[2],2)
  print(summary(mod2))
  
  pvalue2 = 'p < 1e-6'
  p2 <-  ggplot(tdf[propagule_pressure==200,],aes(x=as.numeric(Propagule_Richness), y=Frequency40)) +
    geom_jitter(size= 4,width=0.02,height=0,colour='grey') + 
    geom_smooth(method=lm,colour='#F8766D',fill='#F8766D') + theme_classic() + #scale_x_log10(limits=c(10,20)) +
    labs(x = 'Number of Invading Genotypes') +  
    theme(legend.position = c(0.7,0.6)) +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    labs(y = 'Invasion Success',group='',shape='',col='') +
    annotate(geom='text',x = 12,y=0.25,label=paste('\t\t\t\t Slope=', slope2,', ' , pvalue2),size=5)
  p1 = arrangeGrob(p1,top=textGrob('A)',x=0.05,gp=gpar(fontsize=15)))
  p2 = arrangeGrob(p2,top=textGrob('B)',x=0.05,gp=gpar(fontsize=15)))
  ggsave(file = '../Plots/Supp2.png',grid.arrange(arrangeGrob(p1),arrangeGrob(p2)
                                                    ,layout_matrix = rbind(c(1,2))),width = 8, height = 4)
}  


SuppFig3 <- function(){
  row1 = which(plotdata$selection_strength==s[3] & 
                 plotdata$mutation_rate == u[3] &
                 plotdata$resident_popsize == i[3])
  plotdf1 <- plotdata[row1]
  options("scipen"=1000) 
  tdf = plotdf1[which(plotdf1$invader_generations == g[2] & plotdf1$invader_popsize %in% r[10] & plotdf1$propagule_pressure==p[10]),]
  mod1 = lm(tdf$Relative_Propagule_Fitness~tdf$Propagule_Richness)
  mod2 = lm(tdf$Propagule_Fitness_Variance~tdf$Propagule_Richness)

  slope1 = signif(mod1$coefficients[2],2)
  slope2 = signif(mod2$coefficients[2],2)
  # i do this manually for the sake of formatting.
  pvalue1 = 'p = 0.0357'
  pvalue2 = 'p < 1e-06'
  p1 <-  ggplot(tdf,aes(x=Propagule_Richness, y=Relative_Propagule_Fitness)) +
    geom_jitter(size= 4,width=0.02,height=0,colour='grey') + 
    geom_smooth(method='lm',colour='red') +
    theme_classic() + 
    labs(x = 'Number of Invading Genotypes',y = 'Relative Invader Fitness') +  
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    theme(text = element_text(size=18)) +
    annotate(geom='text',x = 12.5,y=1.002,label=paste('\t\t\t\t Slope =', slope1,', ' , pvalue1),size=5)
  
  p2 <-  ggplot(tdf,aes(x=Propagule_Richness, y=Propagule_Fitness_Variance)) +
    geom_jitter(size= 4,width=0.02,height=0,colour='grey') + 
    geom_smooth(method='lm',colour='red') +
    theme_classic() + 
    labs(x = 'Number of Invading Genotypes',y = 'Invader Fitness Variance') +  
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) +
    theme(text = element_text(size=18)) +
    annotate(geom='text',x = 12.5,y=16,label=paste('\t\t\t\t Slope =', slope2,', ' , pvalue2),size=5)

  sum_df = plotdf1[,list(pval = mean(PR_PMpval),
                       mean=mean(PR_PM),sd=sd(PR_PM),
                         Relative_Propagule_Fitness = mean(Relative_Propagule_Fitness)),
                        by=list(invader_popsize,invader_generations,propagule_pressure)]
  sum_df[is.na(sum_df$pval)]$pval=1
  sum_df[is.na(sum_df$mean)]$mean=0
  sum_df = sum_df[propagule_pressure==200,]
  
  sum_df2 = plotdf1[,list(pval = mean(PR_PVpval),
                         mean=mean(PR_PV),sd=sd(PR_PV),
                         Relative_Propagule_Fitness = mean(Relative_Propagule_Fitness)),
                   by=list(invader_popsize,invader_generations,propagule_pressure)]
  sum_df2 = sum_df2[propagule_pressure==200,]  
  p3 <- ggplot(sum_df, aes(mean)) +  geom_histogram(bins=10) + 
    xlab("Effect on Mean Fitness") + ylab('Number of Paramater Combinations')  + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + geom_vline(xintercept = 0,linetype=2,color='Red') +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    theme(text = element_text(size=18)) + geom_hline(yintercept=0,col='red',linetype=2)
  p4 <- ggplot(sum_df2, aes(mean)) +  geom_histogram(bins=10)  + ylab('Number of Paramater Combinations')  + 
    xlab("Effect on Fitness Variance") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + geom_vline(xintercept = 0,linetype=2,color='Red') +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    theme(text = element_text(size=18)) + geom_hline(yintercept=0,col='red',linetype=2)
  
  sum_df3= plotdf1[,list(pval = mean(PR3pval),
                         mean=mean(PR3slope),
                         Relative_Propagule_Fitness = mean(Relative_Propagule_Fitness)),
                   by=list(invader_popsize,invader_generations,propagule_pressure)]
  sum_df3 = sum_df3[propagule_pressure==200,]

  sum_df4 = plotdf1[,list(pval = mean(PR3pval),
                          mean=mean(PR2slope),
                          Relative_Propagule_Fitness = mean(Relative_Propagule_Fitness)),
                    by=list(invader_popsize,invader_generations,propagule_pressure)]
  sum_df4 = sum_df4[propagule_pressure==200,]  

  
  p5 <- ggplot(sum_df3, aes(x = Relative_Propagule_Fitness,y = mean)) +  geom_point(size=2) + 
    xlab("Relative Invader Fitness") +ylab("Effect of Invading Genotype Number") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + geom_vline(xintercept = 1,linetype=2,color='Red') +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    theme(text = element_text(size=18))  + ggtitle('Success/Mean Invader Fitness')
  
  p6 <- ggplot(sum_df4, aes(x = Relative_Propagule_Fitness,y = mean)) +  geom_point(size=2) + 
    xlab("Relative Invader Fitness") +ylab("Effect of Invading Genotype Number") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + geom_vline(xintercept = 1,linetype=2,color='Red') +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    theme(text = element_text(size=18))  + ggtitle('Success/Invader Fitness Variance')
  
  p1 = arrangeGrob(p1,top=textGrob('A)',x=0.05,gp=gpar(fontsize=15)))
  p2 = arrangeGrob(p2,top=textGrob('B)',x=0.05,gp=gpar(fontsize=15)))
  p3 = arrangeGrob(p3,top=textGrob('C)',x=0.05,gp=gpar(fontsize=15)))
  p4 = arrangeGrob(p4,top=textGrob('D)',x=0.05,gp=gpar(fontsize=15)))
  p5 = arrangeGrob(p5,top=textGrob('E)',x=0.05,gp=gpar(fontsize=15)))
  p6 = arrangeGrob(p6,top=textGrob('F)',x=0.05,gp=gpar(fontsize=15))) 
  ggsave(file = '../Plots/Supp3.png',
         grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2), width = 12, height = 16)
}

SuppFig4 <- function(){
  row1 = which(plotdata$selection_strength==s[3] & 
                 plotdata$mutation_rate == u[3] &
                 plotdata$propagule_pressure==p[10] & 
                 plotdata$resident_popsize <=i[3])
  
  plotdf1 <- plotdata[row1]
  sum_df = plotdf1[,list(pval = mean(IRpval),
                          mean=mean(IR),sd=sd(IR),
                          Relative_Propagule_Fitness = mean(Invader_Mean_End/Resident_Mean_End)),
                    by=list(invader_popsize,invader_generations,
                            resident_popsize,propagule_pressure)]
  #NA have arisen if there is no variation in invasion outcome or invader reichness. By defualt we assign this to have no effect so PR =0 and pval = 1
  sum_df[is.na(sum_df$pval)]$pval=1
  sum_df[is.na(sum_df$mean)]$mean=0
  #correct for multiple comparisons when assigning significance (bonferioni)
  sum_df = sum_df[propagule_pressure==200,]
  sum_df$resident_popsize = paste('R = ',sum_df$resident_popsize,sep='')
  p1 <- ggplot(sum_df, aes(x = Relative_Propagule_Fitness,y = mean)) +  geom_point(size=2) + 
    xlab("Relative Invader Fitness") +ylab("Effect of Resident Community Genotypic Richness") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + geom_vline(xintercept = 1,linetype=2,color='Red') +
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))+
    theme(text = element_text(size=18)) +  #+ labs(col='p-value') +
    # scale_color_gradient(low='red',high='black',trans='log10') + 
    facet_wrap(~resident_popsize, scales='free')
  ggsave(file = '../Plots/Supp4.png',p1,width=15,height=5)
}  


SuppFig5<- function(){
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
  p1 <-ggplot(data=plotdf1,aes(x=as.factor(resident_popsize),y=Frequency40*resident_popsize)) + 
    geom_boxplot() + theme_classic() + 
    theme(text = element_text(size=18)) +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + #scale_y_continuous(limits=c(0,1)) +
    labs(x= 'Resident Community Size',y = 'Invader Abundance') + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))
  p2 <-ggplot(data=plotdf2,aes(x=as.factor(resident_popsize),y=Frequency40*resident_popsize)) + 
    geom_boxplot() + theme_classic() + 
    theme(text = element_text(size=18)) +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())   + #scale_y_continuous(limits=c(0,1)) +
    labs(x= 'Resident Community Size',y = 'Invader Abundance') + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))
  
  
  reg_df = data.frame()
  for(i in unique(plotdf3$resident_popsize)){
    mod = lm(Frequency40*resident_popsize~Resident_Richness_End,data=plotdf3[resident_popsize==i,])
    reg_df = rbind(reg_df,c(signif(mod$coefficients[2],2),signif(summary(mod)$coefficients[2,4],2),i))
  }
  mod = lm(Frequency40*resident_popsize~Resident_Richness_End,data=plotdf3)
  print(summary(mod))
  colnames(reg_df) = c('Slope','P-Value','resident_popsize')
  p3 <-ggplot(data=plotdf3,
              aes(x=Resident_Richness_End,y=Frequency40*resident_popsize,col=as.factor(resident_popsize))) + 
    geom_point() + theme_classic() + 
    theme(text = element_text(size=18)) +
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank())  + 
    geom_smooth(method='lm',se=FALSE) +
    geom_segment(x = 2, xend = 18,y = mod$coefficients[1]  + 2*mod$coefficients[2],
                 yend = mod$coefficients[1]  + 18*mod$coefficients[2],
                 col='Black')+
    annotate(x = 20, y =2800, 
             label=paste('slope =', reg_df$Slope[1],', p = ' , reg_df$`P-Value`[1]),
             geom='text',size=5,colour='#F8766D') +
    annotate(x = 20, y =2700, 
             label=paste('slope =', reg_df$Slope[2],', p = ' , reg_df$`P-Value`[2]),
             geom='text',size=5,colour='#00BA38') +
    annotate(x = 20, y =2600, 
             label=paste('slope =', reg_df$Slope[3],', p = ' , reg_df$`P-Value`[3]),
             geom='text',size=5,colour='#619CFF') +
    annotate(x = 20, y =2500, 
             label=paste('slope =', signif(mod$coefficients[2],2),', p = 0.00001'),
             geom='text',size=5,colour='Black') +
    labs(x= 'Resident Community Genotypic Richness',y = 'Invader Abundance',col='R') + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2))
  
  
  options("scipen"=10)
  p1 = arrangeGrob(p1,top=textGrob('A)',x=0.05,gp=gpar(fontsize=15)))
  p2 = arrangeGrob(p2,top=textGrob('B)',x=0.05,gp=gpar(fontsize=15)))
  p3 = arrangeGrob(p3,top=textGrob('C)',x=0.05,gp=gpar(fontsize=15)))
  ggsave(file = '../Plots/Supp5.png',grid.arrange(arrangeGrob(p1),arrangeGrob(p2),arrangeGrob(p3),
                                                    layout_matrix = rbind(c(1,3,3),c(2,3,3))),
         width = 14, height = 8)
  
  return()
}

SuppFig1()
SuppFig2()
SuppFig3()
SuppFig4()
SuppFig5()