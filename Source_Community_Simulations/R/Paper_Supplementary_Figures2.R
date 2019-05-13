library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library('ggpubr')
rm(list=ls())

Supplementary_Figure_6 <- function(){
  var_df = fread('Variance.csv')
  rich_df = fread('Richness.csv')
  mean_df = fread('Fitness.csv')
  mean_df = mean_df[mean_df$Generation %% 100 ==0,]
  add_df = data.frame('Generation' = mean_df$Generation,'Model' = mean_df$Model)
  add_df$Mean = (mean_df$Mean-c(100,mean_df$Mean[-nrow(mean_df)]))/100
  add_df$Max = (mean_df$Max-c(100,mean_df$Max[-nrow(mean_df)]))/100
  add_df$Min = (mean_df$Min-c(100,mean_df$Min[-nrow(mean_df)]))/100
  hist_df = fread('Distribution.csv')
  p1<- ggplot(hist_df, aes(x=Category, y = Mean, fill= as.factor(Generation))) + 
    geom_bar(stat="identity",width=1, position = position_dodge(width=0.5))  + 
    xlab("Fitnesss Category") + ylab("Abundance") + theme_classic() +
    geom_ribbon(aes(ymin = Min, ymax = Max), alpha = 0.5) +
    scale_fill_discrete(name="Generation") + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) + 
    theme(text = element_text(size=18),legend.text = element_text(size=15))
  p2 <- ggplot(var_df, aes(x = Generation, y = Mean))+ geom_ribbon(aes(ymin = Min, ymax = Max), alpha = 0.5) +
    geom_line(size=1) + ylab("Fitness Variance") + theme_classic()+ 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) + 
    theme(text = element_text(size=18))
  p3 <- ggplot(add_df, aes(x = Generation, y = Mean)) + geom_ribbon(aes(x=Generation,ymin=Min,ymax=Max),fill='grey') +
    geom_line(size=1) + ylab("Adaptation Rate") + theme_classic()+ 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) + 
    theme(text = element_text(size=18))
  p4 <- ggplot(rich_df, aes(x = Generation, y = Mean)) + geom_ribbon(aes(x=Generation,ymin=Min,ymax=Max),fill='grey') +
    geom_line(size=1) + ylab("Genotypic Richness") + theme_classic()+ 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 2), axis.line.y = element_line(color="black", size = 2)) + 
    theme(text = element_text(size=18))
  p1 = arrangeGrob(p1,top=textGrob('A)',x=0.05,gp=gpar(fontsize=15)))
  p2 = arrangeGrob(p2,top=textGrob('B)',x=0.05,gp=gpar(fontsize=15)))
  p3 = arrangeGrob(p3,top=textGrob('C)',x=0.05,gp=gpar(fontsize=15)))
  p4 = arrangeGrob(p4,top=textGrob('D)',x=0.05,gp=gpar(fontsize=15)))
  
  t = grid.arrange(p1,p2,p3,p4,layout_matrix=rbind(c(1,1,1),c(2,3,4)))
  ggsave(file = '../../Plots/Supp6.png',t,width = 16, height = 8)
  
}

Supplementary_Figure_7 <- function(){
  DFE = fread('DFE_Simulations.csv')
  DFE[Model=='Equal']$Model = 'Default'
  mean_ad <- aggregate(adaptation_rate ~ Popsize + Model, data = DFE[Selection_strength ==0.01 & Mutation_rate ==1e-04], FUN = function (x) mean(x))
  min_ad <- aggregate(adaptation_rate ~ Popsize + Model, data = DFE[Selection_strength ==0.01 & Mutation_rate ==1e-04], FUN = function (x) min(x))
  max_ad <- aggregate(adaptation_rate ~ Popsize + Model, data = DFE[Selection_strength ==0.01 & Mutation_rate ==1e-04], FUN = function (x) max(x))
  add_df = data.frame('Mean' = mean_ad$adaptation_rate,'Min' = min_ad$adaptation_rate,'Max' = max_ad$adaptation_rate,'Popsize' = mean_ad$Popsize,'Model'=mean_ad$Model)
  p1 <- ggplot(add_df, aes(x = Popsize, y = Mean,color = Model)) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Source Community Size (R)")+ylab("Adaptation Rate") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) +  
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='Source Community Size (R)') + geom_ribbon(aes(ymin = Min, ymax = Max,fill=Model), alpha = 0.5) + 
    guides(color = FALSE) + scale_fill_discrete(name = "") + scale_x_log10() +
    theme(text = element_text(size=18)) + labs(fill='')
  mean_var <- aggregate(variance ~ Popsize + Model, data = DFE[Selection_strength ==0.01 & Mutation_rate ==1e-04], FUN = function (x) mean(x))
  min_var <- aggregate(variance ~ Popsize + Model, data = DFE[Selection_strength ==0.01 & Mutation_rate ==1e-04], FUN = function (x) min(x))
  max_var <- aggregate(variance ~ Popsize + Model, data = DFE[Selection_strength ==0.01 & Mutation_rate ==1e-04], FUN = function (x) max(x))
  var_df = data.frame('Mean' = mean_var$variance,'Min' = min_var$variance,'Max' = max_var$variance,'Popsize' = mean_var$Popsize,'Model'=mean_var$Model)
  p2 <- ggplot(var_df, aes(x = Popsize, y = Mean,color = Model)) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Source Community Size (R)")+ylab("Fitness Variance") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='Source Community Size (R)') + geom_ribbon(aes(ymin = Min, ymax = Max,fill=Model), alpha = 0.5) + 
    guides(color = FALSE) + scale_fill_discrete(name = "") + scale_x_log10() + scale_y_log10() +
    theme(text = element_text(size=18)) + labs(fill='')
  mean_rich <- aggregate(richness ~ Popsize + Model, data = DFE[Selection_strength ==0.01 & Mutation_rate ==1e-04], FUN = function (x) mean(x))
  min_rich <- aggregate(richness ~ Popsize + Model, data = DFE[Selection_strength ==0.01 & Mutation_rate ==1e-04], FUN = function (x) min(x))
  max_rich <- aggregate(richness ~ Popsize + Model, data = DFE[Selection_strength ==0.01 & Mutation_rate ==1e-04], FUN = function (x) max(x))
  rich_df = data.frame('Mean' = mean_rich$richness,'Min' = min_rich$richness,'Max' = max_rich$richness,'Popsize' = mean_rich$Popsize,'Model'=mean_rich$Model)
  p3 <- ggplot(rich_df, aes(x = Popsize, y = Mean,color = Model)) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Source Community Size (R)")+ylab("Genotypic Richness") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='Source Community Size (R)') + geom_ribbon(aes(ymin = Min, ymax = Max,fill=Model), alpha = 0.5) + 
    guides(color = FALSE) + scale_fill_discrete(name = "") + scale_x_log10() + scale_y_log10() +
    theme(text = element_text(size=18)) + labs(fill='')
  
  mean_ad2 <- aggregate(adaptation_rate ~ Selection_strength + Model, data = DFE[Popsize ==10000 & Mutation_rate ==1e-04], FUN = function (x) mean(x))
  min_ad2 <- aggregate(adaptation_rate ~ Selection_strength + Model, data = DFE[Popsize ==10000  & Mutation_rate ==1e-04], FUN = function (x) min(x))
  max_ad2 <- aggregate(adaptation_rate ~ Selection_strength + Model, data = DFE[Popsize ==10000  & Mutation_rate ==1e-04], FUN = function (x) max(x))
  add_df2 = data.frame('Mean' = mean_ad2$adaptation_rate,'Min' = min_ad2$adaptation_rate,'Max' = max_ad2$adaptation_rate,'Selection_strength' = mean_ad2$Selection_strength,'Model'=mean_ad2$Model)
  p4 <- ggplot(add_df2, aes(x = Selection_strength, y = Mean,color = Model)) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Selection Strength (S)")+ylab("Adaptation Rate") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='Selection Strength (S)') + geom_ribbon(aes(ymin = Min, ymax = Max,fill=Model), alpha = 0.5) + 
    guides(color = FALSE) + scale_fill_discrete(name = "") + scale_x_log10() +
    theme(text = element_text(size=18)) + labs(fill='')
  mean_var2 <- aggregate(variance ~ Selection_strength + Model, data = DFE[Popsize ==10000  & Mutation_rate ==1e-04], FUN = function (x) mean(x))
  min_var2 <- aggregate(variance ~ Selection_strength + Model, data = DFE[Popsize ==10000 & Mutation_rate ==1e-04], FUN = function (x) min(x))
  max_var2 <- aggregate(variance ~ Selection_strength + Model, data = DFE[Popsize ==10000  & Mutation_rate ==1e-04], FUN = function (x) max(x))
  var_df2 = data.frame('Mean' = mean_var2$variance,'Min' = min_var2$variance,'Max' = max_var2$variance,'Selection_strength' = mean_var2$Selection_strength,'Model'=mean_var2$Model)
  p5 <- ggplot(var_df2, aes(x = Selection_strength, y = Mean,color = Model)) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Selection Strength (S)")+ylab("Fitness Variance") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='Selection Strength (S)') + geom_ribbon(aes(ymin = Min, ymax = Max,fill=Model), alpha = 0.5) + 
    guides(color = FALSE) + scale_fill_discrete(name = "") + scale_x_log10() + scale_y_log10() +
    theme(text = element_text(size=18)) + labs(fill='')
  mean_rich2 <- aggregate(richness ~ Selection_strength + Model, data = DFE[Popsize ==10000  & Mutation_rate ==1e-04], FUN = function (x) mean(x))
  min_rich2 <- aggregate(richness ~ Selection_strength + Model, data = DFE[Popsize ==10000  & Mutation_rate ==1e-04], FUN = function (x) min(x))
  max_rich2 <- aggregate(richness ~ Selection_strength + Model, data = DFE[Popsize ==10000& Mutation_rate ==1e-04], FUN = function (x) max(x))
  rich_df2 = data.frame('Mean' = mean_rich2$richness,'Min' = min_rich2$richness,'Max' = max_rich2$richness,'Selection_strength' = mean_rich2$Selection_strength,'Model'=mean_rich2$Model)
  p6 <- ggplot(rich_df2, aes(x = Selection_strength, y = Mean,color = Model)) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Selection Strength (S)")+ylab("Genotypic Richness") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='Selection Strength (S)') + geom_ribbon(aes(ymin = Min, ymax = Max,fill=Model), alpha = 0.5) + 
    guides(color = FALSE) + scale_fill_discrete(name = "") + scale_x_log10() + scale_y_log10() +
    theme(text = element_text(size=18)) + labs(fill='')
  

  mean_ad3 <- aggregate(adaptation_rate ~ Mutation_rate + Model, data = DFE[Selection_strength ==0.01 & Popsize ==10000], FUN = function (x) mean(x))
  min_ad3 <- aggregate(adaptation_rate ~ Mutation_rate + Model, data = DFE[Selection_strength ==0.01 & Popsize ==10000], FUN = function (x) min(x))
  max_ad3 <- aggregate(adaptation_rate ~ Mutation_rate + Model, data = DFE[Selection_strength ==0.01 & Popsize ==10000], FUN = function (x) max(x))
  add_df3 = data.frame('Mean' = mean_ad3$adaptation_rate,'Min' = min_ad3$adaptation_rate,'Max' = max_ad3$adaptation_rate,'Mutation_rate' = mean_ad3$Mutation_rate,'Model'=mean_ad3$Model)
  p7 <- ggplot(add_df3, aes(x = Mutation_rate, y = Mean,color = Model)) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Mutation Rate " * mu ~ "")+ylab("Adaptation Rate") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    geom_ribbon(aes(ymin = Min, ymax = Max,fill=Model), alpha = 0.5) + 
    guides(color = FALSE) + scale_fill_discrete(name = "") + scale_x_log10() +
    theme(text = element_text(size=18)) + labs(fill='')
  mean_var3 <- aggregate(variance ~ Mutation_rate + Model, data = DFE[Selection_strength ==0.01 & Popsize ==10000], FUN = function (x) mean(x))
  min_var3 <- aggregate(variance ~ Mutation_rate + Model, data = DFE[Selection_strength ==0.01 & Popsize ==10000], FUN = function (x) min(x))
  max_var3 <- aggregate(variance ~ Mutation_rate + Model, data = DFE[Selection_strength ==0.01 & Popsize ==10000], FUN = function (x) max(x))
  var_df3 = data.frame('Mean' = mean_var3$variance,'Min' = min_var3$variance,'Max' = max_var3$variance,'Mutation_rate' = mean_var3$Mutation_rate,'Model'=mean_var3$Model)
  p8 <- ggplot(var_df3, aes(x = Mutation_rate, y = Mean,color = Model)) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Mutation Rate " * mu ~ "")+ylab("Fitness Variance") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='"Mutation Rate " * mu ~ ""') + geom_ribbon(aes(ymin = Min, ymax = Max,fill=Model), alpha = 0.5) + 
    guides(color = FALSE) + scale_fill_discrete(name = "") + scale_x_log10() + scale_y_log10() +
    theme(text = element_text(size=18)) + labs(fill='')
  mean_rich3 <- aggregate(richness ~ Mutation_rate + Model, data = DFE[Selection_strength ==0.01 & Popsize ==10000], FUN = function (x) mean(x))
  min_rich3 <- aggregate(richness ~ Mutation_rate + Model, data = DFE[Selection_strength ==0.01 & Popsize ==10000], FUN = function (x) min(x))
  max_rich3 <- aggregate(richness ~ Mutation_rate + Model, data = DFE[Selection_strength ==0.01 & Popsize ==10000], FUN = function (x) max(x))
  rich_df3 = data.frame('Mean' = mean_rich3$richness,'Min' = min_rich3$richness,'Max' = max_rich3$richness,'Mutation_rate' = mean_rich3$Mutation_rate,'Model'=mean_rich3$Model)
  p9 <- ggplot(rich_df3, aes(x = Mutation_rate, y = Mean,color = Model)) +  geom_point(size=2) + geom_line(size=1) + 
    xlab("Mutation Rate " * mu ~ "")+ylab("Genotypic Richness") + theme_classic() + 
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
    theme(panel.border= element_blank()) + 
    theme(axis.line.x = element_line(color="black", size = 1), axis.line.y = element_line(color="black", size = 1))+
    labs(color='"Mutation Rate " * mu ~ ""') + geom_ribbon(aes(ymin = Min, ymax = Max,fill=Model), alpha = 0.5) + 
    guides(color = FALSE) + scale_fill_discrete(name = "") + scale_x_log10() + scale_y_log10() +
    theme(text = element_text(size=18)) + labs(fill='')
  
   ggsave(file = '../../Plots/Supp7.png',ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3,nrow=3,common.legend=TRUE,labels='auto',legend='bottom'),height=10,width=15)
}

Supplementary_Figure_7()
Supplementary_Figure_6()