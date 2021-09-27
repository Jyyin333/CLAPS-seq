#!/usr/bin/env Rscript

# import request libraries 
library(dplyr)
library(ggplot2)


toplot <- function(df, labs, outfig, lg=F){
  # calculate correlation coefficient
  rsquare <- cor(df$V1, df$V2, method='spearman')
  r2 <- sprintf("italic(R^2) == %.3f", rsquare)
  
  # skipzeros
  df <- df[(df$V1!=0 & df$V2!=0), ]
  #cond1 <- df$V1 <= (2.5*quantile(df$V1)[4] - 1.5*quantile(df$V1)[2])
  #cond2 <- df$V2 <= (2.5*quantile(df$V2)[4] - 1.5*quantile(df$V2)[2])
  cond1 <- df$V1 <= 60
  cond2 <- df$V2 <= 60
  cond_df <- data.frame(cd1=cond1, cd2=cond2)
  filter_cond <- cond_df$cd1==T & cond_df$cd2==T
  df <- df[filter_cond, ]
  
  # filter duplicate rows
  df <- distinct(df)
  
  if (lg==T){
    df <- log2(df)
  }

  #xlimit <- max(df$V1)*1.5
  #ylimit <- max(df$V2)*1.5
  xlimit <- 80
  ylimit <- 80
  
  x1 <- as.numeric(quantile(c(0,xlimit),probs = 0.1))
  y1 <- as.numeric(quantile(c(0,ylimit),probs = 0.85))
  p <- ggplot(df, aes(V1, V2)) +
    geom_point(shape=21, size=1, color='black',fill='black') +
    scale_x_continuous(limits=c(0,xlimit), breaks=c(0,20,40,60,80), labels=c(0,20,40,60,80)) +
    scale_y_continuous(limits=c(0,ylimit), breaks=c(0,20,40,60,80), labels=c(0,20,40,60,80)) +
    labs(x=labs[1], y=labs[2], title='') + theme_classic() + # + geom_smooth(method = glm)
    theme(axis.title.x=element_text(size=15,face="bold",vjust=0.5)) +
    theme(axis.title.y=element_text(size=15,face="bold",vjust=0.5)) +
    geom_text(data=df,mapping=aes(x=x1, y=y1, label=r2),
              parse = TRUE,inherit.aes = FALSE,size=4)
  
  ggsave(p,filename=outfig,width=8,height=6,dpi=600)

}



df <- read.table('AllBam.bins.rawCounts.mtx',header=F, sep='\t')
header <- c('chr','start','end','wt1','wt2','wt3','vt1','vt2','vi1','vi2')
colnames(df) <- header

lib.name <- c('wt1','wt2','wt3','vt1','vt2','vi1','vi2')
# library sizes for each sample in order of the header of rawCountsNatrix
lib.sizes <- readLines('AllBam.lib.sizes')

# convert read counts to RPM 
new <- data.frame(wt1=df$wt1/lib.sizes[1]*1e6,
                  wt2=df$wt2/lib.sizes[2]*1e6,
                  wt3=df$wt3/lib.sizes[3]*1e6,
                  vt1=df$vt1/lib.sizes[4]*1e6,
                  vt2=df$vt2/lib.sizes[5]*1e6,
                  vi1=df$vi1/lib.sizes[6]*1e6,
                  vi2=df$vi2/lib.sizes[7]*1e6)

# Vivo
vivo_df <- new[,c(6,7)]
colnames(vivo_df) <- c('V1','V2')
toplot(vivo_df, c('RPM vivo rep1','RPM vivo rep2'), 'vivo_reps.tiff')

# Vitro
vitro_df <- new[,c(4,5)]
colnames(vitro_df) <- c('V1','V2')
toplot(vitro_df, c('RPM vitro rep1','RPM vitro rep2'), 'vitro_reps.tiff')


# WT
wt_df1 <- new[,c(1,2)]
wt_df2 <- new[,c(2,3)]
wt_df3 <- new[,c(1,3)]
colnames(wt_df1) <- c('V1','V2')
colnames(wt_df2) <- c('V1','V2')
colnames(wt_df3) <- c('V1','V2')
toplot(wt_df1, c('RPM WT rep1','RPM WT rep2'), 'WT12_reps.tiff')
toplot(wt_df2, c('RPM WT rep2','RPM WT rep3'), 'WT23_reps.tiff')
toplot(wt_df3, c('RPM WT rep1','RPM WT rep3'), 'WT13_reps.tiff')
