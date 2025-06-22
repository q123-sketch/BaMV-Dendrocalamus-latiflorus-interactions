############################FigueS24, S25 and S26 showing gene sets changes in m6A, PAL, 3'UTR length, expression levels, TPM and protein abundance###############################
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(tidyverse)
library(Rmisc)
###load dataframe
bm_ratio_df<- read.delim("F:\\project\\DRS\\result\\02_m6A\\ratio_df\\bm_ratio_df",sep='\t')
ck_ratio_df<- read.delim("F:\\project\\DRS\\result\\02_m6A\\ratio_df\\ck_ratio_df",sep='\t')

abundance <- read.csv("F:\\project\\DRS\\result\\04_expression\\protemics\\test_DEPs\\pro_abundance.csv",sep= ',',check.names = F);colnames(abundance)[1] <- 'gene_name'
tpm <- read.csv("F:\\project\\DRS\\result\\04_expression\\rnaseq\\deseq2\\deseq2_bamv_ck_res_1.5.csv",sep=',',check.names = F)
rownames(tpm) <- tpm$gene_name
#tpm <- log2(tpm[,-1] + 1)
tpm <- tpm[apply(tpm[,8:13],1,sum)!=0,]
classed_geneset <- read.csv("F:\\project\\DRS\\result\\geneset\\gene_list_2.txt",sep='\t',header = T,check.names = F)
ck_median_pal <- read.csv("F:\\project\\DRS\\result\\01_polyA\\process\\ck_polya_featurecount_inte.csv",sep=',')
bm_median_pal <- read.csv("F:\\project\\DRS\\result\\01_polyA\\process\\bamv_polya_featurecount_inte.csv",sep=',')
ck_median_pal <- na.omit(ck_median_pal)
bm_median_pal <- na.omit(bm_median_pal)


dla_pas_bm_ck <- read.csv("F:\\project\\DRS\\result\\03_APA\\02_pas_cluste\\dla_pas_bm_ck.csv",sep='\t')
dla_pas_bm_ck_use <- dla_pas_bm_ck[,c(1,11)]

#######
i='ABA'
for (i in c('SA','JA','ABA','RNAi','APA')){
  print(i)
  classed_geneset_use <- classed_geneset[classed_geneset$Type == i,]
  n=nrow(classed_geneset_use)
  }

classed_geneset_use <- classed_geneset
ratio_plot <- left_join(classed_geneset_use,bm_ratio_df[,c(1,6)],by='gene_name')
colnames(ratio_plot)[8] <- 'bm_ratio'
ratio_plot <- left_join(ratio_plot,ck_ratio_df[,c(1,6)]);colnames(ratio_plot)[9] <- 'ck_ratio'


ratio_plot  <- ratio_plot %>% 
  tidyr::pivot_longer(cols = c('bm_ratio','ck_ratio',), # 变的列，也就是即将合并为一列的列
                      names_to = "sample", # 表示你分组后的 group 列名，即 cols 合并后的列名
                      values_to = "mod_ratio")
ratio_plot['name'] = ratio_plot['use_name']
summarySE(data = ratio_plot,measurevar = "mod_ratio",groupvars = c("sample","source"))


data_long1 <- tpm_SE_dat %>% 
  tidyr::pivot_longer(cols = c('Ck_1_TPM','Ck_2_TPM','Ck_3_TPM','BaMV_1_TPM','BaMV_2_TPM','BaMV_3_TPM'), # 变的列，也就是即将合并为一列的列
                      names_to = "test_group", # 表示你分组后的 group 列名，即 cols 合并后的列名
                      values_to = "test_value") # 表示对应的值的列名
data_long1$source <- factor(rep(c(rep('CK',3),rep('BaMV',3)),n),levels = c('BaMV','CK'))
tpm_plot_dat <- summarySE(data = data_long1,measurevar = "test_value",groupvars = c("name","source"))
ratio_plot$source <- factor(c(rep(c('BaMV','CK'),nrow(ratio_plot)/2)),levels = c('CK','BaMV'))
ratio_plot$func <- factor(ratio_plot$func,levels =c('synthesis','TF','regulator') )
ratio_plot$use_name <- factor(ratio_plot$use_name,levels = unique(ratio_plot$use_name))
#ratio_plot$plot_name <- paste0()
require(plyr)
data_summary <- function(data,varname,group){
  summary_func <- function(x,col){
    c(median <- median(x[[col]],na.rm = TRUE),sd(x[[col]],na.rm = TRUE))}
    
  data_sum<-ddply(data, group, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("median" = varname))
    return(data_sum)
}

df <- data_summary(ratio_plot,varname='mod_ratio',group=c('source','use_name'))
df <- rename(df,c('V1'='median','V2'='sd'))
df$name <- paste(df$source,df$use_name,sep=' ')
ratio_plot$name <- paste(ratio_plot$source,ratio_plot$use_name,sep =' ')
ratio_plot_dat <- left_join(df,classed_geneset_use)
#####pal_plot_dat
classed_geneset_use
pal_plot <- left_join(classed_geneset_use,bm_median_pal[,c(2,6)],by='gene_name')
colnames(pal_plot)[8] <- 'bm_pal'
pal_plot <- cbind(pal_plot,ck_median_pal[,c(2,6)],by='gene_name');colnames(pal_plot)[9] <- 'ck_pal'

pal_plot  <- pal_plot %>% 
  tidyr::pivot_longer(cols = c('bm_pal','ck_pal',), # 变的列，也就是即将合并为一列的列
                      names_to = "sample", # 表示你分组后的 group 列名，即 cols 合并后的列名
                      values_to = "pal_len")

pal_plot$source <- factor(c(rep(c('BaMV','CK'),nrow(pal_plot)/2)),levels = c('CK','BaMV'))
pal_plot$func <- factor(pal_plot$func,levels =unique(pal_plot$func) )
pal_plot$use_name <- factor(pal_plot$use_name,levels = unique(pal_plot$use_name))

write.csv(ratio_plot_dat,file = 'F:\\project\\DRS\\result\\geneset\\geneset_barplot\\replot\\ratio_plot_dat.csv',row.names = F)

 ggplot(data = plt_dat)+
  geom_boxplot(aes(x=use_name,y=,fill=source))+xlab('')
pdf('F:\\project\\DRS\\result\\geneset\\geneset_barplot\\replot\\test_plot.pdf')  
ggplot(data = ratio_plot_dat,aes(x=use_name,y=median,fill=source))+
  geom_errorbar(aes(ymin=median-sd, ymax=median+sd,color=source),color='black', position = position_dodge(3),size=0.5,width = 2.5)+
  geom_point(aes(size=sd,color=source),position = position_dodge(3))+
  scale_color_manual(values=c("#ffb2cc","#8ce3e3"))+
  theme_minimal()+coord_flip()+
  theme(panel.grid.major = element_line(color='gray',size=1),
        panel.grid.minor = element_blank())
  #scale_fill_manual(values=c('black','black'))
  #
  #stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
   #            width = 0.15,position = position_dodge( .9))
dev.off()


#colnames(classed_m6A) <- c('func','class','gene_name')
#c('SA','JA','ABA','RNAi')
for (i in c('SA','JA','ABA','RNAi','APA')){
  print(i)
  classed_m6A <- classed_geneset[classed_geneset$Type == i,]
  n=nrow(classed_m6A)
  
  ####eraser/reader/writer ratio
  
  #writer_eraser_reader_ratio <- full_join(median_ratio_ck[rownames(median_ratio_ck) %in% writer_eraser_reader_id,],
  #median_ratio_bm[rownames(median_ratio_bm) %in% writer_eraser_reader_id,])
  writer_eraser_reader_ratio <- left_join(classed_m6A,median_ratio_ck,by='gene_name')
  writer_eraser_reader_ratio <- left_join(writer_eraser_reader_ratio,median_ratio_bm,by='gene_name')
  #writer_eraser_reader_ratio[is.na(writer_eraser_reader_ratio)]=0
  
  plot_dat <-writer_eraser_reader_ratio %>% 
    tidyr::pivot_longer(cols = c('ck_median_ratio','bm_median_ratio',), # 变的列，也就是即将合并为一列的列
                        names_to = "test_group", # 表示你分组后的 group 列名，即 cols 合并后的列名
                        values_to = "test_value")
  colnames(plot_dat)[7] <- 'median_ratio'
  plot_dat$source <- factor(c(rep(c('CK','BaMV'),n)),levels = c('BaMV','CK'))
  plot_dat$func<-factor(plot_dat$func,unique(classed_m6A$func))
  plot_dat$new_name <- paste0(plot_dat$name,sep='_',plot_dat$type,sep='_',plot_dat$source)
  plot_dat$new_name <- paste0(plot_dat$name,sep='_',plot_dat$type,sep='_',plot_dat$source)
  plot_dat$new_name <- factor(plot_dat$new_name,rev(unique(plot_dat$new_name)))
  #plot_dat$name<-paste0(plot_dat$name,'_',plot_dat$func)
  #RNAi、NBL name:
  plot_dat$name <- paste0(plot_dat$type,'  ',plot_dat$name)
  plot_dat$name<-factor(plot_dat$name,rev(unique(plot_dat$name)),ordered = T)
  #####设置颜色列：plot_dat$func；yticks:plot_dat$func
  nl =unique(classed_m6A$func)
  nl_order=c()
  for (j in c(length(unique(classed_m6A$func)):1)){
    n_str=paste0('n',j)
    assign(n_str,nrow(classed_m6A[classed_m6A$func==nl[j],]))
    print(c(nl[j],get(n_str)))
    nl_order<-append(nl_order,get(n_str))
  }
  #,'#EEE685'
  color_n=c("#BC8F8F","#6CA6CD","#20B2AA",'#5F9EA0','#EEAD0E') 
  x_cols <- rep(color_n[1:length(nl)],times=nl_order)
  HEIGHT=round(nrow(classed_m6A) / 60,0)*17
  #pdf('F:\\project\\DRS\\result\\geneset\\geneset_barplot\\NBL_m6Aratio2.pdf',height = HEIGHT,width = 20)
  p1 <- ggplot(plot_dat,aes(x =name, y =median_ratio))+
    geom_bar(stat = 'identity',aes(fill = source),position = position_dodge(0.8)) + #使柱子并排放置
    theme_bw()+ theme(axis.text.y= element_text(colour=x_cols))+
    theme(text=element_text(size=20,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
          axis.text.x = element_text(size=20))+ # 设置X轴文字大小
    #scale_x_continuous(expand = c(0, 0)) +
    #scale_x_discrete(limits=c('CK','bamv')) +
    theme(legend.position="bottom")+
    theme(legend.position='none') + 
    ylim(0,1.0)+coord_flip()+#+facet_grid(.~class)
    ylab('median m6A ratio')+xlab("")
  p1
  #dev.off()
  ###PAL
  
  colnames(ck_median_pal)[1]<- 'gene_name';colnames(bm_median_pal)[1]<- 'gene_name'
  writer_eraser_reader_pal <- left_join(classed_m6A,ck_median_pal,by='gene_name')
  writer_eraser_reader_pal <- left_join(writer_eraser_reader_pal,bm_median_pal,by='gene_name')
  colnames(writer_eraser_reader_pal)[6:7] <- c('ck_median_pal','bm_median_pal')
  
  #writer_eraser_reader_pal[is.na(writer_eraser_reader_pal$bm_median_pal),5] = 0
  #writer_eraser_reader_pal[is.na(writer_eraser_reader_pal$ck_median_pal),4] = 0
  
  pal_plot_dat <-writer_eraser_reader_pal %>% 
    tidyr::pivot_longer(cols = c('ck_median_pal','bm_median_pal',), # 变的列，也就是即将合并为一列的列
                        names_to = "test_group", # 表示你分组后的 group 列名，即 cols 合并后的列名
                        values_to = "test_value")
  colnames(pal_plot_dat)[7] <- 'median_pal'
  pal_plot_dat$source <- factor(c(rep(c('CK','BaMV'),n)),levels = c('BaMV','CK'))
  pal_plot_dat$func<-factor(pal_plot_dat$func,unique(classed_m6A$func))
  pal_plot_dat$new_name <- paste0(pal_plot_dat$name,sep='_',pal_plot_dat$type,sep='_',pal_plot_dat$source)
  pal_plot_dat <- pal_plot_dat[match(plot_dat$new_name,pal_plot_dat$new_name),]
  pal_plot_dat$new_name <- factor(pal_plot_dat$new_name,rev(pal_plot_dat$new_name))
  #plot_dat$name<-paste0(plot_dat$gene_name,'_',plot_dat$func)
  pal_plot_dat$name <- paste0(pal_plot_dat$type,'  ',pal_plot_dat$name)
  pal_plot_dat$name<-factor(pal_plot_dat$name,rev(unique(pal_plot_dat$name)),ordered = T)
  
  #pdf('F:\\project\\DRS\\result\\geneset\\geneset_barplot\\NBL_m6Aapl2.pdf',height = HEIGHT,width = 20)
  p2 <- ggplot(pal_plot_dat,aes(x =name, y =median_pal))+
    geom_bar(stat = 'identity',aes(fill = source),position = position_dodge(0.8)) + #使柱子并排放置
    theme_bw()+ theme(axis.text.y= element_blank(),axis.ticks.y= element_blank())+#xlab("median PAL")+
    theme(text=element_text(size=20,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
          axis.text.x = element_text(size=20),axis.title.y =element_blank() )+  # 设置X轴文字大小
    #scale_x_continuous(expand = c(0, 0)) +
    #scale_x_discrete(limits=c('CK','bamv')) +
    theme(legend.position="bottom")+
    theme(legend.position='none') +
    ylim(0,200)+coord_flip()+#+facet_grid(.~class)
    ylab('median PAL')
  p2
  ###pas
  
  writer_eraser_reader_pas <- left_join(classed_m6A,dla_pas_bm_ck_use,by='gene_name')
  pas_plot_dat <- writer_eraser_reader_pas
  pas_plot_dat<-pas_plot_dat[match(unique(plot_dat$gene_name),pas_plot_dat$gene_name),] 
  pas_plot_dat$name <- paste0(pas_plot_dat$type,'  ',pas_plot_dat$name)
  pas_plot_dat$name <-factor(pas_plot_dat$name,rev(pas_plot_dat$name))
  k=0
  if (!all(is.na(pas_plot_dat$RED))){
    k=k+1
    a=min(na.omit(pas_plot_dat$RED))-0.5;b=max(na.omit(pas_plot_dat$RED))
    p3 <- ggplot(pas_plot_dat,aes(x = name, y = RED))+
      geom_bar(stat = 'identity',position = position_dodge(0.8),fill = "#FF6A6A") + #使柱子并排放置
      theme_bw()+ theme(axis.text.y= element_blank(),axis.ticks.y= element_blank())+#xlab("median PAL")+
      theme(text=element_text(size=20,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
            axis.text.x = element_text(size=20) ,axis.title.y =element_blank())+  # 设置X轴文字大小和轴标签
      #scale_x_continuous(expand = c(0, 0)) +
      #scale_x_discrete(limits=c('CK','bamv')) +
      theme(legend.position="bottom")+ 
      theme(legend.position='none') + #标签
      coord_flip()+#+facet_grid(.~class)#坐标转换和分面
      ylab('RED')
    p3
  }
  
  ####TPM
  writer_eraser_reader_tpm <- left_join(classed_m6A,tpm[,c(18:24)],by='gene_name')
  row.names(writer_eraser_reader_tpm) <- writer_eraser_reader_tpm$gene_name
  tpm_SE_dat <- writer_eraser_reader_tpm[,c(2,4:11)]
  #tpm_SE_dat[is.na(tpm_SE_dat[,2]),2:7] = 0
  
  data_long1 <- tpm_SE_dat %>% 
    tidyr::pivot_longer(cols = c('Ck_1_TPM','Ck_2_TPM','Ck_3_TPM','BaMV_1_TPM','BaMV_2_TPM','BaMV_3_TPM'), # 变的列，也就是即将合并为一列的列
                        names_to = "test_group", # 表示你分组后的 group 列名，即 cols 合并后的列名
                        values_to = "test_value") # 表示对应的值的列名
  data_long1$source <- factor(rep(c(rep('CK',3),rep('BaMV',3)),n),levels = c('BaMV','CK'))
  tpm_plot_dat <- summarySE(data = data_long1,measurevar = "test_value",groupvars = c("name","source"))
  tpm_plot_dat <- merge(tpm_plot_dat,writer_eraser_reader_tpm,by='name')
  tpm_plot_dat$new_name <- paste0(tpm_plot_dat$name,'_',tpm_plot_dat$type,'_',tpm_plot_dat$source)
  tpm_plot_dat<-tpm_plot_dat[match(plot_dat$new_name,tpm_plot_dat$new_name),]
  tpm_plot_dat$new_name <- factor(tpm_plot_dat$new_name,levels=rev(unique(plot_dat$new_name)),ordered =F)
  tpm_plot_dat$name <- paste0(tpm_plot_dat$type,'  ',tpm_plot_dat$name)
  tpm_plot_dat$name <- factor(tpm_plot_dat$name,rev(unique(tpm_plot_dat$name)),ordered = F)
  p4 <-ggplot(tpm_plot_dat,aes(x=name,y=test_value,fill=source))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_errorbar(aes(ymin=test_value-se, ymax=test_value+se),width=.1,position=position_dodge(.5))+
    #scale_fill_manual(values=c("gray30","gray60","gray80")) +
    theme_bw()+xlab("")+ylab("")+
    theme(text=element_text(size=20,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
          axis.text.x = element_text(size=20))+ # 设置X轴文字大小
    theme(axis.text = element_blank()) + 
    theme(axis.ticks = element_blank()) + xlab("")+ylab("TPM")+
    theme(legend.position='none') +
    coord_flip()
  p4
  
  ###abundance
  
  writer_eraser_reader_abundance <- left_join(classed_m6A,abundance[,c(1,8:13)],by='gene_name')
  row.names(writer_eraser_reader_abundance) <- writer_eraser_reader_abundance$gene_name
  abundance_SE_dat <- writer_eraser_reader_abundance[,c(2,4:11)]
  data_long2 <- abundance_SE_dat %>% 
    tidyr::pivot_longer(cols = c('bamv_rep1(Normalized)','bamv_rep2(Normalized)','bamv_rep3(Normalized)','ck_rep1(Normalized)','ck_rep2(Normalized)','ck_rep3(Normalized)'), # 变的列，也就是即将合并为一列的列
                        names_to = "test_group", # 表示你分组后的 group 列名，即 cols 合并后的列名
                        values_to = "test_value") # 表示对应的值的列名
  data_long2$source <- factor(rep(c(rep('BaMV',3),rep('CK',3)),n),levels = c('BaMV','CK'))
  abundance_plot_dat <- summarySE(data = data_long2,measurevar = "test_value",groupvars = c("name","source"))
  abundance_plot_dat <- merge(abundance_plot_dat,writer_eraser_reader_abundance,by='name')
  abundance_plot_dat$new_name <- paste0(abundance_plot_dat$name,'_',abundance_plot_dat$type,'_',abundance_plot_dat$source)
  abundance_plot_dat<-abundance_plot_dat[match(plot_dat$new_name,abundance_plot_dat$new_name),]
  abundance_plot_dat$new_name <- factor(abundance_plot_dat$new_name,levels=rev(unique(plot_dat$new_name)),ordered =F)
  abundance_plot_dat$name <- paste0(abundance_plot_dat$type,'  ',abundance_plot_dat$name)
  abundance_plot_dat$name <- factor(abundance_plot_dat$name,rev(unique(abundance_plot_dat$name)),ordered = F)
  #pdf('F:\\project\\DRS\\result\\geneset\\geneset_barplot\\NBL_m6Aratio2.pdf',height = HEIGHT,width = 20)
  p5 <- ggplot(abundance_plot_dat,aes(x=new_name,y=test_value,fill=source))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_errorbar(aes(ymin=test_value-se, ymax=test_value+se),width=.1,position=position_dodge(.5))+
    #scale_fill_manual(values=c("gray30","gray60","gray80")) +
    theme_bw()+xlab("")+ylab("Abundance")+
    theme(text=element_text(size=20,face = "bold"), #设置文字的字体字号（确保汉字可以显示）
          axis.text.x = element_text(size=20))+ # 设置X轴文字大小
    theme(axis.text = element_blank()) + #轴标题文字
    theme(axis.ticks = element_blank()) +#坐标轴刻度文字
    coord_flip()
  
  p5
  
  pdf(paste0('F:\\project\\DRS\\result\\geneset\\geneset_barplot\\',i,'_medianratio_PAL_RED_TPM_Abundance.pdf'),height = HEIGHT,width = 25)
  if (k==1){
    design="123445"
    print(wrap_plots(list(p1,p2,p3,p4,p5),design =design))}
  else{design="12445"
  print(wrap_plots(list(p1,p2,p4,p5),design =design))}
  
  dev.off()
  
}
 
#######

all_gene_median_ratio <- read.csv("F:\\project\\DRS\\result\\02_m6A\\ratio_df\\all_gene_median_ratio.csv",sep='\t')
mergeclassed_geneset
 
