
# Analysis of cdc5as vs WT data
# part of a TMT10plex
# 
# Author: Lori Koch

library(tidyverse)

# functions
strsplit.extract <- function (x, split, index) {
    x <- as.character(x)
    x[x == ""] <- ";"
    x <- sapply(strsplit(x, split, fixed = TRUE), "[[", index)
    return(x)
}

print_plot<-function(plot,filename.pdf,width,height){
    print(plot)
    pdf(filename.pdf, width=width, height=height) 
    print(plot)
    dev.off()
}

# load data tables
protg <- read.delim("cdc5as_proteinGroups.txt")

phosphos <- read.delim("cdc5as_phospho (STY)Sites.txt")

#Remove reverse & contaminant #####

protg <- subset(protg, Reverse != "+" & Potential.contaminant != "+")
phosphos<- subset(phosphos, Reverse != "+" & Potential.contaminant != "+")

#Change 'id' columns to reflect source table #####

colnames(protg)[colnames(protg) == "id"] <- "protg.id"
colnames(phosphos)[colnames(phosphos) == "id"] <- "phosphos.id"

### Make gene name cols ####

### #load yeast gene map table
#I downloaded this file from https://www.uniprot.org/docs/yeast
#I modified the original file so it only includes the table and not the heading part at the top

gene.map.table <-read.csv("Yeast_Uniprot2.csv")

#data prep so there isn't white space and only the first gene name is in the gene column
gene.map.table$Gene<-trimws(gene.map.table$Gene, "right")
gene.map.table$Gene.name <- strsplit.extract(gene.map.table$Gene, ";",1)

#### Create new ID columns ####
#Create a Gene.name column in phosphos table

phosphos$Gene.name <- gene.map.table$Gene.name[match(phosphos$Protein, gene.map.table$OLN)]

#### Create new ID cols in the protg table ####
#Create new Majority.protein.ID col with only FIRST (one) majority protein ID
#for cross referencing prot between tables

protg$Majority.protein.ID <- strsplit.extract(protg$Majority.protein.IDs, ";",1)

#Create Gene.name col in the protg table
protg$Gene.name <- gene.map.table$Gene.name[match(protg$Majority.protein.ID,gene.map.table$OLN)]

# plot raw data

# protg_raw_cols<-grep("Reporter.intensity.corrected.(.+)_N",colnames(protg),value=T)
# 
# phos_raw_cols<-grep("Reporter.intensity.corrected.(.+)_PE___1",colnames(phosphos),value=T)
# 
# protg_raw<-protg %>%
#     select(Gene.name,all_of(protg_raw_cols)) %>%
#     rename_with(.,function(x) gsub("Reporter.intensity.corrected.","",x),.cols=all_of(protg_raw_cols)) %>%
#     pivot_longer(-Gene.name,names_to="sample",values_to="intensity") %>%
#     ggplot(aes(x=fct_inorder(sample),y=log10(intensity),fill=sample))+
#     geom_boxplot()+
#     theme_classic()+
#     theme(axis.text.x = element_text(angle=90),legend.position="none",
#           axis.text=element_text(size=rel(1)))
# 
# print_plot(protg_raw,"protg_raw.pdf",6,4)
# 
# phos_raw<-phosphos %>%
#     select(Gene.name,all_of(phos_raw_cols)) %>%
#     rename_with(.,function(x) gsub("Reporter.intensity.corrected.","",x),.cols=all_of(phos_raw_cols)) %>%
#     pivot_longer(-Gene.name,names_to="sample",values_to="intensity") %>%
#     ggplot(aes(x=fct_inorder(sample),y=log10(intensity),fill=sample))+
#     geom_boxplot()+
#     theme_classic()+
#     theme(axis.text.x = element_text(angle=90),legend.position="none",
#           axis.text=element_text(size=rel(1)))
#
# print_plot(phos_raw,"phos_raw.pdf",6,5)

# analyse WT and cdc5 as samples
protg_raw_cols<-grep("Reporter.intensity.corrected.(.+)_N",colnames(protg),value=T)

phos_raw_cols<-grep("Reporter.intensity.corrected.(.+)_PE___1",colnames(phosphos),value=T)

#  samples 1-2 are WT reps 1,2
#  samples 7-8 are cdc5as reps 1,2

protg_raw_cols<-protg_raw_cols[c(1:2,7:8)]

phos_raw_cols<-phos_raw_cols[c(1:2,7:8)]

# NORMALIZATION
# change zeros to NAs to avoid them during normalization
library(data.table)

# new cols that have zeros changed to NAs
protg.na.cols<-paste0("na.",protg_raw_cols)
phos.na.cols<-paste0("na.",phos_raw_cols)

protg.TMT.DT<-as.data.table(protg)
phos.TMT.DT<-as.data.table(phosphos)

protg.TMT.DT[,c(protg.na.cols) := lapply(.SD,function(x) ifelse(x==0,NA,x)), .SD=protg_raw_cols]
phos.TMT.DT[,c(phos.na.cols) := lapply(.SD,function(x) ifelse(x==0,NA,x)), .SD=phos_raw_cols]

# convert data tables back to data frames
# new columns are called protg.na.cols and phos.na.cols
protg<-as.data.frame(protg.TMT.DT)
phosphos<-as.data.frame(phos.TMT.DT)

# colMedian normalization

# colMedians function
colMedians <- function (df, na.rm = T) {
    apply(df, 2, median, na.rm)
}

# create new column names for normalised intensities

norm.tmt.cols <- paste0("norm_",seq(1:4))

# calculate reporter MEDIAN intensities
protg_tot_med <- colMedians(protg[,protg.na.cols],na.rm=T)
phos_tot_med <- colMedians(phosphos[,phos.na.cols],na.rm=T)

# # Identify the target median, which in this case is the median of the median of each experiment

protg_target_median <- median(protg_tot_med)
phos_target_median <- median(phos_tot_med)

# # # Identify the scaling factors needed to scale each individual experiment to reach the target_med
protg_scaling_factors <- as.numeric( protg_target_median / protg_tot_med )
phos_scaling_factors <- as.numeric( phos_target_median / phos_tot_med )

# # Scale the TMT intensities of each experiment so the median
#  becomes the target median

protg[,norm.tmt.cols] <- sweep(protg[,protg.na.cols], 2, protg_scaling_factors, FUN = "*")
phosphos[,norm.tmt.cols] <- sweep(phosphos[,phos.na.cols], 2, phos_scaling_factors, FUN = "*")

# # check the colMedians
colMedians(protg[,norm.tmt.cols])
colMedians(phosphos[,norm.tmt.cols])

##### New values called "stoic" for "stoichiometric" #####

#create new phos.stoic.cols 
phos.stoic.cols <- paste0("stoic_",seq(1:4))

# data prep
#copy protg IDs in phosphos df to new column called first.protg.id
#this is necessary for matching protg.ids between the phosphos and protg table

phosphos$first.protg.id <- strsplit.extract(phosphos$Protein.group.IDs, ';', 1)

phos.stoic.temp.df <- protg[match(phosphos$first.protg.id, protg$protg.id), norm.tmt.cols]

#this divides the intensity of the phosphos over the protg
#normalized intensity then multiplies by 1000 (to get easier number)
#This creates "stoichiometric" phos. site changes...

phosphos[,phos.stoic.cols] <- mapply(function(x,y) {x/y},phosphos[,norm.tmt.cols],  phos.stoic.temp.df[,norm.tmt.cols]) * 1000

# Plot the phos.stoic.colss (normalized to protg)
boxplot(log2(phosphos[,phos.stoic.cols]),main="Phos/protg norm phosphos",ylab="log2 intensity")

# change NAs back to zeros in the norm/stoic cols
protg.TMT.DT<-as.data.table(protg)
phos.TMT.DT<-as.data.table(phosphos)

protg.TMT.DT[,c(norm.tmt.cols) := lapply(.SD,function(x) ifelse(is.na(x),0,x)), .SD=norm.tmt.cols]
phos.TMT.DT[,c(phos.stoic.cols) := lapply(.SD,function(x) ifelse(is.na(x),0,x)), .SD=phos.stoic.cols]

protg<-as.data.frame(protg.TMT.DT)
phosphos<-as.data.frame(phos.TMT.DT)

# plot of intensity distrubtions after normalization
# protein
protg_norm<-protg[,norm.tmt.cols] %>%
    pivot_longer(cols=norm.tmt.cols,names_to="sample",values_to="intensity",names_prefix="norm_") %>%
    ggplot(aes(x=fct_inorder(sample),y=log2(intensity),fill=sample)) +
    geom_boxplot(aes(fill=sample,alpha=0.2),position="dodge",outlier.size = 0.25, size = 0.25)+
    labs(x="",y="log2(intensity)",title="protein")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90,vjust=0.6,size=rel(1.5)),legend.position="none",axis.text.y=element_text(size=rel(1.5)),axis.title=element_text(size=rel(1.5)),plot.title=element_text(size=rel(2)))+
    coord_flip()

print_plot(protg_norm,"protg_norm.pdf",6,4)


# phospho sites
phos_norm<-phosphos[,phos.stoic.cols] %>%
    pivot_longer(cols=phos.stoic.cols,names_to="sample",values_to="intensity",names_prefix="stoic_") %>%
    ggplot(aes(x=fct_inorder(sample),y=log2(intensity),fill=sample)) +
    geom_boxplot(aes(fill=sample,alpha=0.2),position="dodge",outlier.size = 0.25, size = 0.25)+
    labs(x="",y="log2(intensity)",title="phospho-sites")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90,vjust=0.6,size=rel(1.5)),legend.position="none",axis.text.y=element_text(size=rel(1.5)),axis.title=element_text(size=rel(1.5)),plot.title=element_text(size=rel(2)))+
    coord_flip()

print_plot(phos_norm,"phos_norm.pdf",6,4)


# phospho-site table data prep
# extract first sequence window 
# (some phos-sites have multiple and
# this messes up motif analysis)
phosphos$Sequence.window <- strsplit.extract(phosphos$Sequence.window, ';', 1)

# Add gene name site col
phosphos <- phosphos %>%
    mutate(gene.name_site=paste(Gene.name,Amino.acid,Positions.within.proteins,sep="-"))

# filter to high localiz score
phos_filt<-phosphos %>%
    filter(Localization.prob>0.75)

#### General numbers in each sample

# named cols
named_cols<-c("WT_1","WT_2","cdc5as_1","cdc5as_2")  

# copy norm cols into named cols
protg[,named_cols]<-protg[,norm.tmt.cols]

phos_filt[,named_cols]<-phos_filt[,phos.stoic.cols]

# how many proteins were detected at least once in each sample
# # add new col indicating whether prot was detected in each strain
protg<-protg %>%
    mutate(WT_1rep=apply(.[,named_cols[1:2]],1,function(x) any(x>0))) %>%
    mutate(cdc5as_1rep=apply(.[,named_cols[3:4]],1,function(x) any(x>0))) %>%
    mutate(WT_2rep=apply(.[,named_cols[1:2]],1,function(x) all(x>0))) %>%
    mutate(cdc5as_2rep=apply(.[,named_cols[3:4]],1,function(x) all(x>0))) %>%
    mutate(cdc5as_WT_2rep=apply(.[,named_cols],1,function(x) all(x>0)))

    
phos_filt<-phos_filt %>%
    mutate(WT_1rep=apply(.[,named_cols[1:2]],1,function(x) any(x>0))) %>%
    mutate(cdc5as_1rep=apply(.[,named_cols[3:4]],1,function(x) any(x>0))) %>%
    mutate(WT_2rep=apply(.[,named_cols[1:2]],1,function(x) all(x>0))) %>%
    mutate(cdc5as_2rep=apply(.[,named_cols[3:4]],1,function(x) all(x>0))) %>%
    mutate(cdc5as_WT_2rep=apply(.[,named_cols],1,function(x) all(x>0)))

# variable to hold colnames
cond_cols<-grep("rep",colnames(protg),value=T)
# separate samples not including all 4
cond_cols2<-cond_cols[1:4]

# summarize and plot
# each sample separately
protg_all_numbers<-protg %>%
    summarize(across(named_cols,function(x) sum(x>0))) %>%
    pivot_longer(cols=all_of(named_cols),names_to="condition",values_to="number") %>%
    ggplot(aes(x=fct_inorder(condition),y=number,fill=condition))+
    geom_bar(stat="identity")+
    geom_text(aes(label=number),vjust=0.5,hjust=-0.2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90),legend.position="none",
          axis.text=element_text(size=rel(1)))+
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
    coord_flip()+
    labs(x="",y="",title="number of proteins detected in each sample")

print_plot(protg_all_numbers,"protg_all_numbers.pdf",6,4)

# how many in each strain/condition
protg_cond_numbers<-protg %>%
    summarize(across(cond_cols,function(x) sum(x))) %>%
    pivot_longer(cols=all_of(cond_cols),names_to="condition",values_to="number") %>%
    ggplot(aes(x=fct_inorder(condition),y=number,fill=condition))+
    geom_bar(stat="identity")+
    geom_text(aes(label=number),vjust=0.5,hjust=-0.2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90),legend.position="none",
          axis.text=element_text(size=rel(1)))+
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
    coord_flip()+
    labs(x="",y="",title="cumulative number of proteins")

print_plot(protg_cond_numbers,"protg_cond_numbers.pdf",6,4)


# number of sites in each sample separately
phos_all_numbers<-phos_filt %>%
    summarize(across(named_cols,function(x) sum(x>0))) %>%
    pivot_longer(cols=all_of(named_cols),names_to="condition",values_to="number") %>%
    ggplot(aes(x=fct_inorder(condition),y=number,fill=condition))+
    geom_bar(stat="identity")+
    geom_text(aes(label=number),vjust=0.5,hjust=-0.2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90),legend.position="none",
          axis.text=element_text(size=rel(1)))+
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
    coord_flip()+
    labs(x="",y="",title="number of sites detected in each sample")

print_plot(phos_all_numbers,"phos_all_numbers.pdf",6,4)

# number of sites in each condition/strain
phos_cond_numbers<-phos_filt %>%
    summarize(across(cond_cols,function(x) sum(x))) %>%
    pivot_longer(cols=all_of(cond_cols),names_to="condition",values_to="number") %>%
    ggplot(aes(x=fct_inorder(condition),y=number,fill=condition))+
    geom_bar(stat="identity")+
    geom_text(aes(label=number),vjust=0.5,hjust=-0.2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90),legend.position="none",
          axis.text=element_text(size=rel(1)))+
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
    coord_flip()+
    labs(x="",y="",title="cumulative number of sites")

print_plot(phos_cond_numbers,"phos_cond_numbers.pdf",6,4)

## motifs of all sequences detected in each strain
library(ggseqlogo)
# logo of all sites detected in each condition
short_logo<-function(df,meth,title){
    ggseqlogo(substr(df %>% select(Sequence.window) %>% unlist(),11,21),method=meth)+
        ggtitle(paste0(title,nrow(df)))+
        theme(legend.position="none")+
        scale_x_continuous(breaks=1:11,labels=-5:5)
}

logo1<-short_logo(phos_filt %>% filter(WT_2rep),"prob","WT 2 reps, n=")
logo2<-short_logo(phos_filt %>% filter(cdc5as_2rep),"prob","cdc5as 2 reps, n=")

library(patchwork)
print_plot((logo1+logo2),"all_sites_logos.pdf",10,4)
# no obvious differences

# calculate fold changes between strains
protg<-protg %>%
    mutate(cdc5as_WT_fc=apply(.[,named_cols],1,function(x) {median(c((x[3]/x[1]),(x[4]/x[2])))}))
   
phos_filt<-phos_filt %>%
    mutate(cdc5as_WT_fc=apply(.[,named_cols],1,function(x) {median(c((x[3]/x[1]),(x[4]/x[2])))}))

    
# variables to call the fcs 
fcs<-grep("fc",colnames(protg),value=T)


# distributions of fcs
# # order by channel number differences to see if it trends
fc_dist<-phos_filt %>%
    select(all_of(fcs)) %>% 
    filter(if_all(fcs,is.finite)) %>%
    select(all_of(fcs)) %>%
    pivot_longer(cols=all_of(fcs),names_to="sample",values_to="fc") %>%
    #mutate(sample=factor(sample, levels=c("spo13_WT_fc","spo13cdc5as_WT_fc","cdc5as_WT_fc","spo13m2_WT_fc","spo13_spo13cdc5as_fc","spo13_cdc5as_fc","spo13m2_spo13_fc","spo13m2_spo13cdc5as_fc","spo13m2_cdc5as_fc"))) %>%
    ggplot(aes(x=log10(fc)))+
    geom_histogram(bins=100)+
    facet_wrap(~sample)+
    geom_vline(xintercept=log10(1/2),linetype="dashed",color="red")+
    geom_vline(xintercept=log10(2),linetype="dashed",color="red")

print_plot(fc_dist,"fc_dist.pdf",8,5)

# FC CUTOFFS

# add FC 1.5 sig col
phos_filt<-phos_filt %>%
    mutate(cdc5as_WT_fc_1.5_sig=ifelse(cdc5as_2rep&WT_2rep&
                                           is.finite(cdc5as_WT_fc) & !is.na(cdc5as_WT_fc),
                                       ifelse(cdc5as_WT_fc<(1/1.5),"DECR",
                                              ifelse(cdc5as_WT_fc>1.5,"INCR","NC")),"NC")) %>%
    # make sure factor variables are always the same
    # adjust so factor levels are all the same
    mutate(cdc5as_WT_fc_1.5_sig=factor(cdc5as_WT_fc_1.5_sig,levels=c("NC","INCR","DECR")))
    

# need to update filtered dfs to include new columns
# filtered dfs with only sites found in both strains for comparisons

phos_filt_cdc5as_WT<-phos_filt %>% filter(WT_2rep&cdc5as_2rep)

# pie chart
# pie charts numbers incr/decr in general strain 
library(ggrepel)

pie<-function(filtereddf,sigcat,title){
    filtereddf %>%
        group_by(.data[[sigcat]]) %>%
        summarize(n=n()) %>%
        ggplot(aes(x="",y=n,fill=.data[[sigcat]],label=n))+
        geom_bar(stat="identity")+
        geom_label_repel()+
        coord_polar("y", start=0)+
        theme_void()+
        scale_fill_manual(values=c("#7CAE00", "#00BFC4", "#C77CFF","blue"))+
        ggtitle(title)+
        theme(legend.title=element_blank())
}


p1<-pie(phos_filt_cdc5as_WT,"cdc5as_WT_fc_1.5_sig","cdc5as vs WT FC 1.5")


print_plot(p1,"FC1.5_pie.pdf",4,3)


# motifs 1.5 FC alone

# DECR

fc_logo1<-phos_filt_cdc5as_WT %>%
    filter(cdc5as_WT_fc_1.5_sig=="DECR") %>%
    short_logo(.,"prob","cdc5as/WT < (1/1.5), n=")

# INCR

fc_logo2<-phos_filt_cdc5as_WT %>%
    filter(cdc5as_WT_fc_1.5_sig=="INCR") %>%
    short_logo(.,"prob","cdc5as/WT > 1.5, n=")

print_plot(fc_logo1+fc_logo2,"cdc5as_WT_fc_1.5_logos.pdf",10,3)

# export sequences for icelogo

# excels of seq for icelogos
library(openxlsx)

phos_seq_dfs<-list()

phos_seq_dfs[[1]]<-phos_filt_cdc5as_WT %>%
    mutate(Sequence.window_5=substr(Sequence.window,11,21)) %>%
    arrange(cdc5as_WT_fc_1.5_sig) %>%
    group_by(cdc5as_WT_fc_1.5_sig) %>%
    mutate(id=row_number()) %>%
    pivot_wider(names_from=cdc5as_WT_fc_1.5_sig,values_from=Sequence.window_5,id_cols=id) %>%
    select(-id)

names(phos_seq_dfs)<-c("cdc5as_WT_fc_1.5_sig")

write.xlsx(phos_seq_dfs,file="cdc5as_WT_FC_1.5_sig_cat_seqs.xlsx")

# sig cat IDS for supplement
phos_seq_dfs<-list()

phos_seq_dfs[[1]]<-phos_filt_cdc5as_WT %>%
    #mutate(Sequence.window_5=substr(Sequence.window,11,21)) %>%
    arrange(cdc5as_WT_fc_1.5_sig) %>%
    group_by(cdc5as_WT_fc_1.5_sig) %>%
    mutate(id=row_number()) %>%
    pivot_wider(names_from=cdc5as_WT_fc_1.5_sig,values_from=gene.name_site,id_cols=id) %>%
    select(-id)

names(phos_seq_dfs)<-c("cdc5as_WT_fc_1.5_sig")

write.xlsx(phos_seq_dfs,file="cdc5as_WT_FC_1.5_sig_cat_ids.xlsx")

# what is percent of motif matching sites in each ssample 

# show each sample separately
cond_cols2<-cond_cols[1:4]


phos_filt<-phos_filt %>%
    mutate(STP=grepl("...............[ST]P..............",Sequence.window)) %>%
    mutate(STPxKR=grepl("...............[ST]P.[KR]............",Sequence.window)) %>%
    mutate(DENxST=grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
    mutate(DENxSTF=grepl(".............[DEN].[ST]F..............",Sequence.window)) %>%
    mutate(DENxSTFG=grepl(".............[DEN].[ST][FG]..............",Sequence.window)) %>%
    mutate(DENxSThyd=grepl(".............[DEN].[ST][LIMFY]..............",Sequence.window))
    

match_df<-function(motif){
    rbind(match=phos_filt %>%
              mutate(motif=grepl(motif,Sequence.window)) %>%
              filter(motif) %>%
              summarize(across(cond_cols2,function(x) sum(x))),
          no_match=
              phos_filt %>%
              mutate(motif=!grepl(motif,Sequence.window)) %>%
              filter(motif) %>%
              summarize(across(cond_cols2,function(x) sum(x)))) %>%
        rbind(total=colSums(.)) %>%
        rbind(perc_match=(.[1,]/.[3,])*100)
}

perc_bar<-function(motif,title){
    match_df(motif) %>%
        tail(1) %>%
        pivot_longer(cols=all_of(cond_cols2),names_to="condition",values_to="proportion") %>%
        mutate(proportion=round(proportion,3)) %>%
        ggplot(aes(x=condition,y=proportion,fill=condition))+
        geom_bar(stat="identity")+
        geom_text(aes(label=proportion),vjust=0.5,hjust=-0.2)+
        theme_classic()+
        theme(axis.text.x = element_text(angle=90),legend.position="none",
              axis.text=element_text(size=rel(1)))+
        scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
        coord_flip()+
        labs(x="",y="percent",title=title)
}

# plot of percent matching sites detected in each sample
polo_detect<-perc_bar(".............[DEN].[ST]...............","percent [DEN]x[ST]* sites detected")

poloF_detect<-perc_bar(".............[DEN].[ST]F..............","percent [DEN]x[ST]*F sites detected")

poloFG_detect<-perc_bar(".............[DEN].[ST][FG]..............","percent [DEN]x[ST]*[FG] sites detected")
    
polohyd_detect<-perc_bar(".............[DEN].[ST][LIMFY]..............","percent [DEN]x[ST]*[LIMFY] sites detected")
    
STP_detect<-perc_bar("...............[ST]P..............","percent [ST]*P sites detected")
    
STPxKR_detect<-perc_bar("...............[ST]P.[KR]............","percent [ST]*Px[KR] sites detected")
    

print_plot((polo_detect+poloF_detect+poloFG_detect)/(polohyd_detect+STP_detect+STPxKR_detect),
           "motif_detect_bar.pdf",15,6)

# total number of motif matching sites rather than perc
num_match_bar<-function(motif,title){
    match_df(motif) %>%
        head(1) %>%
        pivot_longer(cols=all_of(cond_cols2),names_to="condition",values_to="proportion") %>%
        mutate(proportion=round(proportion,3)) %>%
        ggplot(aes(x=condition,y=proportion,fill=condition))+
        geom_bar(stat="identity")+
        geom_text(aes(label=proportion),vjust=0.5,hjust=-0.2)+
        theme_classic()+
        theme(axis.text.x = element_text(angle=90),legend.position="none",
              axis.text=element_text(size=rel(1)))+
        scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
        coord_flip()+
        labs(x="",y="number_matches",title=title)
}



# plot of number matching sites detected in each sample
polo_detect<-num_match_bar(".............[DEN].[ST]...............","[DEN]x[ST]* sites detected")

poloF_detect<-num_match_bar(".............[DEN].[ST]F..............","[DEN]x[ST]*F sites detected")

poloFG_detect<-num_match_bar(".............[DEN].[ST][FG]..............","[DEN]x[ST]*[FG] sites detected")

polohyd_detect<-num_match_bar(".............[DEN].[ST][LIMFY]..............","[DEN]x[ST]*[LIMFY] sites detected")

STP_detect<-num_match_bar("...............[ST]P..............","[ST]*P sites detected")

STPxKR_detect<-num_match_bar("...............[ST]P.[KR]............","[ST]*Px[KR] sites detected")


print_plot((polo_detect+poloF_detect+poloFG_detect)/(polohyd_detect+STP_detect+STPxKR_detect),
           "motif_detect_num_bar.pdf",15,6)


# sub-polo motifs


# cdc5as sub-polo motif
logo1<-phos_filt %>%
    filter(cdc5as_WT_fc_1.5_sig=="DECR") %>%
    filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
    short_logo(.,"prob","cdc5as vs WT DECR fc 1.5, n=")

logo2<-phos_filt %>%
    filter(cdc5as_WT_fc_1.5_sig=="INCR") %>%
    filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
    short_logo(.,"prob","cdc5as vs WT INCR fc 1.5, n=")


print_plot((logo1+logo2),"cdc5as_WT_DEN_sub_logo.pdf",10,3)


### Fisher tests for motifs in different sets ####

library(scales)
library(ggsignif)

fisher<-function(df1,df2,motif){
    total.df1 <- df1 %>% nrow() 
    total.df2 <- df2 %>% nrow()
    match.df1 <- df1 %>% filter(grepl(motif,Sequence.window)) %>% nrow()
    match.df2 <- df2  %>% filter(grepl(motif,Sequence.window)) %>% nrow()
    corr.mat<-matrix(c(match.df1,(total.df1-match.df1),match.df2,(total.df2-match.df2)),
                     nrow=2,
                     dimnames= list(c("match","no match"),
                                    c("df1","df2")))
    fisher.result<-fisher.test(corr.mat)
    
    corr.mat<-corr.mat %>% as.data.frame() %>%
        rbind(colSums(.)) %>%
        cbind(rowSums(.)) %>%
        rbind(perc_match=(.[1,]/.[3,])*100)
    
    rownames(corr.mat)[3]<-"col_totals"
    colnames(corr.mat)<-c("df1","df2","row_totals")
    
    corr.mat<-corr.mat %>%
        cbind("pval"=fisher.result$p.value)
    return(corr.mat)}

fisher_barplot<-function(df1,df2,motif,subset1,subset2,title){
    out_df<-fisher(df1,df2,motif) %>%
        tail(1) %>%
        mutate(stars=ifelse(pval<0.001,"***",
                            ifelse(pval<0.01,"**",
                                   ifelse(pval<0.05,"*",NA)))) %>%
        select(df1,df2,pval,stars) %>%
        pivot_longer(-c(pval,stars),values_to="perc",names_to="sample") %>%
        mutate(sample=ifelse(grepl("df1",sample),subset1,subset2)) %>%
        mutate(sample=factor(sample,levels=c(subset1,subset2))) # have to make sure levels are consistent
    
    out_df %>%
        ggplot(aes(x=sample,y=perc,fill=sample))+
        geom_bar(aes(fill=sample),stat="identity")+
        geom_signif(comparisons = list(c(subset1,subset2)), 
                    map_signif_level=TRUE,
                    annotations = out_df$stars[1],
                    size=0.25)+
        scale_y_continuous(expand=expansion(mult=c(0,0.2)),labels = label_number(accuracy=0.1))+
        theme_classic()+
        theme(axis.text.x = element_text(angle=90,vjust=0.6),legend.position="none",axis.text=element_text(size=rel(1)))+
        labs(y="percent matches",x="",title=title,,axis.line=element_line(linewidth=0.1))
}

## CDC5-as v WT comparisons

# FC 1.5
p1<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="DECR"),".............[DEN].[ST]...............",
                   "NC","DECR","[DEN]x[ST]*")
p2<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="DECR"),".............[DEN].[ST]F..............",
                   "NC","DECR","[DEN]x[ST]*F")
p3<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="DECR"),".............[DEN].[ST][G]..............",
                   "NC","DECR","[DEN]x[ST]*G")
p4<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="DECR"),".............[DEN].[ST][LIMFY]..............",
                   "NC","DECR","[DEN]x[ST]*[LIMFY]")
p5<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="DECR"),"...............[ST]P..............",
                   "NC","DECR","[ST]*P")
p6<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="DECR"),"...............[ST]P.[KR]............",
                   "NC","DECR","[ST]*Px[KR]")
p7<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="DECR"),"............R..[ST]...............",
                   "NC","DECR","Rxx[ST]*")
p8<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="DECR"),".............R.[ST]...............",
                   "NC","DECR","Rx[ST]*")

print_plot(((p1+p2+p3+p7)/(p4+p5+p6+p8)),"cdc5as_WT_fc_1.5_DECR_fish.pdf",6,14)


p1<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="INCR"),".............[DEN].[ST]...............",
                   "NC","INCR","[DEN]x[ST]*")
p2<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="INCR"),".............[DEN].[ST]F..............",
                   "NC","INCR","[DEN]x[ST]*F")
p3<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="INCR"),".............[DEN].[ST][G]..............",
                   "NC","INCR","[DEN]x[ST]*G")
p4<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="INCR"),".............[DEN].[ST][LIMFY]..............",
                   "NC","INCR","[DEN]x[ST]*[LIMFY]")
p5<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="INCR"),"...............[ST]P..............",
                   "NC","INCR","[ST]*P")
p6<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="INCR"),"...............[ST]P.[KR]............",
                   "NC","INCR","[ST]*Px[KR]")
p7<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="INCR"),"............R..[ST]...............",
                   "NC","INCR","Rxx[ST]*")
p8<-fisher_barplot(phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="NC"),phos_filt_cdc5as_WT %>% filter(cdc5as_WT_fc_1.5_sig=="INCR"),".............R.[ST]...............",
                   "NC","INCR","Rx[ST]*")

print_plot(((p1+p2+p3+p4)/(p5+p6+p7+p8)),"cdc5as_WT_fc_1.5_INCR_fish.pdf",6,14)
