
# Analysis of cdc28-as vs WT experiment LK29
# part of a TMT16plex

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# function
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


TMT16.protg <- read.delim("cdc28as_proteinGroups.txt")

TMT16.phos <- read.delim("cdc28as_phospho (STY)Sites.txt")

#Remove reverse & contaminant #####

TMT16.protg <- subset(TMT16.protg, Reverse != "+" & Potential.contaminant != "+")
TMT16.phos <- subset(TMT16.phos, Reverse != "+" & Potential.contaminant != "+")

#Change 'id' columns to reflect source table #####

colnames(TMT16.protg)[colnames(TMT16.protg) == "id"] <- "protg.id"
colnames(TMT16.phos)[colnames(TMT16.phos) == "id"] <- "phosphos.id"


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

TMT16.phos$Gene.name <- gene.map.table$Gene.name[match(TMT16.phos$Protein, gene.map.table$OLN)]

#### Create new ID cols in the protg table ####
#Create new Majority.protein.ID col with only FIRST (one) majority protein ID
#for cross referencing prot between tables

TMT16.protg$Majority.protein.ID <- strsplit.extract(TMT16.protg$Majority.protein.IDs, ";",1)

#Create Gene.name col in the protg table
TMT16.protg$Gene.name <- gene.map.table$Gene.name[match(TMT16.protg$Majority.protein.ID,gene.map.table$OLN)]

# variable to hold reporter column names
protg.cols<-grep("Reporter.intensity.corrected.",colnames(TMT16.protg),value=T)

phos.cols<-grep("Reporter.intensity.corrected.(.+).16plex_PE___1",colnames(TMT16.phos),value=T)

### how many sites are 'high localization'

TMT16.phos.filt<-TMT16.phos %>%
    filter(Localization.prob>0.75)

# change zeros to NAs to avoid them during normalization
library(data.table)

# new cols that have zeros changed to NAs
protg.na.cols<-paste0("na.",protg.cols)
phos.na.cols<-paste0("na.",phos.cols)

protg.TMT.DT<-as.data.table(TMT16.protg)
phos.TMT.DT<-as.data.table(TMT16.phos.filt)

protg.TMT.DT[,c(protg.na.cols) := lapply(.SD,function(x) ifelse(x==0,NA,x)), .SD=protg.cols]
phos.TMT.DT[,c(phos.na.cols) := lapply(.SD,function(x) ifelse(x==0,NA,x)), .SD=phos.cols]

TMT16.protg<-as.data.frame(protg.TMT.DT)
TMT16.phos.filt<-as.data.frame(phos.TMT.DT)

# colMedian normalization

# colMedians function
colMedians <- function (df, na.rm = T) {
    apply(df, 2, median, na.rm)
}

# create new column names for normalised intensities
# 
# ONLY ANALYSE WT AND CDC28 samples
# # 1,2=WT_MI
# 3,4=cdc28_MI
# 9,10=WT_MII
# 11,12=cdc28_MII

# so columns 1-4 and 9-12

norm.tmt.cols <- paste0("norm_",c(1:4,9:12))

# only analyse cdc28 samples
protg.na.cols<-protg.na.cols[c(1:4,9:12)]
phos.na.cols<-phos.na.cols[c(1:4,9:12)]

# calculate reporter MEDIAN intensities
protg_tot_med <- colMedians(TMT16.protg[,protg.na.cols],na.rm=T)
phos_tot_med <- colMedians(TMT16.phos.filt[,phos.na.cols],na.rm=T)

# # Identify the target median, which in this case is the median of the median of each experiment

protg_target_median <- median(protg_tot_med)
phos_target_median <- median(phos_tot_med)

# # # Identify the scaling factors needed to scale each individual experiment to reach the target_med
protg_scaling_factors <- as.numeric( protg_target_median / protg_tot_med )
phos_scaling_factors <- as.numeric( phos_target_median / phos_tot_med )

# # Scale the TMT intensities of each experiment so the median
#  becomes the target median

TMT16.protg[,norm.tmt.cols] <- sweep(TMT16.protg[,protg.na.cols], 2, protg_scaling_factors, FUN = "*")
TMT16.phos.filt[,norm.tmt.cols] <- sweep(TMT16.phos.filt[,phos.na.cols], 2, phos_scaling_factors, FUN = "*")

# # check the colMedians
colMedians(TMT16.protg[,norm.tmt.cols])
colMedians(TMT16.phos.filt[,norm.tmt.cols])

##### New values called "stoic" for "stoichiometric" #####

# extract first sequence window 
# (some phos-sites have multiple and
# this messes up motif analysis)
TMT16.phos.filt$Sequence.window <- strsplit.extract(TMT16.phos.filt$Sequence.window, ';', 1)

#create new phos.stoic.cols 
phos.stoic.cols <- paste0("stoic_",(c(1:4,9:12)))

# data prep
#copy protg IDs in phosphos df to new column called first.protg.id
#this is necessary for matching protg.ids between the phosphos and protg table

TMT16.phos.filt$first.protg.id <- strsplit.extract(TMT16.phos.filt$Protein.group.IDs, ';', 1)

phos.stoic.temp.df <- TMT16.protg[match(TMT16.phos.filt$first.protg.id, TMT16.protg$protg.id), norm.tmt.cols]

#this divides the intensity of the phosphos over the protg
#normalized intensity then multiplies by 1000 (to get easier number)
#This creates "stoichiometric" phos. site changes...

TMT16.phos.filt[,phos.stoic.cols] <- mapply(function(x,y) {x/y},TMT16.phos.filt[,norm.tmt.cols],  phos.stoic.temp.df[,norm.tmt.cols]) * 1000

# Plot the phos.stoic.colss (normalized to protg)
boxplot(log2(TMT16.phos.filt[,phos.stoic.cols]),main="Phos/protg norm phosphos",ylab="log2 intensity")

# change NAs back to zeros
protg.TMT.DT<-as.data.table(TMT16.protg)
phos.TMT.DT<-as.data.table(TMT16.phos.filt)

protg.TMT.DT[,c(norm.tmt.cols) := lapply(.SD,function(x) ifelse(is.na(x),0,x)), .SD=norm.tmt.cols]
phos.TMT.DT[,c(phos.stoic.cols) := lapply(.SD,function(x) ifelse(is.na(x),0,x)), .SD=phos.stoic.cols]

TMT16.protg<-as.data.frame(protg.TMT.DT)
TMT16.phos.filt<-as.data.frame(phos.TMT.DT)


# plots of intensity distrubtions after normalization
TMT16.phos.filt[,phos.stoic.cols] %>%
    pivot_longer(cols=phos.stoic.cols,names_to="sample",values_to="intensity",names_prefix="stoic_") %>%
    ggplot(aes(x=fct_inorder(sample),y=intensity,fill=sample)) +
    geom_boxplot()+
    scale_y_log10()+
    labs(x="sample",y="log10(intensity)")+
    theme_classic()+
theme(legend.position = "none")


# Add gene name site col
TMT16.phos.filt <- TMT16.phos.filt %>%
    mutate(gene.name_site=paste(Gene.name,Amino.acid,Positions.within.proteins,sep="-"))

# what are the identities of the columns
# 1,2=WT_MI
# 3,4=cdc28_MI
# 9,10=WT_MII
# 11,12=cdc28_MII

stoic.named.cols<-c("WT_MI_1","WT_MI_2","cdc28_MI_1","cdc28_MI_2",
                    "WT_MII_1","WT_MII_2","cdc28_MII_1","cdc28_MII_2")
    
# rename phos_stoic_cols
# (copies phos.stoic.cols to new columns named stoic.named cols)
TMT16.phos.filt[,stoic.named.cols]<-TMT16.phos.filt[,phos.stoic.cols]

# rename protg norm.tmt.cols
# (copies norm.tmt.colss to new cols where namess are the same as phos.stoic.cols)
protg.named.cols<-stoic.named.cols
TMT16.protg[,protg.named.cols]<-TMT16.protg[,norm.tmt.cols]


# plots of intensity distrubtions after normalization
protg.dist.plot<-TMT16.protg[,protg.named.cols] %>%
    pivot_longer(cols=protg.named.cols,names_to="sample",values_to="intensity",names_prefix="norm_") %>%
    ggplot(aes(x=fct_inorder(sample),y=log2(intensity),fill=sample)) +
    geom_boxplot(aes(fill=sample,alpha=0.2),position="dodge",outlier.size = 0.25, size = 0.25)+
    labs(x="",y="log2(intensity)",title="proteins")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90,vjust=0.6,size=rel(1.5)),legend.position="none",axis.text.y=element_text(size=rel(1.5)),axis.title=element_text(size=rel(1.5)),plot.title=element_text(size=rel(2)))+
    coord_flip()

print_plot(protg.dist.plot,"protg_dist_plot.pdf",6,4)


phos.dist.plot<-TMT16.phos.filt[,stoic.named.cols] %>%
    pivot_longer(cols=stoic.named.cols,names_to="sample",values_to="intensity") %>%
    ggplot(aes(x=fct_inorder(sample),y=log2(intensity),fill=sample)) +
    geom_boxplot(aes(fill=sample,alpha=0.2),position="dodge",outlier.size = 0.25, size = 0.25)+
    labs(x="",y="log2(intensity)",title="phospho-sites")+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90,vjust=0.6,size=rel(1.5)),legend.position="none",axis.text.y=element_text(size=rel(1.5)),axis.title=element_text(size=rel(1.5)),plot.title=element_text(size=rel(2)))+
    coord_flip()

print_plot(phos.dist.plot,"phos_dist_plot.pdf",6,4)

# colnames by strain
WT_cols<-grep("WT",stoic.named.cols,value=T)
cdc28_cols<-grep("cdc28",stoic.named.cols,value=T)

#nrow(TMT16.phos.filt) #2723

# create cols indicating whether site was detected in 1/2 reps
TMT16.protg<-TMT16.protg %>%
    mutate(WT_MI_1rep=apply(.[protg.named.cols[1:2]],1,function(x) any(x>0))) %>%
    mutate(cdc28_MI_1rep=apply(.[protg.named.cols[3:4]],1,function(x) any(x>0))) %>%
    mutate(WT_MII_1rep=apply(.[protg.named.cols[5:6]],1,function(x) any(x>0))) %>%
    mutate(cdc28_MII_1rep=apply(.[protg.named.cols[7:8]],1,function(x) any(x>0)))


TMT16.phos.filt<-TMT16.phos.filt %>%
    mutate(WT_MI_1rep=apply(.[stoic.named.cols[1:2]],1,function(x) any(x>0))) %>%
    mutate(cdc28_MI_1rep=apply(.[stoic.named.cols[3:4]],1,function(x) any(x>0))) %>%
    mutate(WT_MII_1rep=apply(.[stoic.named.cols[5:6]],1,function(x) any(x>0))) %>%
    mutate(cdc28_MII_1rep=apply(.[stoic.named.cols[7:8]],1,function(x) any(x>0)))

# plot number of sites found in at least 1 rep of each strain

# grouped condition cols
cond_cols<-grep("1rep",colnames(TMT16.phos.filt),value=T)

total_protg_bar<-TMT16.protg %>%
    summarize(across(cond_cols,function(x) sum(x))) %>%
    pivot_longer(cols=all_of(cond_cols),names_to="condition",values_to="number") %>%
    ggplot(aes(x=condition,y=number,fill=condition))+
    geom_bar(stat="identity")+
    geom_text(aes(label=number),vjust=0.5,hjust=-0.2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90),legend.position="none",
          axis.text=element_text(size=rel(1)))+
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
    coord_flip()+
    labs(x="",y="",title="number of proteins detected in each condition")

total_sites_bar<-TMT16.phos.filt %>%
    summarize(across(cond_cols,function(x) sum(x))) %>%
    pivot_longer(cols=all_of(cond_cols),names_to="condition",values_to="number") %>%
    ggplot(aes(x=condition,y=number,fill=condition))+
    geom_bar(stat="identity")+
    geom_text(aes(label=number),vjust=0.5,hjust=-0.2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90),legend.position="none",
          axis.text=element_text(size=rel(1)))+
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
    coord_flip()+
    labs(x="",y="",title="number of sites detected in each condition")

print_plot(total_protg_bar,"total_protg_bar.pdf",6,4)
print_plot(total_sites_bar,"total_sites_bar.pdf",6,4)

# calculate number of matches of given motif in each condition

match_df<-function(motif){
    rbind(match=TMT16.phos.filt %>%
              mutate(motif=grepl(motif,Sequence.window)) %>%
              filter(motif) %>%
              summarize(across(cond_cols,function(x) sum(x))),
          no_match=
              TMT16.phos.filt %>%
              mutate(motif=!grepl(motif,Sequence.window)) %>%
              filter(motif) %>%
              summarize(across(cond_cols,function(x) sum(x)))) %>%
        rbind(total=colSums(.)) %>%
        rbind(perc_match=(.[1,]/.[3,])*100)
}


# percent of matching sites in each condition
perc_min_cdk<-match_df("...............[ST]P..............") %>%
    tail(1) %>%
    pivot_longer(cols=all_of(cond_cols),names_to="condition",values_to="proportion") %>%
    mutate(proportion=round(proportion,3)) %>%
    ggplot(aes(x=condition,y=proportion,fill=condition))+
    geom_bar(stat="identity")+
    geom_text(aes(label=proportion),vjust=0.5,hjust=-0.2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90),legend.position="none",
          axis.text=element_text(size=rel(1)))+
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
    coord_flip()+
    labs(x="",y="percent",title="percent [ST]P sites")

perc_strict_cdk<-match_df("...............[ST]P.[KR]............") %>%
    tail(1) %>%
    pivot_longer(cols=all_of(cond_cols),names_to="condition",values_to="proportion") %>%
    mutate(proportion=round(proportion,3)) %>%
    ggplot(aes(x=condition,y=proportion,fill=condition))+
    geom_bar(stat="identity")+
    geom_text(aes(label=proportion),vjust=0.5,hjust=-0.2)+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90),legend.position="none",
          axis.text=element_text(size=rel(1)))+
    scale_y_continuous(expand=expansion(mult=c(0,0.2)))+ # 20% more plot area after max value
    coord_flip()+
    labs(x="",y="percent",title="percent [ST]Px[KR] sites")

print_plot(perc_min_cdk,"perc_min_cdk.pdf",6,4)
print_plot(perc_strict_cdk,"perc_strict_cdk.pdf",6,4)

# select conditions to compare by fisher test

two_fish<-function(df,set1,set2){
    no_match1<-df %>% select(any_of(set1)) %>% slice(2) %>% unlist() %>% unname()
    no_match2<-df %>% select(any_of(set2)) %>% slice(2) %>% unlist() %>% unname()
    match1<-df %>% select(any_of(set1)) %>% slice(1) %>% unlist() %>% unname()
    match2<-df %>% select(any_of(set2)) %>% slice(1) %>% unlist() %>% unname()
    
 
    mat<-matrix(c(match1,(no_match1),match2,(no_match2)),nrow=2,dimnames= list(c("match","no match"),
                               c("set1","set2")))
    
    #return(mat)
    fish_result<-fisher.test(mat)
    
    #return(fish_result)
    mat<-mat %>% as.data.frame() %>%
        rbind("total"=colSums(.)) %>%
        rbind(perc_match=(.[1,]/.[3,])*100) %>%
        cbind("pval"=fish_result$p.value) %>%
    mutate(stars=ifelse(pval<0.001,"***",
                        ifelse(pval<0.01,"**",
                               ifelse(pval<0.05,"*",NA))))

    return(mat)
}

library(ggsignif)
library(scales)
library(patchwork)

two_fish_bar<-function(motif,set1,set2,title){
   out_df<- two_fish(match_df(motif),set1,set2) %>%
        tail(1) %>%
        pivot_longer(-c(pval,stars),values_to="perc",names_to="sample") %>%
        mutate(sample=ifelse(grepl("set1",sample),set1,set2)) %>%
        mutate(sample=factor(sample,levels=c(set1,set2))) #keep levels consistent
   
   #return(out_df)
   out_df %>%
       ggplot(aes(x=sample,y=perc,fill=sample))+
       geom_bar(aes(fill=sample),stat="identity")+
       geom_signif(comparisons = list(c(set1,set2)),
                   map_signif_level=TRUE,
                   annotations = out_df$stars[1],
                   size=0.25)+
       scale_y_continuous(expand=expansion(mult=c(0,0.2)),labels = label_number(accuracy=0.1))+
       theme_classic()+
       theme(axis.text.x = element_text(angle=90,vjust=0.6),legend.position="none",axis.text=element_text(size=rel(1)))+
       labs(y="percent matches",x="",title=title,axis.line=element_line(linewidth=0.1))
}


# looking for motif depletion in all conditions
# cdk motif
# min
p1<-two_fish_bar("...............[ST]P..............","WT_MI_1rep","cdc28_MI_1rep","[ST]P sites detected")
p2<-two_fish_bar("...............[ST]P..............","WT_MII_1rep","cdc28_MII_1rep","[ST]P sites detected")

# strict
p3<-two_fish_bar("...............[ST]P.[KR]............","WT_MI_1rep","cdc28_MI_1rep","[ST]Px[KR] sites detected")
p4<-two_fish_bar("...............[ST]P.[KR]............","WT_MII_1rep","cdc28_MII_1rep","[ST]Px[KR] sites detected")

print_plot((p1+p2),"STP_fish_bar.pdf",8,6)
print_plot((p3+p4),"STPxKR_fish_bar.pdf",8,6)

library(ggseqlogo)
# logo of all sites detected in each condition
short_logo<-function(df,meth,title){
    ggseqlogo(substr(df %>% select(Sequence.window) %>% unlist(),11,21),method=meth)+
        ggtitle(paste0(title,nrow(df)))+
        theme(legend.position="none")+
        scale_x_continuous(breaks=1:11,labels=-5:5)
}

logo1<-short_logo(TMT16.phos.filt %>% filter(WT_MI_1rep),"prob","WT_MI_1rep, n=")
logo2<-short_logo(TMT16.phos.filt %>% filter(WT_MII_1rep),"prob","WT_MII_1rep, n=")
logo3<-short_logo(TMT16.phos.filt %>% filter(cdc28_MI_1rep),"prob","cdc28_MI_1rep, n=")
logo4<-short_logo(TMT16.phos.filt %>% filter(cdc28_MII_1rep),"prob","cdc28_MII_1rep, n=")

print_plot((logo1+logo2+logo3+logo4),"detection_logos.pdf",8,5)

#### COUNTING SITES THAT ARE DETECTED IN 1/2 reps, THOSE THAT change by a given fold change
# calculate fold changes
#stoic.named.cols

TMT16.phos.filt<-TMT16.phos.filt %>%
    mutate(cdc28_WT_MI_fc=apply(.[,stoic.named.cols],1,function(x) median((x[3]/x[1]),(x[4]/x[2]),na.rm=T))) %>%
    mutate(cdc28_WT_MII_fc=apply(.[,stoic.named.cols],1,function(x) median((x[7]/x[5]),(x[8]/x[6]),na.rm=T)))


### MODIFY To only evaluate finite FCs
# need function for number of match/no match in given two groups
match_df_sub<-function(motif,group1,group2,fc_group,cutoff){
    mat<-cbind(
        above=c(match=TMT16.phos.filt %>%
    mutate(motifl=grepl(motif,Sequence.window)) %>%
    filter(.data[[fc_group]]>cutoff) %>%
    filter(motifl&.data[[group1]]&.data[[group2]]&is.finite(.data[[fc_group]])) %>%
    nrow(),
no_match=TMT16.phos.filt %>%
    mutate(motifl=!grepl(motif,Sequence.window)) %>%
    filter(.data[[fc_group]]>cutoff) %>%
    filter(motifl&.data[[group1]]&.data[[group2]]&is.finite(.data[[fc_group]])) %>%
    nrow()),
below=c(match=TMT16.phos.filt %>%
    mutate(motifl=grepl(motif,Sequence.window)) %>%
    filter(.data[[fc_group]]<cutoff) %>%
    filter(motifl&.data[[group1]]&.data[[group2]]&is.finite(.data[[fc_group]])) %>%
    nrow(),
no_match=TMT16.phos.filt %>%
    mutate(motifl=!grepl(motif,Sequence.window)) %>%
    filter(.data[[fc_group]]<cutoff) %>%
    filter(motifl&.data[[group1]]&.data[[group2]]&is.finite(.data[[fc_group]])) %>%
    nrow()))
    
    fish_result<-fisher.test(mat)
    
    mat<-mat %>% as.data.frame() %>%
        rbind("total"=colSums(.)) %>%
        rbind(perc_match=(.[1,]/.[3,])*100) %>%
        cbind("pval"=fish_result$p.value) %>%
        mutate(stars=ifelse(pval<0.001,"***",
                            ifelse(pval<0.01,"**",
                                   ifelse(pval<0.05,"*",NA))))
    
    return(mat)
}

# modified fish bar function

two_fish_bar_sub<-function(motif,group1,group2,fc_group,cutoff,title){
    out_df<- match_df_sub(motif,group1,group2,fc_group,cutoff) %>%
        tail(1) %>%
        pivot_longer(-c(pval,stars),values_to="perc",names_to="sample") %>%
        mutate(sample=ifelse(grepl("above",sample),"above","below")) %>%
        mutate(sample=factor(sample,levels=c("above","below"))) #keep levels consistent
    
   # return(out_df)
    out_df %>%
        ggplot(aes(x=sample,y=perc,fill=sample))+
        geom_bar(aes(fill=sample),stat="identity")+
        geom_signif(comparisons = list(c("above","below")),
                    map_signif_level=TRUE,
                    annotations = out_df$stars[1],
                    size=0.25)+
        scale_y_continuous(expand=expansion(mult=c(0,0.2)),labels = label_number(accuracy=0.1))+
        theme_classic()+
        theme(axis.text.x = element_text(angle=90,vjust=0.6),legend.position="none",axis.text=element_text(size=rel(1)))+
        labs(y="percent matches",x="",title=title,axis.line=element_line(linewidth=0.1))
}


# export sequences for icelogos
library(openxlsx)


## fold change cutoff 1.5

TMT16.phos.filt<-TMT16.phos.filt %>%
    mutate(cdc28_WT_MI_1.5_cat=ifelse(is.finite(cdc28_WT_MI_fc)&WT_MI_1rep&cdc28_MI_1rep,
                                  ifelse(cdc28_WT_MI_fc>1.5,"cdc28_INCR",
                                         ifelse(cdc28_WT_MI_fc<(1/1.5), "cdc28_DECR","NC")),"NC")) %>%
    mutate(cdc28_WT_MII_1.5_cat=ifelse(is.finite(cdc28_WT_MII_fc)&WT_MII_1rep&cdc28_MII_1rep,
                                   ifelse(cdc28_WT_MII_fc>1.5,"cdc28_INCR",
                                          ifelse(cdc28_WT_MII_fc<(1/1.5), "cdc28_DECR","NC")),"NC")) %>%
    # need to make sure factor levels are consistent
    mutate(WT_cdc28_MI_1.5_cat=factor(cdc28_WT_MI_1.5_cat,levels=c("NC","cdc28_DECR","cdc28_INCR"))) %>%
    mutate(WT_cdc28_MII_1.5_cat=factor(cdc28_WT_MII_1.5_cat,levels=c("NC","cdc28_DECR","cdc28_INCR")))
    
    
# Logos FC 1.5 cutoff

p1<-short_logo(TMT16.phos.filt %>% filter(cdc28_WT_MI_1.5_cat=="cdc28_DECR"),"prob","cdc28/WT MI fc < 1/1.5, n=")
p2<-short_logo(TMT16.phos.filt %>% filter(cdc28_WT_MI_1.5_cat=="cdc28_INCR"),"prob","cdc28/WT MI fc > 1.5, n=")
p3<-short_logo(TMT16.phos.filt %>% filter(cdc28_WT_MII_1.5_cat=="cdc28_DECR"),"prob","cdc28/WT MII fc < 1/1.5, n=")
p4<-short_logo(TMT16.phos.filt %>% filter(cdc28_WT_MII_1.5_cat=="cdc28_INCR"),"prob","cdc28/WT MII fc > 1.5, n=")

print_plot((p1+p2),"cdc28_WT_MI_1.5_logos.pdf",8,3)
print_plot((p3+p4),"cdc28_WT_MII_1.5_logos.pdf",8,3)

# Fisher tests

# fold change 1.5
# STP motif
p1<-two_fish_bar_sub("...............[ST]P..............","WT_MI_1rep","cdc28_MI_1rep","cdc28_WT_MI_fc",1.5,"[ST]P MI cdc28/WT fc 1.5")

p2<-two_fish_bar_sub("...............[ST]P..............","WT_MII_1rep","cdc28_MII_1rep","cdc28_WT_MII_fc",1.5,"[ST]P MII cdc28/WT fc 1.5")

print_plot((p1+p2),"STP_fc_1.5.pdf",8,6)

# STPxKR motif
p1<-two_fish_bar_sub("...............[ST]P.[KR]............","WT_MI_1rep","cdc28_MI_1rep","cdc28_WT_MI_fc",1.5,"[ST]Px[KR] MI cdc28/WT fc 1.5")

p2<-two_fish_bar_sub("...............[ST]P.[KR]............","WT_MII_1rep","cdc28_MII_1rep","cdc28_WT_MII_fc",1.5,"[ST]Px[KR] MII cdc28/WT fc 1.5")


print_plot((p1+p2),"STPxKR_fc_1.5.pdf",8,6)


# export sequences for icelogos
library(openxlsx)

# fold change 1.5
phos_seq_dfs<-list()

# cdc28 vs WT MI
phos_seq_dfs[[1]]<-TMT16.phos.filt %>%
    filter(WT_MI_1rep&cdc28_MI_1rep) %>%  # note fc finite-ness already assessed in sigcat calling
    mutate(Sequence.window_5=substr(Sequence.window,11,21)) %>%
    arrange(cdc28_WT_MI_1.5_cat) %>%
    group_by(cdc28_WT_MI_1.5_cat) %>%
    mutate(id=row_number()) %>%
    pivot_wider(names_from=cdc28_WT_MI_1.5_cat,values_from=Sequence.window_5,id_cols=id) %>%
    select(-id)

# cdc28 vs WT MII
phos_seq_dfs[[2]]<-TMT16.phos.filt %>%
    filter(WT_MII_1rep&cdc28_MII_1rep) %>%  # note fc finite-ness already assessed in sigcat calling
    mutate(Sequence.window_5=substr(Sequence.window,11,21)) %>%
    arrange(cdc28_WT_MII_1.5_cat) %>%
    group_by(cdc28_WT_MII_1.5_cat) %>%
    mutate(id=row_number()) %>%
    pivot_wider(names_from=cdc28_WT_MII_1.5_cat,values_from=Sequence.window_5,id_cols=id) %>%
    select(-id)


names(phos_seq_dfs)<-c("cdc28_WT_MI_1.5_cat","cdc28_WT_MII_1.5_cat")

write.xlsx(phos_seq_dfs,file="cdc28_WT_1.5_sig_cat_seqs.xlsx")

# sig cat IDS for supplement
# fold change 1.5
phos_seq_dfs<-list()

# cdc28 vs WT MI
phos_seq_dfs[[1]]<-TMT16.phos.filt %>%
    filter(WT_MI_1rep&cdc28_MI_1rep) %>%  # note fc finite-ness already assessed in sigcat calling
    #mutate(Sequence.window_5=substr(Sequence.window,11,21)) %>%
    arrange(cdc28_WT_MI_1.5_cat) %>%
    group_by(cdc28_WT_MI_1.5_cat) %>%
    mutate(id=row_number()) %>%
    pivot_wider(names_from=cdc28_WT_MI_1.5_cat,values_from=gene.name_site,id_cols=id) %>%
    select(-id)

# cdc28 vs WT MII
phos_seq_dfs[[2]]<-TMT16.phos.filt %>%
    filter(WT_MII_1rep&cdc28_MII_1rep) %>%  # note fc finite-ness already assessed in sigcat calling
    #mutate(Sequence.window_5=substr(Sequence.window,11,21)) %>%
    arrange(cdc28_WT_MII_1.5_cat) %>%
    group_by(cdc28_WT_MII_1.5_cat) %>%
    mutate(id=row_number()) %>%
    pivot_wider(names_from=cdc28_WT_MII_1.5_cat,values_from=gene.name_site,id_cols=id) %>%
    select(-id)

names(phos_seq_dfs)<-c("cdc28_WT_MI_1.5_cat","cdc28_WT_MII_1.5_cat")

write.xlsx(phos_seq_dfs,file="cdc28_WT_1.5_sig_cat_ids.xlsx")

# pie charts
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


pie1<-pie(TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep),"cdc28_WT_MI_1.5_cat","cdc28/WT MI fc 1.5")
pie2<-pie(TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep),"cdc28_WT_MII_1.5_cat","cdc28/WT MII fc 1.5")

print_plot(pie1,"cdc28_WT_MI_fc_1.5_pie.pdf",4,3)

print_plot(pie2,"cdc28_WT_MII_fc_1.5_pie.pdf",4,3)


## New fisher barplot section to match cdc5as analysis

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


#fisher_barplot<-function(df1,df2,motif,subset1,subset2,title){

p1<-fisher_barplot(TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep&is.finite(cdc28_WT_MI_fc)&cdc28_WT_MI_1.5_cat=="NC"),
                   TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep&is.finite(cdc28_WT_MI_fc)&cdc28_WT_MI_1.5_cat=="cdc28_DECR"),
                   "...............[ST]P..............","NC","DECR","cdc28/WT MI [ST]*P")


p2<-fisher_barplot(TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep&is.finite(cdc28_WT_MI_fc)&cdc28_WT_MI_1.5_cat=="NC"),
                   TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep&is.finite(cdc28_WT_MI_fc)&cdc28_WT_MI_1.5_cat=="cdc28_DECR"),
                   "...............[ST]P.[KR]............","NC","DECR","cdc28/WT MI [ST]*Px[KR]")

p3<-fisher_barplot(TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep&is.finite(cdc28_WT_MI_fc)&cdc28_WT_MI_1.5_cat=="NC"),
                   TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep&is.finite(cdc28_WT_MI_fc)&cdc28_WT_MI_1.5_cat=="cdc28_DECR"),
                   ".............[DEN].[ST]...............","NC","DECR","cdc28/WT MI [DEN]x[ST]*")


p4<-fisher_barplot(TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep&is.finite(cdc28_WT_MII_fc)&cdc28_WT_MII_1.5_cat=="NC"),
                   TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep&is.finite(cdc28_WT_MII_fc)&cdc28_WT_MII_1.5_cat=="cdc28_DECR"),
                   "...............[ST]P..............","NC","DECR","cdc28/WT MII [ST]*P")

p5<-fisher_barplot(TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep&is.finite(cdc28_WT_MII_fc)&cdc28_WT_MII_1.5_cat=="NC"),
                   TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep&is.finite(cdc28_WT_MII_fc)&cdc28_WT_MII_1.5_cat=="cdc28_DECR"),
                   "...............[ST]P.[KR]............","NC","DECR","cdc28/WT MII [ST]*Px[KR]")

p6<-fisher_barplot(TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep&is.finite(cdc28_WT_MII_fc)&cdc28_WT_MII_1.5_cat=="NC"),
                   TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep&is.finite(cdc28_WT_MII_fc)&cdc28_WT_MII_1.5_cat=="cdc28_DECR"),
                   ".............[DEN].[ST]P..............","NC","DECR","cdc28/WT MII [DEN]x[ST]*")

print_plot((p1+p2+p3+p4)/(p5+p6+p6+p6),"cdc28_WT_cdk_fish_bar.pdf",6,14)

# Proportion of motif sites respond to Cdc28 inhibition?

# MI proportions
p1<-pie(TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep&is.finite(cdc28_WT_MI_fc)) %>%
        filter(grepl("...............[ST]P..............",Sequence.window)),"cdc28_WT_MI_1.5_cat","cdc28/WT MI [ST]*P fc 1.5")

p2<-pie(TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep&is.finite(cdc28_WT_MI_fc)) %>%
        filter(grepl("...............[ST]P.[KR]............",Sequence.window)),"cdc28_WT_MI_1.5_cat","cdc28/WT MI [ST]*Px[KR] fc 1.5")

p3<-pie(TMT16.phos.filt %>% filter(WT_MI_1rep&cdc28_MI_1rep&is.finite(cdc28_WT_MI_fc)) %>%
        filter(grepl(".............[DEN].[ST]...............",Sequence.window)),"cdc28_WT_MI_1.5_cat","cdc28/WT MI [DEN]x[ST]* fc 1.5")

print_plot((p1+p2+p3),"cdc28_WT_MI_motif_pies.pdf",10,5)

# MII proportions
p1<-pie(TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep&is.finite(cdc28_WT_MI_fc)) %>%
            filter(grepl("...............[ST]P..............",Sequence.window)),"cdc28_WT_MII_1.5_cat","cdc28/WT MII [ST]*P fc 1.5")

p2<-pie(TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep&is.finite(cdc28_WT_MI_fc)) %>%
            filter(grepl("...............[ST]P.[KR]............",Sequence.window)),"cdc28_WT_MII_1.5_cat","cdc28/WT MII [ST]*Px[KR] fc 1.5")

p3<-pie(TMT16.phos.filt %>% filter(WT_MII_1rep&cdc28_MII_1rep&is.finite(cdc28_WT_MI_fc)) %>%
            filter(grepl(".............[DEN].[ST]...............",Sequence.window)),"cdc28_WT_MII_1.5_cat","cdc28/WT MII [DEN]x[ST]* fc 1.5")

print_plot((p1+p2+p3),"cdc28_WT_MII_motif_pies.pdf",10,5)

