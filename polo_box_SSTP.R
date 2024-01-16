# Analysis of number of polo boxes in proteins differentially phos in
# WT vs spo13
# 
# This script takes lists of proteins that were differentially phos in 
# WT vs spo13 
# from the timecourse
# and metaphase I experiments
# and then I retrieved tables from YeastMine
# of those gene lists that had the entire ORF protein sequence
# of those genes/proteins.
# Then I read those tables into R and used regex to find motif matches
# create plots
# and perform fisher tests.

library(dplyr)
library(clipr)

# read in results table from metaphase arrest (exp40e)
# and timecourse experiment 
# these csvs are produced by their respective R scripts
Exp40E_phos<-read.csv("Exp40E_phos.csv")
TC_phos<-read.csv("WT_spo13_timecourse_phos.csv")

# need to edit gene name so it matches SGD
Exp40E_phos<-Exp40E_phos %>%
    mutate(Protein=gsub("YGL251C_1","YGL251C",Protein))

# TC_phos %>%
#     filter(Protein=="YGL251C_1") %>% nrow() #0

# metaphase arrest gene lists
# proteins with decreased phoss in spo13del
# Exp40E_phos %>%
#     filter(spo13_WT_sig=="DECR") %>%
#     #filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() # n=454 genes
# 
# 
# # proteins with decreased phos in spo13del and which match the canonical Polo motif
# Exp40E_phos %>%
#     filter(spo13_WT_sig=="DECR") %>%
#     filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() #n=199
# 
# # proteins with decreased phoss in spo13m2
# Exp40E_phos %>%
#     filter(spo13m2_WT_sig=="DECR") %>%
#     #filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() #n-173
# 
# # proteins with decreased phos in spo13del and which match the canonical Polo motif
# Exp40E_phos %>%
#     filter(spo13m2_WT_sig=="DECR") %>%
#     filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() #n-78
# 
# # NC gene lists
# Exp40E_phos %>%
#     filter(spo13_WT_sig=="NC") %>%
#     #filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() # n=1777 genes
# #     
# # proteins with NC phos in spo13del and which match the canonical Polo motif
# Exp40E_phos %>%
#     filter(spo13_WT_sig=="NC") %>%
#     filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() #n=918
# 
# Exp40E_phos %>%
#     filter(spo13m2_WT_sig=="NC") %>%
#     #filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() # n=1839 genes
# #     
# # proteins with NC phos in spo13del and which match the canonical Polo motif
# Exp40E_phos %>%
#     filter(spo13m2_WT_sig=="NC") %>%
#     filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() #n=984
#  
# # TC gene lists
# # proteins with decreased phos at any 1/10 tiempoints in WT vs spo13
# TC_phos %>%
#     filter(sigdiff.sig.cat=="sigdiff_spo13_DECR") %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() #n=173
# 
# # proteins with decr phos in spo13 and matching polo motif
# TC_phos %>%
#     filter(sigdiff.sig.cat=="sigdiff_spo13_DECR") %>%
#     filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() #n=61
# 
# # NC TC lists
# TC_phos %>%
#     filter(sigdiff.sig.cat=="NC") %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() #941
# 
# TC_phos %>%
#     filter(sigdiff.sig.cat=="NC") %>%
#     filter(grepl(".............[DEN].[ST]...............",Sequence.window)) %>%
#     distinct(Protein) %>%
#     select(Protein) %>%
#     unlist() %>%
#     unname() %>% write_clip() #n=414

# retrieve and read in list of protein sequences from YeastMine
spo13_DECR_df<-read.delim("Exp40_DECR_spo13.tsv")
spo13_polo_DECR_df<-read.delim("Exp40_DECR_polo_spo13.tsv")
spo13m2_DECR_df<-read.delim("Exp40_DECR_spo13m2.tsv")
spo13m2_polo_DECR_df<-read.delim("Exp40_DECR_polo_spo13m2.tsv")
TC_spo13_DECR_df<-read.delim("TC_DECR_spo13.tsv")
TC_spo13_polo_DECR_df<-read.delim("TC_DECR_polo_spo13.tsv")

spo13_NC_df<-read.delim("Exp40_NC_spo13.tsv",quote = "" )
spo13_polo_NC_df<-read.delim("Exp40_NC_polo_spo13.tsv",quote = "" )
spo13m2_NC_df<-read.delim("Exp40_NC_spo13m2.tsv",quote = "" )
spo13m2_polo_NC_df<-read.delim("Exp40_NC_polo_spo13m2.tsv",quote = "" )
TC_spo13_NC_df<-read.delim("TC_NC_spo13.tsv",quote = "" )
TC_spo13_polo_NC_df<-read.delim("TC_NC_polo_spo13.tsv",quote = "" )

# trim long column names
colnames(spo13_DECR_df)<-gsub("Gene...","",colnames(spo13_DECR_df))
colnames(spo13_polo_DECR_df)<-gsub("Gene...","",colnames(spo13_polo_DECR_df))
colnames(spo13m2_DECR_df)<-gsub("Gene...","",colnames(spo13m2_DECR_df))
colnames(spo13m2_polo_DECR_df)<-gsub("Gene...","",colnames(spo13m2_polo_DECR_df))
colnames(TC_spo13_DECR_df)<-gsub("Gene...","",colnames(TC_spo13_DECR_df))
colnames(TC_spo13_polo_DECR_df)<-gsub("Gene...","",colnames(TC_spo13_polo_DECR_df))

colnames(spo13_NC_df)<-gsub("Gene...","",colnames(spo13_NC_df))
colnames(spo13_polo_NC_df)<-gsub("Gene...","",colnames(spo13_polo_NC_df))
colnames(spo13m2_NC_df)<-gsub("Gene...","",colnames(spo13m2_NC_df))
colnames(spo13m2_polo_NC_df)<-gsub("Gene...","",colnames(spo13m2_polo_NC_df))
colnames(TC_spo13_NC_df)<-gsub("Gene...","",colnames(TC_spo13_NC_df))
colnames(TC_spo13_polo_NC_df)<-gsub("Gene...","",colnames(TC_spo13_polo_NC_df))

# replace extra dots in colnames
colnames(spo13_DECR_df)<-gsub("\\.\\.\\.",".",colnames(spo13_DECR_df))
colnames(spo13_polo_DECR_df)<-gsub("\\.\\.\\.",".",colnames(spo13_polo_DECR_df))
colnames(spo13m2_DECR_df)<-gsub("\\.\\.\\.",".",colnames(spo13m2_DECR_df))
colnames(spo13m2_polo_DECR_df)<-gsub("\\.\\.\\.",".",colnames(spo13m2_polo_DECR_df))
colnames(TC_spo13_DECR_df)<-gsub("\\.\\.\\.",".",colnames(TC_spo13_DECR_df))
colnames(TC_spo13_polo_DECR_df)<-gsub("\\.\\.\\.",".",colnames(TC_spo13_polo_DECR_df))

colnames(spo13_NC_df)<-gsub("\\.\\.\\.",".",colnames(spo13_NC_df))
colnames(spo13_polo_NC_df)<-gsub("\\.\\.\\.",".",colnames(spo13_polo_NC_df))
colnames(spo13m2_NC_df)<-gsub("\\.\\.\\.",".",colnames(spo13m2_NC_df))
colnames(spo13m2_polo_NC_df)<-gsub("\\.\\.\\.",".",colnames(spo13m2_polo_NC_df))
colnames(TC_spo13_NC_df)<-gsub("\\.\\.\\.",".",colnames(TC_spo13_NC_df))
colnames(TC_spo13_polo_NC_df)<-gsub("\\.\\.\\.",".",colnames(TC_spo13_polo_NC_df))

# check list lengths
nrow(spo13_DECR_df) #454
nrow(spo13_polo_DECR_df) #199
nrow(spo13m2_DECR_df) #173
nrow(spo13m2_polo_DECR_df) #78
nrow(TC_spo13_DECR_df) #172
nrow(TC_spo13_polo_DECR_df) #61
nrow(spo13_NC_df) #1777
nrow(spo13_polo_NC_df) #918
nrow(spo13m2_NC_df) #1839
nrow(spo13m2_polo_NC_df) #984
nrow(TC_spo13_NC_df) #941
nrow(TC_spo13_polo_NC_df) #414

# add new col that lists all positions of motif in the protein

spo13_DECR_df<-spo13_DECR_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

spo13_polo_DECR_df<-spo13_polo_DECR_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

spo13m2_DECR_df<-spo13m2_DECR_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

spo13m2_polo_DECR_df<-spo13m2_polo_DECR_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

TC_spo13_DECR_df<-TC_spo13_DECR_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

TC_spo13_polo_DECR_df<-TC_spo13_polo_DECR_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

#NC
spo13_NC_df<-spo13_NC_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

spo13_polo_NC_df<-spo13_polo_NC_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

spo13m2_NC_df<-spo13m2_NC_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

spo13m2_polo_NC_df<-spo13m2_polo_NC_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos)%>% unlist())

TC_spo13_NC_df<-TC_spo13_NC_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number =ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos) %>% unlist())

TC_spo13_polo_NC_df<-TC_spo13_polo_NC_df %>%
    mutate(polo_box_pos=sapply(gregexpr("S[ST]P", Proteins.Sequence.Residues), function(x) ifelse(x == -1, 0, x))) %>%
    # logical column to indicate if any matches exist for subset later
    mutate(polo_box_motif=sapply(polo_box_pos, function(x) ifelse(any(x==0), F, T))) %>%
    mutate(polo_motif_number = ifelse(polo_box_motif==T,lengths(polo_box_pos),polo_box_pos) %>% unlist())

# Plot distribution of motifs in DECR v NC gene sets
library(ggplot2)
library(dplyr)
library(tidyr)
    
# spo13 v WT metaphase I
spo13_WT_box<-merge(spo13_NC_df %>%
          select(Proteins.Standard.Name,polo_motif_number),
      spo13_DECR_df %>%
          select(Proteins.Standard.Name,polo_motif_number), 
      by="Proteins.Standard.Name",all=T) %>%
    replace(.=="NULL", NA) %>%# replace NULL with NA 
    rename_with(~ gsub(".x", "_NC", .x, fixed = TRUE)) %>%
    rename_with(~ gsub(".y", "_DECR", .x, fixed = TRUE)) %>%
    rename_with(~gsub("polo_motif_number_","",.x)) %>%
    pivot_longer(-Proteins.Standard.Name,
                 values_to="number",names_to="group") %>%
    filter(!is.na(number)) %>%
    mutate(number=as.integer(number)) %>%
    group_by(group) %>%
    mutate(group_n=n()) %>% 
    mutate(group_n=paste0(group," (",group_n,")")) %>%
    ggplot(aes(y=number,x=group,fill=group_n))+
    geom_boxplot(alpha=0.5)+
    scale_y_continuous(breaks=seq(0,10,2))+
    coord_flip()+
    theme_classic()+
    labs(x="",y="number",title="spo13 v WT, number of S[ST]P motifs")

# spo13 v WT polo sites only metaphase I
spo13_WT_polo_box<-merge(spo13_polo_NC_df %>%
          select(Proteins.Standard.Name,polo_motif_number),
      spo13_polo_DECR_df %>%
          select(Proteins.Standard.Name,polo_motif_number), 
      by="Proteins.Standard.Name",all=T) %>%
    replace(.=="NULL", NA) %>%# replace NULL with NA 
    rename_with(~ gsub(".x", "_NC", .x, fixed = TRUE)) %>%
    rename_with(~ gsub(".y", "_DECR", .x, fixed = TRUE)) %>%
    rename_with(~gsub("polo_motif_number_","",.x)) %>%
    pivot_longer(-Proteins.Standard.Name,
                 values_to="number",names_to="group") %>%
    filter(!is.na(number)) %>%
    mutate(number=as.integer(number)) %>%
    group_by(group) %>%
    mutate(group_n=n()) %>% 
    mutate(group_n=paste0(group," (",group_n,")")) %>%
    ggplot(aes(y=number,x=group,fill=group_n))+
    geom_boxplot(alpha=0.5)+
    scale_y_continuous(breaks=seq(0,10,2))+
    coord_flip()+
    theme_classic()+
    labs(x="",y="number",title="spo13 v WT, number of S[ST]P motifs, polo sites only")

# spo13m2 v WT metaphase I
spo13m2_WT_box<-merge(spo13m2_NC_df %>%
          select(Proteins.Standard.Name,polo_motif_number),
      spo13m2_DECR_df %>%
          select(Proteins.Standard.Name,polo_motif_number), 
      by="Proteins.Standard.Name",all=T) %>%
    replace(.=="NULL", NA) %>%# replace NULL with NA 
    rename_with(~ gsub(".x", "_NC", .x, fixed = TRUE)) %>%
    rename_with(~ gsub(".y", "_DECR", .x, fixed = TRUE)) %>%
    rename_with(~gsub("polo_motif_number_","",.x)) %>%
    pivot_longer(-Proteins.Standard.Name,
                 values_to="number",names_to="group") %>%
    filter(!is.na(number)) %>%
    mutate(number=as.integer(number)) %>%
    group_by(group) %>%
    mutate(group_n=n()) %>% 
    mutate(group_n=paste0(group," (",group_n,")")) %>%
    ggplot(aes(y=number,x=group,fill=group_n))+
    geom_boxplot(alpha=0.5)+
    scale_y_continuous(breaks=seq(0,10,2))+
    coord_flip()+
    theme_classic()+
    labs(x="",y="number",title="spo13m2 v WT, number of S[ST]P motifs")

# spo13m2 v WT polo sites only metaphase I
spo13m2_WT_polo_box<-merge(spo13m2_polo_NC_df %>%
          select(Proteins.Standard.Name,polo_motif_number),
      spo13m2_polo_DECR_df %>%
          select(Proteins.Standard.Name,polo_motif_number), 
      by="Proteins.Standard.Name",all=T) %>%
    replace(.=="NULL", NA) %>%# replace NULL with NA 
    rename_with(~ gsub(".x", "_NC", .x, fixed = TRUE)) %>%
    rename_with(~ gsub(".y", "_DECR", .x, fixed = TRUE)) %>%
    rename_with(~gsub("polo_motif_number_","",.x)) %>%
    pivot_longer(-Proteins.Standard.Name,
                 values_to="number",names_to="group") %>%
    filter(!is.na(number)) %>%
    mutate(number=as.integer(number)) %>%
    group_by(group) %>%
    mutate(group_n=n()) %>% 
    mutate(group_n=paste0(group," (",group_n,")")) %>%
    ggplot(aes(y=number,x=group,fill=group_n))+
    geom_boxplot(alpha=0.5)+
    scale_y_continuous(breaks=seq(0,10,2))+
    coord_flip()+
    theme_classic()+
    labs(x="",y="number",title="spo13m2 v WT, number of S[ST]P motifs, polo sites only")

# TC boxplots
TC_spo13_WT_box<-merge(TC_spo13_NC_df %>%
          select(Proteins.Standard.Name,polo_motif_number),
      TC_spo13_DECR_df %>%
          select(Proteins.Standard.Name,polo_motif_number), 
      by="Proteins.Standard.Name",all=T) %>%
    replace(.=="NULL", NA) %>%# replace NULL with NA 
    rename_with(~ gsub(".x", "_NC", .x, fixed = TRUE)) %>%
    rename_with(~ gsub(".y", "_DECR", .x, fixed = TRUE)) %>%
    rename_with(~gsub("polo_motif_number_","",.x)) %>%
    pivot_longer(-Proteins.Standard.Name,
                 values_to="number",names_to="group") %>%
    filter(!is.na(number)) %>%
    mutate(number=as.integer(number)) %>%
    group_by(group) %>%
    mutate(group_n=n()) %>% 
    mutate(group_n=paste0(group," (",group_n,")")) %>%
    ggplot(aes(y=number,x=group,fill=group_n))+
    geom_boxplot(alpha=0.5)+
    scale_y_continuous(breaks=seq(0,10,2))+
    coord_flip()+
    theme_classic()+
    labs(x="",y="number",title="TC spo13 v WT, number of S[ST]P motifs")

# spo13 v WT polo sites only metaphase I
TC_spo13_WT_polo_box<-merge(TC_spo13_polo_NC_df %>%
          select(Proteins.Standard.Name,polo_motif_number),
      TC_spo13_polo_DECR_df %>%
          select(Proteins.Standard.Name,polo_motif_number), 
      by="Proteins.Standard.Name",all=T) %>%
    replace(.=="NULL", NA) %>%# replace NULL with NA 
    rename_with(~ gsub(".x", "_NC", .x, fixed = TRUE)) %>%
    rename_with(~ gsub(".y", "_DECR", .x, fixed = TRUE)) %>%
    rename_with(~gsub("polo_motif_number_","",.x)) %>%
    pivot_longer(-Proteins.Standard.Name,
                 values_to="number",names_to="group") %>%
    filter(!is.na(number)) %>%
    mutate(number=as.integer(number)) %>%
    group_by(group) %>%
    mutate(group_n=n()) %>% 
    mutate(group_n=paste0(group," (",group_n,")")) %>%
    ggplot(aes(y=number,x=group,fill=group_n))+
    geom_boxplot(alpha=0.5)+
    scale_y_continuous(breaks=seq(0,10,2))+
    coord_flip()+
    theme_classic()+
    labs(x="",y="number",title="TC spo13 v WT, number of S[ST]P motifs, polo sites only")

print_plot<-function(plot,filename.pdf,width,height){
    print(plot)
    pdf(filename.pdf, width=width, height=height) 
    print(plot)
    dev.off()
}

print_plot(spo13_WT_box,"spo13_WT_box.pdf",6,4)
print_plot(spo13_WT_polo_box,"spo13_WT_polo_box.pdf",6,4)
print_plot(spo13m2_WT_box,"spo13m2_WT_box.pdf",6,4)
print_plot(spo13m2_WT_polo_box,"spo13m2_WT_polo_box.pdf",6,4)
print_plot(TC_spo13_WT_box,"TC_spo13_WT_box.pdf",6,4)
print_plot(TC_spo13_WT_polo_box,"TC_spo13_WT_polo_box.pdf",6,4)


library(forcats)


# fisher tests
two_fish<-function(df1,df2,set1,set2){
    match1<-df1 %>% filter(polo_box_motif==TRUE) %>% nrow()
    no_match1<-df1 %>% filter(polo_box_motif==FALSE) %>% nrow()
    match2<-df2 %>% filter(polo_box_motif==TRUE) %>% nrow()
    no_match2<-df2 %>% filter(polo_box_motif==FALSE) %>% nrow()
    
    
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

two_fish(spo13_NC_df,spo13_DECR_df)

library(ggsignif)
library(scales)

two_fish_bar<-function(two_fish_df,set1,set2,title){
    out_df<- two_fish_df %>%
        tail(1) %>%
        pivot_longer(-c(pval,stars),values_to="perc",names_to="sample") %>%
        mutate(sample=ifelse(grepl("set1",sample),set1,set2)) %>%
        mutate(sample=factor(sample,levels=c(set1,set2))) #keep levels consistent
    
    
    ratios<-two_fish_df %>%
        rbind(ratio=paste0(.[1,],"/",.[3,])) %>%
        tail(1) %>%
        select(set1,set2) %>%
        pivot_longer(cols=c(set1,set2),names_to="set",values_to="ratio") %>%
        select(ratio)
    
    out_df<-out_df %>%
        cbind(ratios)
    
    
    #return(out_df)
    out_df %>%
        ggplot(aes(x=sample,y=perc,fill=sample))+
        geom_bar(aes(fill=sample),stat="identity")+
        geom_text(aes(label=ratio),vjust=2)+
        geom_signif(comparisons = list(c(set1,set2)),
                    map_signif_level=TRUE,
                    annotations = out_df$stars[1],
                    size=0.25)+
        scale_y_continuous(expand=expansion(mult=c(0,0.2)),labels = label_number(accuracy=0.1))+
        theme_classic()+
        theme(axis.text.x = element_text(angle=90,vjust=0.6),legend.position="none",axis.text=element_text(size=rel(1)))+
        labs(y="percent containing S[ST]P",x="",title=title,axis.line=element_line(linewidth=0.1))
}
    
p1<-two_fish_bar(two_fish(spo13_NC_df,spo13_DECR_df),"NC","DECR","spo13 v WT phospho-proteins")
p2<-two_fish_bar(two_fish(spo13m2_NC_df,spo13m2_DECR_df),"NC","DECR","spo13m2 v WT phospho-proteins")
p3<-two_fish_bar(two_fish(spo13_polo_NC_df,spo13_polo_DECR_df),"NC","DECR","spo13 v WT polo phospho-proteins")
p4<-two_fish_bar(two_fish(spo13m2_polo_NC_df,spo13m2_polo_DECR_df),"NC","DECR","spo13m2 v WT polo phospho-proteins")
p5<-two_fish_bar(two_fish(TC_spo13_NC_df,TC_spo13_DECR_df),"NC","DECR","spo13 v WT TC phospho-proteins")
p6<-two_fish_bar(two_fish(TC_spo13_polo_NC_df,TC_spo13_polo_DECR_df),"NC","DECR","spo13 v WT TC polo phospho-proteins")


library(patchwork)
print_plot((p1+p2+p3)/(p4+p5+p6),"motif_fish_bar.pdf",11,10)