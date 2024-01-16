
# Each time the DEP impute function is run, slightly different values
# are imputed.
# Therefore I run the script once and saved the results as an R data object
# which should be loaded back in to revise any plots

# Loading packages
library("dplyr")
library("DEP")

# 
# # # load the tables
# protg<-read.delim("/Users/lkoch/Lori MS RStudio/2401_Cdc5_IP/proteinGroups.txt")
# 
# phos<-read.delim("/Users/lkoch/Lori MS RStudio/2401_Cdc5_IP/Phospho (STY)Sites.txt")
# 
# # remove rev/con
# 
# phos <- subset(phos, Reverse != "+" & Potential.contaminant != "+")
# protg <- subset(protg, Reverse != "+" & Potential.contaminant != "+")
# 
# # Add Gene.name cols #
# gene.map.table <-read.csv("/Users/lkoch/Lori MS RStudio/Reference dataframes/Yeast_Uniprot2.csv")
# 
# #data prep so there isn't white space and only the first gene name is in the gene column
# 
# strsplit.extract <- function (x, split, index) {
#     x <- as.character(x)
#     x[x == ""] <- ";"
#     x <- sapply(strsplit(x, split, fixed = TRUE), "[[", index)
#     return(x)
# }
# 
# gene.map.table$Gene<-trimws(gene.map.table$Gene, "right")
# gene.map.table$Gene.name <- strsplit.extract(gene.map.table$Gene, ";",1)
# 
# #### Create new ID columns ####
# #Create a Gene.name column in phos table
# 
# phos$Gene.name <- gene.map.table$Gene.name[match(phos$Protein, gene.map.table$OLN)]
# 
# #### Create new ID cols in the protg table ####
# #Create new Majority.protein.ID col with only FIRST (one) majority protein ID
# #for cross referencing prot between tables
# 
# protg$Majority.protein.ID <- strsplit.extract(protg$Majority.protein.IDs, ";",1)
# 
# #Create Gene.name col in the protg table
# protg$Gene.name <- gene.map.table$Gene.name[match(protg$Majority.protein.ID,gene.map.table$OLN)]


# # use DEP function make_unique
# 
# data_unique<-make_unique(protg,"Gene.name","id",delim = ";")
# 
# # needs to be gene.name_site for phos....
# 
# # create gene.name.site col
# # do not include S/T/Y central residue for space reasons on plot
# # labels will be written as MPS1-185 for example
# 
# phos<-phos %>%
#     mutate(gene.name.site=paste0(Gene.name,"_",Positions.within.proteins))
# 
# 
# phos_unique<-make_unique(phos,"gene.name.site","id",delim = ";")
# 
# 
# # Generate a SummarizedExperiment object
# 
# # # protg
# LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
# exp_design<-data.frame(label=c("N1","N2","N3","N4","N5","N6","N7","N8","N9"),
#                        condition=c(rep(c("no_tag","WT","spo13"),3)),
#                        replicate=c(rep("1",3),rep("2",3),rep("3",3)))
# data_se <- make_se(data_unique, LFQ_columns, exp_design)
# 
# #phos
# phos_LFQ_columns <- grep("Intensity.PE.___1", colnames(phos_unique)) # get LFQ column numbers
# PE_exp_design<-data.frame(label=c("PE1___1","PE2___1","PE3___1","PE4___1","PE5___1","PE6___1","PE7___1","PE8___1","PE9___1"),
#                        condition=c(rep(c("no_tag","WT","spo13"),3)),
#                        replicate=c(rep("1",3),rep("2",3),rep("3",3)))
# phos_se <- make_se(phos_unique, phos_LFQ_columns, PE_exp_design)
# #
# # # Plot a barplot of the protein identification overlap between samples
# plot_frequency(data_se)
# 
# plot_frequency(phos_se)
# 
# plot_numbers(data_se)
# 
# plot_numbers(phos_se)
# 
# # # Normalize the data
# data_norm <- normalize_vsn(data_se)
# 
# phos_norm <- normalize_vsn(phos_se)
# #
# # # Visualize normalization by boxplots for all samples before and after normalization
# plot_normalization(data_se, data_norm)
# 
# plot_normalization(phos_se, phos_norm)
# 
# # Plot a heatmap of proteins with missing values
# plot_missval(data_se)
# 
# plot_missval(phos_se)
# # most missval are in no tag so that's good
# 
# # Plot intensity distributions and cumulative fraction of proteins with and without missing values
# plot_detect(data_se)
# 
# plot_detect(phos_se)
# 
# # DO NOT RUN AGAIN
# # Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
# # NOTE THIS ACTUALLY RETURNS DIFFERENT VALUES EVERY TIME IT IS RUN
# # TO MAKE SURE YOU CAN MODIFY INITIAL OUTPUT FROM THIS,
# # SAVE OUTPUT as R object and then load that in, rather than re-running
# # entire script
# data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
# 
# phos_imp <- impute(phos_norm, fun = "MinProb", q = 0.01)
# 
# 
# #Plot intensity distributions before and after imputation
# plot_imputation(data_norm, data_imp)
# 
# plot_imputation(phos_norm, phos_imp)
# 
# 
# # Test all possible comparisons of samples
# data_diff <- test_diff(data_imp, type = "all")
# 
# phos_diff <- test_diff(phos_imp, type = "all")
# #Tested contrasts: no_tag_vs_WT, no_tag_vs_spo13, WT_vs_spo13
# 
# # Denote significant proteins based on user defined cutoffs
# dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))
# 
# phos_dep  <- add_rejections(phos_diff, alpha = 0.05, lfc = log2(1.5))
# 
# # Generate a results table that is a dataframe
# # easier to manipulate in R
# data_results <- get_results(dep)
# 
# phos_results <-get_results(phos_dep)
# 
# # Number of significant proteins
# data_results %>% filter(significant) %>% nrow()#88 significantly varying prots
# 
# phos_results %>% filter(significant) %>% nrow()#58 significantly varying sites
# 
# # run once after initial run
# saveRDS(data_results,"240111_protein.RDS") # save small results table
# saveRDS(phos_results,"240111_phospho.RDS") # save small results table
# saveRDS(dep,"240111_DEP_protein.RDS") # save large results table ## IMPORTANT!! for plotting.
# saveRDS(phos_dep,"240111_DEP_phospho.RDS") # save large results table

# DO NOT RUN AGAIN #####

# LOAD PREVIOUS RESULTS OF THE SCRIPT
data_results<-readRDS("240111_protein.RDS")

phos_results<-readRDS("240111_phospho.RDS")

dep<-readRDS("240111_DEP_protein.RDS")

phos_dep<-readRDS("240111_DEP_phospho.RDS")

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()#50 significantly varying prots

phos_results %>% filter(significant) %>% nrow()#34 significantly varying sites

# Generate a wide data.frame of all of the data
df_wide <- get_df_wide(dep)

phos_df <- get_df_wide(phos_dep)

# volcano plots

# load plotting packages
library(ggplot2)
library(plotly)
library(ggrepel)

# Plot unadjusted p-values 
# label signficance that was called based on adjusted p-values

# create an unlabeled PDF and plotly HTML volcano with unadjusted pvals and cutoffs
get_plotly_volcano_adj <- function(df,contrast,filename.html,filename.pdf) {
    contrast_cols<-grep(contrast,names(df),value=T)
    to_plot_cols<-c("name",contrast_cols)
    #this df has just the necessary columns
    to_plot<-df[,to_plot_cols]
    #this will be x and y of the plot
    x<-to_plot[,grep("diff",names(to_plot),value=T)]
    y<-to_plot[,grep("p.val",names(to_plot),value=T)] #plot UNadjusted p-value
    sig_annot<-to_plot[,grep("significant",names(to_plot),value=T)]
    # labeled as sig based on ADJUSTED p-value 
    
    # add x and y as columns to the subset df
    to_plot<-cbind(to_plot,x,y,sig_annot)
    
    contrast_volcano<-ggplot(to_plot, aes(x = x, y = -log10(y),color=sig_annot,text= paste("log2(fc): ", x, "<br>", 
                                                                                           "p-value: ", y, "<br>",
                                                                                           "protein: ", name, "<br>",
                                                                                           sep = ""))) + 
        geom_point(aes(fill=sig_annot)) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
        theme_classic() +
        theme(legend.position="none")+
        labs(title = contrast) +
        xlab("log2(fold change)") +
        ylab("-log10(p-value)")
    
    #print to the console and save as pdf
    print(contrast_volcano)
    
    pdf(filename.pdf, width=5, height=8) 
    print(contrast_volcano)
    dev.off()
    
    
    # #now make it an interactive html file.
    # 
    contrast_volcano_plotly<-ggplotly(contrast_volcano,tooltip=c("text"))
    htmlwidgets::saveWidget(contrast_volcano_plotly, filename.html)
    
    plots<-list(contrast_volcano,contrast_volcano_plotly)
    return(plots)
}

# 3 PROTEIN comparisons
get_plotly_volcano_adj(df_wide,"no_tag_vs_WT",filename.html = "no_tag_vs_WT_Cdc5_prot.html",filename.pdf="no_tag_vs_WT_Cdc5_prot.pdf")
get_plotly_volcano_adj(df_wide,"no_tag_vs_spo13",filename.html = "no_tag_vs_spo13_Cdc5_prot.html",filename.pdf="no_tag_vs_spo13_Cdc5_prot.pdf")
get_plotly_volcano_adj(df_wide,"WT_vs_spo13",filename.html = "WT_vs_spo13_Cdc5_prot.html",filename.pdf="WT_vs_spo13_Cdc5_prot.pdf")

# 3 PHOSPHO comparisons
get_plotly_volcano_adj(phos_df,"no_tag_vs_WT",filename.html = "no_tag_vs_WT_Cdc5_phos.html",filename.pdf="no_tag_vs_WT_Cdc5_phos.pdf")
get_plotly_volcano_adj(phos_df,"no_tag_vs_spo13",filename.html = "no_tag_vs_spo13_Cdc5_phos.html",filename.pdf="no_tag_vs_spo13_Cdc5_phos.pdf")
get_plotly_volcano_adj(phos_df,"WT_vs_spo13",filename.html = "WT_vs_spo13_Cdc5_phos.html",filename.pdf="WT_vs_spo13_Cdc5_phos.pdf")


# create a LABELED PDF of the volcano 
# where proteins or phospho called "SIGNFICANT" are labeled with the gene name
# also "SELECT" proteins are also labeled, even if they are not significant

get_volcano_labeled_select <- function(df,contrast,select,filename.pdf) {
    contrast_cols<-grep(contrast,names(df),value=T)
    to_plot_cols<-c("name",contrast_cols)
    #this df has just the necessary columns
    to_plot<-df[,to_plot_cols]
    #this will be x and y of the plot
    x<-to_plot[,grep("diff",names(to_plot),value=T)]
    y<-to_plot[,grep("p.val",names(to_plot),value=T)] #unadjusted p-value
    sig_annot<-to_plot[,grep("significant",names(to_plot),value=T)] 
    
    # add x and y as columns to the subset df
    to_plot<-cbind(to_plot,x,y,sig_annot)
    
    to_plot<-to_plot %>%
        mutate(sig_label=ifelse(sig_annot,name,NA)) %>%
        mutate(sel_label=ifelse(grepl(select,name),name,NA)) %>%
        mutate(sel_point=ifelse(name %in% select,TRUE,FALSE)) %>%
        mutate(sel_point_label=ifelse(sel_point,name,NA))
    
    contrast_volcano<-ggplot(to_plot, aes(x = x, y = -log10(y),label=sig_label,color=sig_annot)) + 
        geom_point(aes(color=sig_annot)) +
        geom_text_repel(color="black",max.overlaps = 50)+
        geom_text_repel(data=to_plot %>% filter(sel_point),aes(label=sel_point_label),color="black")+
        geom_point(data=to_plot %>% filter(sel_point),color="black")+
        scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
        theme_classic() +
        theme(legend.position="none")+
        labs(title = contrast) +
        xlab("log2(fold change)") +
        ylab("-log10(p-value)")
    
    #print to the console and save as pdf
    print(contrast_volcano)
    
    pdf(filename.pdf, width=5, height=8) 
    print(contrast_volcano)
    dev.off()
    
    return(contrast_volcano)
}

# PROTEIN

# no tag v WT
get_volcano_labeled_select(df_wide,"no_tag_vs_WT",c("CLB1","HRR25"),filename.pdf="no_tag_vs_WT_Cdc5_protg_label_CLB1.pdf")

# WT v spo13
get_volcano_labeled_select(df_wide,"WT_vs_spo13",c("CLB1","HRR25"),filename.pdf="WT_vs_spo13_Cdc5_protg_label_CLB1.pdf")

# WT v spo13 also label APC and DDK 
get_volcano_labeled_select(df_wide,"WT_vs_spo13",c("CLB1","HRR25","APC1","DBF4","CDC7"),filename.pdf="WT_vs_spo13_Cdc5_protg_label_CLB1_APC.pdf")

# PHOSPHO
# what are all the Clb1 sites by name
#phos_df %>% filter(grepl("CLB1",name)) %>% select(name) 
# CLB1_84, CLB1_24, CLB1_137, CLB1_15

# CLB1_137 and CLB1_15 are already labeled signficant
get_volcano_labeled_select(phos_df,"WT_vs_spo13",c("CLB1_84", "CLB1_24"),filename.pdf="WT_vs_spo13_Cdc5_phos_label_CLB1.pdf")


# PLOT BARPLOTS OF CLB1 DATA
plot_single(dep,"CLB1",type="centered") # saved as pdf portrait 4x6

# note the panels don't come out in right order because of mixed char/numbers
# re-arrange them in Illustrator
plot_single(phos_dep,c("CLB1_15","CLB1_24","CLB1_84","CLB1_137"),type="centered") #saved as pdf portrait 5x8


