# > sessionInfo()
# R version 4.2.3 (2023-03-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] umap_0.2.10.0               org.Mm.eg.db_3.16.0         AnnotationDbi_1.60.2        ComplexHeatmap_2.14.0      
# [5] ggsci_3.0.3                 ggpubr_0.6.0                gridExtra_2.3               ggplot2_3.5.0              
# [9] dplyr_1.1.4                 stringr_1.5.1               edgeR_3.40.2                DESeq2_1.38.3              
# [13] SummarizedExperiment_1.28.0 Biobase_2.58.0              MatrixGenerics_1.10.0       matrixStats_1.2.0          
# [17] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2           
# [21] BiocGenerics_0.44.0         limma_3.54.2                circlize_0.4.16             clusterProfiler_4.6.2      
# 
# loaded via a namespace (and not attached):
#   [1] shadowtext_0.1.3       backports_1.4.1        fastmatch_1.1-4        plyr_1.8.9             igraph_1.5.1          
# [6] lazyeval_0.2.2         splines_4.2.3          BiocParallel_1.32.6    digest_0.6.35          foreach_1.5.2         
# [11] yulab.utils_0.1.4      GOSemSim_2.24.0        viridis_0.6.5          GO.db_3.16.0           fansi_1.0.6           
# [16] magrittr_2.0.3         memoise_2.0.1          cluster_2.1.6          doParallel_1.0.17      Biostrings_2.66.0     
# [21] annotate_1.76.0        graphlayouts_1.1.1     askpass_1.2.0          enrichplot_1.18.4      colorspace_2.1-0      
# [26] blob_1.2.4             ggrepel_0.9.5          crayon_1.5.2           RCurl_1.98-1.14        jsonlite_1.8.8        
# [31] scatterpie_0.2.2       iterators_1.0.14       ape_5.7-1              glue_1.7.0             polyclip_1.10-6       
# [36] gtable_0.3.4           zlibbioc_1.44.0        XVector_0.38.0         GetoptLong_1.0.5       DelayedArray_0.24.0   
# [41] car_3.1-2              shape_1.4.6.1          abind_1.4-5            scales_1.3.0           DOSE_3.24.2           
# [46] DBI_1.2.2              rstatix_0.7.2          Rcpp_1.0.12            viridisLite_0.4.2      xtable_1.8-4          
# [51] clue_0.3-65            reticulate_1.35.0      gridGraphics_0.5-1     tidytree_0.4.6         bit_4.0.5             
# [56] httr_1.4.7             fgsea_1.24.0           RColorBrewer_1.1-3     pkgconfig_2.0.3        XML_3.99-0.16.1       
# [61] farver_2.1.1           locfit_1.5-9.9         utf8_1.2.4             ggplotify_0.1.2        tidyselect_1.2.1      
# [66] rlang_1.1.3            reshape2_1.4.4         munsell_0.5.1          tools_4.2.3            cachem_1.0.8          
# [71] downloader_0.4         cli_3.6.2              generics_0.1.3         RSQLite_2.3.6          gson_0.1.0            
# [76] broom_1.0.5            fastmap_1.1.1          ggtree_3.6.2           bit64_4.0.5            fs_1.6.3              
# [81] tidygraph_1.3.0        purrr_1.0.2            KEGGREST_1.38.0        ggraph_2.2.1           nlme_3.1-164          
# [86] aplot_0.2.2            compiler_4.2.3         rstudioapi_0.16.0      png_0.1-8              ggsignif_0.6.4        
# [91] treeio_1.22.0          tibble_3.2.1           tweenr_2.0.3           geneplotter_1.76.0     stringi_1.8.3         
# [96] RSpectra_0.16-1        lattice_0.22-6         Matrix_1.6-4           vctrs_0.6.5            pillar_1.9.0          
# [101] lifecycle_1.0.4        GlobalOptions_0.1.2    data.table_1.15.4      cowplot_1.1.3          bitops_1.0-7          
# [106] patchwork_1.2.0        qvalue_2.30.0          R6_2.5.1               codetools_0.2-20       MASS_7.3-60.0.1       
# [111] openssl_2.1.1          rjson_0.2.21           withr_3.0.0            GenomeInfoDbData_1.2.9 parallel_4.2.3        
# [116] ggfun_0.1.4            tidyr_1.3.1            HDO.db_0.99.1          carData_3.0-5          ggforce_0.4.2   

##
library(clusterProfiler)
library(circlize)
library(limma)
library(DESeq2)
library(edgeR)
library(stringr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(org.Mm.eg.db)
library(umap)

#####create dataframe#####
#take direct output from nf-core pipeline. Unscaled.
counts <- read.table("~/Desktop/data_jy/salmon.merged.gene_counts.jy.tsv", header=T, sep="\t", row.names = 1, check.name=FALSE,
                     stringsAsFactor=FALSE)
counts <- counts[2:26]
construct <- c("A","A","B","B","C","C","C","D","D","D","E","E","F","F","G","G","H","H","I","I","I","J","J","K","K")
replicate <- c("1","2","1","2","1","2","3","1","2","3","1","2","1","2","1","2","1","2","1","2","3","1","2","1","2")
colnames(counts) <- c(paste0(construct,replicate))
annot <- cbind(construct, replicate)
annot <- as.data.frame(annot)
rownames(annot) <- c(paste0(construct,replicate))

####
colnames(counts) <- c(rownames(annot))

y <- DGEList(counts=counts,samples=annot)
samplenames <- colnames(y) 
samplenames

construct <- as.factor(y$samples$construct)
#time <- as.factor(y$samples$time)
replicate <- as.factor(y$samples$replicate)

y$geneid <- rownames(y)
lib.size <- as.factor(y$samples$lib.size)

# Subset samples:
x <- y[,,keep.lib.sizes=TRUE] ##subset by what?
dim(x)

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

table(rowSums(x$counts==0)==25)

keep.exprs <- rowSums(cpm>1)>=6 ### changed from at least 1 cpm in 2 samples to
x <- x[keep.exprs,, keep.lib.sizes=FALSE] # (=F) recounts lib size after trim

dim(x)

## Plot Log-CPM RAW and Filtered data
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(8, "Paired") # find pallete to accommodate >12 samples if necessary
par(mfrow=c(1,2))

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")



lcpm_x <- cpm(x, log=TRUE)
x_norm <- calcNormFactors(x, method = "TMM")
x_norm$samples$norm.factors

lcpm_x_norm <- cpm(x_norm, log=TRUE)
x_norm <- x
x_norm$samples$norm.factors <- 1

par(mfrow=c(1,2))
lcpm <- cpm(x_norm, log=TRUE)

## Boxplots of Data pre and post-normalization
boxplot(lcpm_x, las=2, col=col, main="")
title(main="Unnormalised data",ylab="Log-cpm")
x_norm <- calcNormFactors(x_norm)
x_norm$samples$norm.factors

lcpm <- cpm(x, log=TRUE)
boxplot(lcpm_x_norm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")

df_lcpm <- as.data.frame(lcpm)
df_lcpm2 <- df_lcpm

batch <- c("A","A","A","A","C","C","C","C","C","C","A","A","A","A","A","A","A","A","A","A","A","A","A","A", "A")
df_lcpmbatch <- as.data.frame(removeBatchEffect(df_lcpm, batch=batch))

#colnames(df_lcpmbatch) <- colnames(counts[2:26])
#write.table(df_lcpmbatch, "~/Desktop/Geneexpression_countmatrix.csv", sep =',', quote = FALSE) 

#####
#####Perform UMAP####

counts.umap <- umap(t(df_lcpmbatch), n_neighbors = 5)
layoutdf <- as.data.frame(counts.umap$layout)
colnames(layoutdf) <- c("one", "two")
layoutdf$samplename <- construct
layoutdf$samplename <- c("ESC 290KO", "ESC 290KO", "ESC 302KO", "ESC 302KO", "ESC DGCR8KO", "ESC DGCR8KO", "ESC DGCR8KO",
                         "ESC WT", "ESC WT", "ESC WT","iPS 290/302KO C1", "iPS 290/302KO C1", "iPS 290/302KO C2", "iPS 290/302KO C2",
                         "iPS 290KO C1", "iPS 290KO C1", "iPS 290KO C2", "iPS 290KO C2", "iPS WT", "iPS WT", "iPS WT",
                         "iPS 302KO C1", "iPS 302KO C1", "iPS 302KO C2", "iPS 302KO C2")

# layoutdf$samplename <- factor(layoutdf$samplename, 
#                                        levels = c("N WT", "N MLL3KO", "N DKO", "N dCD", "F WT", "F MLL3KO", "F DKO", "F dCD" ))  
#layoutdf$samplename <- factor(layoutdf$samplename, 
                            #  levels = c("N WT", "N MLL3KO", "N DKO", "F WT", "F MLL3KO", "F DKO"))  
ggplot(layoutdf, mapping = aes(x=one , y=two,color = samplename)) +
  geom_point(size=2) +
  #scale_fill_manual(values = c("N WT" = "deepskyblue3", "F WT" = "magenta3", "N DKO" = "deepskyblue3", "F DKO" = "magenta3", "N MLL3KO" = "deepskyblue3", "F MLL3KO" = "magenta3")) +
  #  scale_fill_manual(values = c("N WT" = "deepskyblue3", "F WT" = "magenta3", "N DKO" = "deepskyblue3", "F DKO" = "magenta3", "N MLL3KO" = "deepskyblue3", "F MLL3KO" = "magenta3", "N dCD" = "deepskyblue3", "F dCD" = "magenta3" )) +
  #scale_shape_manual(values = c(21, 22, 23, 24, 21, 22, 23, 24)) +
  #scale_shape_manual(values = c(21, 22, 23, 21, 22, 23)) +
  # #geom_text_repel(aes(label = samplename),
  #   box.padding = 0.1, 
  #   point.padding = 0.1,
  theme_bw() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(aspect.ratio = 1/1.5, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2))
# xlim(c(-10.5, 10.5)) +
# ylim(c(-3.5,3.5)) 


#####
#####Marker heatmap#####

df_lcpmbatch <- as.data.frame(removeBatchEffect(df_lcpm, batch=batch))
counts <- read.table("~/Desktop/data_jy/salmon.merged.gene_counts.jy.tsv", header=T, sep="\t", row.names = 1, check.name=FALSE,
                     stringsAsFactor=FALSE)
counts <- counts[2:26]
counts <- cbind(counts[1:4], counts[11:25])
df_lcpmprimarydata <- cbind(df_lcpm[1:4], df_lcpm[11:25])
colnames(df_lcpmprimarydata) <- colnames(counts)
df_lcpmprimarydata$gene.names <- rownames(df_lcpmprimarydata)

df_lcpmbatch<- df_lcpmprimarydata
esc290ave <- (df_lcpmbatch$esc290_REP1 + df_lcpmbatch$esc290_REP1)/2
esc302ave <- (df_lcpmbatch$esc302_REP1 + df_lcpmbatch$esc302_REP1)/2
ipswtave <- (df_lcpmbatch$ips_wt_REP1 + df_lcpmbatch$ips_wt_REP2 + df_lcpmbatch$ips_wt_REP3)/3
ips290b3ave <- (df_lcpmbatch$ips_b3_REP1 + df_lcpmbatch$ips_b3_REP2)/2
ips290b10ave <- (df_lcpmbatch$ips_b10_REP1 + df_lcpmbatch$ips_b10_REP2)/2
ips302c1ave <- (df_lcpmbatch$ips302_c1_REP1 + df_lcpmbatch$ips302_c1_REP2)/2
ips302c2ave <- (df_lcpmbatch$ips302_c2_REP1 + df_lcpmbatch$ips302_c2_REP2)/2
ips290302a6ave <- (df_lcpmbatch$ips_a6_REP1 + df_lcpmbatch$ips_a6_REP2)/2
ips290302a9ave <- (df_lcpmbatch$ips_a9_REP1 + df_lcpmbatch$ips_a9_REP2)/2

df_ave <- as.data.frame(cbind(esc290ave, esc302ave, ipswtave, ips290b3ave, ips290b10ave,
  ips302c1ave, ips302c2ave, ips290302a6ave, ips290302a9ave))
rownames(df_ave) <- df_lcpmbatch$gene.names
df_ave$gene.names <- df_lcpmbatch$gene.names

genesofinterest <- c("Pou5f1", "Sox2", "Nanog", "Klf4", "Esrrb", "Zfp42", "Otx2", 
                    "Dnmt3b", "Grhl2", "Hand1", "T", "Gata4", "Pax6", "Sox1")

df_lcpmbatchselect <- df_ave %>% filter(gene.names %in% genesofinterest)

col_fun1 = colorRamp2(c(0, 15), c("white","red"))
ha <- rowAnnotation(labels = anno_text(rownames(df_lcpmbatchselect), which = "row"), 
                    width = max(grobWidth(textGrob(rownames(df_lcpmbatchselect)))))

ht1 <- Heatmap(df_lcpmbatchselect[1:9], col = col_fun1,
               show_row_names = TRUE,
               row_order = genesofinterest,
               row_names_gp = gpar(col = "black", fontsize = 10),
               border_gp = gpar(col = "black", lty = 1))
               #left_annotation = ha
               #width = ncol(df_lcpmbatchselect[1:25])*unit(10, "mm"))

ht1


xx <- df_ave[1:9] - ipswtave
xx = subset(xx, select = -c(ipswtave))
xx$gene.names <- rownames(xx)
df_lcpmbatchselectnorm <- xx %>% filter(gene.names %in% genesofinterest)

col_fun2 = colorRamp2(c(-6, 0, 6), c("blue", "white","red"))
ha <- rowAnnotation(labels = anno_text(rownames(df_lcpmbatchselectnorm), which = "row"), 
                    width = max(grobWidth(textGrob(rownames(df_lcpmbatchselectnorm)))))

ht2 <- Heatmap(df_lcpmbatchselectnorm[1:8], col = col_fun2,
               show_row_names = TRUE,
               row_names_gp = gpar(col = "black", fontsize = 10),
               border_gp = gpar(col = "black", lty = 1))

ht1 + ht2

#####
#####DESeq2: differential analysis#####

# read data
counts <- read.table("~/Desktop/data_jy/salmon.merged.gene_counts.jy.tsv", header=T, sep="\t", row.names = 1, check.name=FALSE,
                     stringsAsFactor=FALSE)
counts <- counts[2:26]

ips290302 <- counts[ ,11:14]
ips290 <- counts[ ,15:18]
ipswt <- counts[ ,19:21]
ips302 <- counts[ ,22:25]

# perform desired comparison
counts <- data.frame(ipswt, ips302)
counts2 <- mutate_all(counts, function(x) as.integer(as.numeric(x)))
counts <- counts2
output <- "~/Desktop/data_jy/Diffgenes290302vsWTfull-1.07.24.anno.csv"
#output <- "~/Desktop/data_jy/Diffgenes/DiffgenesFCKOvsFDKO-11.20.22.anno.csv"

construct <- c("A", "A", "A", "B", "B", "B", "B")
replicate <- c("1", "2", "3","1", "2", "3", "4" )
annot <- cbind(construct, replicate)
annot <- as.data.frame(annot)
# rownames(annot) <- c("A1","A2","A3","B1","B2","B3")

colnames(counts) <- c(rownames(annot))

library('DESeq2')
# counts matrix: genes as rows and samples as columns; annotation samples as rows and info in columns
DESeqOurcts <- as.matrix(counts)
DESeqOurcoldata <- as.matrix(annot)

all(rownames(DESeqOurcoldata) %in% colnames(DESeqOurcts))
all(rownames(DESeqOurcoldata) == colnames(DESeqOurcts))

##This is where the DESeq2 vignette starts
# ensure 'control level' is the 1st level of annotated input
DESeqOurdds <- DESeqDataSetFromMatrix(countData = DESeqOurcts,
                                      colData = DESeqOurcoldata,
                                      design = ~ construct)
DESeqOurdds

# add additional metadata
DESeqOurfeatureData <- data.frame(gene=rownames(DESeqOurcts))
mcols(DESeqOurdds) <- DataFrame(mcols(DESeqOurdds), DESeqOurfeatureData)
mcols(DESeqOurdds)

# prefiltering
DESeqOurdds <- DESeqOurdds[ rowSums(counts(DESeqOurdds)) > 1, ]

# explicitly set reference (default is 1st by alphabetic); not necessary if using contrasts
##DESeqOurdds$condition <- relevel(DESeqOurdds$group, ref="e75") ##obvs 7.5 is not in this but what would be good ref?

### Differential Expression
OurddsMF <- DESeqOurdds

design(OurddsMF) <- formula(~ construct)
OurddsMF <- DESeq(OurddsMF)
OurresMF <- results(OurddsMF)
head(OurresMF)


# MA plot
plotMA(OurresMF, alpha = 0.05, ylim=c(-8,8))

resMFOrdered <- OurresMF[order(OurresMF$padj),] # reordered by padj
summary(resMFOrdered)
sum(OurresMF$padj < 0.05, na.rm=TRUE) # count sig genes with padj < 0.05

#filter reads for padj and log2FC and generate .csv file named by output above

resMFOrderedfiltered <-  subset(resMFOrdered, padj < 0.05 & (log2FoldChange < -1 | log2FoldChange > 1))
resMFOrderedfiltered <-  subset(resMFOrdered, padj < 0.05)

write.csv(resMFOrderedfiltered, output)
sum(OurresMF$padj < 0.05, na.rm=TRUE) # count sig genes with padj < 0.05
# converted <- convertensemblgenetosymbol(resMFOrderedfiltered, rownames(resMFOrderedfiltered))
# write.csv(converted, output)

# sigup <- as.data.frame(subset(resMFOrdered, padj < 0.05 & (log2FoldChange > 1)))
# sigdown  <- as.data.frame(subset(resMFOrdered, padj < 0.05 & (log2FoldChange < -1)))

#####
#####Clusterprofiler on differential genes#####

species_database <- "org.Mm.eg.db"
library(species_database, character.only = TRUE)

raw_df1 <- read.csv("~/Desktop/data_jy/Diffgenes290vsWTfull-1.07.24.anno.csv")
raw_df1$comparison <- "290KO"
raw_df2 <- read.csv("~/Desktop/data_jy/Diffgenes302vsWTfull-1.07.24.anno.csv")
raw_df2$comparison <- "302KO"
raw_df3 <- read.csv("~/Desktop/data_jy/Diffgenes290302vsWTfull-1.07.24.anno.csv")
raw_df3$comparison <- "290302KO"

input_deseq_array <- rbind(raw_df1, raw_df2, raw_df3)
#comparison_group <- c("nwt_fwt", "fwt_fdko", "nwt_ndko", "nwt_ndcd", "fwt_fdcd","nwt_ncko","fwt_fcko")
#comparison_group <- c("nwt_fwt", "nwt_ndko", "fwt_fdko")
comparison_group <- c("290KO", "302KO", "290302KO")
comparison_group <- c("290302KO")
comparison_group <- c("290KO")


setwd("~")
plotarray_loss = list()
plotarray_gain = list()
counter=0
for(i in 1:length(comparison_group)) {
  counter = counter + 1
  input_df <- comparison_group[[i]]
  
  deseq_loss <- input_deseq_array %>% filter(comparison == comparison_group[[i]] & log2FoldChange < -1)
  deseq_gain <- input_deseq_array %>% filter(comparison == comparison_group[[i]] & log2FoldChange > 1)
  
  gene_vector <- as.character(deseq_loss$X)
  gene_vector2 <- bitr(gene_vector, fromType = "SYMBOL",
                       toType = c("ENSEMBL", "ENTREZID"),
                       OrgDb = species_database)
  gene_vector3 <- as.character(gene_vector2$ENTREZID)
  
  ego_bp <- enrichGO(gene          = gene_vector3,
                     OrgDb         = species_database,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     #pvalueCutoff  = 0.01,
                     qvalueCutoff  = 1,
                     readable      = TRUE)
  
  p1 <- dotplot(ego_bp, showCategory = 10, color = "qvalue") + ggplot2::theme_bw() + 
    scale_y_discrete(labels = scales::wrap_format(40)) + 
    ggplot2::theme(aspect.ratio =1.5, 
                   text = element_text(size = 18),
                   panel.border = element_rect(fill=NA, colour = "black", size=3))
  plotarray_loss[[counter]] <- print(p1)
  #ggsave(paste0(comparison_group[[i]], "_loss.go_analysis.svg"))
  
  gene_vector <- as.character(deseq_gain$X)
  gene_vector2 <- bitr(gene_vector, fromType = "SYMBOL",
                       toType = c("ENSEMBL", "ENTREZID"),
                       OrgDb = species_database)
  gene_vector3 <- as.character(gene_vector2$ENTREZID)
  
  ego_bp <- enrichGO(gene          = gene_vector3,
                     OrgDb         = species_database,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     #pvalueCutoff  = 0.01,
                     qvalueCutoff  = 1,
                     readable      = TRUE)
  
  p2 <- dotplot(ego_bp, showCategory = 10, color = "qvalue") + 
    ggplot2::theme_bw() +
    scale_y_discrete(labels = scales::wrap_format(40)) + 
    ggplot2::theme(aspect.ratio =1.5, 
                   text = element_text(size = 18),
                   panel.border = element_rect(fill=NA, colour = "black", size=3))#+ ggplot2::theme(aspect.ratio =4)
  
  plotarray_gain[[counter]] <- print(p2)
}

plotarray_loss[[1]]
plotarray_gain[[1]]

#####
#####CDF plots#####
counts <- read.table("~/Desktop/data_jy/salmon.merged.gene_counts.jy.tsv", header=T, sep="\t", row.names = 1, check.name=FALSE,
                     stringsAsFactor=FALSE)
counts <- counts[2:26]
colnames(df_lcpm2) <- colnames(counts)
df_lcpm2$gene.names <- rownames(df_lcpm2)


esc290ave <- (df_lcpm2$esc290_REP1 + df_lcpm2$esc290_REP1)/2
esc302ave <- (df_lcpm2$esc302_REP1 + df_lcpm2$esc302_REP1)/2
ipswtave <- (df_lcpm2$ips_wt_REP1 + df_lcpm2$ips_wt_REP2 + df_lcpm2$ips_wt_REP3)/3
ips290b3ave <- (df_lcpm2$ips_b3_REP1 + df_lcpm2$ips_b3_REP2)/2
ips290b10ave <- (df_lcpm2$ips_b10_REP1 + df_lcpm2$ips_b10_REP2)/2
ips302c1ave <- (df_lcpm2$ips302_c1_REP1 + df_lcpm2$ips302_c1_REP2)/2
ips302c2ave <- (df_lcpm2$ips302_c2_REP1 + df_lcpm2$ips302_c2_REP2)/2
ips290302a6ave <- (df_lcpm2$ips_a6_REP1 + df_lcpm2$ips_a6_REP2)/2
ips290302a9ave <- (df_lcpm2$ips_a9_REP1 + df_lcpm2$ips_a9_REP2)/2

df_ave <- as.data.frame(cbind(esc290ave, esc302ave, ipswtave, ips290b3ave, ips290b10ave,
                              ips302c1ave, ips302c2ave, ips290302a6ave, ips290302a9ave))
rownames(df_ave) <- df_lcpm2$gene.names
df_ave$gene.names <- df_lcpm2$gene.names
df_lcpm2$gene.name <- rownames(df_lcpm2)

df_avefold <- df_ave[1:9] - ipswtave
df_avefold$gene.names <- rownames(df_avefold)
xxtarget <- df_avefold %>% filter(gene.names %in% targetscan$Target.gene)
xxother <- df_avefold %>% filter(!gene.names %in% targetscan$Target.gene)

comparison_group <- c("esc290ave", "esc302ave", "ips290b3ave", "ips290b10ave", 
                      "ips302c1ave", "ips302c2ave", "ips290302a6ave", "ips290302a9ave" )
setwd("~")
plotlist=list()
combinedplotlist=list()
counter=0
for(i in 1:length(comparison_group)) {
  counter = counter + 1
  coi <- comparison_group[[i]]
  xxtarget$kind <- "290target"
  xxxtarget <- xxtarget %>% pivot_longer(cols = coi)
  xxother$kind <- "allother"
  xxxother <- xxother %>% pivot_longer(cols = coi)
  
  zz <- rbind(xxxtarget, xxxother)
  
  plotlist[[paste0("a", counter)]] <- print(ggplot(zz, aes(value, col=kind)) +
                                                  stat_ecdf() +
                                                  xlim(c(-2.5,2.5)) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none"))
  combinedplotlist[[counter]] <- plotlist[[counter]]}
wrap_plots(combinedplotlist)
       
#####