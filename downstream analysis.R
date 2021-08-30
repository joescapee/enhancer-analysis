#set working directory
setwd("~/dev/repos/Rothstein_2019/downsteam_R_analysis")

#read bedtools files for HH6 & HH9
bedtools_output_HH6<-read.table("../analysis/bedtools_output_HH6",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
bedtools_output_HH9<-read.table("../analysis/bedtools_output_HH9",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")

#read homer files for HH6 & HH9
homer_output_HH6<-read.table("../analysis/homer_output_HH6.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
homer_output_HH9<-read.table("../analysis/homer_output_HH9.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")


# generate annotations dataframe
annotations_list<-homer_output_HH6[,c('Entrez.ID', 'Gene.Name')]

# remove duplicated gene ID rows
annotations_list_clean<-annotations_list[!duplicated(annotations_list$Entrez.ID),]


#renaming column names for homer output
colnames(homer_output_HH6)[1] <- "peakname"
colnames(homer_output_HH9)[1] <- "peakname"


#############################################################################

# frequency plots

#install.packages("ggplot2")
#install.packages("extrafont")
#font_import()

library(ggplot2)
library(extrafont)
windowsFonts()
loadfonts(device = "win")


output_path = 'R output'


# extract and clean up annotations category for both stages

annotation_peaks_HH6<-as.factor(sub("[[:space:]].*","",homer_output_HH6[,"Annotation"]))
annotation_peaks_HH9<-as.factor(sub("[[:space:]].*","",homer_output_HH9[,"Annotation"]))

# tables of absolute and relative frequencies
homer_output_HH6$Annotation <- sub("\\(.*", "", homer_output_HH6$Annotation)
homer_output_HH9$Annotation <- sub("\\(.*", "", homer_output_HH9$Annotation)

table_HH6 <- as.data.frame(table(homer_output_HH6$Annotation))
table_HH9 <- as.data.frame(table(homer_output_HH9$Annotation))

relative_table_HH6 <- as.data.frame((table(homer_output_HH6$Annotation)/length(homer_output_HH6$Annotation)))
relative_table_HH9 <- as.data.frame((table(homer_output_HH9$Annotation)/length(homer_output_HH9$Annotation)))

# extract and simplify annotations for categorisation

annotation_peaks_HH6 <- as.factor(sub('.*', "",homer_output_HH6[,"Annotation"]))
annotation_peaks_HH9 <- as.factor(sub('.*', "",homer_output_HH9[,"Annotation"]))

# order frequency

freq_data_HH6 <- as.data.frame(prop.table(table(annotation_peaks_HH6))[order(prop.table(table(annotation_peaks_HH6)))])
colnames(freq_data_HH6) = c('peaks', 'Frequency')

freq_data_HH9 <- as.data.frame(prop.table(table(annotation_peaks_HH9))[order(prop.table(table(annotation_peaks_HH9)))])
colnames(freq_data_HH9) = c('peaks', 'Frequency')

# frequency plot of peak annotations HH6

png(paste0(output_path, "peak_annotation_frequency_HH6.png"), height = 10, width = 10, family = 'Arial', units = 'cm', res = 400)
ggplot(freq_data_HH6, aes(x = peaks, y = Frequency*100)) +
  geom_bar(stat='identity', fill='steelblue') +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5))
graphics.off()

# frequency plot of peak annotations HH9

png(paste0(output_path, "peak_annotation_frequency_HH9.png"), height = 10, width = 10, family = 'Arial', units = 'cm', res = 400)
ggplot(freq_data_HH9, aes(x = peaks, y = Frequency)) +
  geom_bar(stat='identity', fill='steelblue') +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5))
graphics.off()


############################################################################

# filter for intergenic and introns (removing TSS, exon and promoter-TSS) for later analysis

library(dplyr)
library(stringr)
library(utils)


### using grep 

# HH6 & HH9

enhancer_HH6<-homer_output_HH6[grep('Intergenic|intron',homer_output_HH6$Annotation),]
enhancer_HH9<-homer_output_HH9[grep('Intergenic|intron',homer_output_HH9$Annotation),]

### using filter

# HH6 & HH9

check_enhancer_HH6<-homer_output_HH6 %>% filter(str_detect(str_to_lower(Annotation), "int"))
check_enhancer_HH9<-homer_output_HH9 %>% filter(str_detect(str_to_lower(Annotation), "int"))

#check that no TTS, exons or p-TSS are included
#why doesn't the wildcard-argument ("intron*) work?!?!

any('Intergenic'%in%enhancer_HH6$Annotation)
any("intron.*"%in%enhancer_HH6$Annotation, dotall = TRUE)
any('TTS'%in%enhancer_HH6$Annotation)
any('exon'%in%enhancer_HH6$Annotation)
any('promoter-TSS'%in%enhancer_HH6$Annotation)


############################################################################

# plotting putative enhancer profiles and heatmaps for HH6/9

setRepositories() 

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install(version='devel')
#BiocManager::install("S4Vectors")

#install.packages("ChIPpeakAnno")
library("ChIPpeakAnno")
library("rtracklayer")
library("extrafont")
library("Rsamtools")

# troubleshoot "package is not available"
# <- available.packages()
#view(ap) 
#"ChIPpeakAnno" %in% rownames(ap)


#output_path = "/R output/"


##### using filtered homer file 

# peaks HH6

peaks_HH6 <- GRanges(seqnames=enhancer_HH6[,2],
                 ranges=IRanges(start=enhancer_HH6[,3],
                                end=enhancer_HH6[,4],
                                names=enhancer_HH6[,1]))

# peaks HH9

peaks_HH9 <- GRanges(seqnames=enhancer_HH9[,2],
                     ranges=IRanges(start=enhancer_HH9[,3],
                                    end=enhancer_HH9[,4],
                                    names=enhancer_HH9[,1]))

##### using original homer file 

# peaks HH6

#peaks_HH6 <- GRanges(seqnames=homer_output_HH6[,2],
#ranges=IRanges(start=homer_output_HH6[,3],
# end=homer_output_HH6[,4],
# names=homer_output_HH6[,1]))


# peaks HH9

#peaks_HH9 <- GRanges(seqnames=homer_output_HH9[,2],
#ranges=IRanges(start=homer_output_HH9[,3],
#end=homer_output_HH9[,4],
# names=homer_output_HH9[,1]))


# find centre of putative enhancers and get coordinates for +/-2kb (HH6)

peaks.recentered_HH6 <- peaks.center_HH6 <- peaks_HH6
start(peaks.center_HH6) <- start(peaks_HH6) + floor(width(peaks_HH6)/2)
width(peaks.center_HH6) <- 1
start(peaks.recentered_HH6) <- start(peaks.center_HH6) - 2000
end(peaks.recentered_HH6) <- end(peaks.center_HH6) + 2000



# find centre of putative enhancers and get coordinates for +/-2kb (HH9)

peaks.recentered_HH9 <- peaks.center_HH9 <- peaks_HH9
start(peaks.center_HH9) <- start(peaks_HH9) + floor(width(peaks_HH9)/2)
width(peaks.center_HH9) <- 1
start(peaks.recentered_HH9) <- start(peaks.center_HH9) - 2000
end(peaks.recentered_HH9) <- end(peaks.center_HH9) + 2000


### IMPORT of bigWigs ###


# import bigwig files and select regions corresponding to putative enhancer peaks (HH6 & HH9)

bigwig_files_ATAC <- list.files("./bigwigs/atac_bigwig", pattern = "bigWig", full.names = T)
bigwig_files_CUTandRUN <- list.files("./bigwigs/cutandrun_igv_bigwig", pattern = "bigWig", full.names = T)


# import of ATAC bigWig files for HH6

ATAC.bw_HH6_1 <- import(bigwig_files_ATAC[grepl("ATAC_HH6_R1.mLb.clN", bigwig_files_ATAC, ignore.case = T)], format="BigWig", which=peaks.recentered_HH6, as="RleList")
#ATAC.bw_HH6_2 <- import(bigwig_files_ATAC[grepl("ATAC_HH6_R2.mLb.clN", bigwig_files_ATAC, ignore.case = T)], format="BigWig", which=peaks.recentered_HH6, as="RleList")

# import of ATAC bigWig files for HH9

ATAC.bw_HH9_1 <- import(bigwig_files_ATAC[grepl("ATAC_HH9_R1.mLb.clN", bigwig_files_ATAC, ignore.case = T)], format="BigWig", which=peaks.recentered_HH9, as="RleList")
#ATAC.bw_HH9_2 <- import(bigwig_files_ATAC[grepl("ATAC_HH9_R2.mLb.clN", bigwig_files_ATAC, ignore.case = T)], format="BigWig", which=peaks.recentered_HH9, as="RleList")

# import of CUT&RUN bigWig files for HH6

CUTandRUN.bw_HH6_1 <- import(bigwig_files_CUTandRUN[grepl("H3K27ac_HH6_R1", bigwig_files_CUTandRUN, ignore.case = T)], format="BigWig", which=peaks.recentered_HH6, as="RleList")
#CUTandRUN.bw_HH6_2 <- import(bigwig_files_CUTandRUN[grepl("H3K27ac_HH6_R2", bigwig_files_CUTandRUN, ignore.case = T)], format="BigWig", which=peaks.recentered_HH6, as="RleList")

# import of CUT&RUN bigWig files for HH9

CUTandRUN.bw_HH9_1 <- import(bigwig_files_CUTandRUN[grepl("H3K27ac_HH9_R1", bigwig_files_CUTandRUN, ignore.case = T)], format="BigWig", which=peaks.recentered_HH9, as="RleList")
#CUTandRUN.bw_HH9_2 <- import(bigwig_files_CUTandRUN[grepl("H3K27ac_HH9_R3", bigwig_files_CUTandRUN, ignore.case = T)], format="BigWig", which=peaks.recentered_HH9, as="RleList")

# "Input" serves as control data for HH6 and HH9

igg.bw_HH6 <- import(bigwig_files_CUTandRUN[grepl("igg_R1", bigwig_files_CUTandRUN, ignore.case = T)], format="BigWig", which=peaks.recentered_HH6, as="RleList")
igg.bw_HH9 <- import(bigwig_files_CUTandRUN[grepl("igg_R1", bigwig_files_CUTandRUN, ignore.case = T)], format="BigWig", which=peaks.recentered_HH9, as="RleList")


### create list of bigWig files for each stage (HH6 & 9) ###

 # ATAC bw lists HH6 & HH9 

bw_ATAC_HH6 <- list(ATAC = ATAC.bw_HH6_1)
bw_ATAC_HH9 <- list(ATAC = ATAC.bw_HH9_1)

 # CUT & RUN + igg  bw lists HH6 & HH9

bw_CRigg_HH6 <- list(H3K27ac = CUTandRUN.bw_HH6_1, IGG = igg.bw_HH6)
bw_CRigg_HH9 <- list(H3K27ac = CUTandRUN.bw_HH9_1, IGG = igg.bw_HH9)


### extract signal for +/-2kb around enhancer peak for visualisation for HH6 & 9###

 # ATAC signal HH6 & HH9 

sig_ATAC_HH6 <- featureAlignedSignal(bw_ATAC_HH6, peaks.recentered_HH6,
                            upstream=2000, downstream=2000)

sig_ATAC_HH9 <- featureAlignedSignal(bw_ATAC_HH9, peaks.recentered_HH9,
                                     upstream=2000, downstream=2000)

# CUT & RUN + igg signal HH6 & HH9 

sig_CRigg_HH6 <- featureAlignedSignal(bw_CRigg_HH6, peaks.recentered_HH6,
                                      upstream=2000, downstream=2000)

sig_CRigg_HH9 <- featureAlignedSignal(bw_CRigg_HH9, peaks.recentered_HH9,
                                upstream=2000, downstream=2000)


# plot profile around ATAC peaks for HH6 & HH9

 # ATAC plot profile HH6 & HH9

png(paste0("metaprofile_ATAC_HH6.png"), width=20, height=17, family = 'Arial', units = 'cm', res = 400)
featureAlignedDistribution(sig_ATAC_HH6, peaks.recentered_HH6, upstream=2000, downstream=2000, type="l")
graphics.off()

png(paste0("metaprofile_ATAC_HH9.png"), width=20, height=17, family = 'Arial', units = 'cm', res = 400)
featureAlignedDistribution(sig_ATAC_HH9, peaks.recentered_HH9, upstream=2000, downstream=2000, type="l")
graphics.off()

# C&R + igg plot profile HH6 & HH9

png(paste0("metaprofile_CRigg_HH6.png"), width=20, height=17, family = 'Arial', units = 'cm', res = 400)
featureAlignedDistribution(sig_CRigg_HH6, peaks.recentered_HH6, upstream=2000, downstream=2000, type="l")
graphics.off()

png(paste0("metaprofile_CRigg_HH9.png"), width=20, height=17, family = 'Arial', units = 'cm', res = 400)
featureAlignedDistribution(sig_CRigg_HH9, peaks.recentered_HH9, upstream=2000, downstream=2000, type="l")
graphics.off()


# plot heatmaps for HH6 & 9

 # ATAC heatmaps HH6 & HH9

png(paste0("heatmap_ATAC_HH6.png"), width=15, height=15, family = 'Arial', units = 'cm', res = 400)
featureAlignedHeatmap(sig_ATAC_HH6, peaks.recentered_HH6, upstream=2000, downstream=2000, upper.extreme=2.5)
graphics.off()

png(paste0("heatmap_ATAC_HH9.png"), width=15, height=15, family = 'Arial', units = 'cm', res = 400)
featureAlignedHeatmap(sig_ATAC_HH9, peaks.recentered_HH9, upstream=2000, downstream=2000, upper.extreme=2.5)
graphics.off()

# CUT & RUN + igg heatmaps HH6 & HH9

png(paste0("heatmap_CRigg_HH6.png"), width=15, height=15, family = 'Arial', units = 'cm', res = 400)
featureAlignedHeatmap(sig_CRigg_HH6, peaks.recentered_HH6, upstream=2000, downstream=2000, upper.extreme=500)
graphics.off()

png(paste0("heatmap_CRigg_HH9.png"), width=15, height=15, family = 'Arial', units = 'cm', res = 400)
featureAlignedHeatmap(sig_CRigg_HH9, peaks.recentered_HH9, upstream=2000, downstream=2000, upper.extreme=500)
graphics.off()

############################################################################

# Motif enrichment analysis from Homer output

#install.packages("ggseqlogo")
#install.packages("gridExtra")
#install.packages("cowplot")
#install.packages("gridGraphics")

library("grid")
library("gridGraphics")
library("ggseqlogo")
library("gridExtra")
library("cowplot")
library("ggplot2")
library("extrafont")


### create bedfiles from filter homer output

 # HH6 & HH9

enhancer_bedfile_HH6 <- select(enhancer_HH6, c(2, 3, 4, 1))
colnames(enhancer_bedfile_HH6) <- c('chrom', 'chromStart', 'chromEnd', 'name')
write.table(enhancer_bedfile_HH6, "enhancer_bedfile_HH6.bed", sep = '\t', quote=F, row.names=F, col.names = F)


enhancer_bedfile_HH9 <- select(enhancer_HH9, c(2, 3, 4, 1))
colnames(enhancer_bedfile_HH9) <- c('chrom', 'chromStart', 'chromEnd', 'name')
write.table(enhancer_bedfile_HH6, "enhancer_bedfile_HH9.bed", sep = '\t', quote=F, row.names=F, col.names = F)

######### intermediate steps performed in Ubuntu terminal using Homer package ######


 ### read in data (HH6 & HH9)


 #### read in logo data

motif_logos_HH6 = list()
for(i in 1:20){
  motif_logos_HH6[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/global/HH_6/knownResults/known', i, '.motif'))[1:4])
  rownames(motif_logos_HH6[[paste(i)]]) = c('A', 'C', 'G', 'T')}

motif_logos_HH9 = list()
for(i in 1:20){
  motif_logos_HH9[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/global/HH_9/knownResults/known', i, '.motif'))[1:4])
  rownames(motif_logos_HH9[[paste(i)]]) = c('A', 'C', 'G', 'T')}


# read in motif info

motif_meta_HH6 = read.delim(paste0('Homer_motifs/global/HH_6/knownResults.txt'))[1:20,c(1,5)]

motif_meta_HH9 = read.delim(paste0('Homer_motifs/global/HH_9/knownResults.txt'))[1:20,c(1,5)]

# create list of all identified motifs - clean names - and cutoff at q-value <= 0.05
global_list_HH6 = read.delim(paste0('Homer_motifs/global/HH_6/knownResults.txt'))[,c(1,5)]
global_list_HH9 = read.delim(paste0('Homer_motifs/global/HH_9/knownResults.txt'))[,c(1,5)]

global_list_HH6$Motif.Name <-sub ("\\(.*", "", global_list_HH6$Motif.Name)
global_list_HH9$Motif.Name <-sub ("\\(.*", "", global_list_HH9$Motif.Name)

global_list_HH6 <- filter(global_list_HH6, q.value..Benjamini. <= 0.05)
global_list_HH9 <- filter(global_list_HH9, q.value..Benjamini. <= 0.05)

# strip name

motif_meta_HH6[,1] <- sub("\\(.*", "", motif_meta_HH6[,1])

motif_meta_HH9[,1] <- sub("\\(.*", "", motif_meta_HH9[,1])


####### prepare grobs (HH6 & HH9)

# gene names

gene_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                      text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

gene_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                          text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

motif_names_HH6 <- lapply(motif_meta_HH6[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                             text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

motif_names_HH9 <- lapply(motif_meta_HH9[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                     text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

# motifs

motif_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                       text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

motif_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                           text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))


motif_logos_HH6 = lapply(motif_logos_HH6, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})

motif_logos_HH9 = lapply(motif_logos_HH9, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})


# q-values

qval_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                      text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))

qval_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                          text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))

motif_qval_HH6 <- lapply(motif_meta_HH6[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                            text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

motif_qval_HH9 <- lapply(motif_meta_HH9[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                    text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})


######## plot grobs

png(paste0('top20_motifs_HH6.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(gene_HH6, motif_names_HH6, motif_HH6, motif_logos_HH6, qval_HH6, motif_qval_HH6), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

png(paste0('top20_motifs_HH9.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(gene_HH9, motif_names_HH9, motif_HH9, motif_logos_HH9, qval_HH9, motif_qval_HH9), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

####### q-value fix (shows as 0 because smaller than 2.225074e-308) ###########

#.Machine$double.xmin
###############################################################################


########  plot selected motifs (HH6 & HH9) - tbd!!

#motifs_of_interest_HH6 <- c('CTCF', 'Sox3', 'Sox21', 'Sox2', 'Lhx3', 'DLX2', '4-SOX2-TCF-NAI', 'En1')
#motifs_of_interest_HH6 <- which(motif_meta_HH6$Motif.Name %in% motifs_of_interest_HH6)

#motifs_of_interest_HH9 <- c('Sox3', 'CTCF', 'Sox21', '4-SOX2-TCF-NAI', 'Sox2', 'DLX2', 'Lhx3', 'Sox10')
#motifs_of_interest_HH9 <- which(motif_meta_HH9$Motif.Name %in% motifs_of_interest_HH9)


#png(paste0('selected_motifs_HH6.png'), width = 480, height = 480, family = 'Arial', units = 'cm', res = 400)
#grid.arrange(grobs=c(gene_HH6, motif_names_HH6[motifs_of_interest_HH6], motif_HH6, motif_logos_HH6[motifs_of_interest_HH6], pval_HH6, motif_pval_HH6[motifs_of_interest_HH6]), ncol=3, widths = c(1, 4, 1), as.table=FALSE)
#graphics.off()

#png(paste0('selected_motifs_HH9.png'), width = 150, height = 150, family = 'Arial', units = 'cm', res = 400)
#grid.arrange(grobs=c(gene_HH9, motif_names_HH9[motifs_of_interest_HH9], motif_HH9, motif_logos_HH9[motifs_of_interest_HH9], pval_HH9, motif_pval_HH9[motifs_of_interest_HH9]), ncol=3, widths = c(1, 4, 1), as.table=FALSE)
#graphics.off()

################################################################################
################################################################################
################################################################################

###### clean version below hashed lines ######

######################## start integrating global (above) enhancers with scRNA-seq data #################

  # import gene module (GM) list from RNA-seq

#install.packages('stringi')
library('stringi')
library('data.table')
library(stringr)

# import Gene Module list, transpose and organise it

#transposed_gm_list <- transpose(read.table('gm.txt',fill=T,sep=""))
#colnames(transposed_gm_list) <- paste(transposed_gm_list[1,], transposed_gm_list[2,])
#GM_list <- transposed_gm_list[-c(1,2),]

# rename column names according to tissue

#setnames(GM_list, old = c('GM: 7;','GM: 8;', 'GM: 9;'), new = c('placodal','placodal','placodal'))
#setnames(GM_list, old = c('GM: 40;','GM: 41;'), new = c('NC','NC'))
#setnames(GM_list, old = c('GM: 36;','GM: 37;','GM: 39;','GM: 43;'), new = c('neural','neural','neural','neural'))

# create list reduced to reduce to NC, neural and placodal columns

#tissue_list <- list()
#tissue_list$placodal = unname(unlist(GM_list[grepl('placodal',names(GM_list))]))
#tissue_list$neural = unname(unlist(GM_list[grepl('neural',names(GM_list))]))
#tissue_list$NC = unname(unlist(GM_list[grepl('NC',names(GM_list))]))



######## clean version ##########


# read GM list into R

new_list <- scan('gm.txt',what="", sep="\n")

# extract "GM: XY;" elements as future names for the elements (i.e. modules) of the list

names_list <- gsub(";.*","",new_list)

# remove "GM: XY;" elements just stored to "names_list" from GM list elements

new_list <- gsub(".*; ","",new_list)

# split the elements (i.e. genes) in the list at the commas between them

new_list <- str_split(new_list,pattern = ", ")

# apply names (i.e. GM modules) to GM list

names(new_list) <- names_list

# group tissues (placodal, NC, neural) together

tissue_list <- list('placodal' = new_list[c('GM: 7','GM: 8','GM: 9')], 'NC' = new_list[c('GM: 40','GM: 41')], "neural" = new_list[c('GM: 36','GM: 37','GM: 39','GM: 43')])

abc = list()
abc$placodal = unname(unlist(tissue_list$placodal))
abc$NC = unname(unlist(tissue_list$NC))
abc$neural = unname(unlist(tissue_list$neural))

# fill empty "Gene Name" fields with "Intrez.ID" (i.e. Gene ID) in enhancer_HHx dataframes

filled_enhancer_HH6 <- enhancer_HH6
filled_enhancer_HH9 <- enhancer_HH9

filled_enhancer_HH6$Gene.Name <- ifelse(filled_enhancer_HH6$Gene.Name == "", filled_enhancer_HH6$Entrez.ID, filled_enhancer_HH6$Gene.Name)
filled_enhancer_HH9$Gene.Name <- ifelse(filled_enhancer_HH9$Gene.Name == "", filled_enhancer_HH9$Entrez.ID, filled_enhancer_HH9$Gene.Name)

# subset filled_enhancer dfs based in tissue list (i.e. create placodal, neural and NC subset df)

#intersect(filled_enhancer_HH6$Gene.Name, abc$NC)

# tissue enhancers for HH6

placodal_enhancer_HH6 <- subset(filled_enhancer_HH6, filled_enhancer_HH6$Gene.Name %in% abc$placodal)
NC_enhancer_HH6 <- subset(filled_enhancer_HH6, filled_enhancer_HH6$Gene.Name %in% abc$NC)
neural_enhancer_HH6 <- subset(filled_enhancer_HH6, filled_enhancer_HH6$Gene.Name %in% abc$neural)

# tissue enhancers for HH9

placodal_enhancer_HH9 <- subset(filled_enhancer_HH9, filled_enhancer_HH9$Gene.Name %in% abc$placodal)
NC_enhancer_HH9 <- subset(filled_enhancer_HH9, filled_enhancer_HH9$Gene.Name %in% abc$NC)
neural_enhancer_HH9 <- subset(filled_enhancer_HH9, filled_enhancer_HH9$Gene.Name %in% abc$neural)

# confirm number of unique overlaps per tissue

length(abc$placodal)
# 74

length(abc$NC)
# 18

length(abc$neural)
# 66

length(unique(placodal_enhancer_HH6$Gene.Name))
length(intersect(placodal_enhancer_HH6$Gene.Name, abc$placodal))
# both 59 - ok

length(unique(NC_enhancer_HH6$Gene.Name))
length(intersect(NC_enhancer_HH6$Gene.Name, abc$NC))
#both 16 - ok

length(unique(neural_enhancer_HH6$Gene.Name))
length(intersect(neural_enhancer_HH6$Gene.Name, abc$neural))
#both 51 - ok

length(unique(placodal_enhancer_HH9$Gene.Name))
length(intersect(placodal_enhancer_HH9$Gene.Name, abc$placodal))
# both 53 - ok

length(unique(NC_enhancer_HH9$Gene.Name))
length(intersect(NC_enhancer_HH9$Gene.Name, abc$NC))
#both 16 - ok

length(unique(neural_enhancer_HH9$Gene.Name))
length(intersect(neural_enhancer_HH9$Gene.Name, abc$neural))
#both 48 - ok


### arrange filtered tissue data frames (df) and create bedfiles 

dir.create('~/dev/repos/Rothstein_2019/downsteam_R_analysis/tissue_bedfiles')

# placodal HH6 & HH9

placodal_bedfile_HH6 <- select(placodal_enhancer_HH6, c(2, 3, 4, 1))
colnames(placodal_bedfile_HH6) <- c('chrom', 'chromStart', 'chromEnd', 'name')

placodal_bedfile_HH9 <- select(placodal_enhancer_HH9, c(2, 3, 4, 1))
colnames(placodal_bedfile_HH9) <- c('chrom', 'chromStart', 'chromEnd', 'name')

# NC HH6 & HH9

NC_bedfile_HH6 <- select(NC_enhancer_HH6, c(2, 3, 4, 1))
colnames(NC_bedfile_HH6) <- c('chrom', 'chromStart', 'chromEnd', 'name')

NC_bedfile_HH9 <- select(NC_enhancer_HH9, c(2, 3, 4, 1))
colnames(NC_bedfile_HH9) <- c('chrom', 'chromStart', 'chromEnd', 'name')

# neural HH6 & HH9

neural_bedfile_HH6 <- select(neural_enhancer_HH6, c(2, 3, 4, 1))
colnames(neural_bedfile_HH6) <- c('chrom', 'chromStart', 'chromEnd', 'name')

neural_bedfile_HH9 <- select(neural_enhancer_HH9, c(2, 3, 4, 1))
colnames(neural_bedfile_HH9) <- c('chrom', 'chromStart', 'chromEnd', 'name')

# write bedfiles from arranged dfs

# create list with arranged tissue dfs (selected columns required for bedfile format) 

df_list <- list('placodal_HH6' = placodal_bedfile_HH6, 'placodal_HH9'= placodal_bedfile_HH9, 'NC_HH6' = NC_bedfile_HH6, 'NC_HH9' = NC_bedfile_HH9, 'neural_HH6' =  neural_bedfile_HH6 , 'neural_HH9' = neural_bedfile_HH9)

# create bedfiles based on dfs in df_list

# either
for (i in names(df_list)){
  write.table(df_list[[i]], file = paste0('~/dev/repos/Rothstein_2019/downsteam_R_analysis/tissue_bedfiles/', i, '.bed'), sep = '\t', quote=F, row.names=F, col.names = F)}

# or
lapply(names(df_list), function(x) write.table(df_list[[x]], paste0('~/dev/repos/Rothstein_2019/downsteam_R_analysis/tissue_bedfiles/', x, '.bed'), sep = '\t', quote=F, row.names=F, col.names = F))


######### intermediate steps performed in Ubuntu terminal using Homer package ######

### read in data (for 3 tissues and 2 time points)


### read in logo data

 # placodal 

placodal_motif_logos_HH6 = list()
for(i in 1:20){
  placodal_motif_logos_HH6[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/placodal/HH_6/knownResults/known', i, '.motif'))[1:4])
  rownames(placodal_motif_logos_HH6[[paste(i)]]) = c('A', 'C', 'G', 'T')}

placodal_motif_logos_HH9 = list()
for(i in 1:20){
  placodal_motif_logos_HH9[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/placodal/HH_9/knownResults/known', i, '.motif'))[1:4])
  rownames(placodal_motif_logos_HH9[[paste(i)]]) = c('A', 'C', 'G', 'T')}

 # NC

NC_motif_logos_HH6 = list()
for(i in 1:20){
  NC_motif_logos_HH6[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/NC/HH_6/knownResults/known', i, '.motif'))[1:4])
  rownames(NC_motif_logos_HH6[[paste(i)]]) = c('A', 'C', 'G', 'T')}

# top 18 moitifs significant only
NC_motif_logos_HH9 = list() 
for(i in 1:18){ 
  NC_motif_logos_HH9[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/NC/HH_9/knownResults/known', i, '.motif'))[1:4])
  rownames(NC_motif_logos_HH9[[paste(i)]]) = c('A', 'C', 'G', 'T')}

 # neural

neural_motif_logos_HH6 = list()
for(i in 1:20){
  neural_motif_logos_HH6[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/neural/HH_6/knownResults/known', i, '.motif'))[1:4])
  rownames(neural_motif_logos_HH6[[paste(i)]]) = c('A', 'C', 'G', 'T')}

neural_motif_logos_HH9 = list()
for(i in 1:20){
  neural_motif_logos_HH9[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/neural/HH_9/knownResults/known', i, '.motif'))[1:4])
  rownames(neural_motif_logos_HH9[[paste(i)]]) = c('A', 'C', 'G', 'T')}

### read in motif info

 # placodal

placodal_motif_meta_HH6 = read.delim(paste0('Homer_motifs/placodal/HH_6/knownResults.txt'))[1:20,c(1,5)]

placodal_motif_meta_HH9 = read.delim(paste0('Homer_motifs/placodal/HH_9/knownResults.txt'))[1:20,c(1,5)]

 # NC

NC_motif_meta_HH6 = read.delim(paste0('Homer_motifs/NC/HH_6/knownResults.txt'))[1:20,c(1,5)]

NC_motif_meta_HH9 = read.delim(paste0('Homer_motifs/NC/HH_9/knownResults.txt'))[1:18,c(1,5)] # top 18 motifs significant only

 # neural

neural_motif_meta_HH6 = read.delim(paste0('Homer_motifs/neural/HH_6/knownResults.txt'))[1:20,c(1,5)]

neural_motif_meta_HH9 = read.delim(paste0('Homer_motifs/neural/HH_9/knownResults.txt'))[1:20,c(1,5)]

### strip name

 # placodal

placodal_motif_meta_HH6[,1] <- sub("\\(.*", "", placodal_motif_meta_HH6[,1])

placodal_motif_meta_HH9[,1] <- sub("\\(.*", "", placodal_motif_meta_HH9[,1])

 # NC

NC_motif_meta_HH6[,1] <- sub("\\(.*", "", NC_motif_meta_HH6[,1])

NC_motif_meta_HH9[,1] <- sub("\\(.*", "", NC_motif_meta_HH9[,1])

 # neural

neural_motif_meta_HH6[,1] <- sub("\\(.*", "", neural_motif_meta_HH6[,1])

neural_motif_meta_HH9[,1] <- sub("\\(.*", "", neural_motif_meta_HH9[,1])

####### prepare grobs 

### gene names

 # placodal

placodal_gene_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                          text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

placodal_gene_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

placodal_motif_names_HH6 <- lapply(placodal_motif_meta_HH6[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                     text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

placodal_motif_names_HH9 <- lapply(placodal_motif_meta_HH9[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                       text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})
 # NC

NC_gene_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

NC_gene_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

NC_motif_names_HH6 <- lapply(NC_motif_meta_HH6[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                       text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

NC_motif_names_HH9 <- lapply(NC_motif_meta_HH9[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                       text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})                                                
 # neural

neural_gene_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                             text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

neural_gene_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                             text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

neural_motif_names_HH6 <- lapply(neural_motif_meta_HH6[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                           text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

neural_motif_names_HH9 <- lapply(neural_motif_meta_HH9[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                           text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))}) 

### motifs

 # placodal

placodal_motif_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                           text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

placodal_motif_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                    text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

placodal_motif_logos_HH6 = lapply(placodal_motif_logos_HH6, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})

placodal_motif_logos_HH9 = lapply(placodal_motif_logos_HH9, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})

 # NC

NC_motif_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                    text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

NC_motif_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                    text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

NC_motif_logos_HH6 = lapply(NC_motif_logos_HH6, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})

NC_motif_logos_HH9 = lapply(NC_motif_logos_HH9, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})

 # neural

neural_motif_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                              text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

neural_motif_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                              text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

neural_motif_logos_HH6 = lapply(neural_motif_logos_HH6, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})

neural_motif_logos_HH9 = lapply(neural_motif_logos_HH9, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})

### q-values

 # placodal

placodal_qval_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                          text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))
placodal_qval_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))

placodal_motif_qval_HH6 <- lapply(placodal_motif_meta_HH6[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                    text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})
placodal_motif_qval_HH9 <- lapply(placodal_motif_meta_HH9[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                      text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

 # NC

NC_qval_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))
NC_qval_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))

NC_motif_qval_HH6 <- lapply(NC_motif_meta_HH6[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                      text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})
NC_motif_qval_HH9 <- lapply(NC_motif_meta_HH9[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                      text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

 # neural

neural_qval_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                             text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))
neural_qval_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                             text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))

neural_motif_qval_HH6 <- lapply(neural_motif_meta_HH6[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                          text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})
neural_motif_qval_HH9 <- lapply(neural_motif_meta_HH9[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                          text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

### plot grobs

 # placodal

png(paste0('placodal_top20_motifs_HH6.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(placodal_gene_HH6, placodal_motif_names_HH6, placodal_motif_HH6, placodal_motif_logos_HH6, placodal_qval_HH6, placodal_motif_qval_HH6), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

png(paste0('placodal_top20_motifs_HH9.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(placodal_gene_HH9, placodal_motif_names_HH9, placodal_motif_HH9, placodal_motif_logos_HH9, placodal_qval_HH9, placodal_motif_qval_HH9), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

 # NC

png(paste0('NC_top20_motifs_HH6.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(NC_gene_HH6, NC_motif_names_HH6, NC_motif_HH6, NC_motif_logos_HH6, NC_qval_HH6, NC_motif_qval_HH6), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

png(paste0('NC_top20_motifs_HH9.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(NC_gene_HH9, NC_motif_names_HH9, NC_motif_HH9, NC_motif_logos_HH9, NC_qval_HH9, NC_motif_qval_HH9), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

 # neural

png(paste0('neural_top20_motifs_HH6.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(neural_gene_HH6, neural_motif_names_HH6, neural_motif_HH6, neural_motif_logos_HH6, neural_qval_HH6, neural_motif_qval_HH6), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

png(paste0('neural_top20_motifs_HH9.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(neural_gene_HH9, neural_motif_names_HH9, neural_motif_HH9, neural_motif_logos_HH9, neural_qval_HH9, neural_motif_qval_HH9), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

##########################################################################################

### create lists of all enriched enhancers identified by homer

#install.packages('gt')
#install.packages('plyr')
library(plyr)
library(gt)
library(dplyr)
library(stringr)
library(utils)

## read in data

complete_placodal_motif_list_HH6 = read.delim(paste0('Homer_motifs/placodal/HH_6/knownResults.txt'))[,c(1,5)]
complete_placodal_motif_list_HH9 = read.delim(paste0('Homer_motifs/placodal/HH_9/knownResults.txt'))[,c(1,5)]

complete_NC_motif_list_HH6 = read.delim(paste0('Homer_motifs/NC/HH_6/knownResults.txt'))[,c(1,5)]
complete_NC_motif_list_HH9 = read.delim(paste0('Homer_motifs/NC/HH_9/knownResults.txt'))[,c(1,5)]

complete_neural_motif_list_HH6 = read.delim(paste0('Homer_motifs/neural/HH_6/knownResults.txt'))[,c(1,5)]
complete_neural_motif_list_HH9 = read.delim(paste0('Homer_motifs/neural/HH_9/knownResults.txt'))[,c(1,5)]

complete_NC_GM_40_motif_list_HH6 = read.delim(paste0('Homer_motifs/NC_GM_40/HH_6/knownResults.txt'))[,c(1,5)]
complete_NC_GM_40_motif_list_HH9 = read.delim(paste0('Homer_motifs/NC_GM_40/HH_9/knownResults.txt'))[,c(1,5)]

complete_NC_GM_41_motif_list_HH6 = read.delim(paste0('Homer_motifs/NC_GM_41/HH_6/knownResults.txt'))[,c(1,5)]
complete_NC_GM_41_motif_list_HH9 = read.delim(paste0('Homer_motifs/NC_GM_41/HH_9/knownResults.txt'))[,c(1,5)]

# cutoff q-values above 0.05

cutoff_placodal_list_HH6 <- filter(complete_placodal_motif_list_HH6, q.value..Benjamini. <= 0.05)
cutoff_placodal_list_HH9 <- filter(complete_placodal_motif_list_HH9, q.value..Benjamini. <= 0.05)

cutoff_NC_list_HH6 <- filter(complete_NC_motif_list_HH6, q.value..Benjamini. <= 0.05)
cutoff_NC_list_HH9 <- filter(complete_NC_motif_list_HH9, q.value..Benjamini. <= 0.05)

cutoff_neural_list_HH6 <- filter(complete_neural_motif_list_HH6, q.value..Benjamini. <= 0.05)
cutoff_neural_list_HH9 <- filter(complete_neural_motif_list_HH9, q.value..Benjamini. <= 0.05)

cutoff_NC_GM_40_list_HH6 <- filter(complete_NC_GM_40_motif_list_HH6, q.value..Benjamini. <= 0.05)
cutoff_NC_GM_40_list_HH9 <- filter(complete_NC_GM_40_motif_list_HH9, q.value..Benjamini. <= 0.05)

cutoff_NC_GM_41_list_HH6 <- filter(complete_NC_GM_41_motif_list_HH6, q.value..Benjamini. <= 0.05)
cutoff_NC_GM_41_list_HH9 <- filter(complete_NC_GM_41_motif_list_HH9, q.value..Benjamini. <= 0.05)

# clean gene names 

cutoff_placodal_list_HH6[,1] <- sub("\\(.*","", cutoff_placodal_list_HH6[,1])
cutoff_placodal_list_HH6[,1] <- sub("/Homer.*","", cutoff_placodal_list_HH6[,1])
cutoff_placodal_list_HH6[,1] <- sub(".Flag.*","", cutoff_placodal_list_HH6[,1])
cutoff_placodal_list_HH9[,1] <- sub("\\(.*","", cutoff_placodal_list_HH9[,1])
cutoff_placodal_list_HH9[,1] <- sub("/Homer.*","", cutoff_placodal_list_HH9[,1])

cutoff_NC_list_HH6[,1] <- sub("\\(.*","", cutoff_NC_list_HH6[,1])
cutoff_NC_list_HH6[,1] <- sub("/Homer.*","", cutoff_NC_list_HH6[,1])
cutoff_NC_list_HH6[,1] <- sub(".Flag.*","", cutoff_NC_list_HH6[,1])
cutoff_NC_list_HH9[,1] <- sub("\\(.*","", cutoff_NC_list_HH9[,1])

cutoff_neural_list_HH6[,1] <- sub("\\(.*","", cutoff_neural_list_HH6[,1])
cutoff_neural_list_HH6[,1] <- sub("/Homer.*","", cutoff_neural_list_HH6[,1])
cutoff_neural_list_HH9[,1] <- sub("\\(.*","", cutoff_neural_list_HH9[,1])

cutoff_NC_GM_40_list_HH6[,1] <- sub("\\(.*","", cutoff_NC_GM_40_list_HH6[,1])
cutoff_NC_GM_40_list_HH6[,1] <- sub("/Homer.*","", cutoff_NC_GM_40_list_HH6[,1])
cutoff_NC_GM_40_list_HH9[,1] <- sub("\\(.*","", cutoff_NC_GM_40_list_HH9[,1])

cutoff_NC_GM_41_list_HH6[,1] <- sub("\\(.*","", cutoff_NC_GM_41_list_HH6[,1])
cutoff_NC_GM_41_list_HH6[,1] <- sub("/Homer.*","", cutoff_NC_GM_41_list_HH6[,1])
cutoff_NC_GM_41_list_HH9[,1] <- sub("\\(.*","", cutoff_NC_GM_41_list_HH9[,1])


###################### create list of genes/enhancers ###################################

Venn_list <- list('placodal HH6' = cutoff_placodal_list_HH6$Motif.Name,
                  'placodal HH9' = cutoff_placodal_list_HH9$Motif.Name,
                  'NC HH6' = cutoff_NC_list_HH6$Motif.Name,
                  'NC HH9' = cutoff_NC_list_HH9$Motif.Name,
                  'neural HH6' = cutoff_neural_list_HH6$Motif.Name,
                  'neural HH9' = cutoff_neural_list_HH9$Motif.Name,
                  'NC_GM_40 HH6' = cutoff_NC_GM_40_list_HH6$Motif.Name,
                  'NC_GM_40 HH9' = cutoff_NC_GM_40_list_HH9$Motif.Name,
                  'NC_GM_41 HH6' = cutoff_NC_GM_41_list_HH6$Motif.Name,
                  'NC_GM_41 HH9' = cutoff_NC_GM_41_list_HH9$Motif.Name)
capture.output(Venn_list, file = "Venn_list.txt", quote = F)

######################## perform comparison of NC modules (received from Alex) enhancer enrichment (i.e. do separated enhancer plotting for the two modules)

# read GM list into R
new_list <- scan('gm.txt',what="", sep="\n")

# extract names from GM list (i.e. GM names)
names_list <- gsub(";.*","", new_list)

# remove names saved in names_list from initial list
new_list <- gsub(".*;","", new_list)

# split list
split_list <- str_split(new_list, ", ")

# apply names from names_list
names(split_list) <- names_list

# fill empty cells in enhancer gene name column with entrez-ID
enhancer_HH6$Gene.Name <- ifelse(enhancer_HH6$Gene.Name == "", enhancer_HH6$Entrez.ID, enhancer_HH6$Gene.Name)
enhancer_HH9$Gene.Name <- ifelse(enhancer_HH9$Gene.Name == "", enhancer_HH9$Entrez.ID, enhancer_HH9$Gene.Name)

# subset enhancer df's based on NC modules (GM 40 & GM 41) for HH6 & HH9
NC_GM_40_HH6 <- subset(enhancer_HH6, enhancer_HH6$Gene.Name %in% split_list$`GM: 40`)
NC_GM_40_HH9 <- subset(enhancer_HH9, enhancer_HH9$Gene.Name %in% split_list$`GM: 40`)

NC_GM_41_HH6 <- subset(enhancer_HH6, enhancer_HH6$Gene.Name %in% split_list$`GM: 41`)
NC_GM_41_HH9 <- subset(enhancer_HH9, enhancer_HH9$Gene.Name %in% split_list$`GM: 41`)

### create (and name elements of) list containing dfs for subsequent looping of arrangement and file export
NC_GM_list <- list(NC_GM_40_HH6, NC_GM_40_HH9, NC_GM_41_HH6, NC_GM_41_HH9)
names(NC_GM_list) <- c("NC_GM_40_HH6", "NC_GM_40_HH9", "NC_GM_41_HH6", "NC_GM_41_HH9")

#NC_GM_bed <- lapply(names(NC_GM_list), function(x) NC_GM_list[[x]] %>% select(c(2, 3, 4, 1)))
#names(NC_GM_bed) <- c("NC_GM_40_HH6", "NC_GM_40_HH9", "NC_GM_41_HH6", "NC_GM_41_HH9")

NC_GM_list <- lapply(NC_GM_list, function(x) x %>% select(c(2, 3, 4, 1)))
#NC_GM_list <- lapply(NC_GM_list, function(x) colnames(x) <- c('chrom', 'chromStart', 'chromEnd', 'name'))
#############
#remove(NC_GM_list)

#colnames(NC_GM_list$NC_GM_41_HH9) <-c('chrom', 'chromStart', 'chromEnd', 'name')

lapply(names(NC_GM_list), function(x) write.table(NC_GM_list[[x]], paste0('~/dev/repos/Rothstein_2019/downsteam_R_analysis/tissue_bedfiles/', x, '.bed'), sep ='\t', quote = F, row.names = F, col.names = F ))

lapply(names(df_list), function(x) write.table(df_list[[x]], paste0('~/dev/repos/Rothstein_2019/downsteam_R_analysis/tissue_bedfiles/', x, '.bed'), sep = '\t', quote=F, row.names=F, col.names = F))

################## intermediate steps in Ubuntu with Homer ##############################

### read in data (for 2 NC GMs at two time points)

library("grid")
library("gridGraphics")
library("ggseqlogo")
library("gridExtra")
library("cowplot")
library("ggplot2")
library("extrafont")

### read in logo data

# NC GM 40

# top 18 motifs significant only
NC_GM_40_motif_logos_HH6 = list()
for(i in 1:18){
  NC_GM_40_motif_logos_HH6[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/NC_GM_40/HH_6/knownResults/known', i, '.motif'))[1:4])
  rownames(NC_GM_40_motif_logos_HH6[[paste(i)]]) = c('A', 'C', 'G', 'T')}

# for GM 40 at HH9 no significant motifs were identified by Homer!
#NC_GM_40_motif_logos_HH9 = list()
#for(i in 1:14){
  #NC_GM_40_motif_logos_HH9[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/NC_GM_40/HH_9/knownResults/known', i, '.motif'))[1:4])
  #rownames(NC_GM_40_motif_logos_HH9[[paste(i)]]) = c('A', 'C', 'G', 'T')}

# NC GM 41

NC_GM_41_motif_logos_HH6 = list()
for(i in 1:20){
  NC_GM_41_motif_logos_HH6[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/NC_GM_41/HH_6/knownResults/known', i, '.motif'))[1:4])
  rownames(NC_GM_41_motif_logos_HH6[[paste(i)]]) = c('A', 'C', 'G', 'T')}

# top 13 motifs significant only
NC_GM_41_motif_logos_HH9 = list()
for(i in 1:13){
  NC_GM_41_motif_logos_HH9[[paste(i)]] <- t(read.delim(paste0('Homer_motifs/NC_GM_41/HH_9/knownResults/known', i, '.motif'))[1:4])
  rownames(NC_GM_41_motif_logos_HH9[[paste(i)]]) = c('A', 'C', 'G', 'T')}


### read in motif info

NC_GM_40_motif_meta_HH6 = read.delim(paste0('Homer_motifs/NC_GM_40/HH_6/knownResults.txt'))[1:18,c(1,5)] # top 18 motifs significant only
#NC_GM_40_motif_meta_HH9 = read.delim(paste0('Homer_motifs/NC_GM_40/HH_9/knownResults.txt'))[1:14,c(1,5)] # no significant motifs

NC_GM_41_motif_meta_HH6 = read.delim(paste0('Homer_motifs/NC_GM_41/HH_6/knownResults.txt'))[1:20,c(1,5)]
NC_GM_41_motif_meta_HH9 = read.delim(paste0('Homer_motifs/NC_GM_41/HH_9/knownResults.txt'))[1:13,c(1,5)] # top 13 motifs significant only


### strip name

NC_GM_40_motif_meta_HH6[,1] <- sub("\\(.*", "", NC_GM_40_motif_meta_HH6[,1])
#NC_GM_40_motif_meta_HH9[,1] <- sub("\\(.*", "", NC_GM_40_motif_meta_HH9[,1])

NC_GM_41_motif_meta_HH6[,1] <- sub("\\(.*", "", NC_GM_41_motif_meta_HH6[,1])
NC_GM_41_motif_meta_HH9[,1] <- sub("\\(.*", "", NC_GM_41_motif_meta_HH9[,1])


####### prepare grobs 

### gene names

NC_GM_40_gene_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

#NC_GM_40_gene_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   #text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

NC_GM_40_motif_names_HH6 <- lapply(NC_GM_40_motif_meta_HH6[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                       text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

#NC_GM_40_motif_names_HH9 <- lapply(NC_GM_40_motif_meta_HH9[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                       #text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

NC_GM_41_gene_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                              text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

NC_GM_41_gene_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "Gene", cex = 15, col = "black", font=2)))

NC_GM_41_motif_names_HH6 <- lapply(NC_GM_41_motif_meta_HH6[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                       text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

NC_GM_41_motif_names_HH9 <- lapply(NC_GM_41_motif_meta_HH9[,1], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                       text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

### motifs

NC_GM_40_motif_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                    text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

#NC_GM_40_motif_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                    #text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

NC_GM_40_motif_logos_HH6 = lapply(NC_GM_40_motif_logos_HH6, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})

#NC_GM_40_motif_logos_HH9 = lapply(NC_GM_40_motif_logos_HH9, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})


NC_GM_41_motif_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                    text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

NC_GM_41_motif_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                    text(x = 0.5, y = 0.5, "Motif", cex = 15, col = "black", font=2)))

NC_GM_41_motif_logos_HH6 = lapply(NC_GM_41_motif_logos_HH6, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})

NC_GM_41_motif_logos_HH9 = lapply(NC_GM_41_motif_logos_HH9, function(x) {ggseqlogo(x, method = 'prob') + theme_void() + theme(plot.margin = unit(c(2,0,2,0), "cm"))})


### q-values

NC_GM_40_qval_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))
#NC_GM_40_qval_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                  # text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))

NC_GM_40_motif_qval_HH6 <- lapply(NC_GM_40_motif_meta_HH6[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                      text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})
#NC_GM_40_motif_qval_HH9 <- lapply(NC_GM_40_motif_meta_HH9[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                  #    text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

NC_GM_41_qval_HH6 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))
NC_GM_41_qval_HH9 = list(as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                   text(x = 0.5, y = 0.5, "q-value", cex = 15, col = "black", font=2)))

NC_GM_41_motif_qval_HH6 <- lapply(NC_GM_41_motif_meta_HH6[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                      text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})
NC_GM_41_motif_qval_HH9 <- lapply(NC_GM_41_motif_meta_HH9[,2], function(x) {as_grob(~plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n') +
                                                                                      text(x = 0.5, y = 0.5, x, cex = 10, col = "black"))})

### plot grobs

png(paste0('NC_GM_40_top20_motifs_HH6.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(NC_GM_40_gene_HH6, NC_GM_40_motif_names_HH6, NC_GM_40_motif_HH6, NC_GM_40_motif_logos_HH6, NC_GM_40_qval_HH6, NC_GM_40_motif_qval_HH6), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

#png(paste0('NC_GM_40_top20_motifs_HH9.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
#grid.arrange(grobs=c(NC_GM_40_gene_HH9, NC_GM_40_motif_names_HH9, NC_GM_40_motif_HH9, NC_GM_40_motif_logos_HH9, NC_GM_40_qval_HH9, NC_GM_40_motif_qval_HH9), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
#graphics.off()

png(paste0('NC_GM_41_top20_motifs_HH6.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(NC_GM_41_gene_HH6, NC_GM_41_motif_names_HH6, NC_GM_41_motif_HH6, NC_GM_41_motif_logos_HH6, NC_GM_41_qval_HH6, NC_GM_41_motif_qval_HH6), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

png(paste0('NC_GM_41_top20_motifs_HH9.png'), width = 180, height = 300, family = 'Arial', units = 'cm', res = 200)
grid.arrange(grobs=c(NC_GM_41_gene_HH9, NC_GM_41_motif_names_HH9, NC_GM_41_motif_HH9, NC_GM_41_motif_logos_HH9, NC_GM_41_qval_HH9, NC_GM_41_motif_qval_HH9), ncol=3, widths = c(1, 3, 1), as.table=FALSE)
graphics.off()

###############################################################

# plotting Venn diagrams/arranging tables comparing differences and overlaps within tissues (e.g. placodal HH6 vs HH9)
  # see "create lists of all enriched enhancers identified by homer" above for data reading and preparation (line 861ff)

#install.packages('VennDiagram')
#install.packages('tidyverse')
#install.packages('ggforce')
#install.packages('RVenn')
library(RVenn)
library(VennDiagram)
library(tidyverse)
library(ggforce)
library(ggplot2)
library(tidyr)

### plotting per tissue

# placodal
venn.diagram(list("placodal HH6" = Venn_list$`placodal HH6`, "placodal HH9" = Venn_list$`placodal HH9`),
                                 fill = c("green", "blue"), alpha = c(0.5, 0.5), cat.cex = 1.5, cex = 1.0,
                                 margin = 0.15, filename = "Venn_placodal.png", imagetype = "png",
                                 scaled = T)

placodal_list <- list('HH6' = setdiff(Venn_list$`placodal HH6`, Venn_list$`placodal HH9`),
                   'overlap' = intersect(Venn_list$`placodal HH6`,Venn_list$`placodal HH9`),
                   'HH9'= setdiff(Venn_list$`placodal HH9`, Venn_list$`placodal HH6`))
placodal_table <- t(plyr::ldply(placodal_list, rbind))
colnames(placodal_table) <- placodal_table[1,] # applying first row as column name (e.g. "NC HH6")
placodal_table <- placodal_table[-1,] # removed first row containing data label (e.g. "NC HH6")
placodal_table[is.na(placodal_table)] = "" # replace "na" with empty fields
write.csv(placodal_table, file = "placodal.table.csv", quote = F)


# NC
venn.diagram(list("NC HH6" = Venn_list$`NC HH6`, "NC HH9" = Venn_list$`NC HH9`),
             fill = c("green", "blue"), alpha = c(0.5, 0.5), cat.cex = 1.5, cex = 1.0,
             margin = 0.15, filename = "Venn_NC.png", imagetype = "png",
             scaled = T)
NC_list <- list('HH6' = setdiff(Venn_list$`NC HH6`, Venn_list$`NC HH9`),
                      'overlap' = intersect(Venn_list$`NC HH6`,Venn_list$`NC HH9`),
                      'HH9'= setdiff(Venn_list$`NC HH9`, Venn_list$`NC HH6`))
NC_table <- t(plyr::ldply(NC_list, rbind))
colnames(NC_table) <- NC_table[1,] # applying first row as column name (e.g. "NC HH6")
NC_table <- NC_table[-1,] # removed first row containing data label (e.g. "NC HH6")
NC_table[is.na(NC_table)] = "" # replace "na" with empty fields
write.csv(NC_table, file = "NC.table.csv", quote = F)


# neural
venn.diagram(list("neural HH6" = Venn_list$`neural HH6`, "neural HH9" = Venn_list$`neural HH9`),
             fill = c("green", "blue"), alpha = c(0.5, 0.5), cat.cex = 1.5, cex = 1.0,
             margin = 0.15, filename = "Venn_neural.png", imagetype = "png",
             scaled = T)
neural_list <- list('HH6' = setdiff(Venn_list$`neural HH6`, Venn_list$`neural HH9`),
                'overlap' = intersect(Venn_list$`neural HH6`,Venn_list$`neural HH9`),
                'HH9'= setdiff(Venn_list$`neural HH9`, Venn_list$`neural HH6`))
neural_table <- t(plyr::ldply(neural_list, rbind))
colnames(neural_table) <- neural_table[1,] # applying first row as column name (e.g. "NC HH6")
neural_table <- neural_table[-1,] # removed first row containing data label (e.g. "NC HH6")
neural_table[is.na(neural_table)] = "" # replace "na" with empty fields
write.csv(neural_table, file = "neural.table.csv", quote = F)

# NC modules
venn.diagram(list("NC GM 40 HH6" = Venn_list$`NC_GM_40 HH6`, "NC GM 40 HH9" = Venn_list$`NC_GM_40 HH9`),
             fill = c("green", "blue"), alpha = c(0.5, 0.5), cat.cex = 1.5, cex = 1.0,
             margin = 0.15, filename = "Venn_NC_GM_40.png", imagetype = "png",
             scaled = T)
NC_GM_40_list <- list('HH6' = setdiff(Venn_list$`NC_GM_40 HH6`, Venn_list$`NC_GM_40 HH9`),
                    'overlap' = intersect(Venn_list$`NC_GM_40 HH6`,Venn_list$`NC_GM_40 HH9`),
                    'HH9'= setdiff(Venn_list$`NC_GM_40 HH9`, Venn_list$`NC_GM_40 HH6`))
NC_GM_40_table <- t(plyr::ldply(NC_GM_40_list, rbind))
colnames(NC_GM_40_table) <- NC_GM_40_table[1,] # applying first row as column name (e.g. "NC HH6")
NC_GM_40_table <- NC_GM_40_table[-1,] # removed first row containing data label (e.g. "NC HH6")
NC_GM_40_table[is.na(NC_GM_40_table)] = "" # replace "na" with empty fields
write.csv(NC_GM_40_table, file = "NC_GM_40.table.csv", quote = F)

venn.diagram(list("NC GM 41 HH6" = Venn_list$`NC_GM_41 HH6`, "NC GM 41 HH9" = Venn_list$`NC_GM_41 HH9`),
             fill = c("green", "blue"), alpha = c(0.5, 0.5), cat.cex = 1.5, cex = 1.0,
             margin = 0.15, filename = "Venn_NC_GM_41.png", imagetype = "png",
             scaled = T)
NC_GM_41_list <- list('HH6' = setdiff(Venn_list$`NC_GM_41 HH6`, Venn_list$`NC_GM_41 HH9`),
                      'overlap' = intersect(Venn_list$`NC_GM_41 HH6`,Venn_list$`NC_GM_41 HH9`),
                      'HH9'= setdiff(Venn_list$`NC_GM_41 HH9`, Venn_list$`NC_GM_41 HH6`))
NC_GM_41_table <- t(plyr::ldply(NC_GM_41_list, rbind))
colnames(NC_GM_41_table) <- NC_GM_41_table[1,] # applying first row as column name (e.g. "NC HH6")
NC_GM_41_table <- NC_GM_41_table[-1,] # removed first row containing data label (e.g. "NC HH6")
NC_GM_41_table[is.na(NC_GM_41_table)] = "" # replace "na" with empty fields
write.csv(NC_GM_41_table, file = "NC_GM_41.table.csv", quote = F)


### plotting per time point

# HH6
  # 3 tissues (placodal - NC - neural)
venn.diagram(list("placodal" = Venn_list$`placodal HH6`, "NC" = Venn_list$`NC HH6`, "neural" = Venn_list$`neural HH6`),
             fill = c("green", "blue", "orange"), alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex = 1.0,
             margin = 0.15, filename = "Venn_HH6.png", imagetype = "png",
             scaled = T)
HH6_list <- list("placodal" = setdiff(Venn_list$`placodal HH6`, union(Venn_list$`NC HH6`, Venn_list$`neural HH6`)),
                 "NC" = setdiff(Venn_list$`NC HH6`, union(Venn_list$`placodal HH6`, Venn_list$`neural HH6`)),
                   "neural" = setdiff(Venn_list$`neural HH6`, union(Venn_list$`placodal HH6`, Venn_list$`NC HH6`)),
                 "placodal + NC" = setdiff(intersect(Venn_list$`placodal HH6`, Venn_list$`NC HH6`), Venn_list$`neural HH6`),
                 "placodal + neural" = setdiff(intersect(Venn_list$`placodal HH6`, Venn_list$`neural HH6`), Venn_list$`NC HH6`),
                 "NC + neural" = setdiff(intersect(Venn_list$`NC HH6`, Venn_list$`neural HH6`), Venn_list$`placodal HH6`),
                 "shared by all" = intersect(intersect(Venn_list$`placodal HH6`,Venn_list$`NC HH6`),Venn_list$`neural HH6`))
HH6_table <- t(plyr::ldply(HH6_list, rbind))
colnames(HH6_table) <- HH6_table[1,] # applying first row as column name (e.g. "NC HH6")
HH6_table <- HH6_table[-1,] # removed first row containing data label (e.g. "NC HH6")
HH6_table[is.na(HH6_table)] = "" # replace "na" with empty fields
write.csv(HH6_table, file = "HH6.table.csv", quote = F)


  # NC modules
venn.diagram(list("NC GM 40" = Venn_list$`NC_GM_40 HH6`, "NC GM 41" = Venn_list$`NC_GM_41 HH6`),
             fill = c("green", "blue"), alpha = c(0.5, 0.5), cat.cex = 1.5, cex = 1.0,
             margin = 0.15, filename = "Venn_NC_GM_HH6.png", imagetype = "png",
             scaled = T)
NC_GM_HH6_list <- list('GM 40' = setdiff(Venn_list$`NC_GM_40 HH6`, Venn_list$`NC_GM_41 HH6`),
                      'overlap' = intersect(Venn_list$`NC_GM_41 HH6`,Venn_list$`NC_GM_40 HH6`),
                      'GM 41'= setdiff(Venn_list$`NC_GM_41 HH6`, Venn_list$`NC_GM_40 HH6`))
NC_GM_HH6_table <- t(plyr::ldply(NC_GM_HH6_list, rbind))
colnames(NC_GM_HH6_table) <- NC_GM_HH6_table[1,] # applying first row as column name (e.g. "NC HH6")
NC_GM_HH6_table <- NC_GM_HH6_table[-1,] # removed first row containing data label (e.g. "NC HH6")
NC_GM_HH6_table[is.na(NC_GM_HH6_table)] = "" # replace "na" with empty fields
write.csv(NC_GM_HH6_table, file = "NC_GM_HH6.table.csv", quote = F)

# HH9
  # 3 tissues (placodal - NC - neural)
venn.diagram(list("placodal" = Venn_list$`placodal HH9`, "NC" = Venn_list$`NC HH9`, "neural" = Venn_list$`neural HH9`),
             fill = c("green", "blue", "orange"), alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex = 1.0,
             margin = 0.15, filename = "Venn_HH9.png", imagetype = "png",
             scaled = T)
HH9_list <- list("placodal" = setdiff(Venn_list$`placodal HH9`, union(Venn_list$`NC HH9`, Venn_list$`neural HH9`)),
                 "NC" = setdiff(Venn_list$`NC HH9`, union(Venn_list$`placodal HH9`, Venn_list$`neural HH9`)),
                 "neural" = setdiff(Venn_list$`neural HH9`, union(Venn_list$`placodal HH9`, Venn_list$`NC HH9`)),
                 "placodal + NC" = setdiff(intersect(Venn_list$`placodal HH9`, Venn_list$`NC HH9`), Venn_list$`neural HH9`),
                 "placodal + neural" = setdiff(intersect(Venn_list$`placodal HH9`, Venn_list$`neural HH9`), Venn_list$`NC HH9`),
                 "NC + neural" = setdiff(intersect(Venn_list$`NC HH9`, Venn_list$`neural HH9`), Venn_list$`placodal HH9`),
                 "shared by all" = intersect(intersect(Venn_list$`placodal HH9`,Venn_list$`NC HH9`),Venn_list$`neural HH9`))
HH9_table <- t(plyr::ldply(HH9_list, rbind))
colnames(HH9_table) <- HH9_table[1,] # applying first row as column name (e.g. "NC HH6")
HH9_table <- HH9_table[-1,] # removed first row containing data label (e.g. "NC HH6")
HH9_table[is.na(HH9_table)] = "" # replace "na" with empty fields
write.csv(HH9_table, file = "HH9.table.csv", quote = F)

  # NC modules
venn.diagram(list("NC GM 40" = Venn_list$`NC_GM_40 HH9`, "NC GM 41" = Venn_list$`NC_GM_41 HH9`),
             fill = c("green", "blue"), alpha = c(0.5, 0.5), cat.cex = 1.5, cex = 1.0, 
             margin = 0.15, filename = "Venn_NC_GM_HH9.png", imagetype = "png", 
             scaled = T) 
NC_GM_HH9_list <- list('GM 40' = setdiff(Venn_list$`NC_GM_40 HH9`, Venn_list$`NC_GM_41 HH9`),
                       'overlap' = intersect(Venn_list$`NC_GM_41 HH9`,Venn_list$`NC_GM_40 HH9`),
                       'GM 41'= setdiff(Venn_list$`NC_GM_41 HH9`, Venn_list$`NC_GM_40 HH9`))
NC_GM_HH9_table <- t(plyr::ldply(NC_GM_HH9_list, rbind))
colnames(NC_GM_HH9_table) <- NC_GM_HH9_table[1,] # applying first row as column name (e.g. "NC HH6")
NC_GM_HH9_table <- NC_GM_HH9_table[-1,] # removed first row containing data label (e.g. "NC HH6")
NC_GM_HH9_table[is.na(NC_GM_HH9_table)] = "" # replace "na" with empty fields
write.csv(NC_GM_HH9_table, file = "NC_GM_HH9.table.csv", quote = F)

# all (3 significant) NC modules
venn.diagram(list("NC GM 1 - HH6" = Venn_list$`NC_GM_40 HH6`, "NC GM 2 - HH6" = Venn_list$`NC_GM_41 HH6`,
                  "NC GM 2 - HH9" = Venn_list$`NC_GM_41 HH9`),
             fill = c("green", "blue", "red"), alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex = 1.0, 
             margin = 0.15, cat.default.pos = "outer", cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085), filename = "NC GM diagram.png", imagetype = "png", 
             scaled = T) 

Venn_table <- t(plyr::ldply(Venn_list, rbind)) # Venn list created in line 954 and above
colnames(Venn_table) <- Venn_table[1,] # applying first row as column name (e.g. "NC HH6")
Venn_table <- Venn_table[-1,] # removed first row containing data label (e.g. "NC HH6")
Venn_table[is.na(Venn_table)] = "" # replace "na" with empty fields


# create table showing enhancers in all tissues at time points

#install.packages('data.table')
#install.packages('formattable')
library(data.table)
library(formattable)
library(gt)
library(tidyverse)
library(dplyr)
   

# create and store simple .txt file ###############################
write.table(Venn_table, file = "enhancers_file.txt", sep = " ", 
            row.names = T, col.names = T, quote = F)

# turn Venn_table (matrix) into data frame to create nicer tables
Venn_frame <- as.data.frame(Venn_table, row.names = T)
Venn_frame_1 <- Venn_frame %>% add_column(no = 1:144, .after = 0)

# create a (slightly) nicer looking table and export it as png
formattable(Venn_frame_1)

# a more useful file format though
write.csv(Venn_frame_1, file = "all_tissues.csv", quote = F)


######################################################################################################

### create heatmaps of q-values of motifs per tissue

#install.packages('pheatmap')
#install.packages('purrr')
#install.packages('viridis')
#install.packages('gplots')
#install.packages('dendextend')
library(pheatmap)
library(purrr)
library(viridis)
library(gplots)
library(dendextend)
library(tidyverse)
library(dplyr)

## start with reading data from homer output (see line 861ff)

# create list of dataframes containing motifs and q-values by tissue and time point
the_list <- list(cutoff_placodal_list_HH6,cutoff_NC_list_HH6, cutoff_neural_list_HH6,
                 cutoff_placodal_list_HH9, cutoff_NC_list_HH9, cutoff_neural_list_HH9)
names(the_list) <- c('placodal HH6','NC HH6', 'neural HH6',
                     'placodal HH9', 'NC HH9', 'neural HH9')

# clean column names
colnames <- c("motif", "q-value")
the_list <- lapply(the_list, setNames, colnames)

# merge the lists using motif names as rows and q-values of tissue as columns
joined_list <- the_list %>% reduce(full_join, by = "motif")

# turn motif names into row names
joined_list <- column_to_rownames(joined_list, "motif")

# change n/a values into "1" (or other value?)
joined_list[is.na(joined_list)] = 1

# change column names to tissue
colnames(joined_list) <- names(the_list)

# clean up gene names
rownames(joined_list) <- sub("/Homer", "", row.names(joined_list))
rownames(joined_list) <- sub(".Flag-ChIP-Se", "", row.names(joined_list))

# use dendextend for additional plotting - not for now
#as.dendrogram(joined_list) %>% plot(horiz = T) # "no applicable method for 'as.dendrogram' applied to an object of class "data.frame""

# look at q-value distribution in histogram
hist(-log10(unlist(joined_list)))
hist(unlist(joined_list))

# -log10 the q-values before plotting the heatmap
log_list <- -log10(joined_list)
#log_list[which(!is.finite(log_list))] <- 0 # produces error: "default method not implemented for type 'list'" - a bug b/c it is a df?

# change "inf" into 0 and "na" into 10?
log_list[is.na(log_list)] = 10
#log_list[is.infinite(log_list)] = 0 # produces error "default method not implemented for type 'list'"
log_list <- do.call(data.frame, lapply(log_list, function(x) replace(x, is.infinite(x), 0)))
row.names(log_list) <- row.names(joined_list) # do.call deletes the row names... 
colnames(log_list) <- colnames(joined_list) # ... as well as the col names

# histogram of q-values for log_list
hist(unlist(log_list))

# produce heatmap - no annotations for now
pheatmap::pheatmap(log_list, color = viridis(50), fontsize_row =  2, angle_col = 45 ,filename = "global_heatmap.png")

