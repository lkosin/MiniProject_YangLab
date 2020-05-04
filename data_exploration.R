# Data exploration.

# Load packages.
library(tidyverse)

# Uncomment if not installed.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BSgenome")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("PFAM.db")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("AnnotationHub")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Organism.dplyr")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Hs.eg.db")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("igvR")
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(PFAM.db)
library(AnnotationHub)
library(Organism.dplyr)
library(org.Hs.eg.db)
library(MASS)
library(igvR)
library(lme4)

# Storing the genome.
hg19 <- BSgenome.Hsapiens.UCSC.hg19
hg19

# Getting chromosome sequences.
# chr1 <- hg19$chr1
# chr1
# chr2 <- hg19$chr2
# chr2
# chr3 <- hg19$chr3
# chr3
# chr4 <- hg19$chr4
# chr4
# chr5 <- hg19$chr5
# chr5
# chr6 <- hg19$chr6
# chr6
# chr7 <- hg19$chr7
# chr7
# chr8 <- hg19$chr8
# chr8
# chr9 <- hg19$chr9
# chr9
# chr10 <- hg19$chr10
# chr10
# chr11 <- hg19$chr11
# chr11
# chr12 <- hg19$chr12
# chr12
# chr13 <- hg19$chr13
# chr13
# chr14 <- hg19$chr14
# chr14
# chr15 <- hg19$chr15
# chr15
# chr16 <- hg19$chr16
# chr16
# chr17 <- hg19$chr17
# chr17
# chr18 <- hg19$chr18
# chr18
# chr19 <- hg19$chr19
# chr19
# chr20 <- hg19$chr20
# chr20
# chr21 <- hg19$chr21
# chr21
# chr22 <- hg19$chr22
# chr22
# chrX <- hg19$chrX
# chrX
# chrY <- hg19$chrY
# chrY
# chrM <- hg19$chrM
# chrM


# Load cancer mutation data.
mut.data <- read_tsv(file = "~/not_backed_up/MiniProject/test.assignment.txt",
                     col_names = F, col_types = "cdccc")
mut.data
colnames(mut.data) <- c("chr", "loc", "ref", "alt", "SampleID")
mut.data

# Global variables.
window.updown <- 10000 # How far to go upstream/downstream from the starting position for a sequence window.

# Getting mutation density and GC content for 10kb up and downstream of each mutation.
gc.content <- function(start.loc, stop.loc, chr.seq){
  require(Biostrings)
  gc.window <- chr.seq[start.loc: stop.loc]
  nts.window <- alphabetFrequency(gc.window)[c("A", "C", "G", "T")]
  nts.total <- sum(nts.window)
  nts.gc <- sum(nts.window[c("G", "C")])
  gc.percent <- nts.gc / nts.total
  return(gc.percent)
}
mut.gc.content(start.loc = (mut.data$loc[1] - 10000), stop.loc = (mut.data$loc[1] + 10000), chr.seq = chr1)

mut.density <- function(start.loc, stop.loc, loc.list){
  mut.in.window <- loc.list[loc.list >= start.loc & loc.list <= stop.loc]
  return(length(mut.in.window))
}
mut.density(start.loc = (mut.data$loc[1] - window.updown), stop.loc = (mut.data$loc[1] + window.updown),
            loc.list = mut.data[mut.data$chr == 1 & mut.data$SampleID == "TCGA-AA-3516-01A-02D-1554-10",]$loc)

# Checking how long this takes to run on chromosome 1.
mut.data.copy <- mut.data
mut.data$GC <- NA
mut.data$m.density <- NA

start.time <- Sys.time()
mut.data.sample1.chr1 <- mut.data[mut.data$chr == 1 & mut.data$SampleID == "TCGA-AA-3516-01A-02D-1554-10",]
for (i in 1:length(mut.data.sample1.chr1$chr)) {
  mut.data.sample1.chr1$GC[i] <- gc.content(start.loc = (mut.data.sample1.chr1$loc[i] - window.updown),
                                            stop.loc = (mut.data.sample1.chr1$loc[i] + window.updown),
                                            chr.seq = chr1)
  mut.data.sample1.chr1$m.density[i] <- mut.density(start.loc = (mut.data.sample1.chr1$loc[i] - window.updown),
                                                 stop.loc = (mut.data.sample1.chr1$loc[i] + window.updown),
                                                 loc.list = mut.data.sample1.chr1$loc)
}
stop.time <- Sys.time()
stop.time - start.time
mut.data.sample1.chr1

cor.test(mut.data.sample1.chr1$GC, mut.data.sample1.chr1$m.density, method = "spearman")

# Looking at other genomic features.
supportedUCSCtables(genome = "hg19")
hg19.txdb <- makeTxDbFromUCSC(genome = "hg19", tablename = "knownGene")
seqlevels(hg19.txdb) <- "chr1"
hg19.chr1.transcripts <- transcripts(hg19.txdb)
hg19.chr1.transcripts

hg19.chr1.exons <- exons(hg19.txdb)
hg19.chr1.exons

hg19.chr1.promoters <- promoters(hg19.txdb)
hg19.chr1.promoters

hg19.chr1.cds <- cds(hg19.txdb)
hg19.chr1.cds

hg19.chr1.introns <- intronicParts(hg19.txdb)
hg19.chr1.introns

seqlevels(hg19.txdb) <- seqlevels0(hg19.txdb)

# Looking for other features on annotation hub.
ahub <- AnnotationHub::AnnotationHub()
ahub.gr.hs.ucsc.hg19 <- ahub[ahub$rdataclass == "GRanges" & ahub$species == "Homo sapiens" & ahub$dataprovider == "UCSC" & ahub$genome == "hg19",]
ahub.gr.hs.ucsc.hg19$title # Microsatellite is 114
ahub.gr.hs.ucsc.hg19[[114]]
hg19.microsatellites.gr <- ahub.gr.hs.ucsc.hg19[[114]]

# Checking Organism.dplyr
# src.hg19 <- src_organism("TxDb.Hsapiens.UCSC.hg19.knownGene")
# src.hg19
# colnames(src.hg19)

# Splitting genome into non-overlapping windows.
chr1.length <- length(chr1)
names(chr1.length) <- "chr1"
chr1.length

chr1.bins <- tileGenome(chr1.length, tilewidth = 1e5, cut.last.tile.in.chrom = T)
chr1.bins

start(chr1.bins)
end(chr1.bins)

# Checking features via bins.
bin.features.chr1 <- tibble(
  "id" = c(1:length(chr1.bins)),
  "start" = start(chr1.bins),
  "end" = end(chr1.bins),
  "mut.count" = rep(NA, length(chr1.bins)),
  "GC.content" = rep(NA, length(chr1.bins))
)
bin.features.chr1

start.time <- Sys.time()
for (i in 1:length(bin.features.chr1$id)) {
  bin.features.chr1$GC.content[i] <- gc.content(start.loc = bin.features.chr1$start[i],
                                                stop.loc = bin.features.chr1$end[i],
                                                chr.seq = chr1)
  bin.features.chr1$mut.count[i] <- mut.density(start.loc = bin.features.chr1$start[i],
                                                stop.loc = bin.features.chr1$end[i],
                                                loc.list = mut.data.sample1.chr1$loc)
}
stop.time <- Sys.time()
stop.time - start.time

bin.features.chr1
bin.chr1.lm <- lm(
  data = bin.features.chr1,
  formula = mut.count ~ GC.content
)
summary(bin.chr1.lm)

# Scaling up for all chromosomes.
# Selecting the main chromosomes. Note that I include chrM here, and then exclude it
# later in the analysis.
seqlevels(hg19)
chromosomes.norm <- c("chr1", "chr2", "chr3", "chr4", "chr5",
                      "chr6", "chr7", "chr8", "chr9", "chr10",
                      "chr11", "chr12", "chr13", "chr14", "chr15",
                      "chr16", "chr17", "chr18", "chr19", "chr20",
                      "chr21", "chr22", "chrX", "chrY", "chrM")
hg19.lengths <- seqlengths(hg19)
hg19.bins <- tileGenome(hg19.lengths, tilewidth = 1e5, cut.last.tile.in.chrom = T)
hg19.bins

bin.features.all <- tibble(
  "id" = c(1:length(hg19.bins)),
  "chr" = as.vector(seqnames(hg19.bins)),
  "start" = start(hg19.bins),
  "end" = end(hg19.bins),
  "mut.count" = rep(NA, length(hg19.bins)),
  "GC.content" = rep(NA, length(hg19.bins))
)
bin.features.all
tail(bin.features.all)

mut.data$chr <- paste("chr", mut.data$chr, sep = "")
mut.data[mut.data$chr == "chrMT", ]$chr <- "chrM"
mut.data[mut.data$chr == "chrM",]

#####
start.time <- Sys.time()
for (i in 1:length(bin.features.all$id)) {
  if (bin.features.all$chr[i] %in% chromosomes.norm){
    bin.features.all$mut.count[i] <- mut.density(start.loc = bin.features.all$start[i],
                                                 stop.loc = bin.features.all$end[i],
                                                 loc.list = mut.data[mut.data$chr == bin.features.all$chr[i],]$loc)
    bin.features.all$GC.content[i] <- gc.content(start.loc = bin.features.all$start[i],
                                                 stop.loc = bin.features.all$end[i],
                                                 chr.seq = get(bin.features.all$chr[i]))
  }
}
stop.time <- Sys.time()
stop.time - start.time

bin.features.all


# The relationship between mutation density and GC content is clearly non-monotonic. It looks like the switch is about at 0.4
# What if we divide the data in two?
with(bin.features.all[bin.features.all$GC.content < 0.4,],
     cor.test(mut.count, GC.content, method = "pearson"))
with(bin.features.all[bin.features.all$GC.content >= 0.4,],
     cor.test(mut.count, GC.content, method = "pearson"))

# Organizing the data for other features.
hg19.transcripts <- transcripts(hg19.txdb)
hg19.transcripts
seqnames(hg19.transcripts)

hg19.cds <- cds(hg19.txdb)
hg19.cds

# Making a function to convert an GRanges annotation object to a more easily parseable data frame.
# gr.dataframe <- function(gr.object){
#   require(tidyverse)
#   n.gr <- length(gr.object)
#   gr.df <- tibble(
#     "id" = c(1:n.gr),
#     "chr" = as.vector(seqnames(gr.object)),
#     "start" = start(gr.object),
#     "end" = end(gr.object),
#     "strand" = as.vector(strand(gr.object))
#   )
#   return(gr.df)
# }
# hg19.cds.df <- gr.dataframe(hg19.cds)
# hg19.transcripts.df <- gr.dataframe(hg19.transcripts)
# hg19.cds.df
# hg19.transcripts.df

# Making a function to score as 1 or 0 whether a window contains a feature in an annotation list.
# window.feature <- function(start.loc, stop.loc, feature.start.list, feature.stop.list){
#   feature.start.in.window <- feature.start.list[feature.start.list >= start.loc & feature.start.list <= stop.loc]
#   feature.stop.in.window <- feature.stop.list[feature.stop.list >= stop.loc & feature.stop.list <= stop.loc]
#   in.window <- ifelse(length(feature.start.in.window) >= 1 | length(feature.stop.in.window) >= 1, 1, 0)
#   return(in.window)
# }

# Making a function to take a data frame of annotations and score whether a window has the annotation in question.
# INCOMPLETE: SEEMS UNNECESSARY BECAUSE OVERLAPPING REGIONS FUNCTIONS IN GRANGES ALREADY EXIST.
# window.annotate <- function(window.df, feature.df, feature.col, chr.list, chr.col = "chr", start.col = "start", stop.col = "end"){
#   # window.df is a data frame for the windows in question. Needs the window start and stop location, chr, and feature column.
#   # feature.df is the data frame for the feature in question. Needs start and stop location and chr.
#   # feature.col is the column name of the feature in the window.df data frame.
#   # chr.list is the list of chromosomes to restrict the annotation to.
#   # chr.col is the name of the chromosome column in the window.df data frame AND the feature.df data frame.
#   # NOTE: The column names and chromosome identifications need to match between the two data frames!
#   # NOTE: The column names for the start and end for each feature need to be the same!
#   n.window <- nrow(window.df)
#   for (i in 1:n.window) {
#     if (window.df[i, chr.col] %in% chr.list){
#       window.df[i, feature.col] <- window.feature(
#         start.loc = window.df[i, start.col][[1]],
#         stop.loc = window.df[i, stop.col][[1]],
#         feature.start.list = feature.df
#       )
#     }
#   }
# }

hg19.cds.distance <- distanceToNearest(hg19.bins, hg19.cds, ignore.strand = T)
hg19.transcripts.distance <- distanceToNearest(hg19.bins, hg19.transcripts, ignore.strand = T)

mcols(hg19.cds.distance)[[1]]

tail(bin.features.all[bin.features.all$chr %in% chromosomes.norm,])
length(hg19.bins)

# The lengths don't match the total length in hg19.bins.
# hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm]
# distanceToNearest(hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm], hg19.cds, ignore.strand = T)
# 
# hg19.bins[seqnames(hg19.bins) == "chrM"]
# distanceToNearest(hg19.bins[seqnames(hg19.bins) == "chrM"], hg19.cds, ignore.strand = T)
# unique(seqnames(hg19.cds))
# Some of the seqnames in hg19.bins are not in hg19.cds. OK.
# chromosomes.norm[1:24]
# 
# hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm[1:24]]
# distanceToNearest(hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm[1:24]], hg19.cds, ignore.strand = T)

# Adding in other features.
bin.features.all$cds.distance <- NA
bin.features.all$transcript.distance <- NA

hg19.cds.distance <- distanceToNearest(hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm[1:24]], hg19.cds, ignore.strand = T)
hg19.transcripts.distance <- distanceToNearest(hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm[1:24]], hg19.transcripts, ignore.strand = T)
hg19.microsatellites.distance <- distanceToNearest(hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm[1:24]],
                                                   hg19.microsatellites.gr, ignore.strand = T)

bin.features.normchr <- bin.features.all[bin.features.all$chr %in% chromosomes.norm[1:24],]
bin.features.normchr

bin.features.normchr$cds.distance <- mcols(hg19.cds.distance)[[1]]
bin.features.normchr$transcript.distance <- mcols(hg19.transcripts.distance)[[1]]
bin.features.normchr$micro.distance <- mcols(hg19.microsatellites.distance)[[1]]

bin.features.normchr$cds.binary <- ifelse(bin.features.normchr$cds.distance == 0, 1, 0)
bin.features.normchr$transcript.binary <- ifelse(bin.features.normchr$transcript.distance == 0, 1, 0)
bin.features.normchr$micro.binary <- ifelse(bin.features.normchr$micro.distance == 0, 1, 0)
  
bin.features.normchr

# Checking distributions.
hist(bin.features.all$mut.count)
hist(bin.features.all$GC.content)
boxcox(bin.features.all$mut.count + 0.5 ~ 1, lambda = seq(0.8, 1, by = 0.01))
boxcox(bin.features.all$GC.content ~ 1, lambda = seq(-3, -2, by = 0.01))
hist((bin.features.all$mut.count + 0.5) ^ 0.9)
png(filename = "mut_vs_GC.png", width = 480, height = 480)
ggplot(data = bin.features.normchr,
       aes(x = GC.content * 100, y = mut.count)) +
  geom_point() +
  geom_smooth() +
  ylab("Mutation density") +
  xlab("%GC content") +
  theme_bw(base_size = 28)
dev.off()

with(bin.features.normchr,
     cor.test(mut.count, cds.binary, method = "spearman"))
with(bin.features.normchr,
     cor.test(mut.count, transcript.binary, method = "spearman"))
with(bin.features.normchr,
     cor.test(mut.count, cds.distance, method = "spearman"))
with(bin.features.normchr,
     cor.test(mut.count, transcript.distance, method = "spearman"))

# Getting other features.
# Histones
igv <- igvR()
#setGenome(igv, "hg19")
histone.tracks <- AnnotationHub::query(subset(ahub, species == "Homo sapiens"), c("H3K4me3", "hg19", "Peak", "narrow"))
histone.tracks
length(histone.tracks)
histone.tracks[[c("AH25833")]]
histon.anno <- histone.tracks[[c("AH25833")]]
hist.name <- names(histone.tracks)
hist.name[1]
for (i in 2:length(histone.tracks)) {
  histon.anno <- c(histon.anno, histone.tracks[[hist.name[i]]])
}

hg19.normchrom <- hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm[1:24]]
hist.distance <- distanceToNearest(hg19.normchrom, histon.anno)
mcols(hist.distance)
bin.features.normchr$hist.dist <- mcols(hist.distance)[[1]]
bin.features.normchr$hist.binary <- ifelse(bin.features.normchr$hist.dist == 0, 1, 0)


# Checking lincRNAs.
lincRNA.data <- read_tsv(file = "~/not_backed_up/MiniProject/lincRNAsCTColon.txt", col_names = F)
lincRNA.data
names(lincRNA.data) <- c("bin", "chr", "start", "end", "id", "score", "rawscore", "log2rawscore")
lincRNA.data

lincRNA.gr <- GRanges(
  seqnames = factor(lincRNA.data$chr),
  ranges = IRanges(start = lincRNA.data$start, end = lincRNA.data$end),
  score = lincRNA.data$score
)
lincRNA.gr
lincRNA.gr[21472:21626]

lincRNA.distance <- distanceToNearest(hg19.normchrom, lincRNA.gr)
mcols(lincRNA.distance)

bin.features.normchr$lincRNA.distance <- mcols(lincRNA.distance)[[1]]
bin.features.normchr$lincRNA.binary <- ifelse(bin.features.normchr$lincRNA.distance == 0, 1, 0)

bin.features.normchr

# Building linear models to check associations.
bin.all.lm <- lm(
  data = bin.features.normchr,
  formula = mut.count ~ GC.content + cds.binary + micro.binary + hist.binary + lincRNA.binary
)
summary(bin.all.lm)

bin.lowgc.lm <- lm(
  data = bin.features.normchr[bin.features.normchr$GC.content < 0.4,],
  formula = mut.count ~ GC.content + cds.binary + micro.binary + hist.binary + lincRNA.binary
)
summary(bin.lowgc.lm)

bin.highgc.lm <- lm(
  data = bin.features.normchr[bin.features.normchr$GC.content >= 0.4,],
  formula = mut.count ~ GC.content + cds.binary + micro.binary + hist.binary + lincRNA.binary
)
summary(bin.highgc.lm)

# Adding chromosome as a random effect.
bin.lowgc.lme <- lmer(
  data = bin.features.normchr[bin.features.normchr$GC.content < 0.4,],
  formula = mut.count ~ GC.content + cds.binary + micro.binary + hist.binary + lincRNA.binary + (1|chr)
)
summary(bin.lowgc.lme)

bin.highgc.lme <- lmer(
  data = bin.features.normchr[bin.features.normchr$GC.content >= 0.4,],
  formula = mut.count ~ GC.content + cds.binary + micro.binary + hist.binary + lincRNA.binary + (1|chr)
)
summary(bin.highgc.lme)

# Treating patient sample as a random effect.
bin.sample <- rbind(bin.features.normchr, bin.features.normchr, bin.features.normchr, bin.features.normchr, bin.features.normchr)
bin.sample

bin.sample$SampleID <- NA
bin.sample$mut.count <- NA

##################################################
# WARNING: THE FOLLOWING CODE TAKES HOURS TO RUN #
##################################################
for (i in 1:length(unique(mut.data$SampleID))) {
  for (j in 1:length(bin.features.normchr$id)) {
    bin.sample$SampleID[(length(bin.features.normchr$id) * (i - 1)) + j] <- unique(mut.data$SampleID)[i]
    bin.sample$mut.count[(length(bin.features.normchr$id) * (i - 1)) + j] <- mut.density(
      start.loc = bin.sample$start[(length(bin.features.normchr$id) * (i - 1)) + j],
      stop.loc = bin.sample$end[(length(bin.features.normchr$id) * (i - 1)) + j],
      loc.list = mut.data[mut.data$SampleID == unique(mut.data$SampleID)[i] & 
                            mut.data$chr == bin.sample$chr[(length(bin.features.normchr$id) * (i - 1)) + j],]$loc
    )
  }
}
bin.sample
unique(bin.sample$SampleID)

bin.lowgc.sample.lme <- lmer(
  data = bin.sample[bin.sample$GC.content < 0.4,],
  formula = mut.count ~ GC.content + cds.binary + transcript.binary + hist.binary + lincRNA.binary + (1|chr) + (1|SampleID)
)
summary(bin.lowgc.sample.lme)

bin.highgc.sample.lme <- lmer(
  data = bin.sample[bin.sample$GC.content >= 0.4,],
  formula = mut.count ~ GC.content + cds.binary + transcript.binary + hist.binary + lincRNA.binary + (1|chr) + (1|SampleID)
)
summary(bin.highgc.sample.lme)

# Comparing the mixed models vs the fixed effect only models.
bin.features.normchr$fixed.fit <- ifelse(bin.features.normchr$GC.content < 0.4,
                                         predict(bin.lowgc.lm, newdata = bin.features.normchr,
                                                 type = "response"),
                                         predict(bin.highgc.lm, newdata = bin.features.normchr,
                                                 type = "response"))
bin.features.normchr[, c("id", "mut.count", "fixed.fit")]
bin.features.normchr$mixed.fit <- ifelse(
  bin.features.normchr$GC.content < 0.4,
  predict(bin.lowgc.lme, newdata = bin.features.normchr, type = "response", re.form = NA, random.only = F),
  predict(bin.highgc.lme, newdata = bin.features.normchr, type = "response", re.form = NA, random.only = F)
)

cor.test(bin.features.normchr$fixed.fit, bin.features.normchr$mixed.fit, method = "spearman")
summary(lm(
  data = bin.features.normchr,
  formula = mut.count ~ fixed.fit
))
summary(lm(
  data = bin.features.normchr,
  formula = mut.count ~ mixed.fit
))
AIC(bin.lowgc.lm)
AIC(bin.highgc.lm)
AIC(bin.lowgc.lme)
AIC(bin.highgc.lme)

# The AIC is improved, but there does not appear to be any obvious benefit to the mixed
# models for understanding the relationship between mutation density and the genomic
# features examined here. However, note that the models were not properly evaluated for
# their predictive ability, as prediction is not the goal here.