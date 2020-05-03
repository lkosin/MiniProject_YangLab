# Used code.

# Load packages.
library(tidyverse)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationHub)

# Storing the genome.
hg19 <- BSgenome.Hsapiens.UCSC.hg19
hg19

# Load cancer mutation data.
mut.data <- read_tsv(file = "~/not_backed_up/MiniProject/test.assignment.txt",
                     col_names = F, col_types = "cdccc")
mut.data
colnames(mut.data) <- c("chr", "loc", "ref", "alt", "SampleID")
mut.data

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

mut.density <- function(start.loc, stop.loc, loc.list){
  mut.in.window <- loc.list[loc.list >= start.loc & loc.list <= stop.loc]
  return(length(mut.in.window))
}

# Looking at other genomic features.
supportedUCSCtables(genome = "hg19")
hg19.txdb <- makeTxDbFromUCSC(genome = "hg19", tablename = "knownGene")

# Selecting the main chromosomes.
seqlevels(hg19)
chromosomes.norm <- c("chr1", "chr2", "chr3", "chr4", "chr5",
                      "chr6", "chr7", "chr8", "chr9", "chr10",
                      "chr11", "chr12", "chr13", "chr14", "chr15",
                      "chr16", "chr17", "chr18", "chr19", "chr20",
                      "chr21", "chr22", "chrX", "chrY", "chrM")
hg19.lengths <- seqlengths(hg19)
hg19.bins <- tileGenome(hg19.lengths, tilewidth = 1e5, cut.last.tile.in.chrom = T)
hg19.bins

# Making the "chr" column match the chromosome names from UCSC.
mut.data$chr <- paste("chr", mut.data$chr, sep = "")
mut.data[mut.data$chr == "chrMT", ]$chr <- "chrM"
mut.data[mut.data$chr == "chrM",]

# Organizing the data for other features.
hg19.transcripts <- transcripts(hg19.txdb)
hg19.transcripts
seqnames(hg19.transcripts)

hg19.cds <- cds(hg19.txdb)
hg19.cds

# Getting distance to features.
hg19.cds.distance <- distanceToNearest(hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm[1:24]], hg19.cds, ignore.strand = T)
hg19.transcripts.distance <- distanceToNearest(hg19.bins[seqnames(hg19.bins) %in% chromosomes.norm[1:24]], hg19.transcripts, ignore.strand = T)

# Making a data frame to store all the data and run linear analysis.
bin.features.all <- tibble(
  "id" = c(1:length(hg19.bins)),
  "chr" = as.vector(seqnames(hg19.bins)),
  "start" = start(hg19.bins),
  "end" = end(hg19.bins),
  "mut.count" = rep(NA, length(hg19.bins)),
  "GC.content" = rep(NA, length(hg19.bins))
)
bin.features.all

# Parsing this down to just the main chromosomes.
bin.features.normchr <- bin.features.all[bin.features.all$chr %in% chromosomes.norm[1:24],]
bin.features.normchr

bin.features.normchr$cds.distance <- mcols(hg19.cds.distance)[[1]]
bin.features.normchr$transcript.distance <- mcols(hg19.transcripts.distance)[[1]]

bin.features.normchr$cds.binary <- ifelse(bin.features.normchr$cds.distance == 0, 1, 0)
bin.features.normchr$transcript.binary <- ifelse(bin.features.normchr$transcript.distance == 0, 1, 0)

# Checking distributions.
hist(bin.features.all$mut.count)
hist(bin.features.all$GC.content)
boxcox(bin.features.all$mut.count + 0.5 ~ 1, lambda = seq(0.8, 1, by = 0.01))
boxcox(bin.features.all$GC.content ~ 1, lambda = seq(-3, -2, by = 0.01))
hist((bin.features.all$mut.count + 0.5) ^ 0.9)
ggplot(data = bin.features.all,
       aes(x = GC.content, y = (mut.count + 0.5) ^ 0.9)) +
  geom_point() +
  geom_smooth()

# Getting other H3K4me3 histone markers.
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

# Getting colon lincRNAs.
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

# Building linear models to check associations.
bin.all.lm <- lm(
  data = bin.features.normchr,
  formula = mut.count ~ GC.content + cds.binary + transcript.binary + hist.binary + lincRNA.binary
)
summary(bin.all.lm)

bin.lowgc.lm <- lm(
  data = bin.features.normchr[bin.features.normchr$GC.content < 0.4,],
  formula = mut.count ~ GC.content + cds.binary + transcript.binary + hist.binary + lincRNA.binary
)
summary(bin.lowgc.lm)

bin.highgc.lm <- lm(
  data = bin.features.normchr[bin.features.normchr$GC.content >= 0.4,],
  formula = mut.count ~ GC.content + cds.binary + transcript.binary + hist.binary + lincRNA.binary
)
summary(bin.highgc.lm)

# Adding chromosome as a random effect.
bin.lowgc.lme <- lmer(
  data = bin.features.normchr[bin.features.normchr$GC.content < 0.4,],
  formula = mut.count ~ GC.content + cds.binary + transcript.binary + hist.binary + lincRNA.binary + (1|chr)
)
summary(bin.lowgc.lme)

bin.highgc.lme <- lmer(
  data = bin.features.normchr[bin.features.normchr$GC.content >= 0.4,],
  formula = mut.count ~ GC.content + cds.binary + transcript.binary + hist.binary + lincRNA.binary + (1|chr)
)
summary(bin.highgc.lme)

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

# AIC values are better, but the mixed model does not appear to be a better predictor of mutation density.