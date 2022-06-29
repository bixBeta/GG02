#!/usr/bin/env Rscript

args <-  commandArgs(trailingOnly = T)
# check for required argument

if (length(args)==0) {
  print(" Usage = Rscript atacQC.R < human, mouse, fly or 'path to gtf annotation' >")
  stop("Missing gtf annotation path !!! \n", call.=FALSE)

}


library(progress)
pb <- progress_bar$new(total = 100)

for (i in 1:20) {
  pb$tick()
  Sys.sleep(1 / 10)
}


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# fragSize - Distributions

suppressPackageStartupMessages(library(ATACseqQC))
bamfiles <- list.files(pattern = ".bam$")
bamfile.labels <- gsub(".DEDUP.bam", "", basename(bamfiles))

for (i in 1:length(bamfiles)) {
  png(filename = paste0(strsplit(x = bamfiles[i],
                                 split = "\\.DEDUP.bam"[[1]]),".fragSize" ,".png"),
      width = 1080, height = 1080, res = 155)
  fragSizeDist(bamfiles[i], bamfile.labels[i])
  dev.off()
}

for (i in 1:30) {
  pb$tick()
  Sys.sleep(1 / 100)
}

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# Peak Distributions + Peak Annotation
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(clusterProfiler))

for (i in 1:30) {
  pb$tick()
  Sys.sleep(1 / 100)
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# import peak files as list

files <- list.files(path = "peaks.OUT", pattern = ".narrowPeak$", full.names = T)
files <- as.list(files[1:length(files)])
for (i in 1:length(files)) {
  names(files[[i]]) <- gsub(pattern = ".narrowPeak",
                            replacement = "", basename(files[[i]]))
}

if (args[1]=="mouse") {
  gtf.path <-  "/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.96.gtf"
  org <- "Mus musculus"
} else  if (args[1]=="human"){
  gtf.path <-  "/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.96.gtf"
  org <- "Homo sapiens"
} else  if (args[1]=="fly"){
  gtf.path <-  "/workdir/genomes/Drosophila_melanogaster/dm6/ENSEMBL/Drosophila_melanogaster.BDGP6.32.106.gtf"
  org <- "Drosophila melanogaster"
} else {
  gtf.path <- args[1]
  org <- "NA"
}

for (i in 1:20) {
  pb$tick()
  Sys.sleep(1 / 100)
}

txdb <- makeTxDbFromGFF(gtf.path,
                        format = "gtf", organism = org)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# import peak files as list

files <- list.files(path = "peaks.OUT", pattern = ".narrowPeak$", full.names = T)

files <- as.list(files[1:length(files)])

for (i in 1:length(files)) {
  names(files[[i]]) <- gsub(pattern = ".narrowPeak",
                            replacement = "", basename(files[[i]]))
}

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

for (i in 1:length(files)) {
  print(paste0("Processing -- ", names(files[[i]])))
  peak <- readPeakFile(files[[i]])
  tagMatrix <- getTagMatrix(peak, windows=promoter)
  png(filename = paste0(basename(files[[i]]),".Distribution.png"),
      width = 720, height = 640)
  c1 <- plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
                    xlab="Genomic Region (5'->3')",
                    ylab = paste0("Read Count Frequency - ",
                                  basename(files[[i]])))
  print(c1)
  dev.off()

  peakAnno <- annotatePeak(peak, TxDb=txdb)
  png(filename = paste0(basename(files[[i]]),".Annotated.png"),
      width = 720, height = 320)
  c2 <- plotAnnoBar(peakAnno, title = basename(files[[i]]))
  print(c2)
  dev.off()

  write.table(c2$data, file = paste0(c2$labels$title, "__percent__annot__log.txt"), sep = "\t")
}


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# Calculate TSSE scores


txs <- transcripts(txdb)

bamsList = list()

for (i in 1:length(bamfile.labels)) {
  bamsList[[i]] <-  readBamFile(bamfiles[i], bigFile = T)
  names(bamsList)[[i]] <- bamfile.labels[i]
}



tsseScores = lapply(bamsList, function(x){
  tsse = TSSEscore(obj = x, txs = txs)
  score = tsse$TSSEscore
  return(score)

})


tsseScores.df = as.data.frame(round(unlist(tsseScores), 2))
colnames(tsseScores.df) = "value"
tsseScores.df$labels = rownames(tsseScores.df)
library(ggplot2)

png("TSSE_scores.png", width = 1080, height = 1080, res = 100)
ggplot(data = tsseScores.df, aes(x = labels, y = value )) + geom_bar(stat = "identity") +
  geom_text(aes(label=value), vjust=1.6, color="white", size=3.5) + ggtitle("Transcription Start Site (TSS) Enrichment Scores") +
  xlab("Sample") + ylab("TSS Score") + scale_y_continuous() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()





system("mkdir atacQC.out")
system("mv *.png atacQC.out")
system("mv *_peaks.narrowPeak__percent__annot__log.txt atacQC.out")
