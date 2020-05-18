#Rscript of RNA velocity workflow (last edited 5.15.2020)
#Step 2. Indexing the reference genome

#import packages
suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome)
  library(eisaR)
  library(GenomicFeatures)
  library(SummarizedExperiment)
  library(tximeta)
  library(rjson)
  library(SingleCellExperiment)
  library(scater)
})

#extract gRanges object containing genomic coordinates of each annotated transcript and intron
grl <- eisaR::getFeatureRanges(
  gtf = system.file("extdata/gencode.v34.annotation.gtf.gz", package = "eisaR"),
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 90L, 
  joinOverlappingIntrons = FALSE, 
  verbose = TRUE)

#extract the sequences of all features of interest
genome <- Biostrings::readDNAStringSet("GRCh38.primary_assembly.genome.fa.gz")
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)
Biostrings::writeXStringSet(
  seqs, filepath = "gencode.v34.annotation.expanded.fa"
)

#write expanded annotation to .gtf
eisaR::exportToGtf(
  grl, 
  filepath = "gencode.v34.annotation.expanded.gtf"
)

#write metadata of the gRanges object containing a dataframe consists of spliced and unspliced gene ID
write.table(
  metadata(grl)$corrgene, 
  file = "gencode.v34.annotation.expanded.features.tsv",
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

#map transcripts and introns to the corresponding gene identifiers
df <- eisaR::getTx2Gene(
  grl, filepath = "gencode.v34.annotation.expanded.tx2gene.tsv"
)

#create a linked transcriptome with tximeta
tximeta::makeLinkedTxome(
  indexDir = "gencode.v34.annotation.expanded.sidx", 
  source = "GENCODE", genome = "GRCh38", 
  organism = "Homo sapiens", release = "v34", 
  fasta = "gencode.v34.annotation.expanded.fa", 
  gtf = "gencode.v34.annotation.expanded.gtf", 
  write = TRUE, jsonFile = "gencode.v34.annotation.expanded.json"
)
rjson::fromJSON(file = "gencode.v34.annotation.expanded.json")

