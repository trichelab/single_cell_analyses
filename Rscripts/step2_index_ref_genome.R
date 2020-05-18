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

#finally, this might be useful to check which version of R/packages I'm using
sessionInfo()

#R version 4.0.0 (2020-04-24)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)

#Matrix products: default

#Random number generation:
#  RNG:     Mersenne-Twister 
#Normal:  Inversion 
#Sample:  Rounding 

#locale:
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
#[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

#attached base packages:
#[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] scater_1.17.0               ggplot2_3.3.0               SingleCellExperiment_1.11.1 rjson_0.2.20               
#[5] tximeta_1.7.3               SummarizedExperiment_1.19.0 DelayedArray_0.15.0         matrixStats_0.56.0         
#[9] GenomicFeatures_1.41.0      AnnotationDbi_1.51.0        Biobase_2.49.0              eisaR_1.1.0                
#[13] BSgenome_1.57.0             rtracklayer_1.49.0          GenomicRanges_1.41.0        GenomeInfoDb_1.25.0        
#[17] Biostrings_2.57.0           XVector_0.29.0              IRanges_2.23.0              S4Vectors_0.27.0           
#[21] BiocGenerics_0.35.0        

#loaded via a namespace (and not attached):
#[1] ggbeeswarm_0.6.0              colorspace_1.4-1              ellipsis_0.3.0                depmixS4_1.4-2               
#[5] BiocNeighbors_1.7.0           rstudioapi_0.11               bit64_0.9-7                   interactiveDisplayBase_1.27.0
#[9] fansi_0.4.1                   tximport_1.17.0               jsonlite_1.6.1                Rsamtools_2.5.0              
#[13] dbplyr_1.4.3                  shiny_1.4.0.2                 BiocManager_1.30.10           compiler_4.0.0               
#[17] httr_1.4.1                    assertthat_0.2.1              Matrix_1.2-18                 fastmap_1.0.1                
#[21] lazyeval_0.2.2                limma_3.45.0                  cli_2.0.2                     later_1.0.0                  
#[25] BiocSingular_1.5.0            htmltools_0.4.0               prettyunits_1.1.1             tools_4.0.0                  
#[29] rsvd_1.0.3                    gtable_0.3.0                  glue_1.4.0                    GenomeInfoDbData_1.2.3       
#[33] dplyr_0.8.5                   rappdirs_0.3.1                Rcpp_1.0.4.6                  vctrs_0.2.4                  
#[37] nlme_3.1-147                  DelayedMatrixStats_1.11.0     stringr_1.4.0                 mime_0.9                     
#[41] lifecycle_0.2.0               irlba_2.3.3                   ensembldb_2.13.1              XML_3.99-0.3                 
#[45] AnnotationHub_2.21.0          edgeR_3.31.0                  zlibbioc_1.35.0               MASS_7.3-51.5                
#[49] scales_1.1.1                  hms_0.5.3                     promises_1.1.0                ProtGenerics_1.21.0          
#[53] AnnotationFilter_1.13.0       yaml_2.2.1                    curl_4.3                      gridExtra_2.3                
#[57] memoise_1.1.0                 biomaRt_2.45.0                stringi_1.4.6                 RSQLite_2.2.0                
#[61] BiocVersion_3.12.0            BiocParallel_1.23.0           truncnorm_1.0-8               rlang_0.4.6                  
#[65] pkgconfig_2.0.3               bitops_1.0-6                  Rsolnp_1.16                   lattice_0.20-41              
#[69] purrr_0.3.4                   GenomicAlignments_1.25.0      bit_1.1-15.2                  tidyselect_1.1.0             
#[73] magrittr_1.5                  R6_2.4.1                      DBI_1.1.0                     pillar_1.4.4                 
#[77] withr_2.2.0                   RCurl_1.98-1.2                nnet_7.3-13                   tibble_3.0.0                 
#[81] crayon_1.3.4                  BiocFileCache_1.13.0          viridis_0.5.1                 progress_1.2.2               
#[85] locfit_1.5-9.4                grid_4.0.0                    blob_1.2.1                    digest_0.6.25                
#[89] xtable_1.8-4                  httpuv_1.5.2                  openssl_1.4.1                 munsell_0.5.0                
#[93] viridisLite_0.3.0             beeswarm_0.2.3                vipor_0.4.5                   askpass_1.1