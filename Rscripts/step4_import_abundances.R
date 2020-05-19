#Rscript of RNA velocity workflow (last edited 5.18.2020)
#Step 4. Import abundances into R with tximeta

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
  library(scran)
  library(Rtsne)
})

#first, we load the linked transcriptome we created in step 2
tximeta::loadLinkedTxome("gencode.v34.annotation.expanded.json")

#read alevin output
txi <- tximeta::tximeta(coldata = data.frame(
  names = "BMMC_D1T1",
  files = "//nasgw.hpc.vai.org/projects_secondary/triche/Pamela/RNA_velocity_doc/alevin_out/alevin/quants_mat.gz", 
  stringsAsFactors = FALSE
), type = "alevin")

#split counts (with splitSE) into two matrices, one with spliced and one with unspliced abundances, with corresponding rows
cg <- read.delim("gencode.v34.annotation.expanded.features.tsv",
                 header = TRUE, as.is = TRUE)
# Rename the 'intron' column 'unspliced' to make assay names compatible with scVelo
colnames(cg)[colnames(cg) == "intron"] <- "unspliced"
txis <- tximeta::splitSE(txi, cg, assayName = "counts")

#convert txis to a SingleCellExperiment object 
txis <- as(txis, "SingleCellExperiment")
assays(txis) <- list(
  counts = assay(txis, "spliced"),
  spliced = assay(txis, "spliced"),
  unspliced = assay(txis, "unspliced")
)

#removing cells with low gene counts and removing genes that are low across cells
qcstats <- perCellQCMetrics(txis)
qcfilter <- quickPerCellQC(qcstats)
txis <- txis[,!qcfilter$discard]
summary(qcfilter$discard)

#normalize
clusters <- quickCluster(txis)
txis <- computeSumFactors(txis, clusters = clusters)
txis <- scater::logNormCounts(txis)
txis <- scater::runPCA(txis)
txis <- scater::runTSNE(txis, dimred = "PCA")

#save sce object as RDS
saveRDS(txis, "BMMC_D1T1_txi_alevin_abundance.rds")

#class: SingleCellExperiment 
#dim: 60289 3137 
#metadata(6): tximetaInfo quantInfo ... txomeInfo txdbInfo
#assays(4): counts spliced unspliced logcounts
#rownames(60289): ENSG00000223972.5 ENSG00000243485.5 ... ENSG00000210194.1 ENSG00000210196.2
#rowData names(0):
#colnames(3137): GTCAAACTCATGACAC ATAGAGAGTTTGCAGT ... CAGCAATAGTCACTGT CGACAGCTCTGCGGAC
#colData names(1): sizeFactor
#reducedDimNames(2): PCA TSNE
#altExpNames(0):

print(sum(assay(txis, "spliced"))) #32361675
print(sum(assay(txis, "unspliced"))) #13566617

sessionInfo()
#R version 4.0.0 (2020-04-24)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)

#Matrix products: default

#Random number generation:
#RNG:     Mersenne-Twister 
#Normal:  Inversion 
#Sample:  Rounding 

#locale:
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
#[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

#attached base packages:
#[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] Rtsne_0.15                  scran_1.17.0                scater_1.17.0               ggplot2_3.3.0              
#[5] SingleCellExperiment_1.11.1 rjson_0.2.20                tximeta_1.7.3               SummarizedExperiment_1.19.0
#[9] DelayedArray_0.15.0         matrixStats_0.56.0          GenomicFeatures_1.41.0      AnnotationDbi_1.51.0       
#[13] Biobase_2.49.0              eisaR_1.1.0                 BSgenome_1.57.0             rtracklayer_1.49.0         
#[17] GenomicRanges_1.41.0        GenomeInfoDb_1.25.0         Biostrings_2.57.0           XVector_0.29.0             
#[21] IRanges_2.23.0              S4Vectors_0.27.0            BiocGenerics_0.35.0        

#loaded via a namespace (and not attached):
#[1] ggbeeswarm_0.6.0              colorspace_1.4-1              ellipsis_0.3.0                depmixS4_1.4-2               
#[5] BiocNeighbors_1.7.0           rstudioapi_0.11               bit64_0.9-7                   interactiveDisplayBase_1.27.0
#[9] fansi_0.4.1                   tximport_1.17.0               jsonlite_1.6.1                Rsamtools_2.5.0              
#[13] dbplyr_1.4.3                  shiny_1.4.0.2                 BiocManager_1.30.10           compiler_4.0.0               
#[17] httr_1.4.1                    dqrng_0.2.1                   assertthat_0.2.1              Matrix_1.2-18                
#[21] fastmap_1.0.1                 lazyeval_0.2.2                limma_3.45.0                  cli_2.0.2                    
#[25] later_1.0.0                   BiocSingular_1.5.0            htmltools_0.4.0               prettyunits_1.1.1            
#[29] tools_4.0.0                   igraph_1.2.5                  rsvd_1.0.3                    gtable_0.3.0                 
#[33] glue_1.4.0                    GenomeInfoDbData_1.2.3        dplyr_0.8.5                   rappdirs_0.3.1               
#[37] Rcpp_1.0.4.6                  vctrs_0.2.4                   nlme_3.1-147                  DelayedMatrixStats_1.11.0    
#[41] stringr_1.4.0                 beachmat_2.5.0                mime_0.9                      lifecycle_0.2.0              
#[45] irlba_2.3.3                   ensembldb_2.13.1              statmod_1.4.34                XML_3.99-0.3                 
#[49] AnnotationHub_2.21.0          edgeR_3.31.0                  zlibbioc_1.35.0               MASS_7.3-51.5                
#[53] scales_1.1.1                  hms_0.5.3                     promises_1.1.0                ProtGenerics_1.21.0          
#[57] AnnotationFilter_1.13.0       yaml_2.2.1                    curl_4.3                      gridExtra_2.3                
#[61] memoise_1.1.0                 biomaRt_2.45.0                stringi_1.4.6                 RSQLite_2.2.0                
#[65] BiocVersion_3.12.0            BiocParallel_1.23.0           truncnorm_1.0-8               rlang_0.4.6                  
#[69] pkgconfig_2.0.3               bitops_1.0-6                  Rsolnp_1.16                   lattice_0.20-41              
#[73] purrr_0.3.4                   GenomicAlignments_1.25.0      bit_1.1-15.2                  tidyselect_1.1.0             
#[77] magrittr_1.5                  R6_2.4.1                      DBI_1.1.0                     pillar_1.4.4                 
#[81] withr_2.2.0                   RCurl_1.98-1.2                nnet_7.3-13                   tibble_3.0.0                 
#[85] crayon_1.3.4                  BiocFileCache_1.13.0          viridis_0.5.1                 progress_1.2.2               
#[89] locfit_1.5-9.4                grid_4.0.0                    blob_1.2.1                    digest_0.6.25                
#[93] xtable_1.8-4                  httpuv_1.5.2                  openssl_1.4.1                 munsell_0.5.0                
#[97] viridisLite_0.3.0             beeswarm_0.2.3                vipor_0.4.5                   askpass_1.1  