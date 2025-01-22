library(devtools)
devtools::install_github("jtlovell/GENESPACE", upgrade = F)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))

library(GENESPACE)

genomeRepo <- "~/Documents/Plasticity_eQTL/wrkdir_genespace/genomeRepo/"
wd <- "~/Documents/Plasticity_eQTL/wrkdir_genespace_w_col0/"
path2mcscanx <- "~/Downloads/MCScanx/MCScanX-master/"

library(data.table)
sha_bed<-fread("wrkdir_genespace/bed/sha_proteins.bed")
sha_bed<-sha_bed[V8 == "mRNA", .SD, .SDcols = c("V1","V2","V3","V10")]
sha_bed<-sha_bed[, V10 := sub("ID=|;.*", "", V10)]
fwrite(sha_bed,"wrkdir_genespace/bed/Sha.bed",sep="\t")

bay_bed<-fread("wrkdir_genespace/bed/bay_proteins.bed")
bay_bed<-bay_bed[V8 == "mRNA", .SD, .SDcols = c("V1","V2","V3","V10")]
bay_bed<-bay_bed[, V10 := sub("ID=|;.*", "", V10)]
bay_bed[, V10 := sub(";.*", "", V10)]
fwrite(bay_bed,"wrkdir_genespace/bed/Bay.bed", sep="\t")
###########################################################################################
genomes2run<-c("Bay-0","Sha")
parsedPaths<-parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = genomes2run,
  genomeIDs = genomes2run,
  genespaceWd = wd
)
###########################################################################################
gpar <- init_genespace(
  wd = "~/Documents/Plasticity_eQTL/wrkdir_genespace_w_col0/", 
  path2mcscanx = "~/Downloads/MCScanx/MCScanX-master/",
  rawOrthofinderDir = "/Users/mlensink/Documents/Plasticity_eQTL/wrkdir_genespace_w_col0/orthofinder/Results_Sep18"
)
out <- run_genespace(gsParam = gpar)



