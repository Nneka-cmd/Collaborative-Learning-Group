#setwd("D:/Genotyping/2023")
#Dart QC filtering and convert to hapmap

.libPaths()
#--------------------------------------------------------
# Script parameters
#--------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

BiocManager::install("SNPRelate")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocStyle")

library(devtools)
install_github("cran/dartR")

remove.packages("ade4")
remove.packages("adegenet")

#Libraries
library(dartR)
library(data.table)
library(dplyr)
library(magrittr)


## File paths ----
#first run this with SNP data - Report_DS22-6886_SNP_mapping_2.csv
#the run with P/A data Report_DS21-6352_SilicoDArT_1.csv
marker_type<-"SNP" #SNP
file       <- "Report_DS24-9190_SNP_mapping_2"
snpPath    <- paste(file,".csv",sep="")
hapMapOut  <- paste(file,".hmp.txt",sep="")


## Filtering criteria ----
thresholds <- list(
  "indCallrate" = 0.50,
  "snpCallrate" = 0.50,
  "snpRepAvg"   = 0.9,
  "secMethod"   = "best"
)

## DartR / adegenet GenLight objects ----
# glSilico <- gl.read.silicodart(silicoPath)
glSnp   <- gl.read.dart(snpPath)

## Identify absent marker data ----
chromV  <- glSnp@other$loc.metrics$Chrom_Sorghum_bicolor_v5.1
dropped <- glSnp$loc.names[grep("^Chr", as.character(chromV), invert = TRUE)]


## Filter 
glSnpSub <- glSnp %>% 
  gl.drop.loc(dropped) %>% # drop missing data ---
  gl.filter.callrate(method = "loc", threshold = thresholds$snpCallrate) %>%
  gl.filter.callrate(method = "ind", threshold = thresholds$indCallrate) %>% 
  gl.filter.reproducibility(threshold = thresholds$snpRepAvg) %>% 
  gl.filter.secondaries(method = thresholds$secMethod) %>% 
  gl.filter.monomorphs()


## Get genlight numeric matrix ----
glSnpMat <- glSnpSub %>% as.matrix()

## numeric data formating
numeric<-t(glSnpMat)
ID=data.frame(do.call("rbind", strsplit(rownames(numeric), "-", fixed = TRUE)))
numeric=cbind(ID,numeric)

Map<- data.table()
Map$`rs#`       <- glSnpSub$other$loc.metrics$CloneID
Map$alleles     <- glSnpSub$loc.all
Map$chrom       <- glSnpSub$other$loc.metrics$Chrom_Sorghum_bicolor_v5.1
Map$pos         <- glSnpSub$other$loc.metrics$ChromPosSnp_Sorghum_bicolor_v5.1
Map$Ref         <- data.frame(do.call("rbind", strsplit(Map$alleles, "/", fixed = TRUE)))[,1]
Map$Alt         <- data.frame(do.call("rbind", strsplit(Map$alleles, "/", fixed = TRUE)))[,2]

numeric <- cbind(Map,numeric[,-c(1:3)])
numeric$chrom=as.character(numeric$chrom)
## convert chrom to numeric
for(i in 1:10){
  print(i)
  if(i<10){
   numeric$chrom[   which(numeric$chrom==paste("Chr0",i,sep="")) ] = as.character(i)
  } else {
    numeric$chrom[   which(numeric$chrom==paste("Chr",i,sep="")) ] = as.character(i)
    
  }
}

numeric$chrom=as.numeric(numeric$chrom)
# Sort by chrom and pos
numeric=numeric[
  with(numeric, order(chrom, pos)),
]

#remove one of the duplicates not removed my DartR
dim(numeric)
numeric=cbind(data.frame(tempID=paste(numeric$chrom,numeric$pos,sep="_")),numeric)
duplicates=data.frame(table(numeric$tempID))
numeric=numeric[!duplicated(numeric$tempID), ]
dim(numeric)
length(unique(numeric$tempID))

#save numeric data
save(numeric,file=paste("DropXL_2024_Dart_",marker_type,"data_filtered_IndCall",thresholds$indCallrate,"_SNPcall",thresholds$snpCallrate,"_Rep",thresholds$snpRepAvg,"_numeric.Rdata",sep=""))


######Converstion to HapMap#########

#!!!!!!!!!!!!!!!!!!!!!this has only be coded for the SNP data not indel data!!!!!!!!!!!!!!!1

## Initialize character encoded allele matrix ----
glSnpAllele <- matrix(NA, nrow = nrow(glSnpMat), ncol = ncol(glSnpMat))
colnames(glSnpAllele) <- colnames(glSnpMat)
rownames(glSnpAllele) <- rownames(glSnpMat)

## Get genlight numeric matrix ----
glSnpMat <- glSnpSub %>% as.matrix()


## Initialize character encoded allele matrix ----
glSnpAllele <- matrix(NA, nrow = nrow(glSnpMat), ncol = ncol(glSnpMat))
colnames(glSnpAllele) <- colnames(glSnpMat)
rownames(glSnpAllele) <- rownames(glSnpMat)


## Encode each cell for ref/alt/NA values ----
for (i in glSnpMat %>% nrow() %>% seq_len()) {
  for (j in glSnpMat %>% ncol() %>% seq_len()) {
    ref <- alleles(glSnpSub)[j] %>% strsplit("/") %>% unlist() %>% .[1]
    alt <- alleles(glSnpSub)[j] %>% strsplit("/") %>% unlist() %>% .[2]
    
    cell <- glSnpMat[i, j]
    
    glSnpAllele[i, j] <- switch (
      EXPR = cell %>% as.character(),
      `NA` = "NN",
      "0"  = paste0(ref, ref), # Hom ref (CORRECTED)
      "1"  = paste0(ref, alt), # Het (CORRECTED)
      "2"  = paste0(alt, alt)  # Hom alt (CORRECTED)
    )
  }
}


## Transpose ----
tglSnpAllele <- glSnpAllele %>% t()

## Initialize ----
dartHapMap <- data.table()


## Add required columns ----
dartHapMap$`rs#`       <- glSnpSub$other$loc.metrics$CloneID
dartHapMap$alleles     <- glSnpSub$loc.all
dartHapMap$chrom       <- glSnpSub$other$loc.metrics$Chrom_Sorghum_bicolor_v5.1
dartHapMap$pos         <- glSnpSub$other$loc.metrics$ChromPosSnp_Sorghum_bicolor_v5.1
dartHapMap$strand      <- "*"
dartHapMap$`assembly#` <- NA
dartHapMap$center      <- "DArT"
dartHapMap$protLSID    <- NA
dartHapMap$assayLSID   <- NA
dartHapMap$panelLSID   <- NA
dartHapMap$QCode       <- NA


## Add variant data ----
dartHapMap <- cbind(
  dartHapMap,
  tglSnpAllele
)


## Filter (2) duplicate markers ----

dartHapMap$tempID=paste(dartHapMap$chrom,dartHapMap$pos,sep="_")
dartHapMap=dartHapMap[!duplicated(dartHapMap$tempID), ]
#dartHapMap=dartHapMap[,-5649]
#dartHapMap <- dartHapMap[!(pos %in% sumHapDup$pos)] # dropping duplicate markers ----
dartHapMap <- dartHapMap[order(chrom, pos)] # order data based on marker position

sumHapDup <- dartHapMap %>%
  dplyr::group_by(pos) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n > 1)
sumHapDup %>% nrow() %>% paste0("Number of 'duplicate' markers: ", .) %>% print()

## Final summary ----
dartHapMap %>% nrow() %>% paste0("Number of retained markers: ", .) %>% print()
dartHapMap <- dartHapMap[,-c("tempID")] 

## Write to disk ----
data.table::fwrite(
  x     = dartHapMap, 
  file  = hapMapOut, 
  sep   = "\t", # tab identifier
  na    = "NA", # How do we want to write NA data to disk?
  quote = FALSE # when you write to disk do you want quotes to be outputted?
)


###########convert to VCF and impute##################
install.packages("rJava")
library(rJava)
if (!require("devtools")) install.packages("devtools")
devtools::install_bitbucket(
  repo = "bucklerlab/rtassel",
  ref = "master",
  build_vignettes = TRUE
)

## Load hapmap data ----
library(rTASSEL)
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)

#genoPath <- system.file(
  #"C:/Users/ginar/Documents/DropXL_R/Appended_DropXL_2021_DaRTLD",
  #"hapMapOut.hmp.txt",
  #package = "rTASSEL"
#)

genoPath="Report_DS24-9190_SNP_mapping_2.hmp.txt"
Geno <- readGenotypeTableFromPath(genoPath)
Geno %>% exportGenotypeTable(
  file = hapMapOut,
  format = "vcf"
)
#imputation with beagle
system(command = "java -jar ~/Desktop/Old_TempWorkSpace/beagle.28Jun21.220.jar gt=dart_snp_mapping.hmp.txt.vcf out=dart_snp_mapping.BeagleImputed")
#convert to numeric with simple phenotypes
library(simplePHENOTYPES)

as_numeric("~/Google Drive/My Drive/Brian lab folder/HPB/2021 Breeding Data/GenoData/dart_snp_mapping.BeagleImputed.hmp.txt",
           code_as = "012",
           to_file = T)

#add back in allele information
load("~/Google Drive/My Drive/Brian lab folder/HPB/2021 Breeding Data/GenoData/HBP_2021_Dart_SNPdata_filtered_IndCall0.5_SNPcall0.5_Rep0.9_numeric.Rdata")
alleles=numeric[,c(1:7)]
numeric=fread("~/Google Drive/My Drive/Brian lab folder/HPB/2021 Breeding Data/GenoData/dart_snp_mapping.BeagleImputed_numeric.txt")

numeric=merge(alleles,numeric[,-c(2:11)],by.x="rs#",by.y="rs#")
#sort cause merge messed it up
numeric=numeric[
  with(numeric, order(chrom, pos)),
]

save(numeric,file=paste("~/Google Drive/My Drive/Brian lab folder/HPB/2021 Breeding Data/GenoData/HBP_2021_Dart_SNPdata_filtered_IndCall0.5_SNPcall0.5_Rep0.9_numeric_BeagleImputed.Rdata",sep=""))



#######turning Indel data into hapmap
load("~/Google Drive/My Drive/Brian lab folder/HPB/2021 Breeding Data/GenoData/HBP_2021_Dart_SilicoDominatMarkers_filtered_IndCall0.5_SNPcall0.5_numeric.Rdata")

## Initialize ----
dartHapMap <- data.table()

## Add required columns ----
dartHapMap$`rs#`       <- numeric$`rs#`
dartHapMap$alleles     <- "A/T"
dartHapMap$chrom       <- numeric$chrom
dartHapMap$pos         <- numeric$pos
dartHapMap$strand      <- "*"
dartHapMap$`assembly#` <- NA
dartHapMap$center      <- "DArT"
dartHapMap$protLSID    <- NA
dartHapMap$assayLSID   <- NA
dartHapMap$panelLSID   <- NA
dartHapMap$QCode       <- NA



ref="A"
alt="T"
for(i in 1:nrow(dartHapMap)){
  for(j in 12:ncol(dartHapMap) ){
    cell=dartHapMap[i,..j]
IndelAllele[i, j-11] <- switch (
  EXPR = cell %>% as.character(),
  `NA` = "NN",
  "0"  = paste0(ref, ref), # Hom ref (CORRECTED)
  #"1"  = paste0(ref, alt), # Het (CORRECTED)
  "1"  = paste0(alt, alt)  # Hom alt (CORRECTED)
)
  }
}
Markers=numeric[,8:ncol(numeric)]
IndelAllele <- matrix(NA, nrow = nrow(Markers), ncol = ncol(Markers))

# ## Encode each cell for ref/alt/NA values ----
# for (i in  Markers %>% nrow() %>% seq_len()) {
#   print(i)
#   for (j in Markers %>% ncol() %>% seq_len()) {
#     ref="A"
#     alt="T"
#     
#     cell <- Markers[i, j]
#     
#     IndelAllele[i, j] <- switch (
#       EXPR = cell %>% as.character(),
#       `NA` = "NN",
#       "0"  = paste0(ref, ref), # Hom ref (CORRECTED)
#       #"1"  = paste0(ref, alt), # Het (CORRECTED)
#       "1"  = paste0(alt, alt)  # Hom alt (CORRECTED)
#     )
#   }
# }

## Add variant data ----
dartHapMap <- cbind(
  dartHapMap,
  IndelAllele
)
View(data.frame(dartHapMap[,c(1,12)],numeric[,c(2,8)]))
dim(dartHapMap)
colnames(dartHapMap)[12:ncol(dartHapMap)]=colnames(numeric)[8:ncol(numeric)]
hapMapOut  <- paste("dart_Indel_mapping.hmp.txt",sep="")

## Write to disk ----
data.table::fwrite(
  x     = dartHapMap, 
  file  = hapMapOut, 
  sep   = "\t", # tab identifier
  na    = "NA", # How do we want to write NA data to disk?
  quote = FALSE # when you write to disk do you want quotes to be outputted?
)

###########convert to VCF and impute##################
## Load hapmap data ---
library(rTASSEL)
genoPath="dart_Indel_mapping.hmp.txt"
Geno <- readGenotypeTableFromPath(genoPath)
Geno %>% exportGenotypeTable(
  file = hapMapOut,
  format = "vcf"
)
#imputation with beagle
system(command = "java -jar ~/Desktop/Old_TempWorkSpace/beagle.28Jun21.220.jar gt=dart_Indel_mapping.hmp.txt.vcf out=dart_Indel_mapping.BeagleImputed")

#unzip
system(command = "gzip -d dart_Indel_mapping.BeagleImputed.vcf.gz")

genoPath="dart_Indel_mapping.BeagleImputed.vcf"
hapMapOut="dart_Indel_mapping.BeagleImputed"
Geno <- readGenotypeTableFromPath(genoPath)
Geno %>% exportGenotypeTable(
  file = hapMapOut,
  format = "hapmap"
)
#convert to numeric with simple phenotypes

library(simplePHENOTYPES)
hapmap=fread("dart_Indel_mapping.BeagleImputed.hmp.txt",header=T)
hapmap$alleles="T/A"
numeric_Indel=as_numeric(hapmap,
                         method = "reference",
           code_as = "012",
           to_r=T)

View(data.frame(numeric[,c(2,8)],numeric_Indel[,c(1,12)]))
numeric_Indel=merge(numeric[,c(1:7)],numeric_Indel[,-c(2:11)],by.x="rs#",by.y="rs#")
numeric_Indel$alleles="P/A"

# Sort by chrom and pos
numeric_Indel=numeric_Indel[
  with(numeric_Indel, order(chrom, pos)),
]
save(numeric_Indel,file=paste("~/Google Drive/My Drive/Brian lab folder/HPB/2021 Breeding Data/GenoData/HBP_2021_Dart_Indeldata_filtered_IndCall0.5_SNPcall0.5_Rep0.9_numeric_BeagleImputed.Rdata",sep=""))









