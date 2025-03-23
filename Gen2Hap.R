#Load package
library(tidyverse)
library(data.table)

## Get genlight numeric matrix ----
glMat <- Quinby %>% as.matrix()

## numeric data formating
numeric<-t(glMat)
ID=data.frame(do.call("rbind", strsplit(rownames(numeric), "-", fixed = TRUE)))
numeric=cbind(ID,numeric)

Map<- data.table()
Map$`rs#`       <- Quinby$other$loc.metrics$CloneID
Map$alleles     <- Quinby$loc.all
Map$chrom       <- Quinby$other$loc.metrics$Chrom_Sorghum_bicolor_v5.1
Map$pos         <- Quinby$other$loc.metrics$ChromPosSnp_Sorghum_bicolor_v5.1
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
#save(numeric,file=paste("DropXL_2024_Dart_",marker_type,"data_filtered_IndCall",thresholds$indCallrate,"_SNPcall",thresholds$snpCallrate,"_Rep",thresholds$snpRepAvg,"_numeric.Rdata",sep=""))
save(numeric, file = "numreic.rdata")

######Converstion to HapMap#########

#!!!!!!!!!!!!!!!!!!!!!this has only be coded for the SNP data not indel data!!!!!!!!!!!!!!!1

## Initialize character encoded allele matrix ----
glSnpAllele <- matrix(NA, nrow = nrow(glMat), ncol = ncol(glMat))
colnames(glSnpAllele) <- colnames(glMat)
rownames(glSnpAllele) <- rownames(glMat)

## Get genlight numeric matrix ----
glMat <- Quinby %>% as.matrix()


## Initialize character encoded allele matrix ----
glSnpAllele <- matrix(NA, nrow = nrow(glMat), ncol = ncol(glMat))
colnames(glSnpAllele) <- colnames(glMat)
rownames(glSnpAllele) <- rownames(glMat)


## Encode each cell for ref/alt/NA values ----
for (i in glMat %>% nrow() %>% seq_len()) {
  for (j in glMat %>% ncol() %>% seq_len()) {
    ref <- alleles(Quinby)[j] %>% strsplit("/") %>% unlist() %>% .[1]
    alt <- alleles(Quinby)[j] %>% strsplit("/") %>% unlist() %>% .[2]
    
    cell <- glMat[i, j]
    
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
dartHapMap$`rs#`       <- Quinby$other$loc.metrics$CloneID
dartHapMap$alleles     <- Quinby$loc.all
dartHapMap$chrom       <- Quinby$other$loc.metrics$Chrom_Sorghum_bicolor_v5.1
dartHapMap$pos         <- Quinby$other$loc.metrics$ChromPosSnp_Sorghum_bicolor_v5.1
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
  file  = "Quinby.hapmap2.hmp.txt", 
  sep   = "\t", # tab identifier
  na    = "NA", # How do we want to write NA data to disk?
  quote = FALSE # when you write to disk do you want quotes to be outputted?
)

