
#download TCGAbiolink from bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")

#download TCGAbiolinks from git
if(!require(devtools)) install.packages("devtools") 
library("devtools")
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")

install.packages("tidyverse")
install.packages("DT")
install.packages("regexPipes")
install.packagesd('data.table')
#install curl
install.packages("https://github.com/jeroen/curl/archive/master.tar.gz", repos = NULL)

library(bitops)
library(RCurl)
getURL("https://portal.gdc.cancer.gov")
library(TCGAbiolinks)
library(dplyr)
library(tidyverse)
library(DT)
library("regexPipes")
library(data.table)

#get all TCGA cancer subtype codes
tcgacodes <- TCGAbiolinks:::getGDCprojects()$project_id
tcgacodes <- unlist(lapply(tcgacodes[str_detect(tcgacodes, "TCGA")],
                    function(x){strsplit(x, "-")[[1]][2]}))
                                                                                                                                                                                                                                                                    
#download MAF files for all samples and combine into 1 dataframe
df <- data.frame()
for (t in tcgacodes[21:33]){
  maf <- GDCquery_Maf(t, pipelines = "mutect") %>%
    mutate(cancertype = t)
  df <- rbind(df, maf)
}

#you can then get the VAF with df$VAF = df$t_alt_count/df$t_depth
df <- read.csv("df.csv")
write.csv(df,"df.csv")
write.csv(maf,"maf.csv")
tcga <- data.frame(8)
tcga3 = as.data.frame(df$n_alt_count) 
tcga3$t_alt_count <- df$t_alt_count
tcga3$n_ref_count <-  df$n_ref_count
tcga3$t_ref_count <-  df$t_ref_count
tcga3$chr <- df$Chromosome
tcga3$Start_Position <-  df$Start_Position
tcga3$End_Position <- df$End_Position
tcga3$Tumor_Sample_Barcode <- df$Tumor_Sample_Barcode
write.csv(tcga,"maf_1-33.csv")

tcga <- rbind(tcga,tcga2,tcga3)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                write.csv(tcga,"tcga.csv")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                tcga$Matched_df$Matched_Norm_Sample_Barcode
#if you want copy number information you can use the following:
tcgacodes <- TCGAbiolinks:::getGDCprojects()$project_id

query <- GDCquery(project = tcgacodes[str_detect(tcgacodes, "TCGA")],
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")

cnvsamples <- getResults(query) #records of id
GDCdownload(query) #download copy number segment data
dfhg38cnv <- GDCprepare(query)

write.csv(cnvsamples,'cnvsamples.csv') 

##get a manfest of sottoriva's data 
query_sottoriva <- GDCquery(project = c("TCGA-KIRC"),data.category= "Raw Sequencing Data",barcode="TCGA-A3-3373-01A-02D-1421-08",sample.type = "Primary solid Tumor")
query_s1 <- GDCquery(project = c(id[1,2]),data.category= "Raw Sequencing Data",barcode=id[1,1],sample.type = "Primary solid Tumor")
GDCdownload(query_s1)
query_sottoriva[,1]
##load sottoriva's file
install.packages("gdata")
require(gdata)
sottoriva <- read.csv("Sottoriva_TCGA_abbr.csv",header=TRUE)
id <- read.table("/home/xinyu/Work/Phd/Data/Sottoriva_tcga_819.txt")

for (i in 1:nrow(id)){
  id[i,2] <- paste0("TCGA-",sottoriva[i,2])
}
colnames(id) = c("barcode","project")

write.csv(id,'id.csv') 

id <- read.csv('id.csv')[,-1]
####Download 1-10
for (i in 201:240){
  tryCatch({
    print(paste0("Query No.",i))
    query_s1 <- GDCquery(project = id[i,3],data.category= "Raw Sequencing Data",barcode=id[i,2],sample.type = "Primary solid Tumor")
    print("finish querying..downloading..")
    GDCdownload(query_s1,method = "client",token.file = "gdc-user-token.2018-10-26T16_18_38.485Z.txt")
    id[i,4] = TRUE 
  },error=function(e){id[i,4]=FALSE})
}
GDCquery(project = id[1,3],data.category= "Raw Sequencing Data",barcode=id[1,2],sample.type = "Primary solid Tumor")