  #download TCGAbiolinks from git
  devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks",force=TRUE)
  library(TCGAbiolinks)

  #get all TCGA cancer subtype codes
  tcgacodes <- TCGAbiolinks:::getGDCprojects()$project_id
  tcgacodes <- unlist(lapply(tcgacodes[str_detect(tcgacodes, "TCGA")],
                             function(x){strsplit(x, "-")[[1]][2]}))
  
  ##get a manfest of sottoriva's data 
  query <- GDCquery(project = tcgacodes[str_detect(tcgacodes, "TCGA")],
                    data.category = "Copy Number Variation",
                    data.type = "Copy Number Segment")
  
  cnvsamples <- getResults(query) #records of id
  GDCdownload(query) #download copy number segment data
  dfhg38cnv <- GDCprepare(query)
  
  ##get a manfest of RNA
  query <- GDCquery(project = "TCGA-COAD",
                    data.category= "Transcriptome Profiling",
                    data.type= "Raw sequencing data",
                    sample.type = "Primary Tumor")
  manSelect$samplename
  ##get a manfest of RNA
  query <- GDCquery(project = c("TCGA-COAD","TCGA-READ"),
                    barcode = manSelect$samplename,
                    data.category= "Sequencing Reads",
                    #sample.type = "Primary Tumor",
                    access="controlled",
                    sample.type="")
  
  ## Filter query table
  query$results[[1]] <- query$results[[1]][str_detect(query$results[[1]][,"file_name"], "rehead"),]
  
  ## Check query information
  data.table(getResults(query, rows = 1:100,cols = c("file_name","cases")), 
            filter = 'top',
            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
            rownames = FALSE)
  
  GDCdownload(query,method = "client", files.per.chunk = 10,token.file = "gdc-user-token.2021-07-27T09_41_26.963Z.txt")


 GDCdownload(query, method = "api", files.per.chunk = 10)



install.packages("tidyverse")
install.packages("DT")
install.packages("regexPipes")
install.packages('data.table')
#install curl
install.packages("https://github.com/jeroen/curl/archive/master.tar.gz", repos = NULL)

library(bitops)
library(RCurl)
getURL("https://portal.gdc.cancer.gov")
library(TCGAbiolinks)
library(tidyverse)
library(dplyr)
library(DT)
library("regexPipes")
library(data.table)

id <- read.csv('id.csv')
####Download 1-10
for (i in 41:90){
  tryCatch({
    print(paste0("Query No.",i))
    query_s1 <- GDCquery(project = id[i,3],data.category= "Raw Sequencing Data",barcode=id[i,2],sample.type = "Primary solid Tumor")
    print("finish querying..downloading..")
    GDCdownload(query_s1,method = "client",token.file = "/home/xinyu/Work/Phd/TCGA_BAM/gdc-user-token.2018-10-26T16_18_38.485Z.txt")
    id[i,4] = TRUE 
  },error=function(e){id[i,4]=FALSE})
}