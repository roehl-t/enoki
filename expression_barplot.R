# written by Thomas Roehl 1/21/2022

library(ggplot2)


genes <- c("cell division")                  #list of genes to subset
stages <- c("any")                  # use "any" to select all stages
tissues <- c("any")                  # use "any" to select all tissues
andor <- "and"                  # either "and" or "or" -- match one or both?
xaxis <- "combo"                  # "tissue" "stage" or "combo" (for both)
arrangeby <- "tissue"                  # when using "combo" in plot, sort by tissue or stage


data <- read.csv("/mnt/raid1/Flammulina-velutipes_development/Data/mojo-test/analysis/auto-5000/calculations/proteome/named/mojo_genes_qlt05_results_by_tissue_expr.csv")

pheno <- read.csv("/mnt/raid1/Flammulina-velutipes_development/Data/sample-data-for-pipeline.csv")

# keep only specified tissues/stages
pheno$keepS <- rep(0, nrow(pheno))
pheno$keepT <- rep(0, nrow(pheno))
for(i in 1:nrow(pheno)){
  for(stage in stages){
    if(grepl(stage, pheno[i,"type"], ignore.case = T) || stage == "any"){
      pheno[i,"keepS"] <- 1
    }
  }
  for(tissue in tissues){
    if(grepl(tissue, pheno[i,"tissue"], ignore.case = T) || tissue == "any"){
      pheno[i,"keepT"] <- 1
    }
  }
}
pheno$keepB <- pheno$keepS + pheno$keepT

# set up sample titles
pheno$ids <- gsub("-", ".", pheno$ids, fixed = T)
pheno$titlename <- paste(pheno$ids, "_", substr(pheno$type, 1, 3), "_", substr(pheno$tissue, 1, 3), sep = "")
oldnames <- data.frame(names(data))
names(oldnames) <- c("ids")
phenomerge <- merge(pheno, oldnames, all.x = F, all.y = T)
phenomerge$titlename <- ifelse(is.na(phenomerge$titlename), phenomerge$ids, phenomerge$titlename)

# reorder and rename data titles
data <- data[, phenomerge$ids]
names(data) <- phenomerge$titlename

# replace NA with MSTRG
data$protein_name <- ifelse(is.na(data$protein_name), data$query_id, data$protein_name)

# deal with duplicate row names
data$newname <- paste(data$protein_name, "_", data$X, sep = "")
row.names(data) <- data$newname

# remove unwanted columns
if(andor == "and"){
  keeplist <- pheno[pheno$keepB == 2,]
} else {
  keeplist <- pheno[pheno$keepB >= 1,]
}
datanames <- names(data)
first = T
for(dataname in datanames){
  if(!grepl("ROEHL", dataname, fixed = T)){
    if(first){
      datanames2 <- c(dataname)
      first = F
    } else {
      datanames2[length(datanames2)+1] <- dataname
    }
  }
}
keepnames <- append(datanames2, keeplist$titlename, after = length(datanames2))
data <- data[,names(data) %in% keepnames]

# filter genes
data$keep <- rep(F, nrow(data))
for(i in 1:nrow(data)){
  for(gene in genes){
    if(grepl(gene, data$newname[i], ignore.case = T)){
      data$keep[i] <- T
    }
  }
}
data <- data[data$keep == T,]

# create new dataframe: gene name, tissue, stage, sample, expression
first <- T
for(name in names(data)){
  if(grepl("ROEHL", name, fixed = T)){
    if(first){
      samplenames <- c(name)
      first = F
    } else {
      samplenames[length(samplenames)+1] <- name
    }
  }
}
samples <- data[,names(data) %in% samplenames]
rownames(samples) <- 1:nrow(samples)
info <- data[, !(names(data) %in% samplenames)]
rownames(info) <- 1:nrow(info)
first <- T
for(i in 1:ncol(samples)){
  samplename <- names(samples)[i]
  sample <- rep(samplename, nrow(samples))
  FPKM <- samples[,i]
  nameinfo <- strsplit(samplename, "_", fixed = T)[[1]]
  tissue3 <- nameinfo[4]
  tissuename <- switch(tissue3,
                   "myc" = "mycelium",
                   "sti" = "stipe",
                   "pil" = "pileus",
                   "gil" = "gills")
  tissue <- rep(tissuename, nrow(samples))
  stage3 <- nameinfo[3]
  stagename <- switch(stage3,
                  "pri" = "primordium",
                  "you" = "young",
                  "cul" = "cultivated",
                  "nor" = "normal")
  stage <- rep(stagename, nrow(samples))
  combo <- paste(stage, tissue, sep = "_")
  temp <- data.frame(sample, tissue, stage, combo, FPKM)
  for(j in 1:ncol(info)){
    temp[,j+5] <- info[j]
  }
  if(first){
    result <- temp
    first <- F
  } else {
    result2 <- merge(result, temp, all.x = T, all.y = T)
    result <- result2
  }
}
result$log2 <- log2(result$FPKM + 1)
graphdata <- result[,names(result) %in% c("sample", "tissue", "stage","combo", "newname", "log2")]
names(graphdata) <- c("sample", "tissue", "stage", "combo", "gene", "log2")

# give stages/tissues an order
listtissues <- c("mycelium", "stipe", "pileus", "gills")
graphdata$tissue <- factor(graphdata$tissue, levels = listtissues)
liststages <- c("primordium", "young", "cultivated", "normal")
graphdata$stage <- factor(graphdata$stage, levels = liststages)
if(arrangeby == "stage"){
  list1s <- liststages
  list2s <- listtissues
  reverse <- F
} else {
  list2s <- liststages
  list1s <- listtissues
  reverse <- T
}
first <- T
for(list1 in list1s){
  for(list2 in list2s){
    if(!reverse){
      pastename <- paste(list1, list2, sep = "_")
    } else {
      pastename <- paste(list2, list1, sep = "_")
    }
    if(first){
      listcombos <- pastename
      first <- F
    } else {
      listcombos[length(listcombos)+1] <- pastename
    }
  }
}
graphdata$combo <- factor(graphdata$combo, levels = listcombos)

# make graphs
gg <- ggplot(graphdata, aes_string(x = xaxis, y = "log2", fill = "gene", shape = "gene")) + 
  stat_summary(fun = "mean", geom = "bar", position = position_dodge2(preserve = "single")) + 
  geom_jitter(position = position_dodge2(width = 1, preserve = "single")) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(gg)
