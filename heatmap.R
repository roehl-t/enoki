# Create heatmap from csv of FPKM values extracted from Ballgown
# Created by Thomas Roehl on 12/6/2021


# load data
args <- commandArgs(trailingOnly = T)

#manual controls
#args <- c("~/R/Flammulina-velutipes/ballgown/bg_output_2022-01-01/named/300_enoki_genes_qlt05_results_by_tissue_expr.csv",
#"~/R/Flammulina-velutipes/ballgown/fv_sampling_data.csv", 
#"~/R/Flammulina-velutipes/ballgown/bg_output_2022-01-01/heatmaps/cluster4a_heatmap_avg.pdf", 
#"all",
#"type,cul,nor;tissue,gills",
#"all",
#"F",
#"type,tissue",
#"3",
#"~/R/Flammulina-velutipes/ballgown/bg_output_2022-01-01")

data <- read.csv(args[1]) # FPKM data table
pheno <- read.csv(args[2]) # the pheno data file for Ballgown
output <- args[3] # file name to output (PDF)
genes <- strsplit(args[4], ",", fixed = T)[[1]] # list of genes to subset separated by "," -- input "all" to keep all
samplesets <- strsplit(args[5], ";", fixed = T)[[1]] # list of conditions based on column of pheno: each column separated by a semicolon, then column title followed by a comma and a comma-separated list of conditions for that column
organisms <- args[6] # names of organisms to subset separated by "," -- input "all" to keep all
averageT <- args[7] # containing T or t for TRUE
addtonames <- strsplit(args[8], ",", fixed = T)[[1]] # add pheno data from these columns to sample id names
shortennames <- as.numeric(args[9]) # containing an integer to truncate pheno data strings
if(grepl("t", averageT, ignore.case = T)){
    average <- T
} else {
    average <- F
}
setwd(args[10]) # where to set the working directory

# create an output log
zz <- file("heatmap_output_log.txt", open = "wt")
sink(zz, type = "message")

# set up sample titles
pheno$ids <- gsub("-", ".", pheno$ids, fixed = T)
pheno$titlename <- pheno$ids
if(!(args[8] == "none")){
    for(colname in addtonames){
        if(shortennames > 0){
            buffer <- ""
            for(i in 1:shortennames){
                buffer <- paste(buffer, ".", sep = "")
            }
            pheno[,colname] <- paste(pheno[,colname], buffer, sep = "")
            pheno$titlename <- paste(pheno$titlename, substr(pheno[,colname], 1, shortennames), sep = "_")
        } else {
            pheno$titlename <- paste(pheno$titlename, pheno[,colname], sep = "_")
        }
    }
}
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
# a more useable label
#data$newname2 <- paste(data$protein_name, " (", data$UniProt_id, ") _", data$X, sep = "")

# subset organisms
first <- T
if(!(args[6] == "all")){
    for(org in organisms){
        data$keep1 <- grepl(org, data$organism_name, ignore.case = T)
        if(first){
            data$keep2 <- data$keep1
            first <- F
        } else {
            data$keep2 <- (data$keep1 | data$keep2)
        }
    }
    data <- data[(data$keep2 == T),]
    data <- data[,!(colnames(data) %in% c("keep1", "keep2"))]
}
# remove extra columns
data <- data[,!(colnames(data) %in% c("gene_name", "fc", "gene_id", "organism_name", "UniProt_id", "X", "query_id", "protein_name", "newname"))]

# subset samples
first1 <- T
if(!(args[5] == "all")){
    for(sampleset in samplesets){
        first2 <- T
        setlist <- strsplit(sampleset, ",", fixed = T)[[1]]
        if(length(setlist) > 1){
            colname <- setlist[1]
            conditions <- setlist[-1]
            for(condition in conditions){
                pheno$keep1 <- grepl(condition, pheno[,colname], ignore.case = T)
                if(first2){
                    pheno$keep2 <- pheno$keep1
                    first2 <- F
                } else {
                    pheno$keep2 <- (pheno$keep1 | pheno$keep2)
                }
            }
            if(first1){
                pheno$keep3 <- pheno$keep2
                first1 <- F
            } else {
                pheno$keep3 <- (pheno$keep2 & pheno$keep3)
            }
        }
    }
    if(!first1){
        # only do this if you actually added something to keep3
        pheno <- pheno[(pheno$keep3 == T),]
        keepsamples <- pheno$titlename
        keepsamples <- c(keepsamples, c("gene_name", "organism_name", "UniProt_id", "X", "query_id", "protein_name", "newname"))
        data <- data[,(names(data) %in% keepsamples)]
    }
}

# subset genes
first <- T
if(!(args[4] == "all")){
	for(gene in genes){
		data$keep1 <- grepl(tolower(gene), tolower(rownames(data)), fixed = T)
		if(first){
			data$keep2 <- data$keep1
			first <- F
		} else {
			data$keep2 <- (data$keep1 | data$keep2)
		}
	}
	data <- data[(data$keep2 == T),]
	data <- data[,!(colnames(data) %in% c("keep1", "keep2"))]
}

# a more useable label
#row.names(data) <- data$newname2
#data <- data[,!(colnames(data) %in% c("newname2"))]

# instead of showing raw data on the heatmap, take an average for samples with the same phenotypic conditions
if(average){
  stages <- c("pri", "you", "cul", "nor")
  tissues <- c("myc", "sti", "pil", "gil")
  newdata <- data.frame(rep(T, nrow(data)), row.names = row.names(data))
  names(newdata) <- c("X")
  for(stage in stages){
    for(tissue in tissues){
      combo <- paste(stage, tissue, sep = "_")
      counter <- 0
      newdata$newcol <- rep(0, nrow(newdata))
      for(column in names(data)){
        if(grepl(combo, column, fixed = T)){
          counter <- counter+1
          newdata$newcol <- newdata$newcol + data[,column]
        }
      }
      if(counter > 0){
        newdata[,ncol(newdata)+1] <- newdata$newcol/counter
        names(newdata)[ncol(newdata)] <- paste(combo, "_avg_", counter, sep = "")
      }
    }
  }
  newdata <- newdata[,!(names(newdata) %in% c("X", "newcol"))]
  data <- newdata
}

# transform FPKM by log2(FPKM+1)
for(col in colnames(data)){
  data[,col] <- log2(data[,col] + 1)
}

# convert to matrix
datam <- as.matrix(data)

if(nrow(datam) >= 2 & ncol(datam) >= 2){
  # calculate plot size
  ############## not working correctly
  minRowHeight <- 80/(2494*0.1)
  minFileHeight <- minRowHeight*(nrow(datam)*0.1)
  minColWidth <- 10/(23*1)
  minFileWidth <- minColWidth*(ncol(datam)*1)
  filesize <- ifelse(minFileHeight > minFileWidth, minFileHeight, minFileWidth)
  rowscale <- filesize/(minRowHeight*nrow(datam)*0.1)
  colscale <- filesize/(minColWidth*ncol(datam)*1)
  marginFactor <- 2*12/(0.4*10)
  maxMargin <- filesize*0.4*marginFactor
  marginAdj <- ifelse(5*rowscale > maxMargin, (5*rowscale-maxMargin)/marginFactor, 0)
  filesize <- filesize+marginAdj
  
  # plot heatmap and save as pdf
  #pdf(file = output, width = 80, height = 80) # all samples
  #heatmap(datam, dendrogram = "both", Rowv = T, Colv = T, symbreaks = F, key = T, symkey = F, density.info = "none", trace = "none", margins = c(60, 8), cexCol = 5, cexRow = 0.1)
  #pdf(file = output, width = 10, height = 10) # 23c x 26r
  #heatmap(datam, dendrogram = "both", Rowv = T, Colv = T, symbreaks = F, key = T, symkey = F, density.info = "none", trace = "none", margins = c(12,30), cexCol = 1)
  #pdf(file = output, width = filesize, height = filesize)
  #heatmap(datam, dendrogram = "both", Rowv = T, Colv = T, symbreaks = F, key = T, symkey = F, scale = "none", density.info = "none", trace = "none", margins = c(12*colscale, 5*rowscale), cexCol = 1*colscale, cexRow = 0.1*rowscale)
  pdf(file = output, width = 11, height = 7)
  heatmap(datam, dendrogram = "both", Rowv = T, Colv = T, symbreaks = F, key = T, symkey = F, scale = "none", density.info = "none", trace = "none", margins = c(10, 7), cexCol = 0.4, cexRow = 0.4)
  dev.off()
} else {
  print(paste("No data for gene set ", args[4], ", sample set ", args[5], ", and organism set ", args[6], sep = ""))
}
