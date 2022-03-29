# Create heatmap from csv of FPKM values extracted from Ballgown
# Created by Thomas Roehl on 12/6/2021


# load data
args <- commandArgs(trailingOnly = T)

#manual controls
#args <- c("~/R/Flammulina-velutipes/ballgown/bg_output_2022-01-01/named/300_enoki_genes_qlt05_results_by_tissue_expr.csv",
#"~/R/Flammulina-velutipes/ballgown/fv_sampling_data.csv", 
#"~/R/Flammulina-velutipes/ballgown/bg_output_2022-01-01/heatmaps/cluster4a_heatmap_avg.pdf", 
#"all",
#"all",
#"all",
#"F",
#"~/R/Flammulina-velutipes/ballgown/bg_output_2022-01-01")

data <- read.csv(args[1]) # FPKM data table
pheno <- read.csv(args[2]) # the pheno data file for Ballgown
output <- args[3] # file name to output (PDF)
genes <- strsplit(args[4], ",", fixed = T)[[1]] # list of genes to subset separated by "," -- input "all" to keep all
samples <- strsplit(args[5], ",", fixed = T)[[1]] # list of samples or 3-letter sample groups (ex: sti) to subset separated by "," -- input "all" to keep all
samples <- c(samples, c("gene_name", "organism_name", "UniProt_id", "X", "query_id", "protein_name", "newname"))
organisms <- args[6] # names of organisms to subset separated by "," -- input "all" to keep all
averageT <- args[7] # containing T or t for TRUE
if(grepl("t", averageT, ignore.case = T){
    average <- T
} else {
    average <- F
}
setwd(args[8]) # where to set the working directory

# create an output log
zz <- file("heatmap_output_log.txt", open = "wt")
sink(zz, type = "message")

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
data <- data[,!(colnames(data) %in% c("gene_name", "organism_name", "UniProt_id", "X", "query_id", "protein_name", "newname"))]

# subset samples
first <- T
if(!(args[5] == "all")){
	samples <- gsub("-", ".", samples, fixed = T)
	for(sample in samples){
		pheno$keep1 <- grepl(sample, pheno$titlename, ignore.case = T)
		if(first){
			pheno$keep2 <- pheno$keep1
			first <- F
		} else {
			pheno$keep2 <- (pheno$keep1 | pheno$keep2)
		}
	}
	pheno <- pheno[(pheno$keep2 == T),]
	psams <- pheno$titlename
	data <- data[,names(data) %in% psams]
}

# subset genes
first <- T
if(!(args[4] == "all")){
	for(gene in genes){
		data$keep1 <- grepl(tolower(paste(gene," ",sep = "")), tolower(paste(rownames(data)," ",sep = "")), fixed = T)
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

# calculate plot size
############## not working correctly -- this branch is dedicated to fixing the automatic resizing feature
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
