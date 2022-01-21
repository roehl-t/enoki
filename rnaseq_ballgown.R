# from Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095
# modified by Todd Osmundson and Thomas Roehl 2021

#!/usr/bin/env Rscript
# run this in the output directory for rnaseq_pipeline.sh
# passing the pheno data csv file as the only argument 
args = commandArgs(trailingOnly=TRUE)
#if (length(args)==0) {
# assume no output directory argument was given to rnaseq_pipeline.sh

pheno_data_file <- args[1]
data_dir_loc <- args[2]

# set covariate (covname) and confounding variables (adjnames)
covname <- "tissue"
adjnames <- c("type")

# when comparing exactly 2 stages, set covariate as type to take advantage of getFC
setnamesplit <- strsplit(data_dir_loc, "/")[[1]][]
setname <- setnamesplit[length(setnamesplit)]
if(setname == "culnor" | setname == "priyou"){
    covname <- "type"
    adjnames <- c("tissue")
}

setwd(paste(data_dir_loc,"/bg_output", sep = ""))

# create an output log
zz <- file("ballgown_output_log.txt", open="wt")
sink(zz, type="message")

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(ggplot2)
library(ggrepel)


## manual controls
#pheno_data_file <- "~/R/Flammulina-velutipes/ballgown/fv_sampling_data.csv"
#data_dir_loc <- "C:/Users/mycol/Documents/UW-La Crosse/Flammulina velutipes development/Data/all_ballgown_2022-01-18"
#setwd("~/R/Flammulina-velutipes/ballgown/bg_output_2022-01-18/younor")
#covname <- "type"
#adjnames <- c("tissue")


## Read phenotype sample data
pheno_data <- read.csv(pheno_data_file, header=TRUE, stringsAsFactors=F)
print("Phenotypic data read")


# to re-run transcript cutoff manipulation, first restore hidden files from previous runs
direclist <- list.dirs(paste(data_dir_loc, "/hidden", sep = ""), recursive = F)
for(dir in direclist){
    newpath <- sub("/hidden", "", dir)
    if(dir.exists(newpath)){
        # duplicate files were left in /hidden from an old run
        # delete those extra files
        unlink(dir, recursive = T)
    } else {
        # no duplicates exist, restore the hidden files
        file.rename(from = dir, to = newpath)
    }
}


## count the number of distinct transcripts in each t_data file
dirlist <- list.dirs(data_dir_loc, recursive = F)
first <- T
for(dir in dirlist){
    if(file.exists(paste(dir, "/t_data.ctab", sep = ""))){
        tfile <- read.table(paste(dir, "/t_data.ctab", sep = ""), header = T)
        tfile <- tfile[(tfile$FPKM > 0),]
        sample <- strsplit(dir, "/")[[1]]
        sample <- sample[length(sample)]
        newrow <- data.frame(sample, nrow(tfile))
        colnames(newrow) <- c("ids", "transcript_count")
        if(first){
            countdata <- newrow
            first <- F
        } else {
            countdata <- merge(countdata, newrow, all.x = T, all.y = T)
        }
    }
}
write.csv(countdata, file = "transcript_counts.csv", row.names = F)

## remove all files < a cutoff from the pipeline
##    adjust the cutoff to tailor sensitivity
cutoff <- 300
print(paste("Distinct transcripts counted. Removing samples < ", cutoff, "."), sep = "")

removesamples <- countdata[(countdata$transcript_count < cutoff),]

## to manually remove samples, adjust settings here
##    create a new vector for each column in pheno_data you want to check
##    add all vectors to the remAll list
##    add the name of the column to remColList in the same order they appear in remAll
#remStages <- c("pri", "cul")
#remAll <- list(remStages)
#remColList <- c("type", "tissue")
#
#remManual <- data.frame(pheno_data$ids, rep(F, nrow(pheno_data)))
#names(remManual) <- c("ids", "remove")
#for(i in 1:length(remColList)){
#    for(j in 1:length(remAll[[i]])){
#        for(k in 1:nrow(remManual)){
#            if(grepl(remAll[[i]][j], pheno_data[k,remColList[i]], ignore.case = T)){
#                remManual$remove[k] <- T
#            }
#        }
#    }
#}
#remManual <- remManual[(remManual$remove == T),]
#remSamMerge <- merge(removesamples, remManual, all.x = T, all.y = T)
#for(i in 1:nrow(remSamMerge)){
#    if(is.na(remSamMerge$transcript_count[i])){
#        remSamMerge$transcript_count[i] <- 1
#    }
#}
#removesamples <- remSamMerge[,!(names(remSamMerge) %in% c("remove"))]

## remove those sample rows from the pdata file
##    (assumes sample names are stored in "ids" column)
pdatamerge <- merge(pheno_data, removesamples, all.x = T, all.y = F)
pdatamerge$transcript_count <- is.na(pdatamerge$transcript_count)
pheno_data <- pdatamerge[(pdatamerge$transcript_count == T),]
pheno_data <- pheno_data[,!(names(pheno_data) %in% c("transcript_count"))]

## hide the corresponding files from ballgown
if(!dir.exists(paste(data_dir_loc, "/hidden", sep = ""))){
    dir.create(paste(data_dir_loc, "/hidden", sep = ""))
}
dirlist <- list.dirs(data_dir_loc, recursive = F)
for(dir in dirlist){
    for(sample in removesamples$ids){
        if(grepl(sample, dir, fixed = T)){
            filepath <- strsplit(dir, "/")[[1]]
            lastele <- filepath[length(filepath)]
            filepath[length(filepath)] <- "hidden"
            filepath[length(filepath)+1] <- lastele
            newpath <- paste(filepath, collapse = "/", sep = "")
            file.rename(from = dir, to = newpath)
        }
    }
}

# double check that only samples in the directory are included in the pheno data
dirlist <- list.dirs(data_dir_loc, recursive = F)
samplelist <- dirlist
for( i in 1:length(dirlist)){
    sample <- strsplit(dirlist[i], "/")[[1]]
    samplelist[i] <- sample[length(sample)]
}
keeplist <- rep(T, length(dirlist))
keepframe <- data.frame(samplelist, keeplist)
names(keepframe) <- c("ids", "keep")
checkmerge <- merge(pheno_data, keepframe, all.x = T, all.y = F)
checkmerge$keep <- !is.na(checkmerge$keep)
pheno_data <- checkmerge[(checkmerge$keep == T),]
pheno_data <- pheno_data[,!(names(checkmerge) %in% c("keep"))]

print(paste("Removed samples moved to ", data_dir_loc, "/hidden", sep = ""))


## Read in expression data
bg_enoki <- ballgown(dataDir = data_dir_loc, samplePattern="ROEHL-FV-", pData=pheno_data)
print("Expression data read")

## NOTE: If some samples have zero reads mapping to the reference genome, the mergelist.txt file (in the output/hisat2 directory) will not match the phenodata file (pData object). The following error will occur:
        #Error in ballgown(dataDir = "/Volumes/RAID_5_data_array/Todd/D_polymorpha_HiSeq/output/ballgown",  :
        #first column of pData does not match the names of the folders containing the ballgown data.
        #In addition: Warning message:
        #In ballgown(dataDir = "/Volumes/RAID_5_data_array/Todd/D_polymorpha_HiSeq/output/ballgown",  :
        #Rows of pData did not seem to be in the same order as the columns of the expression data. Attempting to rearrange pData...
        #Execution halted
     #Problem can be corrected by removing these samples from the sample files and the phenodata file prior to running the pipeline.
    # Or, the problem can arise when directories in output/ballgown don't match the names in the phenodata file (first column).
## Filter low abundance genes
bg_enoki_filt <- subset(bg_enoki, "rowVars(texpr(bg_enoki)) > 1", genomesubset=TRUE)
print("Low abundance genes filtered")

## Generate PCA plots
# begin with ballgown object
bg_obj <- bg_enoki_filt

# extract pheno data and move sample IDs ("ids") to row names
pdata <- indexes(bg_obj)$pData
rownames(pdata) <- pdata$ids

# extract expression data
gene <- gexpr(bg_obj)
transcript <- texpr(bg_obj)
expr <- gene

# remove "FPKM." previx from sample names
colnames(expr) <- sub("FPKM.", '', colnames(expr))

# conduct PCA on expression matrix
pca_expr <- prcomp(expr)
sink(paste("./PCA/", cutoff, "_enoki_PCA_results.txt", sep = ""))
summary(pca_expr)
sink()
#biplot(pca_expr)

# extract vectors for each PC and store in new data frame
pca_expr_vectors <- data.frame(pca_expr$rotation[,1])
pca_expr_vectors_names <- colnames(pca_expr$rotation)
names(pca_expr_vectors) <- c(pca_expr_vectors_names[1])
for(i in 2:ncol(pca_expr$rotation)){
    pca_expr_vectors[pca_expr_vectors_names[i]] <- pca_expr$rotation[,i]
}

# add in pheno data
c <- ncol(pca_expr_vectors)
for(i in 1:ncol(pdata)){
    # add new column
    t <- rep(pdata[1,i], nrow(pca_expr_vectors))
    pca_expr_vectors[,c+i] <- t
    colnames(pca_expr_vectors)[c+i] <- colnames(pdata)[i]
    
    # copy over pheno data to new column
    for(j in 1:nrow(pca_expr_vectors)){
        pca_expr_vectors[j,c+i] <- pdata[rownames(pca_expr_vectors)[j],i]
    }
}

# create plots

# batch and jar
pca_expr_vectors$batch <- as.factor(pca_expr_vectors$batch)
pca_expr_vectors$jar <- as.factor(pca_expr_vectors$jar)
pdf(file = paste("./PCA/", cutoff, "_enoki_pca_batch-jar.pdf", sep = ""), width = 6, height = 6)
ggplot(pca_expr_vectors, aes(x = PC1, y = PC2, color = jar, shape = batch)) +
geom_vline(xintercept = 0, color = "grey") +
geom_hline(yintercept = 0, color = "grey") +
geom_point() + theme_classic() + theme(text = element_text(family = "serif")) + 
geom_text_repel(aes(label = ids), family = "serif")
dev.off()

# type and tissue
pdf(file = paste("./PCA/", cutoff, "_enoki_pca_type-tissue.pdf", sep = ""), width = 6, height = 6)
ggplot(pca_expr_vectors, aes(x = PC1, y = PC2, color = type, shape = tissue)) +
geom_vline(xintercept = 0, color = "grey") +
geom_hline(yintercept = 0, color = "grey") +
geom_point() + theme_classic() + theme(text = element_text(family = "serif")) + 
geom_text_repel(aes(label = ids), family = "serif")
dev.off()

# position and mushroom
pdf(file = paste("./PCA/", cutoff, "_enoki_pca_position.pdf", sep = ""), width = 6, height = 6)
ggplot(pca_expr_vectors, aes(x = PC1, y = PC2, color = position)) +
geom_vline(xintercept = 0, color = "grey") +
geom_hline(yintercept = 0, color = "grey") +
geom_point() + theme_classic() + theme(text = element_text(family = "serif")) + 
geom_text_repel(aes(label = ids), family = "serif")
dev.off()

print("PCA plots created")


## DE by transcript
results_transcripts <-  stattest(bg_enoki_filt, feature='transcript', covariate=covname,
         adjustvars=adjnames, getFC=TRUE, meas='FPKM')

## DE by gene
results_genes <-  stattest(bg_enoki_filt, feature='gene', covariate=covname, adjustvars=adjnames, getFC=TRUE, meas='FPKM')


## Add gene names
#results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_enoki_filt),
#          geneIDs=ballgown::geneIDs(bg_enoki_filt), results_transcripts)

## Sort results from smallest p-value
results_transcripts <- arrange(results_transcripts, pval)
results_genes <-  arrange(results_genes, pval)

## Write results to CSV
write.csv(results_transcripts, paste(cutoff, "_enoki_transcripts_results_by_", covname, ".csv", sep = ""), row.names=FALSE)
write.csv(results_genes, paste(cutoff, "_enoki_genes_results_by_", covname, ".csv", sep = ""), row.names=FALSE)

## Filter for genes with q-val <0.05
transcripts_q_lt_05 <- subset(results_transcripts, results_transcripts$qval <=0.05)
genes_q_lt_05 <- subset(results_genes, results_genes$qval <=0.05)

## Format for later gene ID using python
# use first t_data.ctab file to match t_name MSTRG numbers to t_id
finding <- T
dirlist <- list.dirs(data_dir_loc, recursive = F)
i <- 1
while(finding){
    if(file.exists(paste(dirlist[i], "/t_data.ctab", sep = ""))){
        tfile <- read.table(paste(dirlist[i], "/t_data.ctab", sep = ""), header = T)
        finding = F
    }
    i <- i+1
}
tfile <- tfile[,(colnames(tfile) %in% c("t_name", "t_id"))]
colnames(transcripts_q_lt_05)[2] <- "t_id"
transcripts_q_lt_05$t_id <- as.numeric(transcripts_q_lt_05$t_id)
tmerge <- merge(transcripts_q_lt_05, tfile, all.x = T, all.y = F)

genes_q_lt_05$gene_id <- paste(genes_q_lt_05$id, ".1", sep = "")

## Write lists of DEGs with q<0.05 (for gene ID) separated by new line
write(tmerge$t_name, file = paste("./lists/", cutoff, "_enoki_transcripts_qlt05_results_by_", covname, ".txt", sep = ""), sep="\n")
write(genes_q_lt_05$gene_id, file = paste("./lists/", cutoff, "_enoki_genes_qlt05_results_by_", covname, ".txt", sep = ""), sep="\n")

print("DE lists written")

## Write expression data for DEs
gene <- gexpr(bg_enoki_filt)
gene <- data.frame(gene)
colnames(gene) <- sub("FPKM.", "", colnames(gene))
gene$id <- rownames(gene)
genemerge <- merge(genes_q_lt_05, gene, all.x = T, all.y = F)
rownames(genemerge) <- genemerge$gene_id
writegenes <- genemerge[,(names(genemerge) %in% c("pval", "qval"))]
write.csv(writegenes, paste(cutoff, "_enoki_genes_qlt05_results_by_", covname, ".csv", sep = ""), row.names=T)
writegenes <- genemerge[,!(names(genemerge) %in% c("id", "feature", "pval", "qval", "gene_id"))]
write.csv(writegenes, paste(cutoff, "_enoki_genes_qlt05_results_by_", covname, "_expr.csv", sep = ""), row.names=T)

transcript <- texpr(bg_enoki_filt)
transcript <- data.frame(transcript)
colnames(transcript) <- sub("FPKM.", "", colnames(transcript))
transcript$t_id <- rownames(transcript)
transmerge <- merge(transcripts_q_lt_05, transcript, all.x = T, all.y = F)
transmerge2 <- merge(transmerge, tfile, all.x = T, all.y = F)
rownames(transmerge2) <- transmerge2$t_name
writetrans <- transmerge2[,(names(transmerge2) %in% c("pval", "qval"))]
write.csv(writetrans, paste(cutoff, "_enoki_transcripts_qlt05_results_by_", covname, ".csv", sep = ""), row.names=T)
writetrans <- transmerge2[,!(names(transmerge2) %in% c("t_id", "feature", "pval", "qval", "t_name"))]
write.csv(writetrans, paste(cutoff, "_enoki_transcripts_qlt05_results_by_", covname, "_expr.csv", sep = ""), row.names=T)

print("DE data written")

## Plotting setup
#tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
#palette(tropical)

## Plotting gene abundance distribution
#fpkm <- texpr(bg_Dpoly, meas='FPKM')
#fpkm <- log2(fpkm +1)
#boxplot(fpkm, col=as.numeric(pheno_data$sex), las=2,ylab='log2(FPKM+1)')

## Plot individual transcripts
#ballgown::transcriptNames(bg_Dpoly)[12]
#plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),
#     main=paste(ballgown::geneNames(bg_Dpoly)[12], ' : ',ballgown::transcriptNames(bg_Dpoly)[12]),
#     pch=19, xlab="Sex", ylab='log2(FPKM+1)')
#points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)), col=as.numeric(pheno_data$sex))

## Plot gene of transcript 1729
#plotTranscripts(ballgown::geneIDs(bg_Dpoly)[1729], bg_Dpoly,
#                main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

## Plot average expression
#plotMeans(ballgown::geneIDs(bg_Dpoly)[203], bg_Dpoly_filt, groupvar="sex", legend=FALSE)

