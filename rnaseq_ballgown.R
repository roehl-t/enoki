# from Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095
# modified by Todd Osmundson and Thomas Roehl 2021

#!/usr/bin/env Rscript
# run this in the output directory for rnaseq_pipeline.sh
# passing the pheno data csv file as the only argument 
args = commandArgs(trailingOnly=TRUE)

pheno_data_file <- args[1] # where the phenotype data file is
data_dir_loc <- args[2] # where the .ctab files are
setname <- args[3] # what to call this run
covname <- args[4] # phenotype data column to use as the covariate
adjnames <- strsplit(args[5], " ")[[1]] # phenotype data columns to use as the confounding variables
pcapairs <- strsplit(args[6], ";")[[1]] # pairs (or singles) of phenotype data columns to be used in generating pca plots
wrkdir <- args[7] # where to save the files
logloc <- args[8] # where to save the log
cutoff <- 0 # ignore files with less than this number of samples (already done in rnaseq_pipeline.sh, but the option is here in case you want to run just this file -- use 0 to skip the cutoff steps)

setwd(wrkdir)

# create an output log
zz <- file(paste(logloc, "/", setname, "_ballgown_output_log.txt", sep = ""), open="wt")
sink(zz, type="message")

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
library(ggplot2)
library(ggrepel)


## manual controls
#pheno_data_file <- "~/R/Flammulina-velutipes/ballgown/fv_sampling_data.csv"
#data_dir_loc <- "C:/Users/mycol/Documents/UW-La Crosse/Flammulina velutipes development/Data/all_ballgown_2022-01-18"
#setname <- "younor"
#covname <- "type"
#adjnames <- c("tissue")
#pcapairs <- c("batch jar", "position", "tissue type")
#wrkdir <- "~/R/Flammulina-velutipes/ballgown/bg_output_2022-01-18/younor"
#logloc <- "~/R/Flammulina-velutipes/ballgown/bg_output_2022-01-18/younor"
#cutoff <- 300
#setwd(wrkdir)


if(cutoff > 0){
    cutofftext <- paste(cutoff, "_", sep = "")
} else {
    cutofftext <- ""
}

# add ballgown version information to version_info.txt if it hasn't already been logged
versionlogfile <- paste(logloc, "/version_info.txt", sep = "")
if(file.exists(versionlogfile)){
    ballgownversion <- paste("ballgown v", packageVersion("ballgown"), sep = "")
    versionlog <- read.delim(versionlogfile)
    if(! grepl(ballgownversion, versionlog)){
        write(ballgownversion, file = versionlogfile, append = T)
    }
} else {
    print("could not open version info file to write ballgown version information")
}


## Read phenotype sample data
pheno_data <- read.csv(pheno_data_file, header=TRUE, stringsAsFactors=F)
print("Phenotypic data read")


# done in rnaseq_pipeline.sh, but included here for user convenience
if(cutoff > 0){
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

    print(paste("Removed samples moved to ", data_dir_loc, "/hidden", sep = ""))
}


# double check that only samples in the directory are included in the pheno data
    # very important, not done elsewhere
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


## Read in expression data
bg_mojo <- ballgown(dataDir = data_dir_loc, samplePattern="ROEHL-FV-", pData=pheno_data)
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
bg_mojo_filt <- subset(bg_mojo, "rowVars(texpr(bg_mojo)) > 1", genomesubset=TRUE)
print("Low abundance genes filtered")


## Generate PCA plots
# begin with ballgown object
bg_obj <- bg_mojo_filt

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
sink(paste("./PCA/", cutofftext, "mojo_PCA_results.txt", sep = ""))
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

for(pcapair in pcapairs){
    pcanames <- strsplit(pcapair, " ")[[1]]
    pair1 <- pcanames[1]
    pca_expr_vectors[,pair1] <- as.factor(pca_expr_vectors[,pair1])
    if(length(pcanames) > 1){
        pair2 <- pcanames[2]
        joiner <- "-"
        pca_expr_vectors[,pair2] <- as.factor(pca_expr_vectors[,pair2])
    } else {
        pair2 <- ""
        joiner <- ""
    }
    pairname <- paste(pair1, joiner, pair2, sep = "")
    
    pdf(file = paste("./PCA/", cutofftext, "mojo_pca_", pairname, ".pdf", sep = ""), width = 6, height = 6)
    if(length(pcanames) > 1){
        pcaplot <- ggplot(pca_expr_vectors, aes_string(x = "PC1", y = "PC2", color = pair1, shape = pair2)) + geom_vline(xintercept = 0, color = "grey") + geom_hline(yintercept = 0, color = "grey") + geom_point() + theme_classic() + theme(text = element_text(family = "serif")) + geom_text_repel(aes(label = ids), family = "serif")
    } else {
        pcaplot <- ggplot(pca_expr_vectors, aes_string(x = "PC1", y = "PC2", color = pair1)) + geom_vline(xintercept = 0, color = "grey") + geom_hline(yintercept = 0, color = "grey") + geom_point() + theme_classic() + theme(text = element_text(family = "serif")) + geom_text_repel(aes(label = ids), family = "serif")
    }
    print(pcaplot)
    dev.off()
}

print("PCA plots created")


## DE by transcript
results_transcripts <-  stattest(bg_mojo_filt, feature='transcript', covariate=covname,
         adjustvars=adjnames, getFC=TRUE, meas='FPKM')

## DE by gene
results_genes <-  stattest(bg_mojo_filt, feature='gene', covariate=covname, adjustvars=adjnames, getFC=TRUE, meas='FPKM')

## Sort results from smallest p-value
results_transcripts <- arrange(results_transcripts, pval)
results_genes <-  arrange(results_genes, pval)

## Write results to CSV
write.csv(results_transcripts, paste(cutofftext, "mojo_transcripts_results_by_", covname, ".csv", sep = ""), row.names=FALSE)
write.csv(results_genes, paste(cutofftext, "mojo_genes_results_by_", covname, ".csv", sep = ""), row.names=FALSE)

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

# To match up with transcriptome, you need to use an isoform ID
# Every gene has at least one isoform, so add ".1" to each gene ID
genes_q_lt_05$gene_id <- paste(genes_q_lt_05$id, ".1", sep = "")
results_genes$gene_id <- paste(results_genes$id, ".1", sep = "")

## Write lists of DEGs with q<0.05 (for gene ID) separated by new line
write(tmerge$t_name, file = paste("./lists/", cutofftext, "mojo_transcripts_qlt05_results_by_", covname, ".txt", sep = ""), sep="\n")
write(genes_q_lt_05$gene_id, file = paste("./lists/", cutofftext, "mojo_genes_qlt05_results_by_", covname, ".txt", sep = ""), sep="\n")

print("DE lists written")

## Write expression data for DEs
gene <- gexpr(bg_mojo_filt)
gene <- data.frame(gene)
colnames(gene) <- sub("FPKM.", "", colnames(gene))
gene$id <- rownames(gene)
genemerge <- merge(genes_q_lt_05, gene, all.x = T, all.y = F)
rownames(genemerge) <- genemerge$gene_id
writegenes <- genemerge[,(names(genemerge) %in% c("pval", "qval"))]
write.csv(writegenes, paste(cutofftext, "mojo_genes_qlt05_results_by_", covname, ".csv", sep = ""), row.names=T)
writegenes <- genemerge[,!(names(genemerge) %in% c("id", "feature", "pval", "qval", "gene_id"))]
write.csv(writegenes, paste(cutofftext, "mojo_genes_qlt05_results_by_", covname, "_expr.csv", sep = ""), row.names=T)

transcript <- texpr(bg_mojo_filt)
transcript <- data.frame(transcript)
colnames(transcript) <- sub("FPKM.", "", colnames(transcript))
transcript$t_id <- rownames(transcript)
transmerge <- merge(transcripts_q_lt_05, transcript, all.x = T, all.y = F)
transmerge2 <- merge(transmerge, tfile, all.x = T, all.y = F)
rownames(transmerge2) <- transmerge2$t_name
writetrans <- transmerge2[,(names(transmerge2) %in% c("pval", "qval"))]
write.csv(writetrans, paste(cutofftext, "mojo_transcripts_qlt05_results_by_", covname, ".csv", sep = ""), row.names=T)
writetrans <- transmerge2[,!(names(transmerge2) %in% c("t_id", "feature", "pval", "qval", "t_name"))]
write.csv(writetrans, paste(cutofftext, "mojo_transcripts_qlt05_results_by_", covname, "_expr.csv", sep = ""), row.names=T)

print("DE data written")

## Write expression data of entire transcriptome for later GO analysis
gene <- gexpr(bg_mojo)
gene <- data.frame(gene)
colnames(gene) <- sub("FPKM.", "", colnames(gene))
gene$id <- rownames(gene)
genemerge <- merge(results_genes, gene, all.x = T, all.y = F)
rownames(genemerge) <- genemerge$gene_id
writegenes <- genemerge[,!(names(genemerge) %in% c("id", "feature", "pval", "qval"))]
write.csv(writegenes, paste(cutofftext, "mojo_genes_results_by_", covname, "_expr.csv", sep = ""), row.names=T)
write(rownames(writegenes), file = paste("./lists/", cutofftext, "mojo_genes_results_by_", covname, ".txt", sep = ""), sep="\n")

transcript <- texpr(bg_mojo)
transcript <- data.frame(transcript)
colnames(transcript) <- sub("FPKM.", "", colnames(transcript))
transcript$id <- rownames(transcript)
transmerge <- merge(results_transcripts, transcript, all.x = T, all.y = F)
names(transmerge)[names(transmerge) == "id"] <- "t_id"
transmerge2 <- merge(transmerge, tfile, all.x = T, all.y = F)
rownames(transmerge2) <- transmerge2$t_name
writetrans <- transmerge2[,!(names(transmerge2) %in% c("t_id", "feature", "pval", "qval", "t_name"))]
write.csv(writetrans, paste(cutofftext, "mojo_transcripts_results_by_", covname, "_expr.csv", sep = ""), row.names=T)
write(rownames(writetrans), file = paste("./lists/", cutofftext, "mojo_transcripts_results_by_", covname, ".txt", sep = ""), sep="\n")

print("transcriptome expression data written")
