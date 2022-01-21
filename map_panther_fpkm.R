#  Created by Thomas Roehl on 12/10/21.

library(dplyr)

# panther output is tab-delimited: sequence ID, panther acc, panther family/subfamily, HMM e-value, HMM bitscore, alignment range
# fpkm tables

# load command line arguments
args = commandArgs(trailingOnly=TRUE)

# assign arguments to variables
data <- read.csv(args[1]) # .csv file with fpkm data and gene names
panther <- read.csv(args[2], sep = "\t", header = F) # .csv file of PANTHER results
pheno <- read.csv(args[3]) # .csv of the pheno data, for column filtering
outfiles <- args[4] # base name for the new files
workdir <- args[5] # where to set the working directory

setwd(workdir)

# create an output log
zz <- file("map_panther_fpkm_log.txt", open="wt")
sink(zz, type="message", append=T)

# set up sample titles
print("Loading and organizing data")

names(panther) <- c("UniProt_id", "panther_acc", "panther_family.subfamily", "HMM_e-value", "HMM_bitscore", "alignment_range")
panther <- panther[, names(panther) %in% c("UniProt_id", "panther_acc")]

pheno$ids <- gsub("-", ".", pheno$ids, fixed = T)
pheno$titlename <- paste(pheno$ids, "_", substr(pheno$type, 1, 3), "_", substr(pheno$tissue, 1, 3), sep = "")
oldnames <- data.frame(names(data))
names(oldnames) <- c("ids")
phenomerge <- merge(pheno, oldnames, all.x = F, all.y = T)
phenomerge$titlename <- ifelse(is.na(phenomerge$titlename), phenomerge$ids, phenomerge$titlename)

keepcolumns <- phenomerge$ids
keepcolumns[length(keepcolumns)+1] <- "UniProt_id"

newnames <- phenomerge$titlename
newnames[length(newnames)+1] <- "UniProt_id"

# reorder and rename data titles
data <- data[, keepcolumns]
names(data) <- newnames


# map PANTHER acc to fpkm table
print("Mapping PANTHER accession to fpkm table")

mapped <- merge(data, panther, all.x = T, all.y = F)

    
# filter and combine columns of interest (ex: pileus)
    # use R to average expr values for each gene
print("Writing subset tables in PANTHER generic mapping file format")

# list of analyses to perform:
grouplist <- c("ROEHL", "myc", "sti", "pil", "gil", "you", "pri", "cul", "nor")

for(group in grouplist){
    # subset columns of interest
    subset <- select(data, contains(group))
    
    # log2 transform data
    for(i in 1:ncol(subset)){
        subset[,i] <- log2(subset[,i] + 1)
    }
    
    # average log2fpkm values for all remaining columns
    subset$mean <- rowMeans(subset)
    
    ################################## does the unique identifier have to match between the query list and reference list?
    # arrange columns in new dataframe -- first column must be a unique identifier
    newdata <- data.frame(paste(data$UniProt_id, "_", rownames(data), sep = ""), data$panther_acc, subset$mean)
    
    ################################### make sure panther_id that doesn't match HMM database is assigned a value of "NOHIT"
    # write tsv -- must be tab separated for PANTHER upload
    write.csv(newdata, paste(outfiles, "_", group, "_panther.tsv", sep = ""), sep = "\t", header = F, row.names = F)
}
