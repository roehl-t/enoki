#  Created by Thomas Roehl on 12/10/21.

# panther output is tab-delimited: sequence ID, panther acc, panther family/subfamily, HMM e-value, HMM bitscore, alignment range
# fpkm tables

# load command line arguments
args = commandArgs(trailingOnly=TRUE)

# assign arguments to variables
data <- read.csv(args[1], header = T) # .csv file with fpkm data and gene names
panther <- read.csv(args[2], sep = "\t", header = F) # .csv file of PANTHER results
up_ncbi <- read.csv(args[3], sep = ",", header = T) # .csv file of UniProtKB to NCBI mapping
pheno <- read.csv(args[4]) # .csv of the pheno data, for column filtering
groups <- args[5] # comma-separated list of subsets to write based on names of phenotypes or files (uses partial matches)
groups <- gsub("-", ".", groups, fixed = T)
groups <- gsub(" ", ".", groups, fixed = T)
grouplist <- strsplit(groups, ",")[[1]]
outfiles <- args[6] # base name for the new files
workdir <- args[7] # where to set the working directory

# list of analyses to perform:
#grouplist <- c("ROEHL", "myc", "sti", "pil", "gil", "you", "pri", "cul", "nor")

setwd(workdir)

# create an output log
zz <- file("map_panther_fpkm_log.txt", open="wt")
sink(zz, type="message", append=T)

## set up sample titles
print("Loading and organizing data")

names(panther) <- c("NCBI_id", "panther_acc", "panther_family.subfamily", "HMM_e-value", "HMM_bitscore", "alignment_range")
panther <- panther[, names(panther) %in% c("NCBI_id", "panther_acc")]

pheno$ids <- gsub("-", ".", pheno$ids, fixed = T)
pheno$type <- gsub(" ", ".", pheno$type, fixed = T)
pheno$titlename <- paste(pheno$ids, pheno$type, pheno$tissue, sep = "_")

sampletitles <- data.frame(pheno$ids, pheno$titlename)

samplenames <- data.frame(pheno$ids, pheno$titlename)
names(samplenames) <- c("ids", "titlename")
keepcolumns <- append(samplenames$ids, c("UniProt_id", "query_id", "protein_name"))

data <- data[, names(data) %in% keepcolumns]

## map PANTHER acc to fpkm table
print("Mapping PANTHER accession to fpkm table")

# map PANTHER/NCBI to NCBI/UniProtKB
panther_uniprot <- merge(panther, up_ncbi, all.x = T, all.y = F)
panther_uniprot <- panther_uniprot[,names(panther_uniprot) %in% c("UniProt_id", "panther_acc")]

mapped <- merge(data, panther_uniprot, all.x = T, all.y = F)

# assign anything that didn't match the PANTHER database the value "NOHIT"
mapped$panther_acc <- ifelse(is.na(mapped$panther_acc), "NOHIT", mapped$panther_acc)
#mapped$discard <- is.na(mapped$panther_acc)
#mapped <- mapped[mapped$discard == F,]
#head(mapped)

    
# filter and combine columns of interest (ex: pileus)
    # use R to average expr values for each gene
print("Writing subset tables in PANTHER generic mapping file format")

# write complete data file
write.csv(mapped, paste(workdir, "/", outfiles, "_complete_mappings.csv", sep = ""))

for(group in grouplist){
    # subset columns of interest
    samplenames$keep <- rep(T, nrow(samplenames))
    samplenames$keep <- grepl(group, samplenames$titlename, fixed = T)
    subsetnames <- samplenames[samplenames$keep == T,]
    keepnames <- append(subsetnames$ids, c("query_id", "panther_acc"), nrow(subsetnames))
    subset <- mapped[,names(mapped) %in% keepnames]
    mean <- rep(0, nrow(subset))
    
    # log2 transform and average data
    colcounter = 0
    for(i in 1:nrow(subset)){
        for(j in 1:ncol(subset)){
            if(is.numeric(subset[i,j])){
                if(i == 1){
                    colcounter = colcounter + 1
                }
                mean[i] <- mean[i] + log2(subset[i,j] + 1)
            }
        }
    }
    subset$mean <- mean/colcounter
    
    # arrange columns in new dataframe -- first column must be a unique identifier
    newdata <- data.frame(subset$query_id, subset$panther_acc, subset$mean)
    names(newdata) <- c("query_id", "panther_acc", "mean")

    # remove duplicate MSTRGs
    condensed <- newdata[1,]
    if(nrow(newdata) > 1){
        for(i in 2:nrow(newdata)){
            found = F
            for(j in 1:nrow(condensed)){
                if(newdata[i, "query_id"] == condensed[j, "query_id"]){
                    found = T
                }
            }
            if(!found){
                condensed <- rbind(condensed, newdata[i,])
            }
        }
    }
    
    # write csv -- must be tab separated for PANTHER upload
    write.table(condensed, paste(workdir, "/", outfiles, "_", group, "_fpkm_tables.csv", sep = ""), sep = "\t", row.names = F, col.names = F)
}
