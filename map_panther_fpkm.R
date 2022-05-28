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
groupsets <- strsplit(args[5], "|", fixed = T)[[1]] # list of conditions based on column of pheno: each column separated by a semicolon, then column title followed by a comma and a comma-separated list of conditions for that column -- make output files based on multiple groups of column conditions by separating each output set with a pipe (|)
outfiles <- args[6] # base name for the new files
workdir <- args[7] # where to set the working directory
logdir <- args[8] # where to write the log
recordfile <- args[9] # print record of a previous run, or "none" -- used to check for previously completed files

# list of analyses to perform:
#grouplist <- c("ROEHL", "myc", "sti", "pil", "gil", "you", "pri", "cul", "nor")

setwd(workdir)

# create an output log
zz <- file(paste(logdir, "/map_panther_fpkm_log.txt", sep = ""), open="wt")
sink(zz, type="message", append=T)

# load the record data
if(file.exists(recordfile) & recordfile != "none"){
    userecord <- T
    record <- read.csv(recordfile, sep = "\t", header = F) 
} else {
    userecord <- F
}


## set up sample titles
print("Loading and organizing data")

names(panther) <- c("NCBI_id", "panther_acc", "panther_family.subfamily", "HMM_e-value", "HMM_bitscore", "alignment_range")
panther <- panther[, names(panther) %in% c("NCBI_id", "panther_acc")]

keepcolumns <- append(pheno$ids, c("UniProt_id", "query_id", "protein_name"))
keepcolumns <- make.names(keepcolumns)

data <- data[, names(data) %in% keepcolumns]

## map PANTHER acc to fpkm table
print("Mapping PANTHER accession to fpkm table")

# map PANTHER/NCBI to NCBI/UniProtKB
panther_uniprot <- merge(panther, up_ncbi, all.x = T, all.y = F)
panther_uniprot <- panther_uniprot[,names(panther_uniprot) %in% c("UniProt_id", "panther_acc")]

mapped <- merge(data, panther_uniprot, all.x = T, all.y = F)

# assign anything that didn't match the PANTHER database the value "NOHIT"
mapped$panther_acc <- ifelse(is.na(mapped$panther_acc), "NOHIT", mapped$panther_acc)

    
# filter and combine columns of interest (ex: pileus)
    # use R to average expr values for each gene
print("Writing subset tables in PANTHER generic mapping file format")

# write complete data file
filename <- paste(workdir, "/", outfiles, "_complete_mappings.csv", sep = "")
message <- paste(filename, "_complete-PANTHER-mapping_written", sep = "")
skip = F
if(userecord){
    record$match <- grepl(message, record[,1], fixed = T)
    check <- record[record$match == T,]
    skip <- (nrow(check) > 0)
}
if(!skip){
    write.csv(mapped, filename)
    print(message)
}

groupnumber <- 0
for(groupset in groupsets){
    groupnumber <- groupnumber + 1
    print(paste("subset ", groupnumber, ": ", groupset), sep = "")
    
    filename <- paste(workdir, "/", outfiles, "_subset", groupnumber, "_fpkm_tables.csv", sep = "")
    message <- paste(filename, "_PANTHER-generic_written", sep = "")
    skip = F
    
        
    # check to see if this set has been completed
    if(userecord){
        record$match <- grepl(message, record[,1], fixed = T)
        check <- record[record$match == T,]
        skip <- (nrow(check) > 0)
    }
    if(!skip){
        
        # subset samples
        samplesets <- strsplit(groupset, ";", fixed = T)[[1]]
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
                phenosub <- pheno[(pheno$keep3 == T),]
                keepsamples <- phenosub$ids
                keepsamples <- c(keepsamples, c("query_id", "panther_acc"))
                keepsamples <- make.names(keepsamples)
                subset <- mapped[,(names(mapped) %in% keepsamples)]
            }
        }
        
        # log2 transform and average data
        mean <- rep(0, nrow(subset))
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
        write.table(condensed, filename, sep = "\t", row.names = F, col.names = F)
        print(message)
    
    }
}
