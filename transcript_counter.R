# created by Thomas Roehl April 2022

args = commandArgs(trailingOnly=TRUE)
input <- args[1] # where the directories containing .ctab files are - OR - the .csv file from a previous run
name <- args[2] # what to call this run
logloc <- args[3] # where to put the logs
outputloc <- args[4] # where to save the files
cutoff <- as.numeric(args[5]) # add any transcripts below this number to a removelist

setwd(outputloc)

# create log
zz <- file(paste(logloc, "/", name, "_trnascript_counter_log.txt", sep = ""), open="wt")
sink(zz, type="message")

# check whether input is .csv file (skip counting) or directory (do counting)
if(grepl(".csv", input, fixed = T)){
    skipcount <- T
    csvname <- input
} else {
    skipcount <- F
    data_dir_loc <- input
}

# count the number of distinct transcripts in each t_data file (unless a .csv file was specified)
if(skipcount == F){
    dirlist <- list.dirs(data_dir_loc, recursive = F)
    first <- T
    for(dir in dirlist){
        if(file.exists(paste(dir, "/t_data.ctab", sep = ""))){
            tfile <- read.table(paste(dir, "/t_data.ctab", sep = ""), header = T)
            tfile <- tfile[(tfile$FPKM > 0),]
            sample <- strsplit(dir, "/")[[1]]
            sample <- sample[length(sample)]
            tcount <- nrow(tfile)
            newrow <- data.frame(sample, tcount)
            colnames(newrow) <- c("ids", "transcript_count")
            if(first){
                countdata <- newrow
                first <- F
            } else {
                countdata <- merge(countdata, newrow, all.x = T, all.y = T)
            }
        }
    }
    # do not use name here so that the shell script can find the file
    csvname <- "transcript_counts.csv"
    write.csv(countdata, file = csvname, row.names = F)
}

# generate a removelist.txt file for the rnaseq_pipeline.sh script
countdata <- read.csv(csvname, headers = T)
countfilt <- countdata[(countdata$transcript_count > cutoff),]
outlist <- ""
first <- T
for(sample in countfilt$ids){
    if(first){
        outlist <- sample
        first <- F
    } else {
        outlist <- paste(outlist, " ", sample)
    }
}
# do not use name here so that the shell script can find the file
cat(outlist, file = "removelist.txt")
