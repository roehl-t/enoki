# created by Thomas Roehl April 2022


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
