# Created by Thomas Roehl on 11/30/21

# load command line arguments
args = commandArgs(trailingOnly=TRUE)

# assign arguments to variables
blastfile <- args[1] # .csv file of blast results
datafolderloc <- args[2] # where to find the folders of tables
datafolders <- list.dirs(datafolderloc, recursive = F)
workdir <- args[3] # where to set the working directory
runname <- args[4] # name of this run to be used in the log file name

setwd(workdir)

# create an output log
zz <- file(paste("./", runname, "_remove_rRNA_output_log.txt", sep = ""), open="wt")
sink(zz, type="message")

# read in blast results
# table not labeled, results in order specified when BLASTed
blast <- read.csv(blastfile, sep = "\t", header = F)

# label query name column with "t_name" and subject title column with "gene_name"
colnames(blast) <- c("t_name","gene_name","sacc","evalue","pident","bitscore",
                     "length","qstart","qend","sstart","send")

# filter blast columns to keep only important ones
blast <- blast[,(colnames(blast) %in% c("t_name", "gene_name"))]

# for each sample folder...
for(folder in datafolders){
    #check to see if the folder contains ctab files, if not it's the wrong folder
  if(file.exists(paste(folder, "/t_data.ctab", sep = ""))){
      print(paste(folder, "/t_data.ctab", sep = ""))
      # read in t_data.ctab (contains MSTRG numbers and t_id numbers)
      # MSTRNG values contained in "t_name", gene names contained in "gene_name"
      tctab <- read.table(paste(folder, "/t_data.ctab", sep = ""), header=T)
      tctab <- tctab[,(colnames(tctab) %in% c("t_name","t_id"))]
      
      # create new table listing which t_id values to remove
      rrna_tids <- merge(tctab, blast, all.x = F, all.y = F)
      
      # add in keep column, populate with TRUE, and keep only t_id and keep columns
      rrna_tids$keep <- rep(T, nrow(rrna_tids))
      rrna_tids <- rrna_tids[,(colnames(rrna_tids) %in% c("t_id", "keep"))]
      
      # for the three tables with t_id columns...
      tablefiles <- c("e2t.ctab", "i2t.ctab", "t_data.ctab")
      for(tablefile in tablefiles){
        
        name <- paste(folder, "/", tablefile, sep = "")
        
        # read in the table
        table <- read.table(name, header=T)
        
        # merge the table and the rrna_tids table
        merged <- merge(table, rrna_tids, all.x = T, all.y = F)
        
        # use is.na to fill in all rows of the keep column
        merged[,"keep"] <- is.na(merged[,"keep"])
        
        # remove rows that are false for keep
        merged <- merged[(merged[,"keep"] == T),]
        
        # remove keep column
        merged <- merged[,!(colnames(merged) %in% c("keep"))]
        
        # write modified files
        write.table(merged, name, quote = F, row.names = F, sep = "\t")
      }
      
      # for _data.ctab files...
      transtabs <- c("e2t.ctab", "i2t.ctab")
      eitabs <- c("e_data.ctab", "i_data.ctab")
      
      for(i in 1:length(transtabs)){
        # read in *_data and *2t
        tablename <- paste(folder, "/", eitabs[i], sep = "")
        rrna_tids <- read.table(paste(folder, "/", transtabs[i], sep = ""), header = T)
        table <- read.table(tablename, header = T)
        
        # merge tables, keep only where both have *_ids
        merged <- merge(table, rrna_tids, all.x = F, all.y = F)
        
        # remove extra t_id column
        merged <- merged[,!(colnames(merged) %in% c("t_id"))]
        
        # write table
        write.table(merged, tablename, quote = F, row.names = F, sep = "\t")
      }
  }
}
