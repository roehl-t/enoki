# Created by Thomas Roehl on 12/3/2021

# load command line arguments
args = commandArgs(trailingOnly=TRUE)

# assign arguments to variables
blastfile <- args[1] # .csv file of readable blast results
resultsfile <- args[2] # .csv file of DE data
outfile <- args[3] # name for the new file
workdir <- args[4] # where to set the working directory


setwd(workdir)

# create an output log
zz <- file("name_genes2_output_log.txt", open="wt")
sink(zz, type="message", append=T)

# read in the blast results
blast <- read.csv(blastfile, header = T)

# keep important columns: UniProt_id, protein_name, organism_name, gene_name, query_id
blast <- blast[,(names(blast) %in% c("UniProt_id", "protein_name", "organism_name", "gene_name", "query_id"))]

# read in the DE data
dedata <- read.csv(resultsfile, header = T)

# label first column (formerly row names) as query_id
colnames(dedata)[1] <- "query_id"

# merge tables based on row name, keep all DE data, but not all blast data
demerge <- merge(dedata, blast, all.x = T, all.y = F)

# write the resulting table
write.csv(demerge, outfile)
