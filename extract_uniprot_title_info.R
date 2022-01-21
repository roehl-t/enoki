# convert BLASTX UniProt output to useable information
# Created by Thomas Roehl on 12/3/2021

args = commandArgs(trailingOnly = T)

# .csv file separated by tags containing UniProt BLAST output
blastfile <- args[1]
# where to write the modified file
outfile <- args[2]
# where to set the working directory
workdir <- args[3]

setwd(workdir)

zz <- file("extract_uniprot_title_info_output_log.txt", open="wt")
sink(zz, type="message", append=T)

# pattern db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
# pattern for isoforms sp|IsoID|EntryName Isoform IsoformName of ProteinName OS=OrganismName OX=OrganismIdentifier[ GN=GeneName]

# examples
#s1 <- "sp|P27748|ACOX_CUPNH Acetoin catabolism protein X OS=Cupriavidus necator (strain ATCC 17699 / H16 / DSM 428 / Stainer 337) OX=381666 GN=acoX PE=4 SV=2"
#s2 <- "sp|P04224|HA22_MOUSE H-2 class II histocompatibility antigen, E-K alpha chain OS=Mus musculus OX=10090 PE=1 SV=1"
#s3 <- "sp|Q4R572-2|1433B_MACFA Isoform Short of 14-3-3 protein beta/alpha OS=Macaca fascicularis OX=9541 GN=YWHAB"
#strs <- c(s1, s2, s3)

blast <- read.csv(blastfile, header = F, sep = "\t")
colnames(blast) <- c("query_id", "subject_title", "subject_accession", "e_value", "percent_identity", "bitscore", "length", "query_start", "query_end", "subject_start", "subject_end")

strs <- blast$subject_title

first <- T

for(s in strs){
	subject_title <- s
	# extract database (either sp for SwissProt or tr for TrEMBL)
		# first two letters
	database <- substr(s, 1, 2)
	
	# extract identifier
		# between | and |
	markers <- gregexpr("|", s, fixed = T)
	UniProt_id <- substr(s, markers[[1]][1]+1, markers[[1]][2]-1)
	
	# extract entry name
		# between second | and first " "
	marker1 <- gregexpr("|", s, fixed = T)[[1]][2]
	marker2 <- gregexpr(" ", s, fixed = T)[[1]][1]
	entry_name <- substr(s, marker1 + 1, marker2 - 1)
	
	# extract isoform name (may not be present)
		# between "Isoform " and " of "
	isoform <- grepl("Isoform", s)
	if(isoform){
		marker1 <- gregexpr("Isoform ", s, fixed = T)[[1]][1]
		marker2 <- gregexpr(" of ", s, fixed = T)[[1]][1]
		isoform_name <- substr(s, marker1 + 8, marker2 - 1)
	} else {
		isoform_name <- NA
	}
	
	# extract protein name (may come after isoform name and the word " of ")
		# if an isoform, between " of " and " OS="
		# if not an isoform, between first " " and " OS="
	if(isoform){
		marker1 <- gregexpr(" of ", s, fixed = T)[[1]][1]
		marker2 <- gregexpr(" OS=", s, fixed = T)[[1]][1]
		protein_name <- substr(s, marker1 + 4, marker2 - 1)
	} else {
		marker1 <- gregexpr(" ", s, fixed = T)[[1]][1]
		marker2 <- gregexpr(" OS=", s, fixed = T)[[1]][1]
		protein_name <- substr(s, marker1 + 1, marker2 - 1)
	}
	
	# extract organism name
		# between "OS=" and " OX="
	marker1 <- gregexpr("OS=", s, fixed = T)[[1]][1]
	marker2 <- gregexpr(" OX=", s, fixed = T)[[1]][1]
	organism_name <- substr(s, marker1 + 3, marker2 - 1)
	
	# extract organism identifier
		# between "OX=" and next " " or str end
	marker1 <- gregexpr("OX=", s, fixed = T)[[1]][1]
	ss <- substr(s, marker1, nchar(s))
	if(grepl(" ", ss, fixed = T)){
		marker2 <- gregexpr(" ", ss, fixed = T)[[1]][1]
		organism_id <- substr(ss, 4, marker2 -1)
	} else {
		marker2 <- nchar(s)
		organism_id <- substr(s, marker1 + 3, marker2)
	}
	
	# extract gene name (may not be present, or may be in different locations)
		# between "GN= and next " " or str end
	if(grepl("GN=", s, fixed = T)){
		marker1 <- gregexpr("GN=", s, fixed = T)[[1]][1]
		ss <- substr(s, marker1, nchar(s))
		if(grepl(" ", ss, fixed = T)){
			marker2 <- gregexpr(" ", ss, fixed = T)[[1]][1]
			gene_name <- substr(ss, 4, marker2 -1)
		} else {
			marker2 <- nchar(s)
			gene_name <- substr(s, marker1 + 3, marker2)
		}
	} else {
		gene_name <- NA
	}
	
	# extract extract protein existence number (may not be present)
		# between "PE=" and " SV="
	if(grepl("PE=", s, fixed = T)){
		marker1 <- gregexpr("PE=", s, fixed = T)[[1]][1]
		marker2 <- gregexpr(" SV=", s, fixed = T)[[1]][1]
		protein_existence <- substr(s, marker1 + 3, marker2 -1)
	} else {
		protein_existence <- NA
	}
	
	# extract sequence version (may not be present)
		# after "SV="
	if(grepl("SV=", s, fixed = T)){
		marker1 <- gregexpr("SV=", s, fixed = T)[[1]][1]
		version <- substr(s, marker1 + 3, nchar(s))
	} else {
		version <- NA
	}
	
	# collect values
	values <- data.frame(subject_title, database, UniProt_id, entry_name, isoform_name, protein_name, organism_name, organism_id, gene_name, protein_existence, version)
	head(values)
	if(first){
		dframe <- values
		first = F
	} else {
		mframe <- merge(dframe, values, all.x = T, all.y = T)
		dframe <- mframe
	}
}

rewrite <- merge(dframe, blast, all.x = T, all.y = T)
rewrite <- rewrite[,!(names(rewrite) %in% c("subject_title", "subject_accession"))]
write.csv(rewrite, outfile)
