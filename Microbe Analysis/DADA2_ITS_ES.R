# Raw data ITS: /data/primary/ME/raaijmakers/PROMISE/GetahunB/Amplicons/EthiopianSoils/BaseClear_data_Project_099037/raw_sequences/Fungi_ITS/Ethiopia


library("devtools")
library("dada2")
library("dplyr")
library("tidyverse")

setwd("~/PROMISE/EthiopianSoils")
path <- "/data/primary/ME/raaijmakers/PROMISE/GetahunB/Amplicons/EthiopianSoils/BaseClear_data_Project_099037/raw_sequences/Fungi_ITS/Ethiopia"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_ITS"), `[`, 1)


# Inspect quality profiles ------------------------------------------------

plotQualityProfile(fnFs[1:3]) #Need to keep only the first 250
plotQualityProfile(fnRs[1:3]) #Need to keep only the first 200

#-#-#-# From here, I followed https://github.com/ErnakovichLab/dada2_ernakovichlab/tree/split_for_premise
project.fp <- "~/PROMISE/EthiopianSoils/ITS/"
preprocess.fp <- file.path(project.fp, "01_preprocess")
filtN.fp <- file.path(preprocess.fp, "filtN")
trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter") 
table.fp <- file.path(project.fp, "03_tabletax") 


# Filter and trim ------------------------------------------------------------------
# Name the N-filtered files to put them in filtN/ subdirectory
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs))
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))

# Filter Ns from reads and put them into the filtN directory
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE) 

# Set up the primer sequences to pass along to cutadapt
## These are not the primers because 1) the primers still show up in the quality profile plot, 2) almost all the reads were deleted after (e.g. from 43k to 1.5k)
#FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## ITS1f
#REV <- "GCTGCGTTCTTCATCGATGC"  ## ITS2

# After manually looking for various ITS primers
FWD <- "GCATCGATGAAGAACGCAGC"  ## ITS3
REV <- "TCCTCCGCTTATTGATATGC"  ## ITS4




# Write a function that creates a list of all orientations of the primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Save the primer orientations to pass to cutadapt
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Write a function that counts how many time primers appear in a sequence
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

library("ShortRead")
library("dplyr")
library("tidyr")
library("Hmisc")
library("ggplot2")
library("plotly")

# To check if there are primers and where exactly they are (this will help the decision-making for cutadapt)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[2]]))

## 1st pair of primers tried
#Forward Complement Reverse RevComp
#FWD.ForwardReads       0          0       0       0
#FWD.ReverseReads       0          0       0       0
#REV.ForwardReads       0          0       0   29815
#REV.ReverseReads     221          0       0       0

## 2nd pair (ITS3 & ITS4)
#Forward Complement Reverse RevComp
#FWD.ForwardReads   29815          0       0       0
#FWD.ReverseReads       0          0       0     221
#REV.ForwardReads       0          0       0     478
#REV.ReverseReads   28894          0       0       0

# Remove primers with cutadapt and assess the output ----------------------

# Create directory to hold the output from cutadapt
if (!dir.exists(trimmed.fp)) dir.create(trimmed.fp)
fnFs.cut <- file.path(trimmed.fp, basename(fnFs))
fnRs.cut <- file.path(trimmed.fp, basename(fnRs))

# Save the reverse complements of the primers to variables
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

##  Create the cutadapt flags ##
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 50") 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 50") 

cutadapt <- "/data/shared/conda/envs/QC-NGS/bin/cutadapt"
# Run Cutadapt
for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-j", 24, R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# As a sanity check, we will check for primers in the first cutadapt-ed sample:
## should all be zero! (it is!)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Check the quality profile again
plotQualityProfile(fnFs.cut[1:3])
plotQualityProfile(fnRs.cut[1:3])


## Now start DADA2 pipeline
# Put filtered reads into separate sub-directories for big data workflow
dir.create(filter.fp)
subF.fp <- file.path(filter.fp, "preprocessed_F") 
subR.fp <- file.path(filter.fp, "preprocessed_R") 
dir.create(subF.fp)
dir.create(subR.fp)

# Move R1 and R2 from trimmed to separate forward/reverse sub-directories
fnFs.Q <- file.path(subF.fp,  basename(fnFs)) 
fnRs.Q <- file.path(subR.fp,  basename(fnRs))
file.copy(from = fnFs.cut, to = fnFs.Q)
file.copy(from = fnRs.cut, to = fnRs.Q)

# File parsing; create file names and make sure that forward and reverse files match
filtpathF <- file.path(subF.fp, "filtered") # files go into preprocessed_F/filtered/
filtpathR <- file.path(subR.fp, "filtered") # ...
fastqFs <- sort(list.files(subF.fp, pattern="fastq.gz"))
fastqRs <- sort(list.files(subR.fp, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")


## FILTER AND TRIM FOR QUALITY
# For this, it is important to check the plots first and see which parts need trimming
# If the number of samples is 20 or less, plot them all, otherwise, just plot 20 randomly selected samples
if( length(fastqFs) <= 20) {
  fwd_qual_plots <- plotQualityProfile(paste0(subF.fp, "/", fastqFs))
  rev_qual_plots <- plotQualityProfile(paste0(subR.fp, "/", fastqRs))
} else {
  rand_samples <- sample(size = 20, 1:length(fastqFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(paste0(subF.fp, "/", fastqFs[rand_samples]))
  rev_qual_plots <- plotQualityProfile(paste0(subR.fp, "/", fastqRs[rand_samples]))
}

fwd_qual_plots ## Remove the last 30 (length drops)
rev_qual_plots ## Keep 210 (remove 90) (it may be too much, will see in the merging part)

# write plots to disk
saveRDS(fwd_qual_plots, paste0(filter.fp, "/fwd_qual_plots.rds"))
saveRDS(rev_qual_plots, paste0(filter.fp, "/rev_qual_plots.rds"))

ggsave(plot = fwd_qual_plots, filename = paste0(filter.fp, "/fwd_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = rev_qual_plots, filename = paste0(filter.fp, "/rev_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")

# NOW we can filter
filt_out <- filterAndTrim(fwd=file.path(subF.fp, fastqFs), filt=file.path(filtpathF, fastqFs),
                          rev=file.path(subR.fp, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                          truncLen=c(270,210), maxN=0, rm.phix=TRUE,
                          compress=TRUE, verbose=TRUE, multithread=FALSE) ### If they're still too little, I can increase MaxEE or reduce TruncQ (original=2)

#Read in 42117 paired-sequences, output 41626 (98.8%) filtered paired-sequences.
#Read in 33477 paired-sequences, output 32990 (98.5%) filtered paired-sequences.
#Read in 54950 paired-sequences, output 54633 (99.4%) filtered paired-sequences.
#Read in 49427 paired-sequences, output 49055 (99.2%) filtered paired-sequences.
#Read in 60097 paired-sequences, output 59510 (99%) filtered paired-sequences.
#Read in 29215 paired-sequences, output 28200 (96.5%) filtered paired-sequences.
#Read in 48211 paired-sequences, output 47760 (99.1%) filtered paired-sequences.
#Read in 66114 paired-sequences, output 65787 (99.5%) filtered paired-sequences.
#Read in 52249 paired-sequences, output 52063 (99.6%) filtered paired-sequences.
#Read in 48849 paired-sequences, output 48459 (99.2%) filtered paired-sequences.
#Read in 55961 paired-sequences, output 55508 (99.2%) filtered paired-sequences.
#Read in 56394 paired-sequences, output 56030 (99.4%) filtered paired-sequences.
#Read in 48331 paired-sequences, output 48065 (99.4%) filtered paired-sequences.
#Read in 33970 paired-sequences, output 33807 (99.5%) filtered paired-sequences.
#Read in 29106 paired-sequences, output 28584 (98.2%) filtered paired-sequences.
#Read in 41928 paired-sequences, output 41506 (99%) filtered paired-sequences.
#Read in 39411 paired-sequences, output 39216 (99.5%) filtered paired-sequences.
#Read in 62886 paired-sequences, output 62397 (99.2%) filtered paired-sequences.
#Read in 44035 paired-sequences, output 43593 (99%) filtered paired-sequences.
#Read in 56932 paired-sequences, output 56212 (98.7%) filtered paired-sequences.
#Read in 55221 paired-sequences, output 54895 (99.4%) filtered paired-sequences.
#Read in 38179 paired-sequences, output 37860 (99.2%) filtered paired-sequences.
#Read in 47509 paired-sequences, output 47307 (99.6%) filtered paired-sequences.
#Read in 49338 paired-sequences, output 48466 (98.2%) filtered paired-sequences.
#Read in 49633 paired-sequences, output 48050 (96.8%) filtered paired-sequences.
#Read in 38691 paired-sequences, output 38434 (99.3%) filtered paired-sequences.
#Read in 32462 paired-sequences, output 32330 (99.6%) filtered paired-sequences.
#Read in 35823 paired-sequences, output 35372 (98.7%) filtered paired-sequences.
#Read in 27819 paired-sequences, output 27688 (99.5%) filtered paired-sequences.
#Read in 36856 paired-sequences, output 36071 (97.9%) filtered paired-sequences.
#Read in 21693 paired-sequences, output 20940 (96.5%) filtered paired-sequences.
#Read in 49332 paired-sequences, output 48708 (98.7%) filtered paired-sequences.
#Read in 40339 paired-sequences, output 39017 (96.7%) filtered paired-sequences.
#Read in 51122 paired-sequences, output 50676 (99.1%) filtered paired-sequences.
#Read in 28514 paired-sequences, output 28223 (99%) filtered paired-sequences.
#Read in 41313 paired-sequences, output 40035 (96.9%) filtered paired-sequences.
#Read in 55647 paired-sequences, output 55357 (99.5%) filtered paired-sequences.
#Read in 35930 paired-sequences, output 35503 (98.8%) filtered paired-sequences.
#Read in 35438 paired-sequences, output 35095 (99%) filtered paired-sequences.
#Read in 35664 paired-sequences, output 35345 (99.1%) filtered paired-sequences.
#Read in 33085 paired-sequences, output 32767 (99%) filtered paired-sequences.
#Read in 35493 paired-sequences, output 34867 (98.2%) filtered paired-sequences.
#Read in 41074 paired-sequences, output 40931 (99.7%) filtered paired-sequences.
#Read in 47008 paired-sequences, output 46670 (99.3%) filtered paired-sequences.
#Read in 36734 paired-sequences, output 36428 (99.2%) filtered paired-sequences.
#Read in 46770 paired-sequences, output 46584 (99.6%) filtered paired-sequences.
#Read in 54397 paired-sequences, output 54003 (99.3%) filtered paired-sequences.
#Read in 35254 paired-sequences, output 35171 (99.8%) filtered paired-sequences.
#Read in 53752 paired-sequences, output 53392 (99.3%) filtered paired-sequences.
#Read in 44886 paired-sequences, output 44604 (99.4%) filtered paired-sequences.
#Read in 51820 paired-sequences, output 51572 (99.5%) filtered paired-sequences.
#Read in 25214 paired-sequences, output 24782 (98.3%) filtered paired-sequences.

write.csv(filt_out, "~/PROMISE/EthiopianSoils/ITS/readsInAndOut_DADA2.txt")


# summary of samples in filt_out by percentage
filt_out %>% 
  data.frame() %>% 
  mutate(Samples = rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything()) %>%
  summarise(min_remaining = paste0(round(min(percent_kept), 2), "%"), 
            median_remaining = paste0(round(median(percent_kept), 2), "%"),
            mean_remaining = paste0(round(mean(percent_kept), 2), "%"), 
            max_remaining = paste0(round(max(percent_kept), 2), "%"))
#min_remaining median_remaining mean_remaining max_remaining
#1        96.53%           99.17%         98.89%        99.76%
  
# If the number of samples greater than 20 figure out which samples, if any, have been filtered out
# so we won't try to plot them, otherwise just plot all the samples that remain
if( length(fastqFs) <= 20) {
  remaining_samplesF <-  fastqFs[
    which(fastqFs %in% list.files(filtpathF))] # keep only samples that haven't been filtered out
  remaining_samplesR <-  fastqRs[
    which(fastqRs %in% list.files(filtpathR))] # keep only samples that haven't been filtered out
  
  fwd_qual_plots_filt <- plotQualityProfile(paste0(filtpathF, "/", remaining_samplesF))
  rev_qual_plots_filt <- plotQualityProfile(paste0(filtpathR, "/", remaining_samplesR))
} else {
  remaining_samplesF <-  fastqFs[rand_samples][
    which(fastqFs[rand_samples] %in% list.files(filtpathF))] # keep only samples that haven't been filtered out
  remaining_samplesR <-  fastqRs[rand_samples][
    which(fastqRs[rand_samples] %in% list.files(filtpathR))] # keep only samples that haven't been filtered out
  fwd_qual_plots_filt <- plotQualityProfile(paste0(filtpathF, "/", remaining_samplesF))
  rev_qual_plots_filt <- plotQualityProfile(paste0(filtpathR, "/", remaining_samplesR))
}

fwd_qual_plots_filt
rev_qual_plots_filt

# write plots to disk
saveRDS(fwd_qual_plots_filt, paste0(filter.fp, "/fwd_qual_plots_filt.rds"))
saveRDS(rev_qual_plots_filt, paste0(filter.fp, "/rev_qual_plots_filt.rds"))

ggsave(plot = fwd_qual_plots_filt, filename = paste0(filter.fp, "/fwd_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = rev_qual_plots_filt, filename = paste0(filter.fp, "/rev_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")


## INFER SEQUENCE VARIANTS
# Housekeeping step - set up and verify the file names for the output
# File parsing
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# Sample names in order
sample.names <- basename(filtFs) # doesn't drop fastq.gz
sample.names <- gsub("_R1_001_000000000-BRLJM.filt.fastq.gz", "", sample.names)
sample.namesR <- basename(filtRs) # doesn't drop fastq.gz 
sample.namesR <- gsub("_R2_001_000000000-BRLJM.filt.fastq.gz", "", sample.namesR)

# Double check
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# Learn the error rates ---------------------------------------------------

errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE)

plotErrors(errF, nominalQ=TRUE)


# Sample inference --------------------------------------------------------

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
# 445 sequence variants were inferred

# Merge paired reads ------------------------------------------------------

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# Construct sequence table ------------------------------------------------

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#52, 11419, 52 samples and 11419 ASVs

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


# Remove chimeras ---------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 670 bimeras out of 11419 input sequences.
dim(seqtab.nochim)
# 52  10749
sum(seqtab.nochim)/sum(seqtab)
# 0.9668496
##  Here chimeras make up about 3.33% of the merged sequence variants



# Track reads through the pipeline ----------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(filt_out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
write.table(track, "ITS/Tracking_reads.txt", quote = FALSE, sep = '\t')


# Assign taxonomy ---------------------------------------------------------

# Downloaded train set from https://zenodo.org/record/4587955#.YjMzmHrMKUk
taxa <- assignTaxonomy(seqtab.nochim, "ITS/sh_general_release_dynamic_10.05.2021.fasta", multithread=FALSE)

# Let’s inspect the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Write results to disk
dir.create(table.fp)
saveRDS(seqtab.nochim, paste0(table.fp, "/seqtab_final.rds"))
saveRDS(taxa, paste0(table.fp, "/tax_final.rds"))


## FORMAT OUTPUT
# Flip table
seqtab.t <- as.data.frame(t(seqtab.nochim))

# Pull out ASV repset
rep_set_ASVs <- as.data.frame(rownames(seqtab.t))
rep_set_ASVs <- mutate(rep_set_ASVs, ASV_ID = 1:n())
rep_set_ASVs$ASV_ID <- sub("^", "ASV_", rep_set_ASVs$ASV_ID)
rep_set_ASVs$ASV <- rep_set_ASVs$`rownames(seqtab.t)` 
rep_set_ASVs$`rownames(seqtab.t)` <- NULL

# Add ASV numbers to table
rownames(seqtab.t) <- rep_set_ASVs$ASV_ID

# Add ASV numbers to taxonomy
taxonomy <- as.data.frame(taxa)
taxonomy$ASV <- as.factor(rownames(taxonomy))
taxonomy <- merge(rep_set_ASVs, taxonomy, by = "ASV")
rownames(taxonomy) <- taxonomy$ASV_ID
taxonomy_for_mctoolsr <- unite(taxonomy, "taxonomy", 
                               c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"),
                               sep = ";")

# Write repset to fasta file
# create a function that writes fasta sequences
writeRepSetFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


# Arrange the taxonomy dataframe for the writeRepSetFasta function
taxonomy_for_fasta <- taxonomy %>%
  unite("TaxString", c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"), 
        sep = ";", remove = FALSE) %>%
  unite("name", c("ASV_ID", "TaxString"), 
        sep = " ", remove = TRUE) %>%
  select(ASV, name) %>%
  rename(seq = ASV)

# write fasta file
writeRepSetFasta(taxonomy_for_fasta, paste0(table.fp, "/repset.fasta"))

# Merge taxonomy and table
seqtab_wTax <- merge(seqtab.t, taxonomy_for_mctoolsr, by = 0)
seqtab_wTax$ASV <- NULL 

# Set name of table in mctoolsr format and save
out_fp <- paste0(table.fp, "/seqtab_wTax_mctoolsr.txt")
names(seqtab_wTax)[1] = "#ASV_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(seqtab_wTax, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Also export files as .txt
write.table(seqtab.t, file = paste0(table.fp, "/seqtab_final.txt"),
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(taxa, file = paste0(table.fp, "/tax_final.txt"), 
            sep = "\t", row.names = TRUE, col.names = NA)

# Better output. Without empty columns and repeated ASV name
taxonomy_for_mctoolsr_fix <- unite(taxonomy, "taxonomy", 
                                   c("Kingdom", "Phylum", "Class", "Order","Family", "Genus"),
                                   sep = ";")
taxonomy_for_mctoolsr_fix <- taxonomy_for_mctoolsr_fix[,-4]
# THIS ONE IS DONE! Taxonomy file


# Including taxonomy
# seqtab.t -- there is the asv table
newnames <- paste('E0',seq(1,5), sep='')
newnames <- c(newnames, 'E05-1')
newnames <- c(newnames, paste('E0', seq(6,9), sep=''))
newnames <- c(newnames, paste('E', seq(10,33), sep=''))
newnames <- c(newnames, 'E33-1')
newnames <- c(newnames, paste('E', seq(34,36), sep=''))
newnames <- c(newnames, 'E36-1')
newnames <- c(newnames, paste('E', seq(37,48), sep=''))
newnames <- c(newnames, 'E48-1')

seqtab.t_fix <- seqtab.t
colnames(seqtab.t_fix) <- newnames

# THIS ONE IS DONE! ASV table

# Merged
data_frame_merged <- merge(seqtab.t_fix, taxonomy_for_mctoolsr_fix, 
                           by = 'row.names', all = TRUE) 
data_frame_merged <- data_frame_merged[,-55]
colnames(data_frame_merged[,1])
names(data_frame_merged)[names(data_frame_merged) == 'Row.names'] <- 'ASV_ID'
# THIS ONE IS DONE! All merged


# Exporting the requested output files ------------------------------------

dir.create("ITS/Final_output")

write.table(taxonomy_for_mctoolsr_fix, "ITS/Final_output/taxonomy.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(seqtab.t_fix, "ITS/Final_output/ASV_table.txt", quote = FALSE, row.names = TRUE, sep = "\t")
write.table(data_frame_merged, "ITS/Final_output/ASVTax_merged.txt", quote = FALSE, row.names = FALSE, sep = "\t")

