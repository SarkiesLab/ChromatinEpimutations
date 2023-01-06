
# Manuscript integrated code 


# This is the complete code for analysis and production of figures in: Wilson et al 2022


# 1. Preparation of Alternate Generation RNA data

# This includes the new B_14 RNA file 

# The RNA data has been intersected with gene names on the Imperial high performance cluster and downloaded to a local folder

RNA_list <- list.files(path = "~/Documents/Computational_work /Alternate_Generation_Experiment/Alt_Gen_RNA_bed_files/", pattern = "_Alt_Gen_RNA_")


# Set wd to the same file

setwd("~/Documents/Computational_work /Alternate_Generation_Experiment/Alt_Gen_RNA_bed_files/")



Counts_total <- c()

for(i in 1:length(RNA_list)){
  Tmp <- read.table(RNA_list[i], sep="\t", stringsAsFactors=F)
  counts <- Tmp[,6]
  
  Counts_total <- cbind(Counts_total, counts)
}



# Now read in the C.elegans gene name file 

C_elegans_gene_names2 <- read.table("~/Documents/Computational_work /C_elegans _reference_objects/C.elegans_gene_names2.bed", sep="\t", stringsAsFactors=F)



# The row names for the Counts table are the rows 1-4 of the C_elegans gene name table 

rownames(Counts_total) <- paste(C_elegans_gene_names2[,1], C_elegans_gene_names2[,2], C_elegans_gene_names2[,3], C_elegans_gene_names2[,4], sep = ":")

colnames(Counts_total) <-
 
 c( 
"A_0",    
"A_10",
"A_12",    
"A_14",    
"A_16",    
"A_18",    
"A_2",    
"A_20",    
"A_4",     
"A_6",     
"A_8",     
"B_0",     
"B_10",    
"B_12",    
"B_16",    
"B_18",    
"B_2",     
"B_20",    
"B_6",     
"B_8",     
"C_0",     
"C_10",    
"C_12",    
"C_16",    
"C_18",    
"C_2",     
"C_20",    
"C_4",     
"C_6",     
"C_8",     
"B_14")


#assume you have created a total counts file, with all entries from the bed represented, which is called "Counts_total"

# Remove any 21- URs

M <- table(C_elegans_gene_names2[,4])
M <- M[-grep("21u", names(M))]


# store the longest isoform for each gene

longest<-c()

for(i in 1:length(M)){q<-which(C_elegans_gene_names2[,4]==names(M)[i]); L <- C_elegans_gene_names2[q,3] - C_elegans_gene_names2[q,2]; longest<-c(longest,q[order(L, decreasing=T)][1])}

Counts_total_longest_isoform <- Counts_total[longest,]


# The final step is to select out entries where there were >1 count in at least one column
# This removes lines that do not have substantial reads in at least one sample.

NEW_simplified_RNAseq2 <- Counts_total_longest_isoform[which(apply(Counts_total_longest_isoform,1,max)>1),]





library(DESeq2)


# We normalize all the data points derived from this experiment together within one matrix 

countData <- NEW_simplified_RNAseq2


condition <- factor(c(
  "A_0",    
  "A_10",
  "A_12",    
  "A_14",    
  "A_16",    
  "A_18",    
  "A_2",    
  "A_20",    
  "A_4",     
  "A_6",     
  "A_8",     
  "B_0",     
  "B_10",    
  "B_12",    
  "B_16",    
  "B_18",    
  "B_2",     
  "B_20",    
  "B_6",     
  "B_8",     
  "C_0",     
  "C_10",    
  "C_12",    
  "C_16",    
  "C_18",    
  "C_2",     
  "C_20",    
  "C_4",     
  "C_6",     
  "C_8",     
  "B_14"))

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)


dds <- estimateSizeFactors(dds)

NEW_normalized_table_RNA <- counts(dds, normalized = T)


# Replace any values of 0 with 1 otherwise will get INF results

NEW_normalized_table_RNA[NEW_normalized_table_RNA == 0] <- 1

#-----------------------------------------------------------------------------------------------------------------------------------
# Identify epimutations using a linear model

# Z score cut off is 2.25 according to simulation (see separate folder)



# Find all the residuals, apply z-score cut off of 2.25 and extract the significant residuals

# So now we have a table of RNA reads for A, B and C generations 2 - 20, some higher/lower than the corresponding read for A/B/C 0
# Split this matrix into 3: A, B, C and then find significant residuals from the 0 generation


Lineage_A_RNA_norm <- cbind(NEW_normalized_table_RNA[,1], NEW_normalized_table_RNA[,7], NEW_normalized_table_RNA[,9], NEW_normalized_table_RNA[,10], NEW_normalized_table_RNA[,11], NEW_normalized_table_RNA[,2], NEW_normalized_table_RNA[,3], NEW_normalized_table_RNA[,4], NEW_normalized_table_RNA[,5], NEW_normalized_table_RNA[,6], NEW_normalized_table_RNA[,8])

colnames(Lineage_A_RNA_norm) <- c("A_0", "A_2", "A_4", "A_6", "A_8", "A_10", "A_12", "A_14", "A_16", "A_18", "A_20")



Lineage_B_RNA_norm <- cbind(NEW_normalized_table_RNA[,12], NEW_normalized_table_RNA[,17], NEW_normalized_table_RNA[,19], NEW_normalized_table_RNA[,20], NEW_normalized_table_RNA[,13], NEW_normalized_table_RNA[,14], NEW_normalized_table_RNA[,31], NEW_normalized_table_RNA[,15], NEW_normalized_table_RNA[,16], NEW_normalized_table_RNA[,18])

colnames(Lineage_B_RNA_norm) <- c("B_0", "B_2", "B_6", "B_8", "B_10", "B_12", "B_14", "B_16", "B_18", "B_20")



Lineage_C_RNA_norm <- cbind(NEW_normalized_table_RNA[,21], NEW_normalized_table_RNA[,26], NEW_normalized_table_RNA[,28], NEW_normalized_table_RNA[,29], NEW_normalized_table_RNA[,30], NEW_normalized_table_RNA[,22], NEW_normalized_table_RNA[,23], NEW_normalized_table_RNA[,24], NEW_normalized_table_RNA[,25], NEW_normalized_table_RNA[,27])

colnames(Lineage_C_RNA_norm) <- c("C_0", "C_2", "C_4", "C_6", "C_8", "C_10", "C_12", "C_16", "C_18", "C_20")


#--------------------------------------------------------------------------------------------


# RNA_A

RNA_Z_score_table_A <- matrix(0, ncol = 10, nrow = nrow(Lineage_A_RNA_norm))

normalized_table_A2_to_A20 <- Lineage_A_RNA_norm[, 2:11]

for (i in 1:ncol(normalized_table_A2_to_A20)) {
  LM_temp <-
    lm(log2(normalized_table_A2_to_A20[, i]) ~ log2(Lineage_A_RNA_norm[, 1]))
  Residuals <- residuals(LM_temp)
  RNA_Z_score_table_A[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}


RNA_binarised_table_A_2.25_no_21URs <- matrix(0, ncol = ncol(RNA_Z_score_table_A), nrow= nrow(RNA_Z_score_table_A))

for(i in 1:ncol(RNA_Z_score_table_A)){
  for(j in 1:nrow(RNA_Z_score_table_A)){
    
    if(RNA_Z_score_table_A[j,i] > 2.25){
      RNA_binarised_table_A_2.25_no_21URs[j, i] <- 1
    }
    if(RNA_Z_score_table_A[j,i] < (-2.25)){
      RNA_binarised_table_A_2.25_no_21URs[j, i] <- (-1)
    }
  }
}

rownames(RNA_binarised_table_A_2.25_no_21URs) <- rownames(normalized_table_A2_to_A20)
colnames(RNA_binarised_table_A_2.25_no_21URs) <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)



rownames(RNA_Z_score_table_A) <- rownames(normalized_table_A2_to_A20)
colnames(RNA_Z_score_table_A) <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)



# RNA_B

RNA_Z_score_table_B <- matrix(0, ncol = 9, nrow = nrow(Lineage_B_RNA_norm))

normalized_table_B2_to_B20 <- Lineage_B_RNA_norm[, 2:10]

for (i in 1:ncol(normalized_table_B2_to_B20)) {
  LM_temp <-
    lm(log2(normalized_table_B2_to_B20[, i]) ~ log2(Lineage_B_RNA_norm[, 1]))
  Residuals <- residuals(LM_temp)
  RNA_Z_score_table_B[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}


RNA_binarised_table_B_2.25_no_21URs <- matrix(0, ncol = ncol(RNA_Z_score_table_B), nrow= nrow(RNA_Z_score_table_B))

for(i in 1:ncol(RNA_Z_score_table_B)){
  for(j in 1:nrow(RNA_Z_score_table_B)){
    
    if(RNA_Z_score_table_B[j,i] > 2.25){
      RNA_binarised_table_B_2.25_no_21URs[j, i] <- 1
    }
    if(RNA_Z_score_table_B[j,i] < (-2.25)){
      RNA_binarised_table_B_2.25_no_21URs[j, i] <- (-1)
    }
  }
}

rownames(RNA_binarised_table_B_2.25_no_21URs) <- rownames(normalized_table_B2_to_B20)
colnames(RNA_binarised_table_B_2.25_no_21URs) <- c(2, 6, 8, 10, 12, 14, 16, 18, 20)



rownames(RNA_Z_score_table_B) <- rownames(normalized_table_B2_to_B20)
colnames(RNA_Z_score_table_B) <- c(2, 6, 8, 10, 12, 14, 16, 18, 20)



# RNA_C

RNA_Z_score_table_C <- matrix(0, ncol = 9, nrow = nrow(Lineage_C_RNA_norm))

normalized_table_C2_to_C20 <- Lineage_C_RNA_norm[, 2:10]

for (i in 1:ncol(normalized_table_C2_to_C20)) {
  LM_temp <-
    lm(log2(normalized_table_C2_to_C20[, i]) ~ log2(Lineage_C_RNA_norm[, 1]))
  Residuals <- residuals(LM_temp)
  RNA_Z_score_table_C[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}


RNA_binarised_table_C_2.25_no_21URs <- matrix(0, ncol = ncol(RNA_Z_score_table_C), nrow= nrow(RNA_Z_score_table_C))

for(i in 1:ncol(RNA_Z_score_table_C)){
  for(j in 1:nrow(RNA_Z_score_table_C)){
    
    if(RNA_Z_score_table_C[j,i] > 2.25){
      RNA_binarised_table_C_2.25_no_21URs[j, i] <- 1
    }
    if(RNA_Z_score_table_C[j,i] < (-2.25)){
      RNA_binarised_table_C_2.25_no_21URs[j, i] <- (-1)
    }
  }
}

rownames(RNA_binarised_table_C_2.25_no_21URs) <- rownames(normalized_table_C2_to_C20)
colnames(RNA_binarised_table_C_2.25_no_21URs) <- c(2, 4, 6, 8, 10, 12, 16, 18, 20)


rownames(RNA_Z_score_table_C) <- rownames(normalized_table_C2_to_C20)
colnames(RNA_Z_score_table_C) <- c(2, 4, 6, 8, 10, 12, 16, 18, 20)


#----------------------------------------------------------------------------------------------------------
# ATAC

# The ATAC data is intersected with coordinates from GFF file (regulatory elements) obtained from Julie Ahringer

# These are coordinates of enhancers and promoters

# Analysing the ATAC data using the Ahringer GFF file



# Set the working directory
setwd("~/Desktop/ATAC_enhancer_promoter_bed/")


ATAC_Ahringer_list <- list.files(path = "~/Desktop/ATAC_enhancer_promoter_bed/", pattern = "Alt_Gen")


ATAC_Ahringer_counts <- c()

for(i in 1:length(ATAC_Ahringer_list)){
  Tmp <- read.table(ATAC_Ahringer_list[i], sep="\t", stringsAsFactors=F)
  
  get_Tmp <- Tmp$V14
  
  ATAC_Ahringer_counts <- cbind(ATAC_Ahringer_counts, get_Tmp)
}



colnames(ATAC_Ahringer_counts) <- c("0_A",  "0_B",  "0_C",   "10_A",
                                    "10_B", "10_C",  "12_A", "12_B",
                                    "12_C", "14_A", "14_B", "14_C",
                                    "16_A", "16_B", "16_C", "18_A",
                                    "18_B", "18_C", "2_A",  "2_B", 
                                    "2_C",  "20_A", "20_B", "20_C",
                                    "4_A",  "4_C",  "6_A",  "6_B", 
                                    "6_C",  "8_A", "8_B")


rownames(ATAC_Ahringer_counts) <- paste(Tmp$V1, Tmp$V4, Tmp$V5, sep = ":")




# Remove the rows which correspond to enhancers/promoters that do not map to genes

library(stringr)

Tmp_no_unmatched <- str_remove(Tmp$V10, "Associated-gene-Name=")
test <- cbind(Tmp, Tmp_no_unmatched)

test_get <- test[test$Tmp_no_unmatched %in% ".", ]

remove_rows <- paste(test_get$V1, test_get$V4, test_get$V5, sep = ":")

j <- which(rownames(ATAC_Ahringer_counts) %in% remove_rows)

ATAC_Ahringer_counts <- ATAC_Ahringer_counts[-j, ]


# At this point also produce a reference table where each gene in the Ahringer data set has domain and tissue specificity annotations 

Ahringer_ref_table <-  data.frame(paste(Tmp$V1, Tmp$V4, Tmp$V5, sep=":"), Tmp$V3, str_remove(Tmp$V10, pattern = "Associated-gene-Name="), 
                                  str_remove(Tmp$V12, pattern = "domain="), str_remove(Tmp$V13, pattern = "Tissue-specificity="))


colnames(Ahringer_ref_table) <- c("Locus", "Enhancer_or_Promoter", "Gene_mapped", "Chromatin_domain", "Tissue_specificity")




b <- str_split(Ahringer_ref_table$Gene_mapped, ",")


Ahringer_single_gene_ref_table <- c()

for(y in 1:length(b)){
  save_row <- Ahringer_ref_table[y, c(1, 2, 4, 5)]
  
  row_by_gene <- c()
  for(t in 1:length(b[[y]])){
    
    get <- data.frame(b[[y]][[t]], save_row)
    
    row_by_gene <- rbind(row_by_gene, get)
    
  }
  
  Ahringer_single_gene_ref_table <- rbind(Ahringer_single_gene_ref_table, row_by_gene)
  
  
}


colnames(Ahringer_single_gene_ref_table) <- c("Gene", "ep_locus", "Enhancer_or_Promoter", "Chromatin_domain", "Tissue_specificity")





# Normalisation of all ATAC counts for all lineages is done together

# library(DESeq2)

countData <- as.matrix(ATAC_Ahringer_counts)

condition <- factor(colnames(ATAC_Ahringer_counts))

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)

dds <- estimateSizeFactors(dds)

normalized_table_ATAC_Ahringer_counts <- counts(dds, normalized = T)


# Replace all values of 0 with 1 to prevent NA / INFINITE values

normalized_table_ATAC_Ahringer_counts[normalized_table_ATAC_Ahringer_counts == 0] <- 1



# Produce the binarised tables through using the linear model to identify epimutations according to z score cut off of 2.25

# ATAC_A

normalized_table_Lin_A <- normalized_table_ATAC_Ahringer_counts[, c(1, 19, 25, 27, 30, 4, 7, 10, 13, 16, 22)]

ATAC_Z_score_table_A_enhancer_promoter <- matrix(0, ncol = 10, nrow = nrow(normalized_table_Lin_A))

normalized_table_A2_to_A20 <- normalized_table_Lin_A[, 2:11]

for (i in 1:ncol(normalized_table_A2_to_A20)) {
  LM_temp <-
    lm(log2(normalized_table_A2_to_A20[, i]) ~ log2(normalized_table_Lin_A[, 1]))
  Residuals <- residuals(LM_temp)
  ATAC_Z_score_table_A_enhancer_promoter[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}


ATAC_binarised_table_A_2.25_ehancer_promoter <- matrix(0, ncol = ncol(ATAC_Z_score_table_A_enhancer_promoter), nrow= nrow(ATAC_Z_score_table_A_enhancer_promoter))

for(i in 1:ncol(ATAC_Z_score_table_A_enhancer_promoter)){
  for(j in 1:nrow(ATAC_Z_score_table_A_enhancer_promoter)){
    
    if(ATAC_Z_score_table_A_enhancer_promoter[j,i] > 2.25){
      ATAC_binarised_table_A_2.25_ehancer_promoter[j, i] <- 1
    }
    if(ATAC_Z_score_table_A_enhancer_promoter[j,i] < (-2.25)){
      ATAC_binarised_table_A_2.25_ehancer_promoter[j, i] <- (-1)
    }
  }
}

rownames(ATAC_binarised_table_A_2.25_ehancer_promoter) <- rownames(normalized_table_A2_to_A20)
colnames(ATAC_binarised_table_A_2.25_ehancer_promoter) <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)


rownames(ATAC_Z_score_table_A_enhancer_promoter) <- rownames(normalized_table_A2_to_A20)
colnames(ATAC_Z_score_table_A_enhancer_promoter) <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)







# ATAC_B


normalized_table_Lin_B <- normalized_table_ATAC_Ahringer_counts[, c(2, 20, 28, 31, 5, 8, 11, 14, 17, 23)]

ATAC_Z_score_table_B_enhancer_promoter <- matrix(0, ncol = 9, nrow = nrow(normalized_table_Lin_B))

normalized_table_B2_to_B20 <- normalized_table_Lin_B[, 2:10]

for (i in 1:ncol(normalized_table_B2_to_B20)) {
  LM_temp <-
    lm(log2(normalized_table_B2_to_B20[, i]) ~ log2(normalized_table_Lin_B[, 1]))
  Residuals <- residuals(LM_temp)
  ATAC_Z_score_table_B_enhancer_promoter[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}


ATAC_binarised_table_B_2.25_ehancer_promoter <- matrix(0, ncol = ncol(ATAC_Z_score_table_B_enhancer_promoter), nrow= nrow(ATAC_Z_score_table_B_enhancer_promoter))

for(i in 1:ncol(ATAC_Z_score_table_B_enhancer_promoter)){
  for(j in 1:nrow(ATAC_Z_score_table_B_enhancer_promoter)){
    
    if(ATAC_Z_score_table_B_enhancer_promoter[j,i] > 2.25){
      ATAC_binarised_table_B_2.25_ehancer_promoter[j, i] <- 1
    }
    if(ATAC_Z_score_table_B_enhancer_promoter[j,i] < (-2.25)){
      ATAC_binarised_table_B_2.25_ehancer_promoter[j, i] <- (-1)
    }
  }
}

rownames(ATAC_binarised_table_B_2.25_ehancer_promoter) <- rownames(normalized_table_B2_to_B20)
colnames(ATAC_binarised_table_B_2.25_ehancer_promoter) <- c(2, 6, 8, 10, 12, 14, 16, 18, 20)


rownames(ATAC_Z_score_table_B_enhancer_promoter) <- rownames(normalized_table_B2_to_B20)
colnames(ATAC_Z_score_table_B_enhancer_promoter) <- c(2, 6, 8, 10, 12, 14, 16, 18, 20)




# ATAC_C


normalized_table_Lin_C <- normalized_table_ATAC_Ahringer_counts[, c(3, 21, 26, 29, 6, 9, 12, 15, 18, 24)]

ATAC_Z_score_table_C_enhancer_promoter <- matrix(0, ncol = 9, nrow = nrow(normalized_table_Lin_C))

normalized_table_C2_to_C20 <- normalized_table_Lin_C[, 2:10]

for (i in 1:ncol(normalized_table_C2_to_C20)) {
  LM_temp <-
    lm(log2(normalized_table_C2_to_C20[, i]) ~ log2(normalized_table_Lin_C[, 1]))
  Residuals <- residuals(LM_temp)
  ATAC_Z_score_table_C_enhancer_promoter[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}


ATAC_binarised_table_C_2.25_ehancer_promoter <- matrix(0, ncol = ncol(ATAC_Z_score_table_C_enhancer_promoter), nrow= nrow(ATAC_Z_score_table_C_enhancer_promoter))

for(i in 1:ncol(ATAC_Z_score_table_C_enhancer_promoter)){
  for(j in 1:nrow(ATAC_Z_score_table_C_enhancer_promoter)){
    
    if(ATAC_Z_score_table_C_enhancer_promoter[j,i] > 2.25){
      ATAC_binarised_table_C_2.25_ehancer_promoter[j, i] <- 1
    }
    if(ATAC_Z_score_table_C_enhancer_promoter[j,i] < (-2.25)){
      ATAC_binarised_table_C_2.25_ehancer_promoter[j, i] <- (-1)
    }
  }
}

rownames(ATAC_binarised_table_C_2.25_ehancer_promoter) <- rownames(normalized_table_C2_to_C20)
colnames(ATAC_binarised_table_C_2.25_ehancer_promoter) <- c(2, 4, 6, 10, 12, 14, 16, 18, 20)


rownames(ATAC_Z_score_table_C_enhancer_promoter) <- rownames(normalized_table_C2_to_C20)
colnames(ATAC_Z_score_table_C_enhancer_promoter) <- c(2, 4, 6, 10, 12, 14, 16, 18, 20)


# ------------------------------------------------------------------------------------------

# Small RNA

# The small RNA data has been pre-processed (see separate folder), these are antisense 22Gs, normalized by DESeq and cutting off to remove low expression (important otherwise there is a lot of noise in the data).  

# This is the input for the linear model (no need to normalize the data as this has already been done)

# Small RNA counts are grouped according to the names of genes with sequences that are antisense to the 22G small RNAs

# Load the data 

AltGen_smallRNA <- load("~/Documents/Computational_work /Alternate_Generation_Experiment/Alt_Gen_small_RNA/Alt_Gen_smallRNAs.Rdata")

Alt_Gen_smallRNA <- eval(parse(text=AltGen_smallRNA))


# Separate by lineage 


# small RNA A

normalized_Lin_A_small_RNA <- data.frame()

normalized_Lin_A_small_RNA <- cbind(Alt_Gen_smallRNA[, 1], Alt_Gen_smallRNA[, 18], Alt_Gen_smallRNA[, 24], Alt_Gen_smallRNA[, 26], Alt_Gen_smallRNA[, 29], Alt_Gen_smallRNA[, 4], 
                                    Alt_Gen_smallRNA[, 7], Alt_Gen_smallRNA[, 10], Alt_Gen_smallRNA[, 12], Alt_Gen_smallRNA[, 15], Alt_Gen_smallRNA[, 21])


colnames(normalized_Lin_A_small_RNA) <-
  c("0_A", "2_A", "4_A", "6_A", "8_A", "10_A", "12_A", "14_A", "16_A", "18_A", "20_A")


# small RNA B

normalized_Lin_B_small_RNA <- data.frame()

normalized_Lin_B_small_RNA <- cbind(Alt_Gen_smallRNA[, 2], Alt_Gen_smallRNA[, 19],  Alt_Gen_smallRNA[, 27], Alt_Gen_smallRNA[, 30], Alt_Gen_smallRNA[, 5], 
                                    Alt_Gen_smallRNA[, 8], Alt_Gen_smallRNA[, 11], Alt_Gen_smallRNA[, 13], Alt_Gen_smallRNA[, 16], Alt_Gen_smallRNA[, 22])


colnames(normalized_Lin_B_small_RNA) <-
  c("0_B", "2_B", "6_B", "8_B", "10_B", "12_B", "14_B", "16_B", "18_B", "20_B")



# small RNA C

normalized_Lin_C_small_RNA <- data.frame()

normalized_Lin_C_small_RNA <- cbind(Alt_Gen_smallRNA[, 3], Alt_Gen_smallRNA[, 20],  Alt_Gen_smallRNA[, 25], Alt_Gen_smallRNA[, 28], Alt_Gen_smallRNA[, 31], Alt_Gen_smallRNA[, 6],
                                    Alt_Gen_smallRNA[, 9], Alt_Gen_smallRNA[, 14], Alt_Gen_smallRNA[, 17], Alt_Gen_smallRNA[, 23])


colnames(normalized_Lin_C_small_RNA) <-
  c("0_C", "2_C", "4_C", "6_C", "8_C", "10_C", "12_C", "16_C", "18_C", "20_C")



# -----------------------------------------------------------------------------------


# Now find the Z scores and epimutations using cut off of 2.25


# small RNA A

# Replace all values of 0 with 1

normalized_Lin_A_small_RNA[normalized_Lin_A_small_RNA== 0] <- 1


normalized_Lin_A_small_RNA_BN <- normalized_Lin_A_small_RNA[, c(2:11)]


smallRNA_Z_score_table_A <- matrix(0, ncol = 10, nrow = nrow(normalized_Lin_A_small_RNA))


for (i in 1:ncol(normalized_Lin_A_small_RNA_BN)) {
  LM_temp <-
    lm(log2(normalized_Lin_A_small_RNA_BN[, i]) ~ log2(normalized_Lin_A_small_RNA[, 1]))
  Residuals <- residuals(LM_temp)
  smallRNA_Z_score_table_A[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}



smallRNA_binarised_table_A_2.25 <- matrix(0, ncol = ncol(smallRNA_Z_score_table_A), nrow= nrow(smallRNA_Z_score_table_A))

for(i in 1:ncol(smallRNA_Z_score_table_A)){
  for(j in 1:nrow(smallRNA_Z_score_table_A)){
    
    if(smallRNA_Z_score_table_A[j,i] > 2.25){
      smallRNA_binarised_table_A_2.25[j, i] <- 1
    }
    if(smallRNA_Z_score_table_A[j,i] < (-2.25)){
      smallRNA_binarised_table_A_2.25[j, i] <- (-1)
    }
  }
}

rownames(smallRNA_binarised_table_A_2.25 ) <- rownames(normalized_Lin_A_small_RNA)
colnames(smallRNA_binarised_table_A_2.25 ) <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)


rownames(smallRNA_Z_score_table_A) <- rownames(normalized_Lin_A_small_RNA)
colnames(smallRNA_Z_score_table_A) <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)


# ---------------------------------------

# small RNA B

# Replace all values of 0 with 1

normalized_Lin_B_small_RNA[normalized_Lin_B_small_RNA== 0] <- 1


normalized_Lin_B_small_RNA_BN <- normalized_Lin_B_small_RNA[, c(2:10)]


smallRNA_Z_score_table_B <- matrix(0, ncol = 9, nrow = nrow(normalized_Lin_B_small_RNA))


for (i in 1:ncol(normalized_Lin_B_small_RNA_BN)) {
  LM_temp <-
    lm(log2(normalized_Lin_B_small_RNA_BN[, i]) ~ log2(normalized_Lin_B_small_RNA[, 1]))
  Residuals <- residuals(LM_temp)
  smallRNA_Z_score_table_B[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}



smallRNA_binarised_table_B_2.25 <- matrix(0, ncol = ncol(smallRNA_Z_score_table_B), nrow= nrow(smallRNA_Z_score_table_B))

for(i in 1:ncol(smallRNA_Z_score_table_B)){
  for(j in 1:nrow(smallRNA_Z_score_table_B)){
    
    if(smallRNA_Z_score_table_B[j,i] > 2.25){
      smallRNA_binarised_table_B_2.25[j, i] <- 1
    }
    if(smallRNA_Z_score_table_B[j,i] < (-2.25)){
      smallRNA_binarised_table_B_2.25[j, i] <- (-1)
    }
  }
}

rownames(smallRNA_binarised_table_B_2.25 ) <- rownames(normalized_Lin_B_small_RNA)
colnames(smallRNA_binarised_table_B_2.25 ) <- c(2, 6, 8, 10, 12, 14, 16, 18, 20)


rownames(smallRNA_Z_score_table_B) <- rownames(normalized_Lin_B_small_RNA)
colnames(smallRNA_Z_score_table_B) <- c(2, 6, 8, 10, 12, 14, 16, 18, 20)


# ----------------------------------------------------------------------------------


# small RNA C

# Replace all values of 0 with 1

normalized_Lin_C_small_RNA[normalized_Lin_C_small_RNA== 0] <- 1


normalized_Lin_C_small_RNA_BN <- normalized_Lin_C_small_RNA[, c(2:10)]


smallRNA_Z_score_table_C <- matrix(0, ncol = 9, nrow = nrow(normalized_Lin_C_small_RNA))


for (i in 1:ncol(normalized_Lin_C_small_RNA_BN)) {
  LM_temp <-
    lm(log2(normalized_Lin_C_small_RNA_BN[, i]) ~ log2(normalized_Lin_C_small_RNA[, 1]))
  Residuals <- residuals(LM_temp)
  smallRNA_Z_score_table_C[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}



smallRNA_binarised_table_C_2.25 <- matrix(0, ncol = ncol(smallRNA_Z_score_table_C), nrow= nrow(smallRNA_Z_score_table_C))

for(i in 1:ncol(smallRNA_Z_score_table_C)){
  for(j in 1:nrow(smallRNA_Z_score_table_C)){
    
    if(smallRNA_Z_score_table_C[j,i] > 2.25){
      smallRNA_binarised_table_C_2.25[j, i] <- 1
    }
    if(smallRNA_Z_score_table_C[j,i] < (-2.25)){
      smallRNA_binarised_table_C_2.25[j, i] <- (-1)
    }
  }
}

rownames(smallRNA_binarised_table_C_2.25 ) <- rownames(normalized_Lin_C_small_RNA)
colnames(smallRNA_binarised_table_C_2.25 ) <- c(2, 4, 6, 8, 10, 12, 16, 18, 20)


rownames(smallRNA_Z_score_table_C) <- rownames(normalized_Lin_C_small_RNA)
colnames(smallRNA_Z_score_table_C) <- c(2, 4, 6, 8, 10, 12, 16, 18, 20)


# ---------------------------------------------------------------------------------------------------------------------
# Figure 1

# Figure 1 A Mechanisms of transgenerational epigenetic inheritance - image made in BioRender

# Figure 1 B Schematic of experimental design - image made in BioRender



# Figure 1 C Example plot to show how the epimutations are identified

# This uses a linear model with lm

LM_temp <- lm(log2(normalized_table_C2_to_C20[, 2]) ~ log2(normalized_table_C2_to_C20[, 1]))

# This allows us to calculate the difference of each individual point from the overall regression line:

residuals_line<-residuals(LM_temp)

# so we can just extract directly the residuals from the linear model

# We extract residuals according to a Z score 

zscore<-(residuals_line-mean(residuals_line))/sd(residuals_line)

Zscore<-abs(zscore)

# The abs function makes them all positive as we just want their size, not if they are up or down at this stage.

# The z-score tells you how many standard deviations a point is away from the mean.  So to select the biggest deviations we can use a cutoff of, for example, 2.25, which means 2.25 standard deviations away from the mean.  

# We could use more or less to be more or less stringent- but 2.25 has been validated through simulation, see Fig.1.D

Sig_points <-which(Zscore> 2.25)

plot(log2(normalized_table_C2_to_C20[,1]), log2(normalized_table_C2_to_C20[,2]), col = "grey", 
     xlab="log2(Baseline (PMA) ATACseq counts)", ylab="log2(Bottlenecked generation ATACseq counts)", 
     main = "Figure 1 C")


points(log2(normalized_table_C2_to_C20[Sig_points,1]),log2(normalized_table_C2_to_C20[Sig_points,2]), col="darkblue")


# -------------

# Figure 1 D see separate folder

#--------------


# Defining the number of transitions UP/DOWN per generation for each data type

# This code counts the number of transitions in chromatin/gene expression state in each generation

# UP transition function

#define function
UP_transition_func<-function(vector_in,input_name){
  
  Gen_ON <- vector()
  
  for(i in 2:length(vector_in)){
    
    if(vector_in[i]==1&vector_in[i-1]== 0){
      
      Gen_ON[[i]] <-1 # there is an UP transition to ON at generation[i]
      
    }
    
    if(vector_in[i]== 1&vector_in[i-1]== -1){
      
      Gen_ON[[i]] <-1 # there is an UP transition to ON at generation[i]
      
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== 0){
      
      Gen_ON[[i]] <- 0 # this is a DOWN transition
      
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== 1){
      
      Gen_ON[[i]] <- 0 # this is a DOWN epimutation
      
    }
    
    
    if(vector_in[i]==1&vector_in[i-1]==1){
      
      Gen_ON[[i]] <- 0 # inside a run of 1 1 epimutations at generation[i]
      
    }
    
    if(vector_in[i]==-1&vector_in[i-1]==-1){
      
      Gen_ON[[i]] <- 0 # inside a run of -1 -1 epimutations at generation[i]
      
    }
    
    if(vector_in[i]==0&vector_in[i-1]==0){
      
      Gen_ON[[i]] <- 0 # between epimutations at generation[i]
      
    }
    
    if(vector_in[i]==0&vector_in[i-1]==1){
      
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==-1){
      
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    output <- Gen_ON
    
  }
  
  return(output)
}






# DOWN transition function


#define function
DOWN_transition_func<-function(vector_in,input_name){
  
  Gen_ON <- vector()
  
  for(i in 2:length(vector_in)){
    
    if(vector_in[i]== -1&vector_in[i-1]== 0){
      
      Gen_ON[[i]] <-1 # there is a DOWN transition to ON at generation[i]
      
    }
    
    if(vector_in[i]== -1&vector_in[i-1]== 1){
      
      Gen_ON[[i]] <-1 # there is a DOWN transition to ON at generation[i]
      
    }
    
    if(vector_in[i]== 1&vector_in[i-1]== 0){
      
      Gen_ON[[i]] <- 0 # this is an UP transition
      
    }
    
    if(vector_in[i]== 1&vector_in[i-1]== -1){
      
      Gen_ON[[i]] <- 0 # this is an UP epimutation
      
    }
    
    
    if(vector_in[i]==1&vector_in[i-1]==1){
      
      Gen_ON[[i]] <- 0 # inside a run of 1 1 epimutations at generation[i]
      
    }
    
    if(vector_in[i]==-1&vector_in[i-1]==-1){
      
      Gen_ON[[i]] <- 0 # inside a run of -1 -1 epimutations at generation[i]
      
    }
    
    if(vector_in[i]==0&vector_in[i-1]==0){
      
      Gen_ON[[i]] <- 0 # between epimutations at generation[i]
      
    }
    
    if(vector_in[i]==0&vector_in[i-1]==1){
      
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    if(vector_in[i]==0&vector_in[i-1]==-1){
      
      Gen_ON[[i]] <- 0 # there is a transition to OFF at generation[i]
    }
    
    output <- Gen_ON
    
  }
  
  return(output)
}

#--------------------------------------------------------------------------------------------
# RNA
# UP


# RNA_A 
RNA_Tab_A<-cbind(rep(0, length=nrow(RNA_binarised_table_A_2.25_no_21URs)),RNA_binarised_table_A_2.25_no_21URs)
colnames(RNA_Tab_A)[1]<-"0"
#shows that the "start" has no expression changes. 


UP_output_RNA_A<- c()

for(i in 1:nrow(RNA_Tab_A)){
  
  UP_output_RNA_A <-rbind(UP_output_RNA_A, UP_transition_func(RNA_Tab_A[i,], input_name=row.names(RNA_Tab_A)[i]))}


colnames(UP_output_RNA_A) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(UP_output_RNA_A) <- rownames(RNA_binarised_table_A_2.25_no_21URs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(UP_output_RNA_A[,2])
transitions_at_4 <- sum(UP_output_RNA_A[,3])
transitions_at_6 <- sum(UP_output_RNA_A[,4])
transitions_at_8 <- sum(UP_output_RNA_A[,5])
transitions_at_10 <- sum(UP_output_RNA_A[,6])
transitions_at_12 <- sum(UP_output_RNA_A[,7])
transitions_at_14 <- sum(UP_output_RNA_A[,8])
transitions_at_16 <- sum(UP_output_RNA_A[,9])
transitions_at_18 <- sum(UP_output_RNA_A[,10])
transitions_at_20 <- sum(UP_output_RNA_A[,11])

UP_RNA_A_Table_of_new_transitions <- rbind(transitions_at_2,
                                           transitions_at_4, 
                                           transitions_at_6,
                                           transitions_at_8, 
                                           transitions_at_10, 
                                           transitions_at_12, 
                                           transitions_at_14, 
                                           transitions_at_16, 
                                           transitions_at_18, 
                                           transitions_at_20)

colnames(UP_RNA_A_Table_of_new_transitions) <- c("A")





# RNA_B

# Put a column of zeros to the left of  RNA_B_expression_status

RNA_Tab_B<-cbind(rep(0, length=nrow(RNA_binarised_table_B_2.25_no_21URs)),RNA_binarised_table_B_2.25_no_21URs)
colnames(RNA_Tab_B)[1]<-"0"


#run function over entire tab

UP_output_RNA_B<- c()

for(i in 1:nrow(RNA_Tab_B)){
  
  UP_output_RNA_B <-rbind(UP_output_RNA_B, UP_transition_func(RNA_Tab_B[i,], input_name=row.names(RNA_Tab_B)[i]))}


colnames(UP_output_RNA_B) <- c("0", "2", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(UP_output_RNA_B) <- rownames(RNA_binarised_table_B_2.25_no_21URs)

# Now take the number of NEW events that occur at each generation, 

transitions_at_2 <- sum(UP_output_RNA_B[,2])
transitions_at_6 <- sum(UP_output_RNA_B[,3])
transitions_at_8 <- sum(UP_output_RNA_B[,4])
transitions_at_10 <- sum(UP_output_RNA_B[,5])
transitions_at_12 <- sum(UP_output_RNA_B[,6])
transitions_at_14 <- sum(UP_output_RNA_B[,7])
transitions_at_16 <- sum(UP_output_RNA_B[,8])
transitions_at_18 <- sum(UP_output_RNA_B[,9])
transitions_at_20 <- sum(UP_output_RNA_B[,10])

UP_RNA_B_Table_of_new_transitions <- rbind(transitions_at_2,
                                           transitions_at_6,
                                           transitions_at_8, 
                                           transitions_at_10, 
                                           transitions_at_12, 
                                           transitions_at_14, 
                                           transitions_at_16, 
                                           transitions_at_18, 
                                           transitions_at_20)

colnames(UP_RNA_B_Table_of_new_transitions) <- c("B")




# RNA_C 

# Put a column of zeros to the left of RNA_C_expression_status

RNA_Tab_C<-cbind(rep(0, length=nrow(RNA_binarised_table_C_2.25_no_21URs)),RNA_binarised_table_C_2.25_no_21URs)
colnames(RNA_Tab_C)[1]<-"0"


#run function over entire tab

UP_output_RNA_C<- c()

for(i in 1:nrow(RNA_Tab_C)){
  
  UP_output_RNA_C <-rbind(UP_output_RNA_C, UP_transition_func(RNA_Tab_C[i,], input_name=row.names(RNA_Tab_C)[i]))}


colnames(UP_output_RNA_C) <- c("0", "2", "4", "6", "8", "10", "12", "16", "18", "20")

row.names(UP_output_RNA_C) <- rownames(RNA_binarised_table_C_2.25_no_21URs)



# Now take the number of NEW events that occur at each generation, 

transitions_at_2 <- sum(UP_output_RNA_C[,2])
transitions_at_4 <- sum(UP_output_RNA_C[,3])
transitions_at_6 <- sum(UP_output_RNA_C[,4])
transitions_at_8 <- sum(UP_output_RNA_C[,5])
transitions_at_10 <- sum(UP_output_RNA_C[,6])
transitions_at_12 <- sum(UP_output_RNA_C[,7])
transitions_at_16 <- sum(UP_output_RNA_C[,8])
transitions_at_18 <- sum(UP_output_RNA_C[,9])
transitions_at_20 <- sum(UP_output_RNA_C[,10])

UP_RNA_C_Table_of_new_transitions <- rbind(transitions_at_2,
                                           transitions_at_4,
                                           transitions_at_6, 
                                           transitions_at_8, 
                                           transitions_at_10, 
                                           transitions_at_12, 
                                           transitions_at_16, 
                                           transitions_at_18, 
                                           transitions_at_20)

colnames(UP_RNA_C_Table_of_new_transitions) <- c("C")





# DOWN
# RNA_A

DOWN_output_RNA_A<- c()

for(i in 1:nrow(RNA_Tab_A)){
  
  DOWN_output_RNA_A <-rbind(DOWN_output_RNA_A, DOWN_transition_func(RNA_Tab_A[i,], input_name=row.names(RNA_Tab_A)[i]))}


colnames(DOWN_output_RNA_A) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_A) <- rownames(RNA_binarised_table_A_2.25_no_21URs)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(DOWN_output_RNA_A[,2])
transitions_at_4 <- sum(DOWN_output_RNA_A[,3])
transitions_at_6 <- sum(DOWN_output_RNA_A[,4])
transitions_at_8 <- sum(DOWN_output_RNA_A[,5])
transitions_at_10 <- sum(DOWN_output_RNA_A[,6])
transitions_at_12 <- sum(DOWN_output_RNA_A[,7])
transitions_at_14 <- sum(DOWN_output_RNA_A[,8])
transitions_at_16 <- sum(DOWN_output_RNA_A[,9])
transitions_at_18 <- sum(DOWN_output_RNA_A[,10])
transitions_at_20 <- sum(DOWN_output_RNA_A[,11])

DOWN_RNA_A_Table_of_new_transitions <- rbind(transitions_at_2,
                                             transitions_at_4, 
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

colnames(DOWN_RNA_A_Table_of_new_transitions) <- c("A")





# RNA_B

# run function over entire tab

DOWN_output_RNA_B<- c()

for(i in 1:nrow(RNA_Tab_B)){
  
  DOWN_output_RNA_B <-rbind(DOWN_output_RNA_B, DOWN_transition_func(RNA_Tab_B[i,], input_name=row.names(RNA_Tab_B)[i]))}


colnames(DOWN_output_RNA_B) <- c("0", "2", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_RNA_B) <- rownames(RNA_binarised_table_B_2.25_no_21URs)

# Now take the number of NEW events that occur at each generation, 

transitions_at_2 <- sum(DOWN_output_RNA_B[,2])
transitions_at_6 <- sum(DOWN_output_RNA_B[,3])
transitions_at_8 <- sum(DOWN_output_RNA_B[,4])
transitions_at_10 <- sum(DOWN_output_RNA_B[,5])
transitions_at_12 <- sum(DOWN_output_RNA_B[,6])
transitions_at_14 <- sum(DOWN_output_RNA_B[,7])
transitions_at_16 <- sum(DOWN_output_RNA_B[,8])
transitions_at_18 <- sum(DOWN_output_RNA_B[,9])
transitions_at_20 <- sum(DOWN_output_RNA_B[,10])

DOWN_RNA_B_Table_of_new_transitions <- rbind(transitions_at_2,
                                             transitions_at_6,
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_14, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

colnames(DOWN_RNA_B_Table_of_new_transitions) <- c("B")




# RNA_C

# run function over entire tab

DOWN_output_RNA_C<- c()

for(i in 1:nrow(RNA_Tab_C)){
  
  DOWN_output_RNA_C <-rbind(DOWN_output_RNA_C, DOWN_transition_func(RNA_Tab_C[i,], input_name=row.names(RNA_Tab_C)[i]))}


colnames(DOWN_output_RNA_C) <- c("0", "2", "4", "6", "8", "10", "12", "16", "18", "20")

row.names(DOWN_output_RNA_C) <- rownames(RNA_binarised_table_C_2.25_no_21URs)



# Now take the number of NEW events that occur at each generation, 

transitions_at_2 <- sum(DOWN_output_RNA_C[,2])
transitions_at_4 <- sum(DOWN_output_RNA_C[,3])
transitions_at_6 <- sum(DOWN_output_RNA_C[,4])
transitions_at_8 <- sum(DOWN_output_RNA_C[,5])
transitions_at_10 <- sum(DOWN_output_RNA_C[,6])
transitions_at_12 <- sum(DOWN_output_RNA_C[,7])
transitions_at_16 <- sum(DOWN_output_RNA_C[,8])
transitions_at_18 <- sum(DOWN_output_RNA_C[,9])
transitions_at_20 <- sum(DOWN_output_RNA_C[,10])

DOWN_RNA_C_Table_of_new_transitions <- rbind(transitions_at_2,
                                             transitions_at_4,
                                             transitions_at_6, 
                                             transitions_at_8, 
                                             transitions_at_10, 
                                             transitions_at_12, 
                                             transitions_at_16, 
                                             transitions_at_18, 
                                             transitions_at_20)

colnames(DOWN_RNA_C_Table_of_new_transitions) <- c("C")


#-----------
# ATAC
# UP

# ATAC_A

ATAC_Tab_A<-cbind(rep(0, length=nrow(ATAC_binarised_table_A_2.25_ehancer_promoter)),ATAC_binarised_table_A_2.25_ehancer_promoter)
colnames(ATAC_Tab_A)[1]<-0
#shows that the "start" has no epimutations. 



UP_output_Lin_A<- c()

for(i in 1:nrow(ATAC_Tab_A)){
  
  UP_output_Lin_A <-rbind(UP_output_Lin_A, UP_transition_func(ATAC_Tab_A[i,], input_name=row.names(ATAC_Tab_A)[i]))}

colnames(UP_output_Lin_A) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(UP_output_Lin_A) <- rownames(ATAC_binarised_table_A_2.25_ehancer_promoter)


# Now take the sum of NEW UP events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(UP_output_Lin_A[,2])
transitions_at_4 <- sum(UP_output_Lin_A[,3])
transitions_at_6 <- sum(UP_output_Lin_A[,4])
transitions_at_8 <- sum(UP_output_Lin_A[,5])
transitions_at_10 <- sum(UP_output_Lin_A[,6])
transitions_at_12 <- sum(UP_output_Lin_A[,7])
transitions_at_14 <- sum(UP_output_Lin_A[,8])
transitions_at_16 <- sum(UP_output_Lin_A[,9])
transitions_at_18 <- sum(UP_output_Lin_A[,10])
transitions_at_20 <- sum(UP_output_Lin_A[,11])

UP_ATAC_A_Table_of_new_transitions <- rbind(transitions_at_2,
                                            transitions_at_4, 
                                            transitions_at_6,
                                            transitions_at_8, 
                                            transitions_at_10, 
                                            transitions_at_12, 
                                            transitions_at_14, 
                                            transitions_at_16, 
                                            transitions_at_18, 
                                            transitions_at_20)

colnames(UP_ATAC_A_Table_of_new_transitions) <- c("A")


# ATAC_B

# Put a column of zeros to the left of Tabulated_Lineage_B

ATAC_Tab_B<-cbind(rep(0, length=nrow(ATAC_binarised_table_B_2.25_ehancer_promoter)),ATAC_binarised_table_B_2.25_ehancer_promoter)
colnames(ATAC_Tab_B)[1]<-"0"


#run function over entire tab

UP_output_Lin_B<- c()

for(i in 1:nrow(ATAC_Tab_B)){
  
  UP_output_Lin_B <-rbind(UP_output_Lin_B, UP_transition_func(ATAC_Tab_B[i,], input_name=row.names(ATAC_Tab_B)[i]))}


colnames(UP_output_Lin_B) <- c("0", "2", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(UP_output_Lin_B) <- rownames(ATAC_binarised_table_B_2.25_ehancer_promoter)

# Now take the number of NEW UP events that occur at each generation, 

transitions_at_2 <- sum(UP_output_Lin_B[,2])
transitions_at_6 <- sum(UP_output_Lin_B[,3])
transitions_at_8 <- sum(UP_output_Lin_B[,4])
transitions_at_10 <- sum(UP_output_Lin_B[,5])
transitions_at_12 <- sum(UP_output_Lin_B[,6])
transitions_at_14 <- sum(UP_output_Lin_B[,7])
transitions_at_16 <- sum(UP_output_Lin_B[,8])
transitions_at_18 <- sum(UP_output_Lin_B[,9])
transitions_at_20 <- sum(UP_output_Lin_B[,10])

UP_ATAC_B_Table_of_new_transitions <- rbind(transitions_at_2,
                                            transitions_at_6,
                                            transitions_at_8, 
                                            transitions_at_10, 
                                            transitions_at_12, 
                                            transitions_at_14, 
                                            transitions_at_16, 
                                            transitions_at_18, 
                                            transitions_at_20)

colnames(UP_ATAC_B_Table_of_new_transitions) <- c("B")




# ATAC_C

# Put a column of zeros to the left of Tabulated_Lineage_C

ATAC_Tab_C<-cbind(rep(0, length=nrow(ATAC_binarised_table_C_2.25_ehancer_promoter)),ATAC_binarised_table_C_2.25_ehancer_promoter)
colnames(ATAC_Tab_C)[1]<-"0"


#run function over entire tab

UP_output_Lin_C<- c()

for(i in 1:nrow(ATAC_Tab_C)){
  
  UP_output_Lin_C <-rbind(UP_output_Lin_C, UP_transition_func(ATAC_Tab_C[i,], input_name=row.names(ATAC_Tab_C)[i]))}


colnames(UP_output_Lin_C) <- c("0", "2", "4", "6", "10", "12", "14", "16", "18", "20")

row.names(UP_output_Lin_C) <- rownames(ATAC_binarised_table_C_2.25_ehancer_promoter)



# Now take the number of NEW UP events that occur at each generation, 

transitions_at_2 <- sum(UP_output_Lin_C[,2])
transitions_at_4 <- sum(UP_output_Lin_C[,3])
transitions_at_6 <- sum(UP_output_Lin_C[,4])
transitions_at_10 <- sum(UP_output_Lin_C[,5])
transitions_at_12 <- sum(UP_output_Lin_C[,6])
transitions_at_14 <- sum(UP_output_Lin_C[,7])
transitions_at_16 <- sum(UP_output_Lin_C[,8])
transitions_at_18 <- sum(UP_output_Lin_C[,9])
transitions_at_20 <- sum(UP_output_Lin_C[,10])

UP_ATAC_C_Table_of_new_transitions <- rbind(transitions_at_2,
                                            transitions_at_4,
                                            transitions_at_6,
                                            transitions_at_10, 
                                            transitions_at_12, 
                                            transitions_at_14, 
                                            transitions_at_16, 
                                            transitions_at_18, 
                                            transitions_at_20)

colnames(UP_ATAC_C_Table_of_new_transitions) <- c("C")




# DOWN
# ATAC_A

DOWN_output_Lin_A<- c()

for(i in 1:nrow(ATAC_Tab_A)){
  
  DOWN_output_Lin_A <-rbind(DOWN_output_Lin_A, DOWN_transition_func(ATAC_Tab_A[i,], input_name=row.names(ATAC_Tab_A)[i]))}


colnames(DOWN_output_Lin_A) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_Lin_A) <- rownames(ATAC_binarised_table_A_2.25_ehancer_promoter)

# Now take the sum of NEW DOWN events (i.e. (0 to -1), (1 to -1) transitions) that occur at each generation, 

transitions_at_2 <- sum(DOWN_output_Lin_A[,2])
transitions_at_4 <- sum(DOWN_output_Lin_A[,3])
transitions_at_6 <- sum(DOWN_output_Lin_A[,4])
transitions_at_8 <- sum(DOWN_output_Lin_A[,5])
transitions_at_10 <- sum(DOWN_output_Lin_A[,6])
transitions_at_12 <- sum(DOWN_output_Lin_A[,7])
transitions_at_14 <- sum(DOWN_output_Lin_A[,8])
transitions_at_16 <- sum(DOWN_output_Lin_A[,9])
transitions_at_18 <- sum(DOWN_output_Lin_A[,10])
transitions_at_20 <- sum(DOWN_output_Lin_A[,11])

DOWN_ATAC_A_Table_of_new_transitions <- rbind(transitions_at_2,
                                              transitions_at_4, 
                                              transitions_at_6,
                                              transitions_at_8, 
                                              transitions_at_10, 
                                              transitions_at_12, 
                                              transitions_at_14, 
                                              transitions_at_16, 
                                              transitions_at_18, 
                                              transitions_at_20)

colnames(DOWN_ATAC_A_Table_of_new_transitions) <- c("A")





# ATAC_B


#run function over entire tab

DOWN_output_Lin_B<- c()

for(i in 1:nrow(ATAC_Tab_B)){
  
  DOWN_output_Lin_B <-rbind(DOWN_output_Lin_B, DOWN_transition_func(ATAC_Tab_B[i,], input_name=row.names(ATAC_Tab_B)[i]))}


colnames(DOWN_output_Lin_B) <- c("0", "2", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_Lin_B) <- rownames(ATAC_binarised_table_B_2.25_ehancer_promoter)

# Now take the sum of NEW DOWN events that occur at each generation, 

transitions_at_2 <- sum(DOWN_output_Lin_B[,2])
transitions_at_6 <- sum(DOWN_output_Lin_B[,3])
transitions_at_8 <- sum(DOWN_output_Lin_B[,4])
transitions_at_10 <- sum(DOWN_output_Lin_B[,5])
transitions_at_12 <- sum(DOWN_output_Lin_B[,6])
transitions_at_14 <- sum(DOWN_output_Lin_B[,7])
transitions_at_16 <- sum(DOWN_output_Lin_B[,8])
transitions_at_18 <- sum(DOWN_output_Lin_B[,9])
transitions_at_20 <- sum(DOWN_output_Lin_B[,10])

DOWN_ATAC_B_Table_of_new_transitions <- rbind(transitions_at_2,
                                              transitions_at_6,
                                              transitions_at_8, 
                                              transitions_at_10, 
                                              transitions_at_12, 
                                              transitions_at_14, 
                                              transitions_at_16, 
                                              transitions_at_18, 
                                              transitions_at_20)

colnames(DOWN_ATAC_B_Table_of_new_transitions) <- c("B")




# ATAC_C

#run function over entire tab

DOWN_output_Lin_C<- c()

for(i in 1:nrow(ATAC_Tab_C)){
  
  DOWN_output_Lin_C <-rbind(DOWN_output_Lin_C, DOWN_transition_func(ATAC_Tab_C[i,], input_name=row.names(ATAC_Tab_C)[i]))}


colnames(DOWN_output_Lin_C) <- c("0", "2", "4", "6", "10", "12", "14", "16", "18", "20")

row.names(DOWN_output_Lin_C) <- rownames(ATAC_binarised_table_C_2.25_ehancer_promoter)



# Now take the number of NEW DOWN events that occur at each generation, 

transitions_at_2 <- sum(DOWN_output_Lin_C[,2])
transitions_at_4 <- sum(DOWN_output_Lin_C[,3])
transitions_at_6 <- sum(DOWN_output_Lin_C[,4])
transitions_at_10 <- sum(DOWN_output_Lin_C[,5])
transitions_at_12 <- sum(DOWN_output_Lin_C[,6])
transitions_at_14 <- sum(DOWN_output_Lin_C[,7])
transitions_at_16 <- sum(DOWN_output_Lin_C[,8])
transitions_at_18 <- sum(DOWN_output_Lin_C[,9])
transitions_at_20 <- sum(DOWN_output_Lin_C[,10])

DOWN_ATAC_C_Table_of_new_transitions <- rbind(transitions_at_2,
                                              transitions_at_4,
                                              transitions_at_6, 
                                              transitions_at_10, 
                                              transitions_at_12, 
                                              transitions_at_14, 
                                              transitions_at_16, 
                                              transitions_at_18, 
                                              transitions_at_20)

colnames(DOWN_ATAC_C_Table_of_new_transitions) <- c("C")

# --------------

# small RNA 

# UP

# small_RNA_A

smallRNA_Tab_A<-cbind(rep(0, length=nrow(smallRNA_binarised_table_A_2.25)),smallRNA_binarised_table_A_2.25)
colnames(smallRNA_Tab_A)[1]<-0
#shows that the "start" has no epimutations. 



UP_output_Lin_A<- c()

for(i in 1:nrow(smallRNA_Tab_A)){
  
  UP_output_Lin_A <-rbind(UP_output_Lin_A, UP_transition_func(smallRNA_Tab_A[i,], input_name=row.names(smallRNA_Tab_A)[i]))}

colnames(UP_output_Lin_A) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(UP_output_Lin_A) <- rownames(smallRNA_binarised_table_A_2.25)

# Now take the sum of NEW UP events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(UP_output_Lin_A[,2])
transitions_at_4 <- sum(UP_output_Lin_A[,3])
transitions_at_6 <- sum(UP_output_Lin_A[,4])
transitions_at_8 <- sum(UP_output_Lin_A[,5])
transitions_at_10 <- sum(UP_output_Lin_A[,6])
transitions_at_12 <- sum(UP_output_Lin_A[,7])
transitions_at_14 <- sum(UP_output_Lin_A[,8])
transitions_at_16 <- sum(UP_output_Lin_A[,9])
transitions_at_18 <- sum(UP_output_Lin_A[,10])
transitions_at_20 <- sum(UP_output_Lin_A[,11])

UP_smallRNA_A_Table_of_new_transitions <- rbind(transitions_at_2,
                                                transitions_at_4, 
                                                transitions_at_6,
                                                transitions_at_8, 
                                                transitions_at_10, 
                                                transitions_at_12, 
                                                transitions_at_14, 
                                                transitions_at_16, 
                                                transitions_at_18, 
                                                transitions_at_20)

colnames(UP_smallRNA_A_Table_of_new_transitions) <- c("A")



# small_RNA_B

smallRNA_Tab_B<-cbind(rep(0, length=nrow(smallRNA_binarised_table_B_2.25)),smallRNA_binarised_table_B_2.25)
colnames(smallRNA_Tab_B)[1]<-0
#shows that the "start" has no epimutations. 



UP_output_Lin_B<- c()

for(i in 1:nrow(smallRNA_Tab_B)){
  
  UP_output_Lin_B <-rbind(UP_output_Lin_B, UP_transition_func(smallRNA_Tab_B[i,], input_name=row.names(smallRNA_Tab_B)[i]))}

colnames(UP_output_Lin_B) <- c("0", "2", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(UP_output_Lin_B) <- rownames(smallRNA_binarised_table_B_2.25)

# Now take the sum of NEW UP events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(UP_output_Lin_B[,2])
transitions_at_6 <- sum(UP_output_Lin_B[,3])
transitions_at_8 <- sum(UP_output_Lin_B[,4])
transitions_at_10 <- sum(UP_output_Lin_B[,5])
transitions_at_12 <- sum(UP_output_Lin_B[,6])
transitions_at_14 <- sum(UP_output_Lin_B[,7])
transitions_at_16 <- sum(UP_output_Lin_B[,8])
transitions_at_18 <- sum(UP_output_Lin_B[,9])
transitions_at_20 <- sum(UP_output_Lin_B[,10])

UP_smallRNA_B_Table_of_new_transitions <- rbind(transitions_at_2,
                                                transitions_at_6,
                                                transitions_at_8, 
                                                transitions_at_10, 
                                                transitions_at_12, 
                                                transitions_at_14, 
                                                transitions_at_16, 
                                                transitions_at_18, 
                                                transitions_at_20)

colnames(UP_smallRNA_B_Table_of_new_transitions) <- c("B")



# small_RNA_C

smallRNA_Tab_C<-cbind(rep(0, length=nrow(smallRNA_binarised_table_C_2.25)),smallRNA_binarised_table_C_2.25)
colnames(smallRNA_Tab_C)[1]<-0
#shows that the "start" has no epimutations. 




UP_output_Lin_C<- c()

for(i in 1:nrow(smallRNA_Tab_C)){
  
  UP_output_Lin_C <-rbind(UP_output_Lin_C, UP_transition_func(smallRNA_Tab_C[i,], input_name=row.names(smallRNA_Tab_C)[i]))}

colnames(UP_output_Lin_C) <- c("0", "2", "4", "6", "8", "10", "12", "16", "18", "20")
row.names(UP_output_Lin_C) <- rownames(smallRNA_binarised_table_C_2.25)

# Now take the sum of NEW events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(UP_output_Lin_C[,2])
transitions_at_4 <- sum(UP_output_Lin_C[,3])
transitions_at_6 <- sum(UP_output_Lin_C[,4])
transitions_at_8 <- sum(UP_output_Lin_C[,5])
transitions_at_10 <- sum(UP_output_Lin_C[,6])
transitions_at_12 <- sum(UP_output_Lin_C[,7])
transitions_at_16 <- sum(UP_output_Lin_C[,8])
transitions_at_18 <- sum(UP_output_Lin_C[,9])
transitions_at_20 <- sum(UP_output_Lin_C[,10])

UP_smallRNA_C_Table_of_new_transitions <- rbind(transitions_at_2,
                                                transitions_at_4,
                                                transitions_at_6, 
                                                transitions_at_8, 
                                                transitions_at_10, 
                                                transitions_at_12, 
                                                transitions_at_16, 
                                                transitions_at_18, 
                                                transitions_at_20)

colnames(UP_smallRNA_C_Table_of_new_transitions) <- c("C")



# DOWN

# small_RNA_A

DOWN_output_Lin_A<- c()

for(i in 1:nrow(smallRNA_Tab_A)){
  
  DOWN_output_Lin_A <-rbind(DOWN_output_Lin_A, DOWN_transition_func(smallRNA_Tab_A[i,], input_name=row.names(smallRNA_Tab_A)[i]))}

colnames(DOWN_output_Lin_A) <- c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_Lin_A) <- rownames(smallRNA_binarised_table_A_2.25)

# Now take the sum of NEW DOWN events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(DOWN_output_Lin_A[,2])
transitions_at_4 <- sum(DOWN_output_Lin_A[,3])
transitions_at_6 <- sum(DOWN_output_Lin_A[,4])
transitions_at_8 <- sum(DOWN_output_Lin_A[,5])
transitions_at_10 <- sum(DOWN_output_Lin_A[,6])
transitions_at_12 <- sum(DOWN_output_Lin_A[,7])
transitions_at_14 <- sum(DOWN_output_Lin_A[,8])
transitions_at_16 <- sum(DOWN_output_Lin_A[,9])
transitions_at_18 <- sum(DOWN_output_Lin_A[,10])
transitions_at_20 <- sum(DOWN_output_Lin_A[,11])

DOWN_smallRNA_A_Table_of_new_transitions <- rbind(transitions_at_2,
                                                  transitions_at_4, 
                                                  transitions_at_6,
                                                  transitions_at_8, 
                                                  transitions_at_10, 
                                                  transitions_at_12, 
                                                  transitions_at_14, 
                                                  transitions_at_16, 
                                                  transitions_at_18, 
                                                  transitions_at_20)

colnames(DOWN_smallRNA_A_Table_of_new_transitions) <- c("A")



# small_RNA_B

DOWN_output_Lin_B<- c()

for(i in 1:nrow(smallRNA_Tab_B)){
  
  DOWN_output_Lin_B <-rbind(DOWN_output_Lin_B, DOWN_transition_func(smallRNA_Tab_B[i,], input_name=row.names(smallRNA_Tab_B)[i]))}

colnames(DOWN_output_Lin_B) <- c("0", "2", "6", "8", "10", "12", "14", "16", "18", "20")
row.names(DOWN_output_Lin_B) <- rownames(smallRNA_binarised_table_B_2.25)

# Now take the sum of NEW DOWN events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(DOWN_output_Lin_B[,2])
transitions_at_6 <- sum(DOWN_output_Lin_B[,3])
transitions_at_8 <- sum(DOWN_output_Lin_B[,4])
transitions_at_10 <- sum(DOWN_output_Lin_B[,5])
transitions_at_12 <- sum(DOWN_output_Lin_B[,6])
transitions_at_14 <- sum(DOWN_output_Lin_B[,7])
transitions_at_16 <- sum(DOWN_output_Lin_B[,8])
transitions_at_18 <- sum(DOWN_output_Lin_B[,9])
transitions_at_20 <- sum(DOWN_output_Lin_B[,10])

DOWN_smallRNA_B_Table_of_new_transitions <- rbind(transitions_at_2,
                                                  transitions_at_6,
                                                  transitions_at_8, 
                                                  transitions_at_10, 
                                                  transitions_at_12, 
                                                  transitions_at_14, 
                                                  transitions_at_16, 
                                                  transitions_at_18, 
                                                  transitions_at_20)

colnames(DOWN_smallRNA_B_Table_of_new_transitions) <- c("B")



# small_RNA_C

DOWN_output_Lin_C<- c()

for(i in 1:nrow(smallRNA_Tab_C)){
  
  DOWN_output_Lin_C <-rbind(DOWN_output_Lin_C, DOWN_transition_func(smallRNA_Tab_C[i,], input_name=row.names(smallRNA_Tab_C)[i]))}

colnames(DOWN_output_Lin_C) <- c("0", "2", "4", "6", "8", "10", "12", "16", "18", "20")
row.names(DOWN_output_Lin_C) <- rownames(smallRNA_binarised_table_C_2.25)

# Now take the sum of NEW DOWN events (i.e. (0 to 1), (-1 to 1) transitions) that occur at each generation, 

transitions_at_2 <- sum(DOWN_output_Lin_C[,2])
transitions_at_4 <- sum(DOWN_output_Lin_C[,3])
transitions_at_6 <- sum(DOWN_output_Lin_C[,4])
transitions_at_8 <- sum(DOWN_output_Lin_C[,5])
transitions_at_10 <- sum(DOWN_output_Lin_C[,6])
transitions_at_12 <- sum(DOWN_output_Lin_C[,7])
transitions_at_16 <- sum(DOWN_output_Lin_C[,8])
transitions_at_18 <- sum(DOWN_output_Lin_C[,9])
transitions_at_20 <- sum(DOWN_output_Lin_C[,10])

DOWN_smallRNA_C_Table_of_new_transitions <- rbind(transitions_at_2,
                                                  transitions_at_4,
                                                  transitions_at_6, 
                                                  transitions_at_8, 
                                                  transitions_at_10, 
                                                  transitions_at_12, 
                                                  transitions_at_16, 
                                                  transitions_at_18, 
                                                  transitions_at_20)

colnames(DOWN_smallRNA_C_Table_of_new_transitions) <- c("C")


# ------------------------------------------------------------------------------------------------------------------

# For each UP/DOWN data set, normalise the number of transitions per generation to a % of the total number of genes / regulatory element loci in that data set

# RNA

UP_RNA_A_Table_of_new_transitions_percent <- UP_RNA_A_Table_of_new_transitions/27779*100
UP_RNA_B_Table_of_new_transitions_percent <- UP_RNA_B_Table_of_new_transitions/27779*100
UP_RNA_C_Table_of_new_transitions_percent <- UP_RNA_C_Table_of_new_transitions/27779*100

DOWN_RNA_A_Table_of_new_transitions_percent <- DOWN_RNA_A_Table_of_new_transitions/27779*100
DOWN_RNA_B_Table_of_new_transitions_percent <- DOWN_RNA_B_Table_of_new_transitions/27779*100
DOWN_RNA_C_Table_of_new_transitions_percent <- DOWN_RNA_C_Table_of_new_transitions/27779*100


# ATAC

UP_ATAC_A_Table_of_new_transitions_percent <- UP_ATAC_A_Table_of_new_transitions/31618*100
UP_ATAC_B_Table_of_new_transitions_percent <- UP_ATAC_B_Table_of_new_transitions/31618*100
UP_ATAC_C_Table_of_new_transitions_percent <- UP_ATAC_C_Table_of_new_transitions/31618*100

DOWN_ATAC_A_Table_of_new_transitions_percent <- DOWN_ATAC_A_Table_of_new_transitions/31618*100
DOWN_ATAC_B_Table_of_new_transitions_percent <- DOWN_ATAC_B_Table_of_new_transitions/31618*100
DOWN_ATAC_C_Table_of_new_transitions_percent <- DOWN_ATAC_C_Table_of_new_transitions/31618*100


# small RNA 

UP_smallRNA_A_Table_of_new_transitions_percent <- UP_smallRNA_A_Table_of_new_transitions/8479*100
UP_smallRNA_B_Table_of_new_transitions_percent <- UP_smallRNA_B_Table_of_new_transitions/8479*100
UP_smallRNA_C_Table_of_new_transitions_percent <- UP_smallRNA_C_Table_of_new_transitions/8479*100

DOWN_smallRNA_A_Table_of_new_transitions_percent <- DOWN_smallRNA_A_Table_of_new_transitions/8479*100
DOWN_smallRNA_B_Table_of_new_transitions_percent <- DOWN_smallRNA_B_Table_of_new_transitions/8479*100
DOWN_smallRNA_C_Table_of_new_transitions_percent <- DOWN_smallRNA_C_Table_of_new_transitions/8479*100

# ----------------------------------------------------------------------------------------------------------------------

# Examining concordance between lineages A, B & C for each data set

library(ggplot2)

# Using the Kruskal-Wallis & Wilcoxon tests to compare the distributions 


# Supplementary Figure 2 contains the following plots

# A. RNA_UP_Lineage_compare
# B. RNA_DOWN_Lineage_compare
# C. ATAC_UP_Lineage_compare
# D. ATAC_DOWN_Lineage_compare
# E. smallRNA_UP_Lineage_compare
# F. smallRNA_DOWN_Lineage_compare



# -------------------
# RNA concordancy

# Compare UP transitions RNA ABC 

datk <- c(UP_RNA_A_Table_of_new_transitions_percent, UP_RNA_B_Table_of_new_transitions_percent, UP_RNA_C_Table_of_new_transitions_percent)
namk <- c(rep("A", length(UP_RNA_A_Table_of_new_transitions_percent)), rep("B", length(UP_RNA_B_Table_of_new_transitions_percent)), rep("C", length(UP_RNA_C_Table_of_new_transitions_percent)))

k <- data.frame(data = datk, Lineage = namk)

kruskal.test(data~Lineage, data = k)

pairwise.wilcox.test(k$data, k$Lineage, p.adj='bonferroni', exact=F)

k$Lineage <- factor(k$Lineage, levels = c("A", "B", "C"))


# Supplementary Figure 2 A

RNA_UP_Lineage_compare <- ggplot(k, aes(x=Lineage, y=data, color = Lineage)) + 
  geom_boxplot(fatten = 1, lwd = 1)+
  scale_color_manual(values=c("dodgerblue4", "turquoise4", "mediumaquamarine"))+
  ylim(0, 3.2)+
  labs(y = "% of genes showing increased \nexpression levels\n", x = "\nLineage")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center', dotsize=0.5)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  ggtitle(paste("Supplementary Figure 2 A"))

RNA_UP_Lineage_compare


# Compare DOWN transitions RNA ABC

datk <- c(DOWN_RNA_A_Table_of_new_transitions_percent, DOWN_RNA_B_Table_of_new_transitions_percent, DOWN_RNA_C_Table_of_new_transitions_percent)
namk <- c(rep("A", length(DOWN_RNA_A_Table_of_new_transitions_percent)), rep("B", length(DOWN_RNA_B_Table_of_new_transitions_percent)), rep("C", length(DOWN_RNA_C_Table_of_new_transitions_percent)))

k <- data.frame(data = datk, Lineage = namk)

kruskal.test(data~Lineage, data = k)

pairwise.wilcox.test(k$data, k$Lineage, p.adj='bonferroni', exact=F)

k$Lineage <- factor(k$Lineage, levels = c("A", "B", "C"))


# Supplementary Figure 2 B


RNA_DOWN_Lineage_compare <- ggplot(k, aes(x=Lineage, y=data, color = Lineage)) + 
  geom_boxplot(fatten = 1, lwd = 1)+
  scale_color_manual(values=c("dodgerblue4", "turquoise4", "mediumaquamarine"))+
  ylim(0, 3.2)+
  labs(y = "% of genes showing decreased \nexpression levels\n", x = "\nLineage")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center', dotsize=0.5)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  ggtitle(paste("Supplementary Figure 2 B"))

RNA_DOWN_Lineage_compare

#-------------------------

# ATAC concordancy

# Compare UP transitions ATAC ABC

datk <- c(UP_ATAC_A_Table_of_new_transitions_percent, UP_ATAC_B_Table_of_new_transitions_percent, UP_ATAC_C_Table_of_new_transitions_percent)
namk <- c(rep("A", length(UP_ATAC_A_Table_of_new_transitions_percent)), rep("B", length(UP_ATAC_B_Table_of_new_transitions_percent)), rep("C", length(UP_ATAC_C_Table_of_new_transitions_percent)))

k <- data.frame(data = datk, Lineage = namk)

kruskal.test(data~Lineage, data = k)

pairwise.wilcox.test(k$data, k$Lineage, p.adj='bonferroni', exact=F)

k$Lineage <- factor(k$Lineage, levels = c("A", "B", "C"))


# Supplementary Figure 2 C 


ATAC_UP_Lineage_compare <- ggplot(k, aes(x=Lineage, y=data, color = Lineage)) + 
  geom_boxplot(fatten = 1, lwd = 1)+
  scale_color_manual(values=c("indianred3", "indianred1", "lightpink"))+
  ylim(0, 3.2)+
  labs(y = "% of loci showing increased \nchromatin accessibility\n", x = "\nLineage")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center', dotsize=0.5)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  ggtitle(paste("Supplementary Figure 2 C"))

ATAC_UP_Lineage_compare


# Compare DOWN transitions ATAC ABC

datk <- c(DOWN_ATAC_A_Table_of_new_transitions_percent, DOWN_ATAC_B_Table_of_new_transitions_percent, DOWN_ATAC_C_Table_of_new_transitions_percent)
namk <- c(rep("A", length(DOWN_ATAC_A_Table_of_new_transitions_percent)), rep("B", length(DOWN_ATAC_B_Table_of_new_transitions_percent)), rep("C", length(DOWN_ATAC_C_Table_of_new_transitions_percent)))

k <- data.frame(data = datk, Lineage = namk)

kruskal.test(data~Lineage, data = k)

pairwise.wilcox.test(k$data, k$Lineage, p.adj='bonferroni', exact=F)

k$Lineage <- factor(k$Lineage, levels = c("A", "B", "C"))


# Supplementary Figure 2 D

ATAC_DOWN_Lineage_compare <- ggplot(k, aes(x=Lineage, y=data, color = Lineage)) + 
  geom_boxplot(fatten = 1, lwd = 1)+
  scale_color_manual(values=c("indianred3", "indianred1", "lightpink"))+
  ylim(0, 3.2)+
  labs(y = "% of loci showing decreased \nchromatin accessibility\n", x = "\nLineage")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center', dotsize=0.5)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  ggtitle(paste("Supplementary Figure 2 D"))

# ATAC_DOWN_Lineage_compare



# ----------------

# Compare UP transitions small RNA ABC 

datk <- c(UP_smallRNA_A_Table_of_new_transitions_percent, UP_smallRNA_B_Table_of_new_transitions_percent, UP_smallRNA_C_Table_of_new_transitions_percent)
namk <- c(rep("A", length(UP_smallRNA_A_Table_of_new_transitions_percent)), rep("B", length(UP_smallRNA_B_Table_of_new_transitions_percent)), rep("C", length(UP_smallRNA_C_Table_of_new_transitions_percent)))

k <- data.frame(data = datk, Lineage = namk)

kruskal.test(data~Lineage, data = k)

pairwise.wilcox.test(k$data, k$Lineage, p.adj='bonferroni', exact=F)

k$Lineage <- factor(k$Lineage, levels = c("A", "B", "C"))


# Supplementary Figure 2 E

smallRNA_UP_Lineage_compare <- ggplot(k, aes(x=Lineage, y=data, color = Lineage)) + 
  geom_boxplot(fatten = 1, lwd = 1)+
  scale_color_manual(values=c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4"))+
  ylim(0, 3.2)+
  labs(y = "% of genes showing increased \nexpression levels\n", x = "\nLineage")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center', dotsize=0.5)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  ggtitle(paste("Supplementary Figure 2 E"))

 smallRNA_UP_Lineage_compare

# Compare DOWN transitions small RNA ABC

datk <- c(DOWN_smallRNA_A_Table_of_new_transitions_percent, DOWN_smallRNA_B_Table_of_new_transitions_percent, DOWN_smallRNA_C_Table_of_new_transitions_percent)
namk <- c(rep("A", length(DOWN_smallRNA_A_Table_of_new_transitions_percent)), rep("B", length(DOWN_smallRNA_B_Table_of_new_transitions_percent)), rep("C", length(DOWN_smallRNA_C_Table_of_new_transitions_percent)))

k <- data.frame(data = datk, Lineage = namk)

kruskal.test(data~Lineage, data = k)

pairwise.wilcox.test(k$data, k$Lineage, p.adj='bonferroni', exact=F)

k$Lineage <- factor(k$Lineage, levels = c("A", "B", "C"))


# Supplementary Figure 2 F

smallRNA_DOWN_Lineage_compare <- ggplot(k, aes(x=Lineage, y=data, color = Lineage)) + 
  geom_boxplot(fatten = 1, lwd = 1)+
  scale_color_manual(values=c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4"))+
  ylim(0, 3.2)+
  labs(y = "% of genes showing decreased \nexpression levels\n", x = "\nLineage")+
  geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center', dotsize=0.5)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  ggtitle(paste("Supplementary Figure 2 F"))

smallRNA_DOWN_Lineage_compare 

# --------------------------------------------------------------------------------------------------------------


# Now we have shown that UP and DOWN chromatin are concordant and UP and DOWN RNA are concordant we can combine the UP and DOWNS within the categories and compare again

# Also now include small RNA


ATAC_UP_ABC <- rbind(UP_ATAC_A_Table_of_new_transitions_percent, UP_ATAC_B_Table_of_new_transitions_percent, UP_ATAC_C_Table_of_new_transitions_percent)

ATAC_DOWN_ABC <- rbind(DOWN_ATAC_A_Table_of_new_transitions_percent, DOWN_ATAC_B_Table_of_new_transitions_percent, DOWN_ATAC_C_Table_of_new_transitions_percent)

RNA_UP_ABC <- rbind(UP_RNA_A_Table_of_new_transitions_percent, UP_RNA_B_Table_of_new_transitions_percent, UP_RNA_C_Table_of_new_transitions_percent)

RNA_DOWN_ABC <-  rbind(DOWN_RNA_A_Table_of_new_transitions_percent, DOWN_RNA_B_Table_of_new_transitions_percent, DOWN_RNA_C_Table_of_new_transitions_percent)


smallRNA_UP_ABC <- rbind(UP_smallRNA_A_Table_of_new_transitions_percent, 
                         UP_smallRNA_B_Table_of_new_transitions_percent, 
                         UP_smallRNA_C_Table_of_new_transitions_percent)


smallRNA_DOWN_ABC <- rbind(DOWN_smallRNA_A_Table_of_new_transitions_percent, 
                           DOWN_smallRNA_B_Table_of_new_transitions_percent, 
                           DOWN_smallRNA_C_Table_of_new_transitions_percent)





ATAC_UP_DOWN_ABC <- rbind(ATAC_UP_ABC, ATAC_DOWN_ABC)

RNA_UP_DOWN_ABC <- rbind(RNA_UP_ABC, RNA_DOWN_ABC)

smallRNA_UP_DOWN_ABC <- rbind(smallRNA_UP_ABC, smallRNA_DOWN_ABC)



datk <- c(ATAC_UP_DOWN_ABC, RNA_UP_DOWN_ABC, smallRNA_UP_DOWN_ABC)
namk <- c(rep("Changes in chromatin state", length(ATAC_UP_DOWN_ABC)), 
          rep("Changes in RNA expression", length(RNA_UP_DOWN_ABC)), 
          rep("Changes in small RNA level", length(smallRNA_UP_DOWN_ABC)))

k <- data.frame(data = datk, Category = namk)



kruskal.test(data~Category, data = k)


pairwise.wilcox.test(k$data, k$Category, p.adj='bonferroni', exact=F)


k$Category <- factor(k$Category, levels = c("Changes in RNA expression",
                                            "Changes in chromatin state",
                                            "Changes in small RNA level"))



# Figure 2 A

# Comparing transitions per generation across all data sets 

smallRNA_ATAC_RNA_overall_category_compare <- ggplot(k, aes(x=Category, y=data, color = Category)) + 
  geom_boxplot(fatten = 1, lwd = 1)+
  scale_color_manual(values=c("brown3", "darkslategray4", "darkgoldenrod"))+
  ylim(0, 3.2)+
  labs(y = "% of genes/regulatory elements undergoing \nchanges per generation\n", x = NULL)+
  geom_dotplot(binaxis='y', binwidth = 0.10, stackdir='center', dotsize=0.5)+
  theme_bw()+
  theme_linedraw()+
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 45, hjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+
  ggtitle(paste("Figure 2 A"))

smallRNA_ATAC_RNA_overall_category_compare 






# -----------------------------------------------------------------------------------------------------

# At this point we produce the RNA-centric table 

# Here we show for each gene for which we have RNAseq data, whether they have enhancer/promoter coordinates as per Ahringer data set, if RNA expression changed and additional data 

# Produce a data table of all genes that have enhancer/promoter coordinates and their RNA expression change and chromatin epimutation status 

RNA_list <- list(RNA_binarised_table_A_2.25_no_21URs, RNA_binarised_table_B_2.25_no_21URs, RNA_binarised_table_C_2.25_no_21URs)

ATAC_ep_list <- list(ATAC_binarised_table_A_2.25_ehancer_promoter, ATAC_binarised_table_B_2.25_ehancer_promoter, ATAC_binarised_table_C_2.25_ehancer_promoter)




# We will need piRNA cluster genes for this

#----------
# Coordinates for piRNA cluster genes (5 - 7.5 & 13 - 17 MB from Ruby et al. 2006) have been intersected with RNA coordinates 

RNA_piRNA_domains <- read.table("~/Documents/Computational_work /Alternate_Generation_Experiment/2. RNA/RNA R objects/RNA_Coords_intersected_piRNA_c.bed", header = FALSE)

piRNA_domains <- RNA_piRNA_domains[RNA_piRNA_domains$V5==1, ]


# When we are considering the genes with chromatin domain annotations we can only consider the genes which map to Ahringer as this is where we get the domain annotations from. 

all_Ahr_RNA_genes <- unique(Ahringer_single_gene_ref_table$Gene)

piRNA_cluster_genes <- intersect(unique(piRNA_domains$V4), all_Ahr_RNA_genes)
#-----------



library(stringr)


lin <- c("A", "B", "C")

All_lins_integrated_table <- c()

for(x in 1:length(RNA_list)){
  
  RNA_bin <- RNA_list[[x]]
  ATAC_bin <- ATAC_ep_list[[x]]
  
  Lineage <- lin[[x]]
  
  
  integrated_table <- c()
  
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    
    coord <- rownames(RNA_bin)[[e]]
    gene <- strsplit(coord, ":")[[1]][4]
    
    
    
    AllUP <- 0
    AllDOWN <- 0
    MixedUPDOWN <- 0
    
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      
      
      # Direction of RNA expression changes  
      
      # Determine the original direction of the RNA events    
      
      events <- sum(abs(RNA_bin[coord, ]))
      
      sum <- sum(RNA_bin[coord, ])
      
      
      if(sum == events){
        
        AllUP <- 1
      }     
      
      
      if(sum == -1*(events)){
        
        AllDOWN <- 1
        
      }     
      
      
      if(!abs(sum) == events){
        
        MixedUPDOWN <- 1
      }  
      
      
      
    }
    
    
    RNA_gens <- "0"
    
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      
      RNA_gens <- names(which(abs(RNA_bin[coord, ]) > 0))
      
      if(length(RNA_gens) > 1){
        
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          
          RNA_inherited <- 1
        } 
      }
    }
    
    
    
    maps_to_Ahr <- 0
    
    ep_mut <- 0
    
    ep_inherited <- 0
    
    ep_gens <- 0
    
    multiple_ep_gens  <- 0
    
    piRNA_cluster <- 0
    
    if(gene %in% piRNA_cluster_genes){
      piRNA_cluster <- 1
    }
    
    
    if(gene %in% Ahringer_single_gene_ref_table$Gene){
      
      maps_to_Ahr <- 1 
      
      chrom_loci <- Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene %in% gene, 2]
      
      if(sum(abs(ATAC_bin[chrom_loci, ]))>0){
        
        ep_mut <- 1
        
        
        if(length(chrom_loci) ==1){
          
          ep_gens <- names(which(abs(ATAC_bin[chrom_loci, ]) > 0)) 
          
        }
        
        
        if(length(chrom_loci) > 1){
          
          test <- abs(ATAC_bin[chrom_loci, ])
          
          ep_gens <- colnames(test)[colSums(test)>0]
        }
        
        
        if(length(ep_gens)> 1){
          
          ep_gens <- paste(ep_gens, collapse = "_")
          
          multiple_ep_gens <- 1
          
          
        }
        
        
        inherited <- c()
        for(u in 1:length(chrom_loci)){
          
          inherit <- 0 
          
          each_row <- names(which(abs(ATAC_bin[chrom_loci[[u]], ]) > 0))
          
          if(2 %in% diff(as.numeric(each_row))){
            inherit <- 1}
          
          inherited <- c(inherited, inherit) 
          
          if(sum(inherited > 0)){
            
            ep_inherited <- 1
          }
          
        }
        
      }}
    
    
    
    time_matched <- 0
    
    
    ATAC_UP <- 0
    ATAC_DOWN <- 0
    ATAC_mixed_dir <- 0
    
    
    epimut_is_Active <- 0
    epimut_is_Regulated <- 0
    epimut_is_Border <- 0
    epimut_is_ChrX <- 0
    Mixed_Domain_epimuts <- 0
    
    epimut_is_Promoter <- 0
    epimut_is_Enhancer <- 0
    Mixed_Element_epimuts <- 0  
    
    save_dir <- c()
    
    Domain <- c()
    
    Element <- c()
    
    if(ep_mut == 1){
      
      
      
      for(u in 1:length(chrom_loci)){
        
        chrom_row <- ATAC_bin[chrom_loci[[u]], ]
        
        events <- sum(abs(chrom_row))
        
        if(events > 0){
          
          
          # Domain
          
          
          domain <-  Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$ep_locus %in% chrom_loci[[u]], 4]
          
          Domain <- c(Domain, domain)
          
          # Element
          
          
          element <-  Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$ep_locus %in% chrom_loci[[u]], 3]
          
          Element <- c(Element, element)
          
          
          # Direction
          
          if(!-1 %in% chrom_row){
            
            dir <- 1
          }
          
          
          if(!1 %in% chrom_row){
            
            dir <- -1
          }
          
          if(abs(sum(chrom_row)) < events){
            
            dir <- 0
          }
          
          save_dir <- c(save_dir, dir)  
          
        }
        
      }
      
      
      # Domain of epimutation
      
      
      if("A" %in% unique(Domain)){ 
        
        epimut_is_Active <- 1
      }  
      
      if("R" %in% unique(Domain)){ 
        
        epimut_is_Regulated <- 1
      }  
      
      if("border" %in% unique(Domain)){ 
        
        epimut_is_Border <- 1
      }
      
      if("." %in% unique(Domain)){ 
        
        epimut_is_ChrX <- 1
      }  
      
      if(length(unique(Domain)) > 1){
        
        Mixed_Domain_epimuts <- 1}
      
      
      
      
      # Element type of epimutation
      
      
      if("putative_enhancer" %in% unique(Element)){ 
        
        epimut_is_Enhancer <- 1
      }  
      
      if("coding_promoter" %in% unique(Element)){ 
        
        epimut_is_Promoter <- 1
      }  
      
      if(length(unique(Element)) > 1){
        
        Mixed_Element_epimuts <- 1}
      
      
      
      
      # Direction of epimutation
      
      #ATAC_mixed_dir <- 1
      
      if(sum(save_dir) == length(save_dir) | sum(save_dir) > 0){
        
        ATAC_UP <- 1
        #ATAC_mixed_dir <- 0
      }
      
      if(sum(save_dir)== -1*length(save_dir) | sum(save_dir) < 0){
        
        ATAC_DOWN <- 1
        #ATAC_mixed_dir <- 0
      }
      
      if(sum(save_dir)== 0){
        
        ATAC_mixed_dir <- 1
        
      }
      
      
      
    }
    
    
    
    Direction_match <- 0
    
    if(ATAC_UP==1&AllUP==1){
      
      Direction_match <- 1
      
    }
    
    if(ATAC_DOWN==1&AllDOWN==1){
      
      Direction_match <- 1
      
    } 
    
    
    
    
    
    if(RNA_mut==1&ep_mut==1){
      
      chrom_gens <- unlist(str_split(ep_gens, "_"))
      
      rna_gens <- unlist(str_split(RNA_gens, "_"))
      
      matched_gens <- chrom_gens[chrom_gens %in% rna_gens]
      
      matched <-  length(which(chrom_gens %in% rna_gens))
      
      
      
      if(matched > 0){
        
        time_matched <- 1
        
      }
      
      
      
      
    }
    
    save <- data.frame(Lineage, gene, coord, maps_to_Ahr, RNA_mut, RNA_gens, RNA_inherited, 
                       epimut_is_Active, epimut_is_Regulated, 
                       epimut_is_Border, epimut_is_ChrX, Mixed_Domain_epimuts, 
                       epimut_is_Promoter,  epimut_is_Enhancer, Mixed_Element_epimuts, 
                       ep_mut, ep_gens,  multiple_ep_gens,  ep_inherited, time_matched, 
                       ATAC_UP, ATAC_DOWN, ATAC_mixed_dir, 
                       AllUP, AllDOWN, MixedUPDOWN, Direction_match, piRNA_cluster)
    
    
    integrated_table <- rbind(integrated_table, save)
    
  }
  All_lins_integrated_table <- rbind(All_lins_integrated_table, integrated_table)
}



colnames(All_lins_integrated_table) <- c(
  "Lineage", "gene", "RNA_coord", "gene_maps_to_ep", "is_RNA_exp_change",
  "RNA_mut_gens", "is_RNA_inherited", "epimut_is_Active", "epimut_is_Regulated","epimut_is_Border", "epimut_is_ChrX", "Mixed_Domain_epimuts", "epimut_is_Promoter",  "epimut_is_Enhancer", "Mixed_Element_epimuts","is_ep_epimut",  
  "ep_mut_gens","Multiple_chrom_gens","is_ep_inherited", "has_time_matched", "ATAC_UP", "ATAC_DOWN", "ATAC_mixed_dir",  "RNA_All_UP", "RNA_All_DOWN", "RNA_Mixed_Muts", "Direction_Match", "piRNA_cluster"
)

#-------------------------------------------------------------------------------------------------------


# Make a table to show the Ahringer data and epimutation state for each enhancer/promoter locus 


Ahringer_single_unique_locus <- unique(rownames(ATAC_binarised_table_A_2.25_ehancer_promoter))

lin <- c("A", "B", "C")

All_lins_ATAC_table <- c() 

for(a in 1:length(ATAC_ep_list)){
  
  ATAC_bin <- ATAC_ep_list[[a]]
  
  Lin <- lin[[a]]
  
  ep_table <- c()
  
  for(w in 1:length(Ahringer_single_unique_locus)){
    
    chrom_locus <-  Ahringer_single_unique_locus[[w]]
    
    gene <- Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$ep_locus %in% chrom_locus, 1]
    
    
    piRNA_cluster <- 0
    
    if(TRUE %in% (gene %in% piRNA_cluster_genes)){
      piRNA_cluster <- 1
    }
    
    
    if(length(gene) > 1){
      
      gene <- paste(gene, collapse = "_")
    }
    
    Ahringer_data <- Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$ep_locus == chrom_locus, ][1, 2:5]
    
    ATAC_muts <- sum(abs(ATAC_bin[chrom_locus, ]))
    
    
    inherited <- 0
    
    if(ATAC_muts > 1){
      
      each_row <- names(which(abs(ATAC_bin[chrom_locus, ]) > 0))
      
      if(2 %in% diff(as.numeric(each_row))){
        inherited <- 1}
      
    }
    
    
    
    
    ATAC_muts[ATAC_muts > 1] <- 1 
    
    Ahringer_data_save <- data.frame(Lin, gene, Ahringer_data, piRNA_cluster, ATAC_muts, inherited) 
    
    colnames(Ahringer_data_save) <- c("Lineage", "Gene", "Locus", "Enhancer_or_Promoter", "Chromatin_domain", "Tissue_specificity", "piRNA_cluster", "Epimutated", "Inherited")
    
    ep_table <- rbind(ep_table, Ahringer_data_save)
    
  }
  
  All_lins_ATAC_table <- rbind(All_lins_ATAC_table, ep_table)
}

#--------------------------------------------------------------------------------------------------------

# A small RNA centric table


list_target_genes <- rownames(smallRNA_binarised_table_A_2.25)

ABC_small_RNA_super_table <- list(smallRNA_binarised_table_A_2.25,
                                  smallRNA_binarised_table_B_2.25, 
                                  smallRNA_binarised_table_C_2.25)

names(ABC_small_RNA_super_table) <- c("A", "B", "C")




# We also include details on the target genes obtained from Beltran et al.

Gene_Features <- read.table("~/Desktop/gene_features_cel.txt")

piRNA_csr_hrde <- Gene_Features[, c(1, 8:10)]

piRNA_target_genes <- unique(piRNA_csr_hrde[piRNA_csr_hrde$piRNA_targets==TRUE, 1])

csr_genes <- unique(piRNA_csr_hrde[piRNA_csr_hrde$csr1_targets==TRUE,  1])

hrde_genes <- unique(piRNA_csr_hrde[piRNA_csr_hrde$hrde1_targets ==TRUE,  1])





small_RNA_annotation_table <- c()

for(t in 1:length(names(ABC_small_RNA_super_table))){
  
  Lineage <- names(ABC_small_RNA_super_table)[[t]]
  
  Lineage_table <- c()
  
  for(i in 1:length(list_target_genes)){
    
    Target_gene <- list_target_genes[[i]]
    
    
    # What chromatin domain is the potential target gene in?
    
    Gene_domain <- 0
    maps_to_Ahr <- 0
    
    if(Target_gene %in% Ahringer_single_gene_ref_table$Gene){
      
      maps_to_Ahr <- 1
      
      Gene_domain <- unique(Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == Target_gene, 4])
      
      if(length(Gene_domain)> 1){ 
        
        Gene_domain <- paste(Gene_domain, collapse = "_")
      } 
    }
    
    
    
    
    # What tissue specificity is pertinent to the potential target gene?   
    
    Tissue_specificity <- 0
    
    if( maps_to_Ahr == 1){
      
      Tissue_specificity <- unique(Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == Target_gene, 5])
      
      if(length(Tissue_specificity)> 1){
        
        Tissue_specificity <- paste(Tissue_specificity, collapse = "_")
      } 
    }
    
    
    
    # Is the potential target gene a target for piRNA triggered 22G RNAs?   
    
    piRNA_target <- 0
    
    if(Target_gene %in% piRNA_target_genes){
      
      piRNA_target <- 1
    } 
    
    
    # Is the potential target gene a target for HRDE-1 triggered 22G RNAs?   
    
    hrde_target <- 0
    
    if(Target_gene %in% hrde_genes){
      
      hrde_target <- 1
    } 
    
    # Is the potential target gene a target for CSR-1 triggered 22G RNAs?   
    
    csr_target <- 0
    
    if(Target_gene %in% csr_genes){
      
      csr_target <- 1
    } 
    
    
    # Is the potential target gene within the piRNA cluster?
    
    piRNA_cluster <- 0
    
    if(Target_gene %in% piRNA_cluster_genes){
      piRNA_cluster <- 1
    }
    
    
    
    
    # For the specific lineage, check the epimutation status of the gene 
    
    
    frame <- ABC_small_RNA_super_table[Lineage][[1]]
    
    epimutated <- 0
    
    if(Target_gene %in% rownames(frame)){
      
      epimutated <- 0
      
      if(sum(abs(ABC_small_RNA_super_table[[t]][Target_gene, ])) > 0){
        
        epimutated <- 1}
    }
    
    
    row <- data.frame(Lineage, Target_gene,  maps_to_Ahr, epimutated, Gene_domain, Tissue_specificity, piRNA_target, hrde_target, csr_target, piRNA_cluster)
    
    Lineage_table <- rbind(Lineage_table, row)
  }
  
  small_RNA_annotation_table <- rbind(small_RNA_annotation_table, Lineage_table)
}


#--------------------------------------------------------------------------------------------------------
# An RNA centric table with RNA genes, whether they map to small RNA, RNA expression change and small RNA epimutation status


RNA_list <- list(RNA_binarised_table_A_2.25_no_21URs, RNA_binarised_table_B_2.25_no_21URs, RNA_binarised_table_C_2.25_no_21URs)

small_RNA_list <- list(smallRNA_binarised_table_A_2.25, smallRNA_binarised_table_B_2.25, smallRNA_binarised_table_C_2.25)

smallRNA_genes <- rownames(smallRNA_binarised_table_A_2.25)

library(stringr)


lin <- c("A", "B", "C")

smallRNA_All_lins_integrated_table <- c()

for(x in 1:length(RNA_list)){
  
  RNA_bin <- RNA_list[[x]]
  smallRNA_bin <- small_RNA_list[[x]]
  
  Lineage <- lin[[x]]
  
  
  integrated_table <- c()
  
  ep_table <- c()
  
  for(e in 1:length(rownames(RNA_bin))){
    
    coord <- rownames(RNA_bin)[[e]]
    gene <- strsplit(coord, ":")[[1]][4]
    
    
    RNA_UP <- 0
    RNA_DOWN <- 0
    RNA_MixedUPDOWN <- 0
    
    RNA_mut <- 0
    
    if(sum(abs(RNA_bin[coord, ])) >0){
      RNA_mut <- 1
      
      
      # Direction of RNA expression changes  
      
      # Determine the original direction of the RNA events    
      
      events <- sum(abs(RNA_bin[coord, ]))
      
      sum <- sum(RNA_bin[coord, ])
      
      
      if(sum == events){
        
        RNA_UP <- 1
      }     
      
      
      if(sum == -1*(events)){
        
        RNA_DOWN <- 1
        
      }     
      
      
      if(!abs(sum) == events){
        
        RNA_MixedUPDOWN <- 1
      }  
      
    }
    
    
    RNA_gens <- "0"
    
    RNA_inherited <- 0
    
    if(RNA_mut ==1){
      
      RNA_gens <- names(which(abs(RNA_bin[coord, ]) > 0))
      
      if(length(RNA_gens) > 1){
        
        RNA_gens <- paste(RNA_gens, collapse = "_")
        
        
        if(2 %in% diff(as.numeric(strsplit(RNA_gens, "_")[[1]]))){
          
          RNA_inherited <- 1
        } 
      }
    }
    
    
    
    # The RNA coding locus has to a) be present in the small RNA genes so as to get a small RNA signal data point
    # b) be present in the Ahringer data set so as to get associated annotations on chromatin domain distribution
    
    
    # First assume the coding locus with RNAseq data is not present in small RNA 
    
    maps_to_Ahr <- 0
    
    maps_to_smallRNA <- 0
    
    smallRNA_mut <- 0
    
    smallRNA_inherited <- 0
    
    smallRNA_ep_gens <- 0
    
    multiple_smallRNA_ep_gens  <- 0
    
    smallRNA_UP <- 0
    smallRNA_DOWN <- 0
    smallRNA_mixed_dir <- 0
    
    piRNA_cluster <- 0
    
    
    if(gene %in% piRNA_cluster_genes){
      
      piRNA_cluster <- 1
    }
    
    if(gene %in% Ahringer_single_gene_ref_table$Gene){
      
      maps_to_Ahr <- 1} 
    
    if(gene %in% smallRNA_genes == T){
      
      maps_to_smallRNA <- 1
      
      if(sum(abs(smallRNA_bin[gene, ]))>0){
        
        
        smallRNA_mut <- 1
        
        target_row <- smallRNA_bin[gene, ]
        events <- sum(abs(smallRNA_bin[gene, ]))
        
        # Direction
        
        if(!-1 %in% target_row){
          
          smallRNA_UP <- 1
        }
        
        
        if(!1 %in% target_row){
          
          smallRNA_DOWN <- 1
        }
        
        if(abs(sum(target_row)) < events){
          
          smallRNA_mixed_dir <- 1
        }
      }
      
      if(smallRNA_mut ==1){
        
        
        # Generations
        
        smallRNA_ep_gens <- names(which(abs(smallRNA_bin[gene, ]) > 0)) 
        
        if(length(smallRNA_ep_gens)> 1){
          
          smallRNA_ep_gens <- paste(smallRNA_ep_gens, collapse = "_")
          
          multiple_smallRNA_ep_gens <- 1
        }
        
        
        # Inherited
        inherited <- c()
        
        inherit <- 0 
        
        each_row <- names(which(abs(smallRNA_bin[gene, ]) > 0))
        
        if(2 %in% diff(as.numeric(each_row))){
          inherit <- 1}
        
        inherited <- c(inherited, inherit) 
        
        if(sum(inherited > 0)){
          
          smallRNA_inherited <- 1
        }
        
        
        
        
      }}
    
    
    time_matched_to_RNA <- 0
    
    smallRNA_target_is_Active <- 0
    smallRNA_target_is_Regulated <- 0
    smallRNA_target_is_Border <- 0
    smallRNA_target_is_ChrX <- 0
    Mixed_Domain_smallRNA_targets <- 0
    Chromatin_epimut_reg_el <- 0
    Chromatin_epimut_inherited <- 0
    
    
    
    
    
    # Determine Domain of target gene
    
    # Domain can only be determined if the gene also maps to Ahringer
    # Similarly, regulatory element chromatin epimutation status can be determined
    
    if(maps_to_Ahr == 1){
      
      Domain <-  Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene %in% gene , 4]
      
      if(All_lins_integrated_table[All_lins_integrated_table$gene == gene & All_lins_integrated_table$Lineage== lin[[x]], 16] ==1 ){
        
        Chromatin_epimut_reg_el <- 1  
        
        Chromatin_epimut_inherited <- All_lins_integrated_table[All_lins_integrated_table$gene == gene & All_lins_integrated_table$Lineage== lin[[x]], 19]
        
      }
      
      
      if("A" %in% unique(Domain)){ 
        
        smallRNA_target_is_Active <- 1
      }  
      
      if("R" %in% unique(Domain)){ 
        
        smallRNA_target_is_Regulated <- 1
      }  
      
      if("border" %in% unique(Domain)){ 
        
        smallRNA_target_is_Border <- 1
      }
      
      if("." %in% unique(Domain)){ 
        
        smallRNA_target_is_ChrX <- 1
      }  
      
      if(length(unique(Domain)) > 1){
        
        Mixed_Domain_smallRNA_targets <- 1
      }
    }
    
    
    # Determine Direction matching
    
    
    Direction_match_to_RNA <- 0
    
    if(smallRNA_UP==1&RNA_UP==1){
      
      Direction_match_to_RNA <- 1
      
    }
    
    if(smallRNA_DOWN==1&RNA_DOWN==1){
      
      Direction_match_to_RNA <- 1
      
    } 
    
    
    
    # Determine time matching of small RNA to RNA
    
    if(RNA_mut==1&smallRNA_mut==1){
      
      smallRNA_gens <- unlist(str_split(smallRNA_ep_gens, "_"))
      
      rna_gens <- unlist(str_split(RNA_gens, "_"))
      
      matched_gens <- smallRNA_gens[smallRNA_gens %in% rna_gens]
      
      matched <-  length(which(smallRNA_gens %in% rna_gens))
      
      
      
      if(matched > 0){
        
        time_matched_to_RNA <- 1
        
      }
      
    }
    
    
    time_matched_to_chrom <- 0
    
    
    # Determine time matching of small RNA to ATAC
    
    if(Chromatin_epimut_reg_el==1&smallRNA_mut==1){
      
      smallRNA_gens <- unlist(str_split(smallRNA_ep_gens, "_"))
      
      chromatin_gens <- unlist(str_split(All_lins_integrated_table[All_lins_integrated_table$Lineage== lin[[x]] & All_lins_integrated_table$gene== gene, 17], "_"))
      
      matched_gens <- smallRNA_gens[smallRNA_gens %in% chromatin_gens]
      
      matched <-  length(which(smallRNA_gens %in% chromatin_gens))
      
      if(matched > 0){
        
        time_matched_to_chrom <- 1
        
      }
      
    }
    
    
    
    
    
    
    
    save <- data.frame(Lineage, gene, coord, maps_to_Ahr, maps_to_smallRNA, Chromatin_epimut_reg_el, Chromatin_epimut_inherited, RNA_mut, RNA_inherited, 
                       smallRNA_mut, smallRNA_target_is_Active, smallRNA_target_is_Regulated, 
                       smallRNA_target_is_Border, smallRNA_target_is_ChrX, Mixed_Domain_smallRNA_targets, 
                       smallRNA_inherited, time_matched_to_RNA, time_matched_to_chrom,
                       smallRNA_UP, smallRNA_DOWN, smallRNA_mixed_dir, 
                       RNA_UP, RNA_DOWN, RNA_MixedUPDOWN, Direction_match_to_RNA, 
                       piRNA_cluster)
    
    
    integrated_table <- rbind(integrated_table, save)
    
  }
  smallRNA_All_lins_integrated_table <- rbind(smallRNA_All_lins_integrated_table, integrated_table)
}


# ------------------------------------------------------------------------------------------------------
# Survival plots 

# Now we have shown that the 3 lineages are behaving similarly i.e. concordant 
# The rates of RNA seq change, ATAC change and small RNA change per generation are quite different 
# Now we want to examine the survival of runs of epigenetic changes/epimutations in each of these data sets


# In creating the below data tables, we assume the missing generations within a run of epimutations are also epimutated  

# At the end we will remove the "suspicious" epimutations, i.e. those that have the same onset and are the same length in all 3 lineages could have been environmentally driven



# ----------------------
# RNA: counting lengths of expression changes


# For RNA we will also identfy whether each specific RNA expression change run has an accompanying ATAC or small RNA epimutation



RNA_Tab_A<-cbind(rep(0, length=nrow(RNA_binarised_table_A_2.25_no_21URs)),RNA_binarised_table_A_2.25_no_21URs)
colnames(RNA_Tab_A)[1]<-"0"
RNA_Tab_A <- as.data.frame(RNA_Tab_A)


RNA_Tab_B<-cbind(rep(0, length=nrow(RNA_binarised_table_B_2.25_no_21URs)),RNA_binarised_table_B_2.25_no_21URs)
colnames(RNA_Tab_B)[1]<-"0"
RNA_Tab_B <- as.data.frame(RNA_Tab_B)


RNA_Tab_C<-cbind(rep(0, length=nrow(RNA_binarised_table_C_2.25_no_21URs)),RNA_binarised_table_C_2.25_no_21URs)
colnames(RNA_Tab_C)[1]<-"0"
RNA_Tab_C <- as.data.frame(RNA_Tab_C)



All_main_frames <- list(RNA_Tab_A, RNA_Tab_B, RNA_Tab_C)

names(All_main_frames) <- c("A", "B", "C")

compare_smallRNA <- list(smallRNA_binarised_table_A_2.25, smallRNA_binarised_table_B_2.25, smallRNA_binarised_table_C_2.25)
compare_ATAC <- list(ATAC_binarised_table_A_2.25_ehancer_promoter, ATAC_binarised_table_B_2.25_ehancer_promoter, ATAC_binarised_table_C_2.25_ehancer_promoter)



RNA_All_output_list_with_missing_as_muts <- list()

for(e in 1:length(All_main_frames)){
  
  main_frame <- as.data.frame(All_main_frames[[e]])
  
  smallRNA_compare_frame <- compare_smallRNA[[e]]
  ATAC_compare_frame <- compare_ATAC[[e]]
  
  
  Overall_output <- c()
  
  for(t in 1:nrow(main_frame)){
    
    vector_in <- main_frame[t, ]
    gen_names <- as.numeric(colnames(vector_in))
    
    number_transitions <- 0
    onset_gen <- 0
    tempL <- 0
    output_frame <- c()
    
    is_up <- 0
    is_down <- 0
    complete <- 0
    
    
    smallRNA_mut <- 0
    smallRNA_up <- 0
    smallRNA_down <- 0
    
    ATAC_mut <- 0
    ATAC_up <- 0
    ATAC_down <- 0
    
    # first deal with no epimutations in the lineage for a locus/gene
    
    if(sum(abs(vector_in)) == 0){
      
      length = 0
      is_up <- 0
      is_down <- 0
      gene_name <- rownames(vector_in)
      gene <- strsplit(gene_name, ":")[[1]][4]
      
      number_transitions <- 0 
      onset_gen <- 0
      
      smallRNA_mut <- 0
      smallRNA_up <- 0
      smallRNA_down <- 0
      
      ATAC_mut <- 0
      ATAC_up <- 0
      ATAC_down <- 0
      
      Lineage <- names(All_main_frames[e])
      complete <- 0
      
      
      save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete,  Lineage, is_up, is_down, smallRNA_mut, smallRNA_down, smallRNA_up, ATAC_mut, ATAC_down, ATAC_up)   
      
      output_frame <- rbind(output_frame, save_mut)
      
    }
    
    
    else
      
      if(sum(abs(vector_in)) > 0){
        
        for(i in 2:length(vector_in)){
          
          
          # if it is an UP transition, turn off any down transitions and start new epimutation
          
          if(vector_in[i]==1&vector_in[i-1]== 0|vector_in[i]== 1&vector_in[i-1]== -1) {
            
            
            if(is_down == 1){ # a down epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 1 
              
              length <-tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              
              
              # is there a  small RNA epimutation that corresponds?
              
              if(gene %in% rownames(smallRNA_compare_frame)){
                
                smallRNA_coord <- gene
                
                check <- onset_check : (i -2)
                
                if(sum(abs(smallRNA_compare_frame[smallRNA_coord, check])) > 0){
                  
                  smallRNA_mut <- 1
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) < sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_down <- 1}
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) == sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_up <- 1}
                }
              }
              
              
              # is there an ATAC epimutation that corresponds?
              
              if(gene %in% unique(Ahringer_single_gene_ref_table$Gene)){
                
                ATAC_coord <- Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == gene, 2]
                
                check <- onset_check : (i -2)
                
                if(sum(abs(ATAC_compare_frame[ATAC_coord, check])) > 0){
                  
                  ATAC_mut <- 1
                  
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == -1*sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_down <- 1}
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_up <- 1}
                }
              }
              
              
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, smallRNA_mut, smallRNA_down, smallRNA_up, ATAC_mut, ATAC_down, ATAC_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            # start the UP transition
            
            is_up <- 1
            is_down <- 0
            tempL <- 1
            
            number_transitions <- 1 # there is a transition to an UP epimutation at generation[i]
            onset_gen <- gen_names[i]
            
            onset_check <- i -1
          }
          
          
          
          # if it is a new DOWN transition, turn off any up transitions and start new epimutation
          
          if(vector_in[i]==-1&vector_in[i-1]== 0|vector_in[i]== -1&vector_in[i-1]== 1) {
            
            if(is_up == 1){ # an up epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 1
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              
              
              # is there a  small RNA epimutation that corresponds?
              
              if(gene %in% rownames(smallRNA_binarised_table_A_2.25)){
                
                
                smallRNA_coord <- gene
                
                check <- onset_check : (i -2)
                
                if(sum(abs(smallRNA_compare_frame[smallRNA_coord, check])) > 0){
                  
                  smallRNA_mut <- 1
                  
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) < sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_down <- 1}
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) == sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_up <- 1}
                }
              }
              
              
              # is there an ATAC epimutation that corresponds?
              
              if(gene %in% unique(Ahringer_single_gene_ref_table$Gene)){
                
                ATAC_coord <- Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == gene, 2]
                
                check <- onset_check : (i -2)
                
                if(sum(abs(ATAC_compare_frame[ATAC_coord, check])) > 0){
                  
                  ATAC_mut <- 1
                  
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == -1*sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_down <- 1}
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_up <- 1}
                }
              }
              
              
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, smallRNA_mut, smallRNA_down, smallRNA_up, ATAC_mut, ATAC_down, ATAC_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            
            # start the down transition
            
            is_down <- 1
            is_up <- 0
            tempL <- 1
            number_transitions <- 1 # there is a transition to a DOWN epimutation at generation[i]
            onset_gen <- gen_names[i]
            onset_check <- i-1
          }
          
          
          
          # if it is an UP epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_up==1){
            
            
            if(vector_in[i]==1&vector_in[i-1]==1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}  # the epimutation continues
            
            
            
            if(vector_in[i]==0&vector_in[i-1]==1){
              
              gene_name <- rownames(vector_in)
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from an UP epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              
              
              # is there a  small RNA epimutation that corresponds?
              
              if(gene %in% rownames(smallRNA_binarised_table_A_2.25)){
                
                
                smallRNA_coord <- gene
                
                check <- onset_check : (i -2)
                
                if(sum(abs(smallRNA_compare_frame[smallRNA_coord, check])) > 0){
                  
                  smallRNA_mut <- 1
                  
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) < sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_down <- 1}
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) == sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_up <- 1}
                }
              }
              
              
              # is there an ATAC epimutation that corresponds?
              
              if(gene %in% unique(Ahringer_single_gene_ref_table$Gene)){
                
                ATAC_coord <- Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == gene, 2]
                
                check <- onset_check : (i -2)
                
                if(sum(abs(ATAC_compare_frame[ATAC_coord, check])) > 0){
                  
                  ATAC_mut <- 1
                  
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == -1*sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_down <- 1}
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_up <- 1}
                }
              }            
              
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, smallRNA_mut, smallRNA_down, smallRNA_up, ATAC_mut, ATAC_down, ATAC_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              length <- 0
              complete <- 0
              
              smallRNA_mut <- 0
              smallRNA_down <- 0
              smallRNA_up <- 0
              
              ATAC_mut <- 0
              ATAC_down <- 0
              ATAC_up <- 0
              
              
            }
            
            
            
            # deal with an up epimutation extending to end of lineage
            
            if(vector_in[i]==1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 0
              
              
              
              # is there a  small RNA epimutation that corresponds?
              
              if(gene %in% rownames(smallRNA_binarised_table_A_2.25)){
                
                
                smallRNA_coord <- gene
                
                check <- onset_check : (i -1)
                
                if(sum(abs(smallRNA_compare_frame[smallRNA_coord, check])) > 0){
                  
                  smallRNA_mut <- 1
                  
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) < sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_down <- 1}
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) == sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_up <- 1}
                }
              }
              
              
              # is there an ATAC epimutation that corresponds?
              
              if(gene %in% unique(Ahringer_single_gene_ref_table$Gene)){
                
                ATAC_coord <- Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == gene, 2]
                
                check <- onset_check : (i -1)
                
                if(sum(abs(ATAC_compare_frame[ATAC_coord, check])) > 0){
                  
                  ATAC_mut <- 1
                  
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == -1*sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_down <- 1}
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_up <- 1}
                }
              }              
              
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, smallRNA_mut, smallRNA_down, smallRNA_up, ATAC_mut, ATAC_down, ATAC_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              length <- 0
              tempL <- 0
              complete <- 0
              
              smallRNA_mut <- 0
              smallRNA_down <- 0
              smallRNA_up <- 0
              
              ATAC_mut <- 0
              ATAC_down <- 0
              ATAC_up <- 0
              
            }
            
          }
          
          
          # if it is a DOWN epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_down==1){
            
            if(vector_in[i]==-1&vector_in[i-1]==-1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}
            
            
            if(vector_in[i]==0&vector_in[i-1]==-1){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from a DOWN epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              
              
              # is there a  small RNA epimutation that corresponds?
              
              if(gene %in% rownames(smallRNA_binarised_table_A_2.25)){
                
                
                smallRNA_coord <- gene
                
                check <- onset_check : (i -2)
                
                if(sum(abs(smallRNA_compare_frame[smallRNA_coord, check])) > 0){
                  
                  smallRNA_mut <- 1
                  
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) < sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_down <- 1}
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) == sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_up <- 1}
                }
              }
              
              
              # is there an ATAC epimutation that corresponds?
              
              if(gene %in% unique(Ahringer_single_gene_ref_table$Gene)){
                
                ATAC_coord <- Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == gene, 2]
                
                check <- onset_check : (i -2)
                
                if(sum(abs(ATAC_compare_frame[ATAC_coord, check])) > 0){
                  
                  ATAC_mut <- 1
                  
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == -1*sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_down <- 1}
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_up <- 1}
                }
              }           
              
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, smallRNA_mut, smallRNA_down, smallRNA_up, ATAC_mut, ATAC_down, ATAC_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              
              smallRNA_mut <- 0
              smallRNA_down <- 0
              smallRNA_up <- 0
              
              ATAC_mut <- 0
              ATAC_down <- 0
              ATAC_up <- 0
              
            }
            
            # deal with a down epimutation extending to end of lineage
            
            if(vector_in[i]==-1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              gene <- strsplit(gene_name, ":")[[1]][4]
              
              number_transitions <- number_transitions + 0
              
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 0
              
              
              
              # is there a  small RNA epimutation that corresponds?
              
              if(gene %in% rownames(smallRNA_binarised_table_A_2.25)){
                
                
                smallRNA_coord <- gene
                
                check <- onset_check : (i -1)
                
                if(sum(abs(smallRNA_compare_frame[smallRNA_coord, check])) > 0){
                  
                  smallRNA_mut <- 1
                  
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) < sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_down <- 1}
                  
                  if(sum(smallRNA_compare_frame[smallRNA_coord, check]) == sum(abs(smallRNA_compare_frame[smallRNA_coord, check]))){
                    smallRNA_up <- 1}
                }
              }
              
              
              # is there an ATAC epimutation that corresponds?
              
              if(gene %in% unique(Ahringer_single_gene_ref_table$Gene)){
                
                ATAC_coord <- Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == gene_name, 2]
                
                check <- onset_check : (i -1)
                
                if(sum(abs(ATAC_compare_frame[ATAC_coord, check])) > 0){
                  
                  ATAC_mut <- 1
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == -1*sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_down <- 1}
                  
                  if(sum(ATAC_compare_frame[ATAC_coord, check]) == sum(abs(ATAC_compare_frame[ATAC_coord, check]))){
                    ATAC_up <- 1}
                }
              }           
              
              
              save_mut <- data.frame(gene_name, gene, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, smallRNA_mut, smallRNA_down, smallRNA_up, ATAC_mut, ATAC_down, ATAC_up)   
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              
              smallRNA_mut <- 0
              smallRNA_down <- 0
              smallRNA_up <- 0
              
              ATAC_mut <- 0
              ATAC_down <- 0
              ATAC_up <- 0
              
            }
            
          }    
          
        }} 
    
    Overall_output <- rbind(Overall_output, output_frame)  
    
  }
  
  RNA_All_output_list_with_missing_as_muts[[e]] <- Overall_output
}

names(RNA_All_output_list_with_missing_as_muts) <- names(All_main_frames)      

#-------------------------------------------------------------------------------------------

# Filter out suspicious epimutations

RNA_All <- rbind(RNA_All_output_list_with_missing_as_muts$A, RNA_All_output_list_with_missing_as_muts$B, 
                 RNA_All_output_list_with_missing_as_muts$C)

All_RNA_info <- paste(RNA_All[,1], RNA_All[,3], RNA_All[,4], sep=":")

RNA_epimutinfo <- table(All_RNA_info)

censor <- which(RNA_epimutinfo==3|RNA_epimutinfo==6|RNA_epimutinfo==9|RNA_epimutinfo==12|RNA_epimutinfo==15)

# selecting only the ones that are the same in all lineages, thus 3 entries

Keep <-which(All_RNA_info%in%names(censor)==F)

# Filtering out those with 3 (or if all zeros), only keeping ones without

RNA_minus_suspect <- RNA_All[Keep,]

# ----------------------------------------------------------------------------------------------------------------------------------------------

# ATAC counting lengths of epimutations

ATAC_Tab_A<-cbind(rep(0, length=nrow(ATAC_binarised_table_A_2.25_ehancer_promoter)),ATAC_binarised_table_A_2.25_ehancer_promoter)
colnames(ATAC_Tab_A)[1]<-"0"



ATAC_Tab_B<-cbind(rep(0, length=nrow(ATAC_binarised_table_B_2.25_ehancer_promoter)),ATAC_binarised_table_B_2.25_ehancer_promoter)
colnames(ATAC_Tab_B)[1]<-"0"



ATAC_Tab_C<-cbind(rep(0, length=nrow(ATAC_binarised_table_C_2.25_ehancer_promoter)), ATAC_binarised_table_C_2.25_ehancer_promoter)
colnames(ATAC_Tab_C)[1]<-"0"




All_main_frames <- list(ATAC_Tab_A, ATAC_Tab_B, ATAC_Tab_C)


names(All_main_frames) <- c("A", "B", "C")




All_output_list_with_missing_as_muts <- list()

for(e in 1:length(All_main_frames)){
  
  main_frame <- as.data.frame(All_main_frames[[e]])
  
  Overall_output <- c()
  
  for(t in 1:nrow(main_frame)){
    
    vector_in <- main_frame[t, ]
    gen_names <- as.numeric(colnames(vector_in))
    
    number_transitions <- 0
    onset_gen <- 0
    tempL <- 0
    output_frame <- c()
    
    is_up <- 0
    is_down <- 0
    complete <- 0
    
    # first deal with no epimutations in the lineage for a locus/gene
    
    if(sum(abs(vector_in)) == 0){
      
      length = 0
      is_up <- 0
      is_down <- 0
      names <- rownames(vector_in)
      number_transitions <- 0 
      onset_gen <- 0
      Lineage <- names(All_main_frames[e])
      complete <- 0
      
      gene_name <- unique(unlist(str_split(Ahringer_ref_table[Ahringer_ref_table$Locus == names, 3], pattern = ",")))
      
      for(t in 1:length(gene_name)){
        
        Ahringer_details <- Ahringer_ref_table[Ahringer_ref_table$Locus == names, c(2, 4, 5)]
        
        save_mut <- data.frame(names, number_transitions, onset_gen, length, complete, Lineage, is_up, is_down, gene_name[[t]], Ahringer_details)   
        
        output_frame <- rbind(output_frame, save_mut)
        
      }
      
    }
    
    
    else
      
      if(sum(abs(vector_in)) > 0){
        
        for(i in 2:length(vector_in)){
          
          
          # if it is an UP transition, turn off any down transitions and start new epimutation
          
          if(vector_in[i]==1&vector_in[i-1]== 0|vector_in[i]== 1&vector_in[i-1]== -1) {
            
            
            if(is_down == 1){ # a down epimutation that is now turning off
              
              names <- rownames(vector_in)
              number_transitions <- number_transitions + 1 
              
              length <-tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <-1
              
              gene_name <- unique(unlist(str_split(Ahringer_ref_table[Ahringer_ref_table$Locus == names, 3], pattern = ",")))
              
              for(t in 1:length(gene_name)){
                
                Ahringer_details <- Ahringer_ref_table[Ahringer_ref_table$Locus == names, c(2, 4, 5)]
                
                save_mut <- data.frame(names, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, gene_name[[t]], Ahringer_details)   
                
                output_frame <- rbind(output_frame, save_mut)
                
              }
              
            }
            
            # start the UP transition
            
            is_up <- 1
            is_down <- 0
            tempL <- 1
            
            number_transitions <- 1 # there is a transition to an UP epimutation at generation[i]
            onset_gen <- gen_names[i]
          }
          
          
          
          # if it is a new DOWN transition, turn off any up transitions and start new epimutation
          
          if(vector_in[i]==-1&vector_in[i-1]== 0|vector_in[i]== -1&vector_in[i-1]== 1) {
            
            if(is_up == 1){ # an up epimutation that is now turning off
              
              names <- rownames(vector_in)
              
              number_transitions <- number_transitions + 1
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              gene_name <- unique(unlist(str_split(Ahringer_ref_table[Ahringer_ref_table$Locus == names, 3], pattern = ",")))
              
              for(t in 1:length(gene_name)){
                
                Ahringer_details <- Ahringer_ref_table[Ahringer_ref_table$Locus == names, c(2, 4, 5)]
                
                save_mut <- data.frame(names, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, gene_name[[t]], Ahringer_details)   
                
                output_frame <- rbind(output_frame, save_mut)
                
              }
              
            }
            
            
            # start the down transition
            
            is_down <- 1
            is_up <- 0
            tempL <- 1
            number_transitions <- 1 # there is a transition to a DOWN epimutation at generation[i]
            onset_gen <- gen_names[i]
          }
          
          
          
          # if it is an UP epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_up==1){
            
            
            if(vector_in[i]==1&vector_in[i-1]==1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}  # the epimutation continues
            
            
            
            if(vector_in[i]==0&vector_in[i-1]==1){
              
              names <- rownames(vector_in)
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from an UP epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete  <- 1
              
              gene_name <- unique(unlist(str_split(Ahringer_ref_table[Ahringer_ref_table$Locus == names, 3], pattern = ",")))
              
              for(t in 1:length(gene_name)){
                
                Ahringer_details <- Ahringer_ref_table[Ahringer_ref_table$Locus == names, c(2, 4, 5)]
                
                save_mut <- data.frame(names, number_transitions, onset_gen, length, complete, Lineage, is_up, is_down, gene_name[[t]], Ahringer_details)   
                
                output_frame <- rbind(output_frame, save_mut)
                
              }
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              length <- 0
              complete <- 0
              
            }
            
            
            
            # deal with an up epimutation extending to end of lineage
            
            if(vector_in[i]==1&i==length(vector_in)){
              
              names <- rownames(vector_in)
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 0
              
              gene_name <- unique(unlist(str_split(Ahringer_ref_table[Ahringer_ref_table$Locus == names, 3], pattern = ",")))
              
              for(t in 1:length(gene_name)){
                
                Ahringer_details <- Ahringer_ref_table[Ahringer_ref_table$Locus == names, c(2, 4, 5)]
                
                save_mut <- data.frame(names, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, gene_name[[t]], Ahringer_details)   
                
                output_frame <- rbind(output_frame, save_mut)
                
              }
              
              is_up <- 0
              is_down <- 0
              length <- 0
              tempL <- 0
              complete <- 0
              
            }
            
          }
          
          
          # if it is a DOWN epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_down==1){
            
            if(vector_in[i]==-1&vector_in[i-1]==-1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}
            
            
            
            if(vector_in[i]==0&vector_in[i-1]==-1){
              
              names <- rownames(vector_in)
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from a DOWN epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              gene_name <- unique(unlist(str_split(Ahringer_ref_table[Ahringer_ref_table$Locus == names, 3], pattern = ",")))
              
              for(t in 1:length(gene_name)){
                
                Ahringer_details <- Ahringer_ref_table[Ahringer_ref_table$Locus == names, c(2, 4, 5)]
                
                save_mut <- data.frame(names, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, gene_name[[t]], Ahringer_details)   
                
                output_frame <- rbind(output_frame, save_mut)
                
              }
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              
              
            }
            
            # deal with a down epimutation extending to end of lineage
            
            if(vector_in[i]==-1&i==length(vector_in)){
              
              names <- rownames(vector_in)
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 0
              
              gene_name <- unique(unlist(str_split(Ahringer_ref_table[Ahringer_ref_table$Locus == names, 3], pattern = ",")))
              
              for(t in 1:length(gene_name)){
                
                Ahringer_details <- Ahringer_ref_table[Ahringer_ref_table$Locus == names, c(2, 4, 5)]
                
                save_mut <- data.frame(names, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, gene_name[[t]], Ahringer_details)   
                
                output_frame <- rbind(output_frame, save_mut)
                
              }
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
            }
            
          }    
          
        }} 
    
    Overall_output <- rbind(Overall_output, output_frame)  
    
  }
  
  All_output_list_with_missing_as_muts[[e]] <- Overall_output
}

names(All_output_list_with_missing_as_muts) <- names(All_main_frames)      



#---------------------------

# Filter out suspicious epimutations

ATAC_ep_All <- rbind(All_output_list_with_missing_as_muts$A, All_output_list_with_missing_as_muts$B, 
                     All_output_list_with_missing_as_muts$C)

All_ATAC_ep_info <- paste(ATAC_ep_All[,1], ATAC_ep_All[,3], ATAC_ep_All[,4], sep=":")

ATAC_ep_epimutinfo <- table(All_ATAC_ep_info)

censor <- which(ATAC_ep_epimutinfo==3|ATAC_ep_epimutinfo==6|ATAC_ep_epimutinfo==9|ATAC_ep_epimutinfo==12|ATAC_ep_epimutinfo==15)

# selecting only the ones that are the same in all lineages, thus 3 entries

Keep <-which(All_ATAC_ep_info%in%names(censor)==F)

# Filtering out those with 3 (or if all zeros), only keeping ones without

ATAC_ep_minus_suspect <- ATAC_ep_All[Keep,]

colnames(ATAC_ep_minus_suspect) <- c("names", "number_transitions", "onset_gen", "length", "complete", "Lineage", "is_up", "is_down", "Gene_mapped", "Enhancer_or_Promoter", "Chromatin_domain", "Tissue_specificity")

# ---------------------------------------------------------------------------------

# Small RNA : counting lengths of epimutations

small_RNA_Tab_A<-cbind(rep(0, length=nrow(smallRNA_binarised_table_A_2.25)),smallRNA_binarised_table_A_2.25)
colnames(small_RNA_Tab_A)[1]<-"0"
small_RNA_Tab_A <- as.data.frame(small_RNA_Tab_A)


small_RNA_Tab_B<-cbind(rep(0, length=nrow(smallRNA_binarised_table_B_2.25)),smallRNA_binarised_table_B_2.25)
colnames(small_RNA_Tab_B)[1]<-"0"
small_RNA_Tab_B <- as.data.frame(small_RNA_Tab_B)


small_RNA_Tab_C<-cbind(rep(0, length=nrow(smallRNA_binarised_table_C_2.25)),smallRNA_binarised_table_C_2.25)
colnames(small_RNA_Tab_C)[1]<-"0"
small_RNA_Tab_C <- as.data.frame(small_RNA_Tab_C)



All_main_frames <- list(small_RNA_Tab_A, small_RNA_Tab_B, small_RNA_Tab_C)

compare_RNA <- list(RNA_binarised_table_A_2.25_no_21URs, RNA_binarised_table_B_2.25_no_21URs, RNA_binarised_table_C_2.25_no_21URs)


names(All_main_frames) <- c("A", "B", "C")




smallRNA_All_output_list_with_missing_as_muts <- list()

for(e in 1:length(All_main_frames)){
  
  main_frame <- as.data.frame(All_main_frames[[e]])
  
  compare_frame <- compare_RNA[[e]]
  
  Overall_output <- c()
  
  for(t in 1:nrow(main_frame)){
    
    vector_in <- main_frame[t, ]
    gen_names <- as.numeric(colnames(vector_in))
    
    number_transitions <- 0
    onset_gen <- 0
    tempL <- 0
    output_frame <- c()
    
    is_up <- 0
    is_down <- 0
    complete <- 0
    
    RNA_mut <- 0
    RNA_up <- 0
    RNA_down <- 0
    
    # first deal with no epimutations in the lineage for a locus/gene
    
    if(sum(abs(vector_in)) == 0){
      
      length = 0
      is_up <- 0
      is_down <- 0
      gene_name <- rownames(vector_in)
      number_transitions <- 0 
      onset_gen <- 0
      RNA_mut <- 0
      RNA_up <- 0
      RNA_down <- 0
      
      Lineage <- names(All_main_frames[e])
      complete <- 0
      
      
      save_mut <- data.frame(gene_name, number_transitions, onset_gen, length,  complete,  Lineage, is_up, is_down, RNA_mut, RNA_down, RNA_up)   
      
      output_frame <- rbind(output_frame, save_mut)
      
    }
    
    
    else
      
      if(sum(abs(vector_in)) > 0){
        
        for(i in 2:length(vector_in)){
          
          
          # if it is an UP transition, turn off any down transitions and start new epimutation
          
          if(vector_in[i]==1&vector_in[i-1]== 0|vector_in[i]== 1&vector_in[i-1]== -1) {
            
            
            if(is_down == 1){ # a down epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              number_transitions <- number_transitions + 1 
              
              length <-tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              # is there an RNA epimutation that corresponds?
              
              RNA_coord <- unique(All_lins_integrated_table[All_lins_integrated_table$gene == gene_name, 3])
              
              check <- onset_check : (i -2)
              
              if(sum(abs(compare_frame[RNA_coord, check])) > 0){
                
                RNA_mut <- 1
                
                
                if(sum(compare_frame[RNA_coord, check]) < sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_down <- 1}
                
                if(sum(compare_frame[RNA_coord, check]) == sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_up <- 1}
              }
              
              save_mut <- data.frame(gene_name, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, RNA_mut, RNA_down, RNA_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            # start the UP transition
            
            is_up <- 1
            is_down <- 0
            tempL <- 1
            
            number_transitions <- 1 # there is a transition to an UP epimutation at generation[i]
            onset_gen <- gen_names[i]
            
            onset_check <- i -1
          }
          
          
          
          # if it is a new DOWN transition, turn off any up transitions and start new epimutation
          
          if(vector_in[i]==-1&vector_in[i-1]== 0|vector_in[i]== -1&vector_in[i-1]== 1) {
            
            if(is_up == 1){ # an up epimutation that is now turning off
              
              gene_name <- rownames(vector_in)
              
              number_transitions <- number_transitions + 1
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              # is there an RNA epimutation that corresponds?
              
              RNA_coord <- unique(All_lins_integrated_table[All_lins_integrated_table$gene == gene_name, 3])
              
              check <- onset_check : (i -2)
              
              if(sum(abs(compare_frame[RNA_coord, check])) > 0){
                
                RNA_mut <- 1
                
                
                if(sum(compare_frame[RNA_coord, check]) < sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_down <- 1}
                
                if(sum(compare_frame[RNA_coord, check]) == sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_up <- 1}
              }
              
              save_mut <- data.frame(gene_name, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, RNA_mut, RNA_down, RNA_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
            }
            
            
            # start the down transition
            
            is_down <- 1
            is_up <- 0
            tempL <- 1
            number_transitions <- 1 # there is a transition to a DOWN epimutation at generation[i]
            onset_gen <- gen_names[i]
            onset_check <- i-1
          }
          
          
          
          # if it is an UP epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_up==1){
            
            
            if(vector_in[i]==1&vector_in[i-1]==1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}  # the epimutation continues
            
            
            
            if(vector_in[i]==0&vector_in[i-1]==1){
              
              gene_name <- rownames(vector_in)
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from an UP epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              # is there an RNA epimutation that corresponds?
              
              RNA_coord <- unique(All_lins_integrated_table[All_lins_integrated_table$gene == gene_name, 3])
              
              check <- onset_check : (i -2)
              
              if(sum(abs(compare_frame[RNA_coord, check])) > 0){
                
                RNA_mut <- 1
                
                
                if(sum(compare_frame[RNA_coord, check]) < sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_down <- 1}
                
                if(sum(compare_frame[RNA_coord, check]) == sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_up <- 1}
                
              }
              save_mut <- data.frame(gene_name, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, RNA_mut, RNA_down, RNA_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              length <- 0
              complete <- 0
              RNA_mut <- 0
              RNA_down <- 0
              RNA_up <- 0
            }
            
            
            
            # deal with an up epimutation extending to end of lineage
            
            if(vector_in[i]==1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              number_transitions <- number_transitions + 0
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 0
              
              # is there an RNA epimutation that corresponds?
              
              RNA_coord <- unique(All_lins_integrated_table[All_lins_integrated_table$gene == gene_name, 3])
              
              check <- onset_check : (i-1)
              
              if(sum(abs(compare_frame[RNA_coord, check])) > 0){
                
                RNA_mut <- 1
                
                
                if(sum(compare_frame[RNA_coord, check]) < sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_down <- 1}
                
                if(sum(compare_frame[RNA_coord, check]) == sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_up <- 1}
              }
              
              save_mut <- data.frame(gene_name, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, RNA_mut, RNA_down, RNA_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              length <- 0
              tempL <- 0
              complete <- 0
              RNA_mut <- 0
              RNA_down <- 0
              RNA_up <- 0
              
            }
            
          }
          
          
          # if it is a DOWN epimutation, does it continue, does it turn off, does it extend to last gen?
          
          if(is_down==1){
            
            if(vector_in[i]==-1&vector_in[i-1]==-1)
              
            {tempL<-tempL+(gen_names[i]-gen_names[i-1])}
            
            
            if(vector_in[i]==0&vector_in[i-1]==-1){
              
              gene_name <- rownames(vector_in)
              
              number_transitions <- number_transitions +1  # there is a transition to OFF from a DOWN epimutation at generation[i]
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 1
              
              # is there an RNA epimutation that corresponds?
              
              RNA_coord <- unique(All_lins_integrated_table[All_lins_integrated_table$gene == gene_name, 3])
              
              check <- onset_check : (i -2)
              
              if(sum(abs(compare_frame[RNA_coord, check])) > 0){
                
                RNA_mut <- 1
                
                
                if(sum(compare_frame[RNA_coord, check]) < sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_down <- 1}
                
                if(sum(compare_frame[RNA_coord, check]) == sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_up <- 1}
              }
              
              save_mut <- data.frame(gene_name, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, RNA_mut, RNA_down, RNA_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              RNA_mut <- 0
              RNA_down <- 0
              RNA_up <- 0
              
              
            }
            
            # deal with a down epimutation extending to end of lineage
            
            if(vector_in[i]==-1&i==length(vector_in)){
              
              gene_name <- rownames(vector_in)
              
              number_transitions <- number_transitions + 0
              
              
              length <- tempL
              
              Lineage <- names(All_main_frames[e])
              
              complete <- 0
              
              # is there an RNA epimutation that corresponds?
              
              RNA_coord <- unique(All_lins_integrated_table[All_lins_integrated_table$gene == gene_name, 3])
              
              check <- onset_check : (i-1)
              
              if(sum(abs(compare_frame[RNA_coord, check])) > 0){
                
                RNA_mut <- 1
                
                
                if(sum(compare_frame[RNA_coord, check]) < sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_down <- 1}
                
                if(sum(compare_frame[RNA_coord, check]) == sum(abs(compare_frame[RNA_coord, check]))){
                  RNA_up <- 1}
                
              }
              save_mut <- data.frame(gene_name, number_transitions, onset_gen, length,  complete, Lineage, is_up, is_down, RNA_mut, RNA_down, RNA_up)   
              
              output_frame <- rbind(output_frame, save_mut)
              
              is_up <- 0
              is_down <- 0
              tempL <- 0
              complete <- 0
              RNA_mut <- 0
              RNA_down <- 0
              RNA_up <- 0
              
            }
            
          }    
          
        }} 
    
    Overall_output <- rbind(Overall_output, output_frame)  
    
  }
  
  smallRNA_All_output_list_with_missing_as_muts[[e]] <- Overall_output
}

names(smallRNA_All_output_list_with_missing_as_muts) <- names(All_main_frames)      



#-----------------------------------

# Filter out suspicious epimutations

small_RNA_All <- rbind(smallRNA_All_output_list_with_missing_as_muts$A, smallRNA_All_output_list_with_missing_as_muts$B, 
                       smallRNA_All_output_list_with_missing_as_muts$C)




All_RNA_info <- paste(small_RNA_All[,1], small_RNA_All[,3], small_RNA_All[,4], sep=":")

RNA_epimutinfo <- table(All_RNA_info)

censor <- which(RNA_epimutinfo==3|RNA_epimutinfo==6|RNA_epimutinfo==9|RNA_epimutinfo==12|RNA_epimutinfo==15)

# selecting only the ones that are the same in all lineages, thus 3 entries

Keep <-which(All_RNA_info%in%names(censor)==F)

# Filtering out those with 3 (or if all zeros), only keeping ones without

small_RNA_minus_suspect <- small_RNA_All[Keep,]

# ----------------------------------------------------------------------------------------------------------------------

# Survival of ATAC, RNA and small RNA separately

# Load required packages
library(survival)
library(dplyr)
library(survminer)
library(ggplot2)


#----------------------------------------------------------------------------------------------------------

# Figure 2 B
# combine all data sets and compare

# **Restrict to true epimutations only, i.e. length is > 1**

RNA_minus_suspect_mut <- RNA_minus_suspect[RNA_minus_suspect$length>1, ]


ATAC_minus_suspect_mut <- ATAC_ep_minus_suspect[ATAC_ep_minus_suspect$length>1, ]


smallRNA_minus_suspect_mut <- small_RNA_minus_suspect[small_RNA_minus_suspect$length>1, ]



#------------------------


colnames(smallRNA_minus_suspect_mut)[1] <- "names"

smallRNA_minus_suspect_mut <- cbind(rep("small_RNA", length = nrow(smallRNA_minus_suspect_mut)), smallRNA_minus_suspect_mut)
colnames(smallRNA_minus_suspect_mut)[1] <- "Data_type"


colnames(RNA_minus_suspect_mut)[1] <- "names"

RNA_minus_suspect_mut <- cbind(rep("RNA", length = nrow(RNA_minus_suspect_mut)), RNA_minus_suspect_mut)
colnames(RNA_minus_suspect_mut)[1] <- "Data_type"



ATAC_minus_suspect_mut <- cbind(rep("ATAC", length = nrow(ATAC_minus_suspect_mut)), ATAC_minus_suspect_mut)
colnames(ATAC_minus_suspect_mut)[1] <- "Data_type"




ALL_SETS_minus_suspect_mut <- rbind(RNA_minus_suspect_mut[, c(1, 2, 4:10)], ATAC_minus_suspect_mut[, 1:9], smallRNA_minus_suspect_mut[, 1:9])


All_SETS_surv_object_ABC <- Surv(time = ALL_SETS_minus_suspect_mut$length, event = ALL_SETS_minus_suspect_mut$complete)



ALL_SETS_ABC_fit <- survfit(All_SETS_surv_object_ABC ~ 1, data = ALL_SETS_minus_suspect_mut)

ALL_SETS_ABC_plot <- ggsurvplot(ALL_SETS_ABC_fit, ALL_SETS_minus_suspect_mut, censor = T)+
  ggtitle("All data compare")
xlab("Time (generations)")

ALL_SETS_ABC_plot$plot <- ALL_SETS_ABC_plot$plot + theme(legend.text = element_text(size = 14))


ALL_SETS_ABC_fit_data_type <-  survfit(All_SETS_surv_object_ABC ~ Data_type, data = ALL_SETS_minus_suspect_mut)


ALL_SETS_ABC_plot_data_type <- ggsurvplot(ALL_SETS_ABC_fit_data_type,  ALL_SETS_minus_suspect_mut, censor = T, break.time.by= 1, pval = T,  pval.coord = c(4, 1), 
                                          font.main = c(16, "bold"),
                                          font.x = 14,
                                          font.y = 14,
                                          font.legend = 14,
                                          font.tickslab = 12)+
  # surv.median.line = "hv")+
  ggtitle("Figure 2 B")+
  xlab("Time (generations)")

ALL_SETS_ABC_plot_data_type$plot +theme(plot.title = element_text(hjust = 0.5))


# n.b. p value derived in ggsurv plot through log rank test


#-----------------------------------------------------------------------------------------------------------
# FIGURE 3: OVERLAPS

# Visualising proportions of genes with RNA expression changes which have ATAC epimutations 


# Figure 3 A: Schematic made in BioRender

# Figure 3 B

# RNA with ATAC


All_genes_RNA_muts_inherited_chrom_muts <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$is_ep_epimut==1, ])


RNA_inherit_chrom_mut_tm_inherit <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_ep_inherited==1&All_lins_integrated_table$has_time_matched==1, ])/All_genes_RNA_muts_inherited_chrom_muts*100

RNA_inherit_chrom_mut_tm_non_inherit <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_ep_inherited==0&All_lins_integrated_table$has_time_matched==1, ])/All_genes_RNA_muts_inherited_chrom_muts*100

RNA_inherit_chrom_mut_non_tm <-  nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$has_time_matched==0, ])/All_genes_RNA_muts_inherited_chrom_muts*100




RNA_Inherit_proportions <-  data.frame(c(RNA_inherit_chrom_mut_non_tm, 
                                         RNA_inherit_chrom_mut_tm_non_inherit,
                                         RNA_inherit_chrom_mut_tm_inherit))


Names <- c(
  "Has chromatin epimutations - not time matched", 
  "Has time matched chromatin epimutations - not inherited", 
  "Has time matched chromatin epimutations - inherited")                                           


Species <- rep("Genes with inherited \nRNA expression changes\n& chromatin epimutations", 3)

RNA_Inherit <- cbind(RNA_Inherit_proportions, Names, Species)


colnames(RNA_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")



# ------------

All_genes_RNA_muts_n_inherited_chrom_muts <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==0&All_lins_integrated_table$is_ep_epimut==1, ])


RNA_n_inherit_chrom_mut_tm_inherit <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==0&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_ep_inherited==1&All_lins_integrated_table$has_time_matched==1, ])/All_genes_RNA_muts_n_inherited_chrom_muts*100

RNA_n_inherit_chrom_mut_tm_non_inherit <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==0&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_ep_inherited==0&All_lins_integrated_table$has_time_matched==1, ])/All_genes_RNA_muts_n_inherited_chrom_muts*100

RNA_n_inherit_chrom_mut_non_tm <-  nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==0&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$has_time_matched==0, ])/All_genes_RNA_muts_n_inherited_chrom_muts*100




RNA_n_Inherit_proportions <-  data.frame(c(RNA_n_inherit_chrom_mut_non_tm, 
                                           RNA_n_inherit_chrom_mut_tm_non_inherit,
                                           RNA_n_inherit_chrom_mut_tm_inherit))


Names <- c(
  "Has chromatin epimutations - not time matched", 
  "Has time matched chromatin epimutations - not inherited", 
  "Has time matched chromatin epimutations - inherited")                                           


Species <- rep("Genes with non inherited \nRNA expression changes\n& chromatin epimutations", 3)

RNA_N_Inherit <- cbind(RNA_n_Inherit_proportions, Names, Species)


colnames(RNA_N_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#-------------

RNA_Inherited_features <- rbind(RNA_N_Inherit, RNA_Inherit)


RNA_Inherited_features <- RNA_Inherited_features[order(nrow(RNA_Inherited_features):1), ] 

RNA_Inherited_features$Epigenetic_feature_of_gene <- factor(RNA_Inherited_features$Epigenetic_feature_of_gene, levels = unique(RNA_Inherited_features$Epigenetic_feature_of_gene))


positions <- unique(RNA_Inherited_features$Species)



# Figure 3 B

RNA_w_ATAC_proportion_plot <-
  
  ggplot(RNA_Inherited_features, aes(fill=Epigenetic_feature_of_gene, y=Percentages, x=Species)) + 
  
  geom_bar(position="stack", stat="identity", width = 0.4, color="black")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c("deeppink", 
                             "deepskyblue", 
                             "aquamarine",
                             "deeppink", 
                             "deepskyblue", 
                             "aquamarine"))+
  theme_bw()+
  labs(x= " ", 
       y = "Percentage of genes (%)\n")+
  theme(plot.margin = unit(c(2,1,2,2), "cm"))+
  theme(plot.title = element_text(size = 16))+
  labs(title = "Figure 3 B")+
  theme(text = element_text(size=14))+
  theme(axis.text.x=element_text(colour="black", size = 16))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name = "RNA expression change status",
                    values=c("deeppink",
                             "deepskyblue",
                             "aquamarine"),
                    labels=c("Simultaneous inherited chromatin epimutation", 
                             "Simultaneous non inherited chromatin epimutation", 
                             "Non simultaneous chromatin epimutation"))

RNA_w_ATAC_proportion_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Co-existing chromatin epimutation status"))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))



#-----------------------------------------------

# RNA with small RNA

All_genes_RNA_muts_inherited_smallRNA_muts <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1, ])

RNA_inherit_smallRNA_mut_tm_inherit <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_inherited==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1, ])/All_genes_RNA_muts_inherited_smallRNA_muts*100

RNA_inherit_smallRNA_mut_tm_non_inherit <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_inherited==0&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1, ])/All_genes_RNA_muts_inherited_smallRNA_muts*100

RNA_inherit_smallRNA_mut_non_tm <-  nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==0, ])/All_genes_RNA_muts_inherited_smallRNA_muts*100


RNA_Inherit_proportions <-  data.frame(c(RNA_inherit_smallRNA_mut_non_tm, 
                                         RNA_inherit_smallRNA_mut_tm_non_inherit,
                                         RNA_inherit_smallRNA_mut_tm_inherit))
Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                           


Species <- rep("Genes with inherited \nRNA expression changes\ntargeted by small RNA epimutations", 3)

RNA_Inherit <- cbind(RNA_Inherit_proportions, Names, Species)


colnames(RNA_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")


# ------------

# Non Inherited RNA seq changes 

All_genes_RNA_muts_n_inherited_smallRNA_muts <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==0&smallRNA_All_lins_integrated_table$smallRNA_mut==1, ])

RNA_n_inherit_smallRNA_mut_tm_inherit <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==0&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_inherited==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1, ])/All_genes_RNA_muts_n_inherited_smallRNA_muts*100

RNA_n_inherit_smallRNA_mut_tm_non_inherit <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==0&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_inherited==0&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1, ])/All_genes_RNA_muts_n_inherited_smallRNA_muts*100

RNA_n_inherit_smallRNA_mut_non_tm <-  nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==0&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==0, ])/All_genes_RNA_muts_n_inherited_smallRNA_muts*100




RNA_n_Inherit_proportions <-  data.frame(c(RNA_n_inherit_smallRNA_mut_non_tm, 
                                           RNA_n_inherit_smallRNA_mut_tm_non_inherit,
                                           RNA_n_inherit_smallRNA_mut_tm_inherit))


Names <- c(
  "Has small RNA epimutations - not time matched", 
  "Has time matched small RNA epimutations - not inherited", 
  "Has time matched small RNA epimutations - inherited")                                             


Species <- rep("Genes with non inherited \nRNA expression changes\ntargeted by small RNA epimutations", 3)

RNA_N_Inherit <- cbind(RNA_n_Inherit_proportions, Names, Species)


colnames(RNA_N_Inherit) <- c("Percentages", "Epigenetic_feature_of_gene", "Species")

#-------------
# Figure 3C

# Plot for RNA seq changes 

RNA_Inherited_features <- rbind(RNA_N_Inherit, RNA_Inherit)


RNA_Inherited_features <- RNA_Inherited_features[order(nrow(RNA_Inherited_features):1), ] 

RNA_Inherited_features$Epigenetic_feature_of_gene <- factor(RNA_Inherited_features$Epigenetic_feature_of_gene, levels = unique(RNA_Inherited_features$Epigenetic_feature_of_gene))


positions <- unique(RNA_Inherited_features$Species)



# Figure 3 C  

# Proportions of genes with inherited and non inherited RNAseq changes that have small RNA changes

RNA_w_smallRNA_proportion_plot <-
  
  ggplot(RNA_Inherited_features, aes(fill=Epigenetic_feature_of_gene, y=Percentages, x=Species)) + 
  
  geom_bar(position="stack", stat="identity", width = 0.4, color="black")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral",
                             "cornflowerblue", 
                             "cornsilk3", 
                             "coral"))+
  theme_bw()+
  labs(x= " ", 
       y = "Percentage of genes (%)\n")+
  theme(plot.margin = unit(c(2,1,2,2), "cm"))+
  theme(plot.title = element_text(size = 16))+
  labs(title = "Figure 3 C")+
  theme(text = element_text(size=14))+
  theme(axis.text.x=element_text(colour="black", size = 12))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name = "RNA expression change status",
                    values=c("cornflowerblue", 
                             "cornsilk3", 
                             "coral"),
                    labels=c("Simultaneous inherited small RNA epimutation", 
                             "Simultaneous non inherited small RNA epimutation", 
                             "Non simultaneous small RNA epimutation"))

RNA_w_smallRNA_proportion_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill=guide_legend(title="Associated small RNA epimutation status"))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))


#----------------------------------------------------------------------------------------------
# Supplementary Figure 4

# Supplementary figure showing the percentage overlap of gene expression/ATAC change for piRNA regions relative to rest of genome].  


# Chromatin epimutations

# Background is all Regulated loci genome wide


Reg_loci <- nrow(All_lins_ATAC_table[All_lins_ATAC_table$Chromatin_domain == "R", ])

# Of which are epimutated 

Epimut_Reg <- nrow(All_lins_ATAC_table[All_lins_ATAC_table$Chromatin_domain == "R" & All_lins_ATAC_table$Epimutated==1, ])


# Of which are piRNA cluster genes 

piRNA_yes <- nrow(All_lins_ATAC_table[All_lins_ATAC_table$Chromatin_domain == "R" & All_lins_ATAC_table$piRNA_cluster==1, ])

# Of which are both


piRNA_yes_epimut <- nrow(All_lins_ATAC_table[All_lins_ATAC_table$Chromatin_domain == "R" & All_lins_ATAC_table$piRNA_cluster==1 & All_lins_ATAC_table$Epimutated==1, ])


contingency_table <- rbind(c(piRNA_yes_epimut, (piRNA_yes - piRNA_yes_epimut)), 
                           c((Epimut_Reg - piRNA_yes_epimut), Reg_loci - (piRNA_yes + Epimut_Reg - piRNA_yes_epimut)))


output <- fisher.test(contingency_table)


#---------------------------



# Of which are NOT piRNA cluster genes 

piRNA_no <- nrow(All_lins_ATAC_table[All_lins_ATAC_table$Chromatin_domain == "R" & All_lins_ATAC_table$piRNA_cluster==0, ])


# Of which are both: not piRNA cluster genes and epimutated

piRNA_no_epimut <- nrow(All_lins_ATAC_table[All_lins_ATAC_table$Chromatin_domain == "R" & All_lins_ATAC_table$piRNA_cluster==0 & All_lins_ATAC_table$Epimutated==1, ])


contingency_table <- rbind(c(piRNA_no_epimut, (piRNA_no - piRNA_no_epimut)), 
                           c((Epimut_Reg - piRNA_no_epimut), Reg_loci - (piRNA_no + Epimut_Reg - piRNA_no_epimut)))


output <- fisher.test(contingency_table)




# Gene expression changes overlap with piRNA regions from all lineages

# Number of genes with RNA expression changes in piRNA cluster

RNA_piRNA <- length(All_lins_integrated_table[All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$piRNA_cluster==1, 2])

# Number of genes with RNA expression changes genome wide   

RNA_genome <- length(All_lins_integrated_table[All_lins_integrated_table$is_RNA_exp_change==1, 2])

# % of genes with RNA expression changes which are in piRNA cluster 

RNA_percent_piRNA <- RNA_piRNA/RNA_genome*100




# ATAC changes overlap with piRNA regions from all lineages

# Number of loci with epimutations in piRNA cluster

ATAC_piRNA <- length(All_lins_ATAC_table[All_lins_ATAC_table$Epimutated==1 & All_lins_ATAC_table$piRNA_cluster==1, 3])

# Number of loci with epimutations genome wide
ATAC_genome <- length(All_lins_ATAC_table[All_lins_ATAC_table$Epimutated==1, 3])

# % of loci with ATAC changes which are in piRNA cluster 

ATAC_percent_piRNA <- ATAC_piRNA/ATAC_genome*100




# small RNA changes overlap with piRNA regions from all lineages

# Number of genes with small RNA epimutations in piRNA cluster

smallRNA_piRNA <- length(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$smallRNA_mut==1 & smallRNA_All_lins_integrated_table$piRNA_cluster==1, 2])

# Number of genes with epimutations genome wide

smallRNA_genome <- length(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$smallRNA_mut==1, 2])


# % of genes with ATAC changes which are in piRNA cluster 

smallRNA_percent_piRNA <- smallRNA_piRNA/smallRNA_genome*100




# Supplementary Figure 4


data_overlap_piRNA <- c(RNA_percent_piRNA, ATAC_percent_piRNA, smallRNA_percent_piRNA)

names <- c("RNA", 
           "ATAC", 
           "small_RNA")

df <- data.frame(data_overlap_piRNA, names)

df$names <- factor(df$names, levels = names )

ggplot(data=df, aes(x=names, y= data_overlap_piRNA)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle(paste("Supplementary Figure 4"))




#-----------------------------------------------------------------------------------------------------------------
# Supplementary Figure 6

# Revision for PLOS genetics 

# Their previous publication (Beltran 2020) states that canonical silencing pathways (including piRNA, HRDE-1, WAGO-1, ERGO-1, and NRDE-3) 
# were enriched for 22G-RNA epimutations but that the canonical licensing pathway CSR-1 is not. If this remains true, 
# it should be briefly discussed here as well as it is important for interpreting the functional relevance of concordant and discordant 
# changes in gene expression vs 22G-RNA epimutations.



# The next step is to look at how the different types of small RNA pathways affect epimutation durations 

# We know our long-lived, short-lived and non-inherited 22G-RNAs

# these gene lists are unique 



list_target_genes <- rownames(smallRNA_binarised_table_A_2.25)

ABC_small_RNA_super_table <- list(smallRNA_binarised_table_A_2.25,
                                  smallRNA_binarised_table_B_2.25, 
                                  smallRNA_binarised_table_C_2.25)

names(ABC_small_RNA_super_table) <- c("A", "B", "C")




gene_features_subset <- Gene_Features[, c(1, 8:14)]

piRNA_genes <- unique(gene_features_subset[gene_features_subset$piRNA_targets==TRUE, 1])

csr1_genes <- unique(gene_features_subset[gene_features_subset$csr1_targets==TRUE,  1])

hrde1_genes <- unique(gene_features_subset[gene_features_subset$hrde1_targets ==TRUE,  1])

wago1_genes <- unique(gene_features_subset[gene_features_subset$wago1_targets ==TRUE,  1])

nrde3_genes <- unique(gene_features_subset[gene_features_subset$nrde3_targets ==TRUE,  1])

ergo1_genes <- unique(gene_features_subset[gene_features_subset$ergo1_targets ==TRUE,  1])

alg34_genes <- unique(gene_features_subset[gene_features_subset$alg34_targets ==TRUE,  1])



all_small_RNA_annotation_table <- c()

for(t in 1:length(names(ABC_small_RNA_super_table))){
  
  Lineage <- names(ABC_small_RNA_super_table)[[t]]
  
  Lineage_table <- c()
  
  for(i in 1:length(list_target_genes)){
    
    Target_gene <- list_target_genes[[i]]
    
    
    # What chromatin domain is the potential target gene in?
    
    Gene_domain <- 0
    maps_to_Ahr <- 0
    
    if(Target_gene %in% Ahringer_single_gene_ref_table$Gene){
      
      maps_to_Ahr <- 1
      
      Gene_domain <- unique(Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == Target_gene, 4])
      
      if(length(Gene_domain)> 1){ 
        
        Gene_domain <- paste(Gene_domain, collapse = "_")
      } 
    }
    
    
    
    
    # What tissue specificity is pertinent to the potential target gene?   
    
    Tissue_specificity <- 0
    
    if(maps_to_Ahr == 1){
      
      Tissue_specificity <- unique(Ahringer_single_gene_ref_table[Ahringer_single_gene_ref_table$Gene == Target_gene, 5])
      
      if(length(Tissue_specificity)> 1){
        
        Tissue_specificity <- paste(Tissue_specificity, collapse = "_")
      } 
    }
    
    
    
    # Is the potential target gene a target for piRNA triggered 22G RNAs?   
    
    piRNA_target <- 0
    
    if(Target_gene %in% piRNA_genes){
      
      piRNA_target <- 1
    } 
    
    
    # Is the potential target gene a target for HRDE-1 triggered 22G RNAs?   
    
    hrde1_target <- 0
    
    if(Target_gene %in% hrde1_genes){
      
      hrde1_target <- 1
    } 
    
    
    # Is the potential target gene a target for WAGO-1 triggered 22G RNAs?   
    
    wago1_target <- 0
    
    if(Target_gene %in% wago1_genes){
      
      wago1_target <- 1
    } 
    
    
    # Is the potential target gene a target for NRDE-3 triggered 22G RNAs?   
    
    nrde3_target <- 0
    
    if(Target_gene %in% nrde3_genes){
      
      nrde3_target <- 1
    } 
    
    
    # Is the potential target gene a target for ERGO-1 triggered 22G RNAs?   
    
    ergo1_target <- 0
    
    if(Target_gene %in% ergo1_genes){
      
      ergo1_target <- 1
    } 
    
    
    
    # Is the potential target gene a target for ALG_3/4 triggered 22G RNAs?   
    
    alg34_target <- 0
    
    if(Target_gene %in% alg34_genes){
      
      alg34_target <- 1
    } 
    
    
    # Is the potential target gene a target for CSR-1 triggered 22G RNAs?   
    
    csr1_target <- 0
    
    if(Target_gene %in% csr1_genes){
      
      csr1_target <- 1
    } 
    
    
    
    
    # Is the potential target gene within the piRNA cluster?
    
    piRNA_cluster <- 0
    
    if(Target_gene %in% piRNA_cluster_genes){
      piRNA_cluster <- 1
    }
    
    
    
    # For the specific lineage, check the epimutation status of the gene 
    
    
    frame <- ABC_small_RNA_super_table[Lineage][[1]]
    
    epimutated <- 0
    
    if(Target_gene %in% rownames(frame)){
      
      epimutated <- 0
      
      if(sum(abs(frame[Target_gene, ])) > 0){
        
        epimutated <- 1}
    }
    
    
    # determine whether the gene is targeted at all by any small RNA pathway 
    
    small_RNA_targeted <- 0
    
    if(sum(piRNA_target, hrde1_target, wago1_target, nrde3_target, ergo1_target, alg34_target, csr1_target) > 0){
      
      small_RNA_targeted <- 1
    }
    
    
    
    
    # if epimutated, determine the duration
    
    long_lived <- 0
    short_lived <- 0
    non_inherited <- 0
    
    
    if(epimutated == 1){
      
      if(Target_gene %in% small_RNA_Top_genes){
        long_lived <- 1
        short_lived <- 0
        non_inherited <- 0
      }
      
      if(Target_gene %in% small_RNA_Muts_Inherited_not_long){
        long_lived <- 0
        short_lived <- 1
        non_inherited <- 0
      }
      
      
      if(Target_gene %in% small_RNA_Muts_Not_Inherited){
        long_lived <- 0
        short_lived <- 0
        non_inherited <- 1
      }
      
    }
    
    row <- data.frame(Lineage, Target_gene,  maps_to_Ahr, epimutated, long_lived, short_lived, non_inherited, small_RNA_targeted, Gene_domain, Tissue_specificity, piRNA_target, hrde1_target, wago1_target, nrde3_target, ergo1_target, alg34_target, csr1_target, piRNA_cluster)
    
    Lineage_table <- rbind(Lineage_table, row)
  }
  
  all_small_RNA_annotation_table <- rbind(all_small_RNA_annotation_table, Lineage_table)
}



# Now we have to determine the relative enrichments of 22G-RNA based epimutations according to the small RNA pathway targets


# We don't care at this stage whether it is Ahringer mapped

# but we should take as our background, all genes which are targeted by at least 1 small RNA pathway (so exlcude genes which are not targeted by any small RNA pathways)

# of all the possible target genes 
# how many are epimutated
# how many are targeted by specific pathway
# how many are both


Background_targeted <- all_small_RNA_annotation_table[all_small_RNA_annotation_table$small_RNA_targeted==1, ]

background_genes <- Background_targeted$Target_gene

epimut_genes <- Background_targeted[Background_targeted$epimutated==1, 2]

Target_Genes <- list(Background_targeted[Background_targeted$piRNA_target==1, 2], 
                     Background_targeted[Background_targeted$hrde1_target==1, 2], 
                     Background_targeted[Background_targeted$wago1_target==1, 2], 
                     Background_targeted[Background_targeted$nrde3_target==1, 2], 
                     Background_targeted[Background_targeted$ergo1_target==1, 2], 
                     Background_targeted[Background_targeted$alg34_target==1, 2], 
                     Background_targeted[Background_targeted$csr1_target==1, 2])

small_RNA_pathways <- c("piRNA_target", "hrde1_target", "wago1_target", "nrde3_target", "ergo1_target", "alg34_target", "csr1_target")


small_RNA_pathway_frame <-c()


for(i in 1:length(Target_Genes)){
  
  target_genes <- Target_Genes[[i]]
  pathway <- small_RNA_pathways[[i]]
  
  j <- which(colnames(Background_targeted)== pathway)
  
  epimut_target <- length(Background_targeted[Background_targeted$epimutated ==1 & Background_targeted[, j]==1, 2])
  epimut_not_target <- length(epimut_genes) - epimut_target
  target_not_epimut <- length(target_genes) - epimut_target
  not_target_not_epimut <- length(background_genes) - (length(epimut_genes) + length(target_genes) - epimut_target)
  contingency_table <-
    rbind(c(epimut_target, epimut_not_target),
          c(target_not_epimut, not_target_not_epimut))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OddsRatio <- FT_out$estim
  output <- data.frame(pathway, OddsRatio, pval)
  
  small_RNA_pathway_frame <-rbind(small_RNA_pathway_frame, output)
  
}




# The next step is to look at how the different types of small RNA pathways affect epimutation duration

# For this, the background is all epimutated genes which are small RNA targeted

Background_targeted_epimut <- Background_targeted[Background_targeted$epimutated==1,]

# Make a matrix to write the results into 

small_RNA_OR_matrix <- matrix(data = NA, nrow = 7, ncol = 3)
small_RNA_pval_matrix <- matrix(data = NA, nrow = 7, ncol = 3)

colnames(small_RNA_OR_matrix) <- c("non_inherited", "short_lived", "long_lived")
rownames(small_RNA_OR_matrix) <- c("piRNA_target", "hrde1_target", "wago1_target", "nrde3_target", "ergo1_target", "alg34_target", "csr1_target")


colnames(small_RNA_pval_matrix) <- c("non_inherited", "short_lived", "long_lived")
rownames(small_RNA_pval_matrix) <- c("piRNA_target", "hrde1_target", "wago1_target", "nrde3_target", "ergo1_target", "alg34_target", "csr1_target")


duration_results <- list()

for(i in 1:length(colnames(small_RNA_OR_matrix))){
  
  duration <-  colnames(small_RNA_OR_matrix)[[i]]
  
  pos_dur <- which(colnames(Background_targeted_epimut) == duration)
  
  
  pathway_results <- c()
  
  
  for(j in 1:length(rownames(small_RNA_OR_matrix))){
    
    pathway <- rownames(small_RNA_OR_matrix)[[j]]
    
    pos_path <- which(colnames(Background_targeted_epimut) == pathway)
    
    
    in_pathway_and_duration <- nrow(Background_targeted_epimut[Background_targeted_epimut[, pos_dur]==1 & Background_targeted_epimut[, pos_path]==1, ])
    
    in_pathway_not_duration <- nrow(Background_targeted_epimut[Background_targeted_epimut[, pos_dur]==0 & Background_targeted_epimut[, pos_path]==1, ])
    
    not_pathway_but_duration <- nrow(Background_targeted_epimut[Background_targeted_epimut[, pos_dur]==1 & Background_targeted_epimut[, pos_path]==0, ])
    
    not_pathway_not_duration <- nrow(Background_targeted_epimut[Background_targeted_epimut[, pos_dur]==0 & Background_targeted_epimut[, pos_path]==0, ])
    
    
    result <- data.frame(in_pathway_and_duration,  in_pathway_not_duration, not_pathway_but_duration, not_pathway_not_duration)
    
    pathway_results <- rbind(pathway_results, result)
    
  }
  rownames(pathway_results) <- rownames(small_RNA_OR_matrix)
  
  duration_results[[i]] <- pathway_results
  
  names(duration_results)[[i]] <- duration
  
}



# Because CSR1 target in long lived categiory has 0 genes in the overlap this will produce a -Inf result for log2(odds ratio)
# This is difficult to plot so transform all the data by adding 1 to every value



for(i in 1:length(duration_results)){
  duration_results[[i]] <- duration_results[[i]] + 1
  
}


# Now proceed to calculate odds ratios


for(i in 1:length(duration_results)){
  
  for(j in 1:nrow(duration_results[[i]])){
    
    result_line <- duration_results[[i]][j, ]
    
    contingency_table <-
      rbind(c(result_line$in_pathway_and_duration,  result_line$in_pathway_not_duration),
            c(result_line$not_pathway_but_duration, result_line$not_pathway_not_duration))
    
    
    FT_out <- fisher.test(contingency_table)
    pval <- FT_out$p.value
    OddsRatio <- FT_out$estimate
    
    small_RNA_OR_matrix[j, i] <- OddsRatio
    small_RNA_pval_matrix[j, i] <- p.adjust(pval, method = "bonferroni", n= 21)
    
  }}


# Make a bubble plot



small_RNA_OR <- as.data.frame(small_RNA_OR_matrix)

small_RNA_pval <- as.data.frame(small_RNA_pval_matrix)



P_1 <- as.numeric(small_RNA_pval$non_inherited)

log_Odds_1 <- log2(as.numeric(small_RNA_OR$non_inherited))

neg_log_p_1 <- -log2(P_1)

bubble_table_1 <- cbind(rownames(small_RNA_OR), small_RNA_pval$non_inherited, neg_log_p_1, log_Odds_1)

bubble_table_1 <- data.frame(rep("This is the analysis- NON INHERIT"), bubble_table_1)

bubble_table_1[, 2] <- factor(bubble_table_1[, 2], levels = bubble_table_1[, 2])

colnames(bubble_table_1) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")




P_2 <- as.numeric(small_RNA_pval$short_lived)

log_Odds_2 <- log2(as.numeric(small_RNA_OR$short_lived))

neg_log_p_2 <- -log2(P_2)

bubble_table_2 <- cbind(rownames(small_RNA_OR), small_RNA_pval$short_lived, neg_log_p_2, log_Odds_2)

bubble_table_2 <- data.frame(rep("This is the analysis- SHORT LIVED"), bubble_table_2)

bubble_table_2[, 2] <- factor(bubble_table_2[, 2], levels = bubble_table_2[, 2])

colnames(bubble_table_2) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_3 <- as.numeric(small_RNA_pval$long_lived)

log_Odds_3 <- log2(as.numeric(small_RNA_OR$long_lived))

neg_log_p_3 <- -log2(P_3)

bubble_table_3 <- cbind(rownames(small_RNA_OR), small_RNA_pval$long_lived, neg_log_p_3, log_Odds_3)

bubble_table_3 <- data.frame(rep("This is the analysis- LONG LIVED"), bubble_table_3)

bubble_table_3[, 2] <- factor(bubble_table_3[, 2], levels = bubble_table_3[, 2])

colnames(bubble_table_3) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3)




find_size <- c()


for(i in 1:nrow(c_tests)){
  
  
  size_point <- 14
  
  if(as.numeric(c_tests[i, 3]) < 0.1){
    
    size_point <- 50
    
  }
  
  find_size <- c(find_size, size_point)
}


smallRNA_c_combined_tests <- cbind(c_tests, find_size)



guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")

smallRNA_c_combined_tests$Analysis <- factor(smallRNA_c_combined_tests$Analysis, levels = unique(smallRNA_c_combined_tests$Analysis))



# Supplementary Figure 6

# Bubble plot to show distribution of RNAs in chromatin domains according to length


ggplot(smallRNA_c_combined_tests, aes(y=Description, x=as.numeric(log_odds)))+
  geom_point(aes(color=Analysis, size=find_size))+
  expand_limits(x=c(-5, 5))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue", "aquamarine2", "darkcyan"))+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  
  
  scale_alpha(name = paste("Transparency of bubble represents \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value < 0.1", "p value > 0.1"))+
  
  
  labs(y="small RNA pathways\n", x=paste('log2(Odds Ratio) for enrichment'))+
  theme(axis.title=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))







#----------------------------------------------------------------------------------------------
# STEPWISE ENRICHMENT ANALYSIS: ATAC

# 1. What are the odds that genes have both chromatin epimutations and RNA changes ?

# background is all genes which have enhancer/promoter coordinates 

background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1, ])

# of which have RNA changes

with_RNA_changes <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1, ])

# of which have ATAC changes 

with_ATAC_changes <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_ep_epimut==1, ])


# of which have both 

with_both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_RNA_exp_change==1, ])




ATAC_RNA <-  with_both

RNA_not_ATAC <- with_RNA_changes - ATAC_RNA

ATAC_not_RNA <- with_ATAC_changes  -  ATAC_RNA

not_ATAC_not_RNA <- background - (with_RNA_changes + with_ATAC_changes -  ATAC_RNA)


contingency_table <-
  rbind(c(ATAC_RNA, RNA_not_ATAC),
        c(ATAC_not_RNA, not_ATAC_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate





# 2. What are the odds that genes with RNA expression changes have both INHERITED RNA expression changes and chromatin epimutations?


# background is all genes which have enhancer/promoter coordinates AND RNA expression changes 

background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1, ])

# of which have ATAC changes of any kind

with_ATAC_changes <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1, ])


# of which have both inherited RNA changes and ATAC changes

with_both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1, ])



ATAC_RNA <-  with_both

RNA_not_ATAC <- with_RNA_changes - ATAC_RNA

ATAC_not_RNA <- with_ATAC_changes  -  ATAC_RNA

not_ATAC_not_RNA <- background - (with_RNA_changes + with_ATAC_changes -  ATAC_RNA)


contingency_table <-
  rbind(c(ATAC_RNA, RNA_not_ATAC),
        c(ATAC_not_RNA, not_ATAC_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_2 <- FT_out$p.value
OddsRatio_2 <- FT_out$estimate






# Time matching

# 3. What are the odds that Genes with RNA expression changes have both inherited RNA expression changes and time matched chromatin epimutations?


# background is all genes which have enhancer/promoter coordinates AND RNA expression changes

background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1 & All_lins_integrated_table$is_RNA_inherited==1, ])

# of which have time matched ATAC changes 

with_ATAC_changes <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$has_time_matched==1, ])


# of which have both inherited RNA changes and time matched ATAC changes

with_both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$has_time_matched==1, ])



ATAC_RNA <-  with_both

RNA_not_ATAC <- with_RNA_changes - ATAC_RNA

ATAC_not_RNA <- with_ATAC_changes  -  ATAC_RNA

not_ATAC_not_RNA <- background - (with_RNA_changes + with_ATAC_changes -  ATAC_RNA)


contingency_table <-
  rbind(c(ATAC_RNA, RNA_not_ATAC),
        c(ATAC_not_RNA, not_ATAC_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_3 <- FT_out$p.value
OddsRatio_3 <- FT_out$estimate






# 4. What are the odds that in Genes with RNA expression changes and chromatin epimutations, the chromatin epimutations are time matched and the RNA expression changes are inherited?



# background is all genes which have enhancer/promoter coordinates AND RNA expression changes and epimutations

background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1, ])

# of which have time matched ATAC changes

with_ATAC_changes <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$has_time_matched==1, ])


# of which have both inherited RNA changes and time matched ATAC changes

with_both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$has_time_matched==1, ])



ATAC_RNA <-  with_both

RNA_not_ATAC <- with_RNA_changes - ATAC_RNA

ATAC_not_RNA <- with_ATAC_changes  -  ATAC_RNA

not_ATAC_not_RNA <- background - (with_RNA_changes + with_ATAC_changes -  ATAC_RNA)


contingency_table <-
  rbind(c(ATAC_RNA, RNA_not_ATAC),
        c(ATAC_not_RNA, not_ATAC_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate



# This becomes Table 1

# Create a figure to show all the above results 1-4:


OR_frame <- as.data.frame(c(OddsRatio_1, OddsRatio_2, OddsRatio_3, OddsRatio_4))
Pval_frame <- as.data.frame(c(pval_1, pval_2, pval_3, pval_4))  


names <- c("Out of all genes, \nlikelihood of genes having both \nRNA epression changes and chromatin epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and chromatin epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and simultaneous chromatin epimutations",
           "Out of genes with RNAseq changes and chromatin state changes, \nlikelihood of genes having both \n inherited RNA expression changes and simultaneous chromatin epimutations")



Table_1 <- cbind(names, OR_frame, Pval_frame)

colnames(Table_1) <- c("Association_tested", "Odds Ratio", "p-value")

#-------------------------------------------------

# STEPWISE ENRICHMENT ANALYSIS: SMALL RNA

# Table 2

# 1. What are the odds that genes have both small RNA epimutations and RNA changes ?

# background is all genes which have small RNA signal, so have to map to small RNA

background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1, ])

# of which have RNA changes

with_RNA_changes <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1, ])

# of which have small RNA changes 

with_smallRNA_changes <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1, ])


# of which have both 

with_both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_mut==1, ])




smallRNA_RNA <-  with_both

RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA

smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA

not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)


contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_1 <- FT_out$p.value
OddsRatio_1 <- FT_out$estimate





# 2. What are the odds that genes with RNA expression changes have both INHERITED RNA expression changes and small RNA epimutations?


# background is all genes which map to small RNA and have AND RNA expression changes 

background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])

# of which have small RNA changes of any kind

with_smallRNA_changes <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1, ])


# of which have both inherited RNA changes and small RNA changes

with_both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])



smallRNA_RNA <-  with_both

RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA

smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA

not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)


contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_2 <- FT_out$p.value
OddsRatio_2 <- FT_out$estimate






# Time matching

# 3. What are the odds that Genes with RNA expression changes have both inherited RNA expression changes and time matched small RNA epimutations?


# background is all genes which map to small RNA and have RNA expression changes

background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])

# of which have time matched small RNA changes 

with_smallRNA_changes <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA ==1, ])


# of which have both inherited RNA changes and time matched small RNA changes

with_both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1, ])



smallRNA__RNA <-  with_both

RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA

smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA

not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)


contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_3 <- FT_out$p.value
OddsRatio_3 <- FT_out$estimate






# 4. What are the odds that in Genes with RNA expression changes and small RNA epimutations, the small RNA epimutations are time matched and the RNA expression changes are inherited?



# background is all genes which have small RNA data points AND RNA expression changes and small RNA epimutations

background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1, ])

# of which have inherited RNA changes

with_RNA_changes <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])

# of which have time matched small RNA changes

with_smallRNA_changes <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1, ])


# of which have both inherited RNA changes and time matched smallRNA changes

with_both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1, ])



smallRNA_RNA <-  with_both

RNA_not_smallRNA <- with_RNA_changes - smallRNA_RNA

smallRNA_not_RNA <- with_smallRNA_changes  -  smallRNA_RNA

not_smallRNA_not_RNA <- background - (with_RNA_changes + with_smallRNA_changes -  smallRNA_RNA)


contingency_table <-
  rbind(c(smallRNA_RNA, RNA_not_smallRNA),
        c(smallRNA_not_RNA, not_smallRNA_not_RNA))
FT_out <- fisher.test(contingency_table)
pval_4 <- FT_out$p.value
OddsRatio_4 <- FT_out$estimate



# Table 2

# Create a table to show all the above results 1-4:


OR_frame <- as.data.frame(c(OddsRatio_1, OddsRatio_2, OddsRatio_3, OddsRatio_4))
Pval_frame <- as.data.frame(c(pval_1, pval_2, pval_3, pval_4))  


names <- c("Out of all genes, \nlikelihood of genes having both \nRNA epression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and small RNA epimutations", 
           "Out of genes with RNAseq changes, \nlikelihood of genes having both \ninherited RNA expression changes and simultaneous small RNA epimutations",
           "Out of genes with RNAseq changes and small RNA changes, \nlikelihood of genes having both \n inherited RNA expression changes and simultaneous small RNA epimutations")


Table_2 <- cbind(names, OR_frame, Pval_frame)

colnames(Table_2) <- c("Association_tested", "Odds Ratio", "p-value")

#--------------------------------------------------
# Figure 3 D

# RNA exp change  with ATAC concordancy

# Now look within the category of epimutated genes only, i.e. remove genes which may have inherited RNA changes but no chromatin epimutation
# Do we see any bias for directionality in terms of enrichment for inherited RNA changes?   


# Direction matched 

# background is all genes which map to ep, have RNA expression changes and chromatin epimutations  

background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1,])

# of which is inherited RNA expression changes 

inh <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_RNA_inherited==1,])

# of which is concordant chromatin epimutations

concordant <-  nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==1, ])

# of which is both inherited RNA and concordant chromatin epimutations

both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==1, ])


contingency_table <- rbind(c(both, (concordant- both)), 
                           c((inh - both), background - (inh + concordant - both)))



FT_out <- fisher.test(contingency_table)
epimut_concord_pval <- FT_out$p.value
epimut_concord_OddsRatio <- FT_out$estimate



# Direction not matched 

# background is all genes which map to ep, have RNA expression changes and chromatin epimutations  


background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1,])

# of which is inherited RNA expression changes 


inh <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_RNA_inherited==1,])

# of which is discordant chromatin epimutations


discordant <-  nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==0, ])

# of which is both inherited RNA and discordant chromatin epimutations


both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==0, ])


contingency_table <- rbind(c(both, (discordant- both)), 
                           c((inh - both), background - (inh + discordant - both)))



FT_out <- fisher.test(contingency_table)
epimut_discord_pval <- FT_out$p.value
epimut_discord_OddsRatio <- FT_out$estimate






bar_plot <-  as.data.frame(c(log2(epimut_concord_OddsRatio), log2(epimut_discord_OddsRatio)))

names <- c("All genes concordant epimutations", 
           "All genes discordant epimutations")

Pvalues <- c(epimut_concord_pval, epimut_discord_pval)

Histo_Bars_1 <- cbind(bar_plot, names, Pvalues)




colnames(Histo_Bars_1) <- c("Log_OR", "Names", "Pvalues")



# When we split into UP and DOWN and look within the category of epimutated genes only, i.e. remove genes which may have inherited RNA changes but no chromatin epimutation
# Do we see any bias for directionality in terms of enrichment for inherited RNA changes?   


# UP

# Direction matched 

background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_UP==1,])

inh <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_UP==1,])

concordant <-  nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_UP==1, ])

both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_UP==1, ])


contingency_table <- rbind(c(both, (concordant- both)), 
                           c((inh - both), background - (inh + concordant - both)))




FT_out <- fisher.test(contingency_table)
epimut_UP_concord_pval <- FT_out$p.value
epimut_UP_concord_OddsRatio <- FT_out$estimate




# Direction not matched 

background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_UP==1,])

inh <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_UP==1,])

discordant <-  nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==0&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_UP==1, ])

both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==0&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_UP==1, ])


contingency_table <- rbind(c(both, (discordant- both)), 
                           c((inh - both), background - (inh + discordant - both)))



FT_out <- fisher.test(contingency_table)
epimut_UP_discord_pval <- FT_out$p.value
epimut_UP_discord_OddsRatio <- FT_out$estimate




# DOWN

# Direction matched 

background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_DOWN==1,])

inh <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_DOWN==1,])

concordant <-  nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_DOWN==1, ])

both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_DOWN==1, ])


contingency_table <- rbind(c(both, (concordant- both)), 
                           c((inh - both), background - (inh + concordant - both)))



FT_out <- fisher.test(contingency_table)
epimut_DOWN_concord_pval <- FT_out$p.value
epimut_DOWN_concord_OddsRatio <- FT_out$estimate




# Direction not matched 

background <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_DOWN==1,])

inh <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_DOWN==1,])

discordant <-  nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==0&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_DOWN==1, ])

both <- nrow(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1&All_lins_integrated_table$is_RNA_inherited==1&All_lins_integrated_table$has_time_matched==1&All_lins_integrated_table$Direction_Match==0&All_lins_integrated_table$is_ep_epimut==1&All_lins_integrated_table$RNA_All_DOWN==1, ])


contingency_table <- rbind(c(both, (discordant- both)), 
                           c((inh - both), background - (inh + discordant - both)))



FT_out <- fisher.test(contingency_table)
epimut_DOWN_discord_pval <- FT_out$p.value
epimut_DOWN_discord_OddsRatio <- FT_out$estimate





bar_plot <-  as.data.frame(c(log2(epimut_UP_concord_OddsRatio), log2(epimut_DOWN_concord_OddsRatio), 
                             log2(epimut_UP_discord_OddsRatio), log2(epimut_DOWN_discord_OddsRatio)))

names <- c("Concordant UP",
           "Concordant DOWN",
           "Discordant UP", 
           "Discordant DOWN")


Pvalues <- c(epimut_UP_concord_pval, epimut_DOWN_concord_pval,
             epimut_UP_discord_pval, epimut_DOWN_discord_pval)

Histo_Bars_2 <- cbind(bar_plot, names, Pvalues)


colnames(Histo_Bars_2) <- c("Log_OR", "Names", "Pvalues")


D_Histo_Bars <- rbind(Histo_Bars_1, Histo_Bars_2)


D_Histo_Bars$Names <- factor(D_Histo_Bars$Names, levels = D_Histo_Bars$Names)


neg_log_P <- -log2(D_Histo_Bars$Pvalues)


D_Histo_Bars <- cbind(D_Histo_Bars, neg_log_P)


find_size <- c()
find_alpha <- c()

for(i in 1:nrow(D_Histo_Bars)){
  
  alpha <- 0.3
  size_point <- 14
  
  if(as.numeric(D_Histo_Bars[i, 3]) < 0.1){
    alpha <- 1
    size_point <- 50
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}


D_Histo_Bars <- cbind(D_Histo_Bars, find_alpha, find_size)


colnames(D_Histo_Bars) <- c("log_OR", "Names", "Pvalues", "neg_log_P", "find_alpha", "find_size")





D_Histo_Bars$Names <- factor(D_Histo_Bars$Names, levels = rev(D_Histo_Bars$Names))


guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")






RNA_ATAC_UP_DOWN <-
  
  ggplot(D_Histo_Bars, aes(y=Names, x=as.numeric(log_OR)))+
  geom_point(aes(color=Names, size=as.numeric(find_size), alpha=find_alpha))+
  expand_limits(x=c(-1, 3))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.1", "p value not significant > 0.1"))+
  
  
  ggtitle(paste("Figure 3 D"))




# Figure 3 D

RNA_ATAC_UP_DOWN <- 
  RNA_ATAC_UP_DOWN + labs(y="Concordance subset\n", x=bquote(paste('log2(Odds Ratio) for enrichment for inherited RNAseq changes')))+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(D_Histo_Bars$Names))+
  theme(axis.text.x = element_text(color="#000000", size=12))+
  theme(axis.text.y = element_text(color="#000000", size = 16))+
  theme(legend.text=element_text(color="#000000", size=12))

RNA_ATAC_UP_DOWN 


#-------------------------------
# Figure 3 E

# RNA exp change and small RNA concordancy


# RNA exp change  with small RNA concordancy

# Now look within the category of small RNA epimutated genes only, i.e. remove genes which may have inherited RNA changes but no small RNA epimutation
# Do we see any bias for directionality in terms of enrichment for inherited RNA changes?   


# Direction matched 

# background is all genes which map to small RNA, have RNA expression changes and small RNA epimutations  

background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1, ])

# of which is inherited RNA expression changes 

inh <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])

# of which is concordant small RNA epimutations, i.e. time matched and direction matched

concordant <-  nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==1,])

# of which is both inherited RNA and concordant small RNA epimutations

both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==1,])


contingency_table <- rbind(c(both, (concordant- both)), 
                           c((inh - both), background - (inh + concordant - both)))



FT_out <- fisher.test(contingency_table)
smallRNA_concord_pval <- FT_out$p.value
smallRNA_concord_OddsRatio <- FT_out$estimate



# Direction not matched 

# background is all genes which map to small RNA, have RNA expression changes and small RNA epimutations  

background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1,])

# of which is inherited RNA expression changes 

inh <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])


# of which is discordant small RNA  epimutations, i.e. time matched but direction not matched

discordant <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==0,])


# of which is both inherited RNA and discordant small RNA  epimutations

both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==0,])


contingency_table <- rbind(c(both, (discordant- both)), 
                           c((inh - both), background - (inh + discordant - both)))



FT_out <- fisher.test(contingency_table)
smallRNA_discord_pval <- FT_out$p.value
smallRNA_discord_OddsRatio <- FT_out$estimate






bar_plot <-  as.data.frame(c(log2(smallRNA_concord_OddsRatio), log2(smallRNA_discord_OddsRatio)))

names <- c("All genes concordant epimutations", 
           "All genes discordant epimutations")

Pvalues <- c(smallRNA_concord_pval, smallRNA_discord_pval)

Histo_Bars_1 <- cbind(bar_plot, names, Pvalues)




colnames(Histo_Bars_1) <- c("Log_OR", "Names", "Pvalues")



# When we split into UP and DOWN and look within the category of epimutated genes only, i.e. remove genes which may have inherited RNA changes but no small RNA epimutation
# Do we see any bias for directionality in terms of enrichment for inherited RNA changes?   


# UP

# Direction matched 

background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1,])

inh <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])

concordant <-  nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==1&smallRNA_All_lins_integrated_table$RNA_UP==1,])

both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==1&smallRNA_All_lins_integrated_table$RNA_UP==1&smallRNA_All_lins_integrated_table$RNA_inherited==1,])


contingency_table <- rbind(c(both, (concordant- both)), 
                           c((inh - both), background - (inh + concordant - both)))




FT_out <- fisher.test(contingency_table)
smallRNA_UP_concord_pval <- FT_out$p.value
smallRNA_UP_concord_OddsRatio <- FT_out$estimate




# Direction not matched 

background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1,])

inh <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])

discordant <-  nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==0&smallRNA_All_lins_integrated_table$RNA_UP==1,])

both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==0&smallRNA_All_lins_integrated_table$RNA_UP==1&smallRNA_All_lins_integrated_table$RNA_inherited==1,])


contingency_table <- rbind(c(both, (discordant- both)), 
                           c((inh - both), background - (inh + discordant - both)))



FT_out <- fisher.test(contingency_table)
smallRNA_UP_discord_pval <- FT_out$p.value
smallRNA_UP_discord_OddsRatio <- FT_out$estimate




# DOWN

# Direction matched 


background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1,])

inh <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])

concordant <-  nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==1&smallRNA_All_lins_integrated_table$RNA_DOWN==1,])

both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==1&smallRNA_All_lins_integrated_table$RNA_DOWN==1&smallRNA_All_lins_integrated_table$RNA_inherited==1,])


contingency_table <- rbind(c(both, (concordant- both)), 
                           c((inh - both), background - (inh + concordant - both)))




FT_out <- fisher.test(contingency_table)
smallRNA_DOWN_concord_pval <- FT_out$p.value
smallRNA_DOWN_concord_OddsRatio <- FT_out$estimate





# Direction not matched 



background <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1,])

inh <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$RNA_inherited==1, ])

discordant <-  nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==0&smallRNA_All_lins_integrated_table$RNA_DOWN==1,])

both <- nrow(smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$maps_to_smallRNA ==1&smallRNA_All_lins_integrated_table$RNA_mut==1&smallRNA_All_lins_integrated_table$smallRNA_mut==1&smallRNA_All_lins_integrated_table$time_matched_to_RNA==1&smallRNA_All_lins_integrated_table$Direction_match_to_RNA==0&smallRNA_All_lins_integrated_table$RNA_DOWN==1&smallRNA_All_lins_integrated_table$RNA_inherited==1,])


contingency_table <- rbind(c(both, (discordant- both)), 
                           c((inh - both), background - (inh + discordant - both)))




FT_out <- fisher.test(contingency_table)
smallRNA_DOWN_discord_pval <- FT_out$p.value
smallRNA_DOWN_discord_OddsRatio <- FT_out$estimate






bar_plot <-  as.data.frame(c(log2(smallRNA_UP_concord_OddsRatio), log2(smallRNA_DOWN_concord_OddsRatio), 
                             log2(smallRNA_UP_discord_OddsRatio), log2(smallRNA_DOWN_discord_OddsRatio)))

names <- c("Concordant UP",
           "Concordant DOWN",
           "Discordant UP", 
           "Discordant DOWN")


Pvalues <- c(smallRNA_UP_concord_pval, smallRNA_DOWN_concord_pval,
             smallRNA_UP_discord_pval, smallRNA_DOWN_discord_pval)

Histo_Bars_2 <- cbind(bar_plot, names, Pvalues)


colnames(Histo_Bars_2) <- c("Log_OR", "Names", "Pvalues")


E_Histo_Bars <- rbind(Histo_Bars_1, Histo_Bars_2)


E_Histo_Bars$Names <- factor(E_Histo_Bars$Names, levels = E_Histo_Bars$Names)


neg_log_P <- -log2(E_Histo_Bars$Pvalues)


E_Histo_Bars <- cbind(E_Histo_Bars, neg_log_P)


find_size <- c()
find_alpha <- c()

for(i in 1:nrow(E_Histo_Bars)){
  
  alpha <- 0.3
  size_point <- 14
  
  if(as.numeric(E_Histo_Bars[i, 3]) < 0.1){
    alpha <- 1
    size_point <- 50
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}


E_Histo_Bars <- cbind(E_Histo_Bars, find_alpha, find_size)


colnames(E_Histo_Bars) <- c("log_OR", "Names", "Pvalues", "neg_log_P", "find_alpha", "find_size")




E_Histo_Bars$Names <- factor(E_Histo_Bars$Names, levels = rev(E_Histo_Bars$Names))


guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")






RNA_smallRNA_UP_DOWN <-
  
  ggplot(E_Histo_Bars, aes(y=Names, x=as.numeric(log_OR)))+
  geom_point(aes(color=Names, size=find_size, alpha=find_alpha))+
  expand_limits(x=c(-1, 3))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  

  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.1", "p value not significant > 0.1"))+
  
  
  ggtitle(paste("Figure 3 E"))




# Figure 3 E

RNA_smallRNA_UP_DOWN <- RNA_smallRNA_UP_DOWN + labs(y="Concordance subset\n", x=bquote(paste('log2(Odds Ratio) for enrichment for inherited RNAseq changes')))+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(E_Histo_Bars$Names))+
  theme(axis.text.x = element_text(color="#000000", size=12))+
  theme(axis.text.y = element_text(color="#000000", size = 16))+
  theme(legend.text=element_text(color="#000000", size=12))


RNA_smallRNA_UP_DOWN 




#-------------------------------------------------------------------------------------------


# Finding long-lived, short-lived and non-inherited epimutations


# Only consider genes with epimutations

RNA_minus_suspect_muts <- RNA_minus_suspect[RNA_minus_suspect$length >0, ]
ATAC_minus_suspect_muts <- ATAC_ep_minus_suspect[ATAC_ep_minus_suspect$length >0, ]
small_RNA_minus_suspect_muts <- small_RNA_minus_suspect[small_RNA_minus_suspect$length >0, ]



# Define long lasting, short lasting and non inherited epimutations for RNA, ATAC and small RNA 

library(lme4)

library(mixtools)

library(ggVennDiagram)

#----------------------------
# RNA length categorisation


RNA_ratio <- (as.numeric(RNA_minus_suspect$length)/as.numeric(RNA_minus_suspect$number_transitions))

RNA_ratio[RNA_ratio %in% NaN] <- 0

RNA_tab <- data.frame(RNA_minus_suspect$gene, RNA_ratio)

colnames(RNA_tab) <- c("ID","meanLength")


#filter so that only epimutations are considered.

RNA_InputTab<-RNA_tab[which(RNA_tab[,2]>0),]



RNAModel<-lmer(meanLength~(1|ID),data=RNA_InputTab)

lengthEst <-ranef(RNAModel)[[1]]




# K means clustering

Clusters <- kmeans(lengthEst, 2)

RNA_KMeans_Cluster_1 <- names(Clusters$cluster[Clusters$cluster == 1])

RNA_KMeans_Cluster_2 <- names(Clusters$cluster[Clusters$cluster == 2])



if(length(RNA_KMeans_Cluster_1) < length(RNA_KMeans_Cluster_2)){
  
  RNA_longest_genes <- RNA_KMeans_Cluster_1
}else{RNA_longest_genes <- RNA_KMeans_Cluster_2}




RNA_Top_genes <- unique(All_lins_integrated_table[All_lins_integrated_table$gene %in% RNA_longest_genes, 3])





# Get the different gene categories with RNA expression changes


never_epimutated <- unique(RNA_All[RNA_All$length ==0, 1])

# Now find the epimutations which have a length of 1

epimutation_length_one <- unique(RNA_minus_suspect[RNA_minus_suspect$length == 1, 1])

# Now find the epimutations which have a length > 1

epimutation_length_two_plus <- unique(RNA_minus_suspect[RNA_minus_suspect$length > 1, 1])

# now set difference between inherited and long lasting

RNA_Muts_Inherited_not_long <- setdiff(epimutation_length_two_plus, RNA_Top_genes)

# Now set difference between non inherited and inherited not long

RNA_Muts_Not_Inherited <- setdiff(epimutation_length_one, RNA_Muts_Inherited_not_long)

# Now set difference between non inherited and long lasting

RNA_Muts_Not_Inherited <- setdiff(RNA_Muts_Not_Inherited, RNA_Top_genes)

# Now remove any that are epimutated from the never epimutated genes

RNA_No_Muts <- setdiff(never_epimutated, c(RNA_Top_genes, RNA_Muts_Inherited_not_long, RNA_Muts_Not_Inherited))





# Get the gene names

RNA_Muts_Inherited_not_long_genes <- c()

for(i in 1:length(RNA_Muts_Inherited_not_long)){
  
  name <-  strsplit(RNA_Muts_Inherited_not_long[[i]], ":")[[1]][4]
  RNA_Muts_Inherited_not_long_genes <- c(RNA_Muts_Inherited_not_long_genes, name)
  
}



RNA_Muts_Not_Inherited_genes <- c()

for(i in 1:length(RNA_Muts_Not_Inherited)){
  
  name <-  strsplit(RNA_Muts_Not_Inherited[[i]], ":")[[1]][4]
  RNA_Muts_Not_Inherited_genes <- c(RNA_Muts_Not_Inherited_genes, name)
  
}


RNA_No_Muts_genes <- c()

for(i in 1:length(RNA_No_Muts)){
  
  name <-  strsplit(RNA_No_Muts[[i]], ":")[[1]][4]
  RNA_No_Muts_genes <- c(RNA_No_Muts_genes, name)
  
}



RNA_longest_genes <- c()

for(i in 1:length(RNA_Top_genes)){
  
  name <-  strsplit(RNA_Top_genes[[i]], ":")[[1]][4]
  RNA_longest_genes <- c(RNA_longest_genes, name)
  
}





# Supplementary Figure 3 A

# Make a density plot 

tab <- data.frame(lengthEst, Clusters$cluster)

colnames(tab) <- c("lmer_est", "cluster")

tab$cluster <- factor(tab$cluster, levels = c(2, 1))

ggplot(tab, aes(x=lmer_est)) + geom_density(aes(group=cluster, colour=cluster, fill=cluster), alpha=0.3)+
  theme_bw()+
  scale_x_continuous("lmer_est", breaks= seq(-1, 10, 0.5), limits = c(-1, 10))+
  scale_y_continuous("density",  breaks = seq(0, 35, 5), limits = c(0, 35))+
  theme_classic()+
  ggtitle("Supplementary Figure 3 A")





# Verify that K means produces reasonable number of genes with Normal Mix Em 



# Normal Mix EM



Normal_MixEM<- normalmixEM(lengthEst,2)


#plot to check

# plot(Normal_MixEM, which=2, xlim= c(-2, 6), xlab2= "lmer estimate", main2 ="RNA Normal MixEM", border = "azure3")


RNA_Normal_MixEM_model <- data.frame(lengthEst, Normal_MixEM$posterior)

# It appears that genes with very low p values in column 1 (or p = 1 in column 2) are in the second long lasting distribution

# To obtain a similar number of genes as in the two methods above I set the p value to 1e-30

MixEM_RNA_longest_genes <- rownames(RNA_Normal_MixEM_model[RNA_Normal_MixEM_model$comp.1 < 1e-30, ])



# Supplementary Figure 3 D

list_2_methods <- list(RNA_longest_genes, MixEM_RNA_longest_genes)

RNA_Venn <- ggVennDiagram(list_2_methods, label = "percent", category.names = c("K-means clustering", "MixEM method"))

RNA_Venn + scale_fill_distiller(palette = "Reds", direction = 1) + ggtitle(paste("Supplementary Figure 3 D")) 
  

                                                                                                                                                 





#---------------------------
# ATAC length categorisation 

ATAC_ABC_ep_ratio <- (as.numeric(ATAC_ep_minus_suspect$length)/as.numeric(ATAC_ep_minus_suspect$number_transitions))

ATAC_ABC_ep_ratio[ATAC_ABC_ep_ratio  %in% NaN] <- 0

ATAC_ABC_ep_tab <- data.frame(ATAC_ep_minus_suspect$names, ATAC_ABC_ep_ratio)

colnames(ATAC_ABC_ep_tab) <- c("ID","meanLength")



# Just select the coordinates which have epimutations

ATAC_ABC_ep_tab <- ATAC_ABC_ep_tab[which(ATAC_ABC_ep_tab[,2]>0),]

# the intercept (shown by '1' can vary by ID)

ATACModel<-lmer(meanLength~(1|ID),data=ATAC_ABC_ep_tab) 

#this makes the mean length depend on the particular gene as a random factor.  

# IE does the identity of the gene influence how long it lasts for on average?

#for each gene an estimate of the influence that it makes on the mean length is generated.  This information is what we want.
#to extract this:

lengthEst <-ranef(ATACModel)[[1]]




# K means clustering

Clusters<-kmeans(lengthEst,2)

ATAC_KMeans_Cluster_1 <- names(Clusters$cluster[Clusters$cluster == 1])

ATAC_KMeans_Cluster_2 <- names(Clusters$cluster[Clusters$cluster == 2])


if(length(ATAC_KMeans_Cluster_1) < length(ATAC_KMeans_Cluster_2)){
  
  ATAC_ep_ABC_Top <- ATAC_KMeans_Cluster_1
}else{ATAC_ep_ABC_Top <- ATAC_KMeans_Cluster_2}






# Get the categories of genes with chromatin epimutations

# find never epimutated genes

never_epimutated <- unique(ATAC_ep_All[ATAC_ep_All$length ==0, 1])


# Now find the epimutations which have a length of 1

epimutation_length_one <- unique(ATAC_ep_minus_suspect[ATAC_ep_minus_suspect$length == 1, 1])


# Now find the epimutations which have a length > 1

epimutation_length_two_plus <- unique(ATAC_ep_minus_suspect[ATAC_ep_minus_suspect$length > 1, 1])


# now set difference between inherited and long lasting

ATAC_ep_Muts_Inherited_not_long <- setdiff(epimutation_length_two_plus, ATAC_ep_ABC_Top)


# Now set difference between non inherited and inherited not long & long

ATAC_ep_Muts_Not_Inherited <- setdiff(epimutation_length_one, c(ATAC_ep_Muts_Inherited_not_long, ATAC_ep_ABC_Top))


# Now remove any that are epimutated from the never epimutated genes

ATAC_ep_No_Muts <- setdiff(never_epimutated, c(ATAC_ep_ABC_Top, ATAC_ep_Muts_Inherited_not_long, ATAC_ep_Muts_Not_Inherited))





# Now get the actual gene names

library(stringr)

# split any gene names at commas and make the overall lists unique


ATAC_ep_No_Muts_genes <- unique(unlist(str_split((Ahringer_ref_table[Ahringer_ref_table$Locus %in% ATAC_ep_No_Muts, 3]), pattern = ",")))

ATAC_ep_Muts_Not_Inherited_genes <- unique(unlist(str_split((Ahringer_ref_table[Ahringer_ref_table$Locus %in% ATAC_ep_Muts_Not_Inherited, 3]), pattern = ",")))

ATAC_ep_Muts_Inherited_not_long_genes <- unique(unlist(str_split((Ahringer_ref_table[Ahringer_ref_table$Locus %in% ATAC_ep_Muts_Inherited_not_long, 3]), pattern = ",")))

ATAC_ep_ABC_Top_genes <- unique(unlist(str_split((Ahringer_ref_table[Ahringer_ref_table$Locus %in% ATAC_ep_ABC_Top, 3]), pattern = ",")))


list_ATAC_ep_genes <- list(ATAC_ep_No_Muts_genes, 
                           ATAC_ep_Muts_Not_Inherited_genes,
                           ATAC_ep_Muts_Inherited_not_long_genes,
                           ATAC_ep_ABC_Top_genes)

names(list_ATAC_ep_genes) <- c("Never_mut", "Single", "Short", "Long") 





# Supplementary Figure 3 B

# Make a density plot 

tab <- data.frame(lengthEst, Clusters$cluster)

colnames(tab) <- c("lmer_est", "cluster")

tab$cluster <- factor(tab$cluster, levels = c(1, 2))

ggplot(tab, aes(x=lmer_est)) + geom_density(aes(group=cluster, colour=cluster, fill=cluster), alpha=0.3)+
  theme_bw()+
  scale_x_continuous("lmer_est", breaks= seq(-1, 10, 0.5), limits = c(-1, 10))+
  scale_y_continuous("density",  breaks = seq(0, 35, 5), limits = c(0, 35))+
  theme_classic()+
  ggtitle(paste("Supplementary Figure 3 B"))






# Mix EM method

lengthEst <- as.matrix(lengthEst)

Normal_MixEM<- normalmixEM(lengthEst,2)

ATAC_Normal_MixEM_model <- data.frame(lengthEst, Normal_MixEM$posterior)

# It appears that genes with very low p values in column 2 (or p = 1 in column 1) are in the long lasting distribution

# To obtain a similar number of genes as in the two methods above I set the p value to 1e-30

# You have to check to see if it is comp.1 or comp.2

MixEM_ATAC_longest_genes <- rownames(ATAC_Normal_MixEM_model[ATAC_Normal_MixEM_model$comp.1 < 1e-30, ])

MixEM_ATAC_longest_genes <- unique(unlist(str_split((Ahringer_ref_table[Ahringer_ref_table$Locus %in% MixEM_ATAC_longest_genes, 3]), pattern = ",")))




# Supplementary Figure 3 E


list_2_methods <- list(ATAC_ep_ABC_Top_genes, MixEM_ATAC_longest_genes)

ATAC_Venn <- ggVennDiagram(list_2_methods, label = "percent", category.names = c("K-means clustering", "MixEM method"))

ATAC_Venn + scale_fill_distiller(palette = "Blues", direction = 1) + ggtitle(paste("Supplementary Figure 3 E")) 





#----------------------------
# small RNA length categorisation 

small_RNA_ratio <- (as.numeric(small_RNA_minus_suspect_muts$length)/as.numeric(small_RNA_minus_suspect_muts$number_transitions))

small_RNA_tab <- data.frame(small_RNA_minus_suspect_muts$gene_name, small_RNA_ratio)

colnames(small_RNA_tab) <- c("ID","meanLength")



#filter so that only epimutations are considered.
small_RNA_InputTab<-small_RNA_tab[which(small_RNA_tab[,2]>0),]

smallRNAModel<-lmer(meanLength~(1|ID),data=small_RNA_InputTab)



lengthEst <-ranef(smallRNAModel)[[1]]


# K means clustering

Clusters<-kmeans(lengthEst,2)

smallRNA_KMeans_Cluster_1 <- names(Clusters$cluster[Clusters$cluster == 1])

smallRNA_KMeans_Cluster_2 <- names(Clusters$cluster[Clusters$cluster == 2])


if(length(smallRNA_KMeans_Cluster_1) < length(smallRNA_KMeans_Cluster_2)){

  small_RNA_Top_genes <- smallRNA_KMeans_Cluster_1
    }else{small_RNA_Top_genes <- smallRNA_KMeans_Cluster_2}


# Find the categories of genes with 22G-RNA epimutations

never_epimutated <- unique(small_RNA_All[small_RNA_All$length ==0, 1])

# Now find the epimutations which have a length of 1
epimutation_length_one <- unique(small_RNA_minus_suspect[small_RNA_minus_suspect$length == 1, 1])

# Now find the epimutations which have a length > 1
epimutation_length_two_plus <- unique(small_RNA_minus_suspect[small_RNA_minus_suspect$length > 1, 1])

# now set difference between inherited and long lasting
small_RNA_Muts_Inherited_not_long <- setdiff(epimutation_length_two_plus, small_RNA_Top_genes)

# Now set difference between non inherited and inherited not long
small_RNA_Muts_Not_Inherited <- setdiff(epimutation_length_one, small_RNA_Muts_Inherited_not_long)

# Now set difference between non inherited and long lasting
small_RNA_Muts_Not_Inherited <- setdiff(small_RNA_Muts_Not_Inherited, small_RNA_Top_genes)

# Now remove any that are epimutated from the never epimutated genes
small_RNA_No_Muts <- setdiff(never_epimutated, c(small_RNA_Top_genes, small_RNA_Muts_Inherited_not_long, small_RNA_Muts_Not_Inherited))




# Supplementary Figure 3 C

# Make a density plot 

tab <- data.frame(lengthEst, Clusters$cluster)

colnames(tab) <- c("lmer_est", "cluster")

tab$cluster <- factor(tab$cluster, levels = c(1, 2))

ggplot(tab, aes(x=lmer_est)) + geom_density(aes(group=cluster, colour=cluster, fill=cluster), alpha=0.3)+
  theme_bw()+
  scale_x_continuous("lmer_est", breaks= seq(-1, 10, 0.5), limits = c(-1, 10))+
  scale_y_continuous("density",  breaks = seq(0, 35, 5), limits = c(0, 35))+
  theme_classic()+
  ggtitle(paste("Supplementary Figure 3 C"))




# Mix EM method

lengthEst <- as.matrix(lengthEst)

Normal_MixEM<- normalmixEM(lengthEst,2)

smallRNA_Normal_MixEM_model <- data.frame(lengthEst, Normal_MixEM$posterior)

# It appears that genes with very low p values in column 1 (or p = 1 in column 2) are in the second long lasting distribution

# To obtain a similar number of genes as in the two methods above I set the p value to 1e-30

MixEM_smallRNA_longest_genes <- rownames(smallRNA_Normal_MixEM_model[smallRNA_Normal_MixEM_model$comp.1 < 1e-30, ])



# Supplementary Figure 3 F

list_2_methods <- list(small_RNA_Top_genes, MixEM_smallRNA_longest_genes)

brewer.pal(n = 3, name = "YlOrRd")

smallRNA_Venn <- ggVennDiagram(list_2_methods, label = "percent", category.names = c("K-means clustering", "MixEM method"))

smallRNA_Venn + scale_fill_distiller(palette = "Greens", direction = 1) + ggtitle(paste("Supplementary Figure 3 F")) 




#--------------------------------------------
# Figure 4 

# Bubble plots for distribution of length categorised epimutations in different chromatin domains

# K means clustering used to obtain long-lived genes

# RNA Active and Regulated domain IDs

# These are RNA_Active_domains and RNA_Regulated_domains


# When we are considering the genes with chromatin domain annotations we can only consider the genes which map to Ahringer as this is where we get the domain annotations from. 

all_Ahr_RNA_genes <- unique(RNA_centric_All_table[, 2])


# The lists of length categorised expression changes must be restricted to genes which have Ahringer annotations 

Ahr_RNA_not_inherit <- intersect(RNA_Muts_Not_Inherited_genes, all_Ahr_RNA_genes)
Ahr_RNA_short <- intersect(RNA_Muts_Inherited_not_long_genes, all_Ahr_RNA_genes)
Ahr_RNA_long <- intersect(RNA_longest_genes, all_Ahr_RNA_genes)


# small_RNA_background are genes which map to Ahringer

Active_genes <- unique(RNA_centric_All_table[RNA_centric_All_table$Domain == "A", 2])
Regulated_genes <- unique(RNA_centric_All_table[RNA_centric_All_table$Domain == "R", 2])
Chr_X_genes <- unique(RNA_centric_All_table[RNA_centric_All_table$Domain == ".", 2])



# the background data set therefore is length(unique(small_RNA_background$Target_gene))


RNA_test_list <- list(Ahr_RNA_not_inherit, Ahr_RNA_short, Ahr_RNA_long)

RNA_Features <- list(Regulated_genes, Active_genes, piRNA_cluster_genes, Chr_X_genes)








# Define the function

epiIntersect <-
  function(common,
           Epimutation_total,
           Feature_total,
           overall_total) {
    contingency_table <-
      rbind(c(common, Epimutation_total - common),
            c(
              Feature_total - common,
              overall_total - (Epimutation_total + Feature_total - common)
            ))
    FT_out <- fisher.test(contingency_table)
    pval <- FT_out$p.value
    OddsRatio <- FT_out$estimate
    return(c(pval, OddsRatio))
  }





Non_Nested_RNA_OR <- c()

Non_Nested_RNA_pval <- c()


for(d in 1:length(RNA_test_list)){
  
  test <- RNA_test_list[[d]]
  
  temp_pval <- matrix(0, ncol = 1, nrow = 4)
  temp_oR <- matrix(0, ncol = 1, nrow = 4)
  
  for (i in 1:length(RNA_Features)) {
    common_in <- length(intersect(RNA_Features[[i]], test))
    Epimutation_total_in <- length(test)
    Feature_total_in <- length(RNA_Features[[i]])
    
    output <- epiIntersect(
      common = common_in,
      Epimutation_total = Epimutation_total_in,
      Feature_total = Feature_total_in,
      overall_total = 12617  # Total number of RNA genes mapping to Ahringer
    )
    temp_pval[i, ] <- output[1]
    temp_oR[i, ] <- output[2]
  }
  
  
  
  temp_pval <- as.matrix(p.adjust(temp_pval, "bonferroni", 4))
  
  
  colnames(temp_pval) <-
    c("pval")
  rownames(temp_pval) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  
  colnames(temp_oR) <-
    c("oR")
  rownames(temp_oR) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  Non_Nested_RNA_OR <- cbind(Non_Nested_RNA_OR, temp_oR)
  
  Non_Nested_RNA_pval <- cbind(Non_Nested_RNA_pval, temp_pval)
  
}


colnames(Non_Nested_RNA_OR) <- c("Not_inherited", "Short", "Long")
colnames(Non_Nested_RNA_pval) <- c("Not_inherited", "Short", "Long")



#------------------
# RNA expression changes distributions in distinct chromatin domains


ordered_sequence <- c("piRNA", "ChrX", "Regulated", "Active")

colnames(Non_Nested_RNA_OR) <- c("Test_1", "Test_2", "Test_3")

colnames(Non_Nested_RNA_pval) <- c("Test_1", "Test_2", "Test_3")


Non_Nested_RNA_OR <- as.data.frame(Non_Nested_RNA_OR[ordered_sequence, ])

Non_Nested_RNA_pval <- as.data.frame(Non_Nested_RNA_pval[ordered_sequence, ])



P_1 <- as.numeric(Non_Nested_RNA_pval$Test_1)

log_Odds_1 <- log2(as.numeric(Non_Nested_RNA_OR$Test_1))

neg_log_p_1 <- -log2(P_1)

bubble_table_1 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_1, neg_log_p_1, log_Odds_1)

bubble_table_1 <- data.frame(rep("This is the analysis- NON INHERIT"), bubble_table_1)

bubble_table_1[, 2] <- factor(bubble_table_1[, 2], levels = bubble_table_1[, 2])

colnames(bubble_table_1) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_2 <- as.numeric(Non_Nested_RNA_pval$Test_2)

log_Odds_2 <- log2(as.numeric(Non_Nested_RNA_OR$Test_2))

neg_log_p_2 <- -log2(P_2)

bubble_table_2 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_2, neg_log_p_2, log_Odds_2)

bubble_table_2 <- data.frame(rep("This is the analysis- SHORT"), bubble_table_2)

bubble_table_2[, 2] <- factor(bubble_table_2[, 2], levels = bubble_table_2[, 2])

colnames(bubble_table_2) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_3 <- as.numeric(Non_Nested_RNA_pval$Test_3)

log_Odds_3 <- log2(as.numeric(Non_Nested_RNA_OR$Test_3))

neg_log_p_3 <- -log2(P_3)

bubble_table_3 <- cbind(rownames(Non_Nested_RNA_OR), Non_Nested_RNA_pval$Test_3, neg_log_p_3, log_Odds_3)

bubble_table_3 <- data.frame(rep("This is the analysis- LONG"), bubble_table_3)

bubble_table_3[, 2] <- factor(bubble_table_3[, 2], levels = bubble_table_3[, 2])

colnames(bubble_table_3) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3)




find_size <- c()


for(i in 1:nrow(c_tests)){
  

  size_point <- 14
  
  if(as.numeric(c_tests[i, 3]) < 0.1){
   
    size_point <- 50
    
  }

  find_size <- c(find_size, size_point)
}


RNA_c_combined_tests <- cbind(c_tests, find_size)






guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")


RNA_c_combined_tests$Analysis <- factor(RNA_c_combined_tests$Analysis, levels = unique(RNA_c_combined_tests$Analysis))




# Figure 4 A 

# Bubble plot to show distribution of RNAs in chromatin domains according to length



RNA_Lengths_Domain <-
  
  ggplot(RNA_c_combined_tests, aes(y=Description, x=as.numeric(log_odds)))+
  geom_point(aes(color=Analysis, size=find_size))+
  expand_limits(x=c(-3, 3))+
  scale_x_continuous(breaks = seq(-3, 3, by = 1))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue", "aquamarine2", "darkcyan"))+
  
  

  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  

  scale_alpha(name = paste("Transparency of bubble represents \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value < 0.1", "p value > 0.1"))+
  
 
  
  ggtitle(paste("Figure 4 A"))



RNA_Lengths_Domain <- 
  
  RNA_Lengths_Domain + labs(y="Chromatin Domains\n", x=paste('log2(Odds Ratio) for enrichment'))+
  theme(axis.title=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

RNA_Lengths_Domain 





#------------------------------------------------------------
# ATAC seq

# Analysis done according to loci
# The background is all Ahringer mapped loci


# Regulatory element Active and Regulated loci

ep_Active_domains <- unique(All_lins_ATAC_table[All_lins_ATAC_table$Chromatin_domain == "A", 3])

ep_Regulated_domains <- unique(All_lins_ATAC_table[All_lins_ATAC_table$Chromatin_domain == "R", 3])

ep_X_domains <- unique(All_lins_ATAC_table[grep("chrX:", All_lins_ATAC_table$Locus), 3])

ep_piRNA_domains <- unique(All_lins_ATAC_table[All_lins_ATAC_table$piRNA_cluster==1, 3])





ATAC_test_list <- list(ATAC_ep_Muts_Not_Inherited, ATAC_ep_Muts_Inherited_not_long, ATAC_ep_ABC_Top)

ep_Features <- list(ep_Active_domains, ep_Regulated_domains, ep_piRNA_domains, ep_X_domains)


Non_Nested_ATAC_OR <- c()

Non_Nested_ATAC_pval <- c()


# For the background use the number of loci in the binarised ATAC table as these all have ahringer annotations 

# nrow(ATAC_binarised_table_A_2.25_ehancer_promoter)



for(d in 1:length(ATAC_test_list)){
  
  test <- ATAC_test_list[[d]]
  
  temp_pval <- matrix(0, ncol = 1, nrow = 4)
  temp_oR <- matrix(0, ncol = 1, nrow = 4)
  
  for (i in 1:length(ep_Features)) {
    common_in <- length(intersect(ep_Features[[i]], test))
    Epimutation_total_in <- length(test)
    Feature_total_in <- length(ep_Features[[i]])
    output <- epiIntersect(
      common = common_in,
      Epimutation_total = Epimutation_total_in,
      Feature_total = Feature_total_in,
      overall_total = 31618
    )
    temp_pval[i, ] <- output[1]
    temp_oR[i, ] <- output[2]
  }
  
  # Bonferroni adjust:
  
  temp_pval <- as.matrix(p.adjust(temp_pval, "bonferroni", 4))
  
  
  colnames(temp_pval) <-
    c("pval")
  rownames(temp_pval) <-
    c("Active", "Regulated","piRNA", "ChrX")
  
  
  colnames(temp_oR) <-
    c("oR")
  rownames(temp_oR) <-
    c("Active", "Regulated", "piRNA", "ChrX")
  
  Non_Nested_ATAC_OR <- cbind(Non_Nested_ATAC_OR, temp_oR)
  
  Non_Nested_ATAC_pval <- cbind(Non_Nested_ATAC_pval, temp_pval)
  
}


colnames(Non_Nested_ATAC_OR) <- c("Not_inherited", "Short", "Long")

colnames(Non_Nested_ATAC_pval) <- c("Not_inherited", "Short", "Long")


# -----------------

# ATAC epimutations distributions in distinct chromatin domains


ordered_sequence <- c("piRNA", "ChrX", "Regulated", "Active")


# ATAC plot


colnames(Non_Nested_ATAC_OR) <- c("Test_1", "Test_2", "Test_3")

colnames(Non_Nested_ATAC_pval) <- c("Test_1", "Test_2", "Test_3")

rownames(Non_Nested_ATAC_OR) <- c("Active", "Regulated", "piRNA", "ChrX")

rownames(Non_Nested_ATAC_pval) <- c("Active", "Regulated", "piRNA", "ChrX")


Non_Nested_ATAC_OR <- as.data.frame(Non_Nested_ATAC_OR[ordered_sequence, ])

Non_Nested_ATAC_pval <- as.data.frame(Non_Nested_ATAC_pval[ordered_sequence, ])




P_1 <- as.numeric(Non_Nested_ATAC_pval$Test_1)

log_Odds_1 <- log2(as.numeric(Non_Nested_ATAC_OR$Test_1))

neg_log_p_1 <- -log2(P_1)

bubble_table_1 <- cbind(rownames(Non_Nested_ATAC_OR), Non_Nested_ATAC_pval$Test_1, neg_log_p_1, log_Odds_1)

bubble_table_1 <- data.frame(rep("This is the analysis- NON INHERIT"), bubble_table_1)

bubble_table_1[, 2] <- factor(bubble_table_1[, 2], levels = bubble_table_1[, 2])

colnames(bubble_table_1) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")





P_2 <- as.numeric(Non_Nested_ATAC_pval$Test_2)

log_Odds_2 <- log2(as.numeric(Non_Nested_ATAC_OR$Test_2))

neg_log_p_2 <- -log2(P_2)

bubble_table_2 <- cbind(rownames(Non_Nested_ATAC_OR), Non_Nested_ATAC_pval$Test_2, neg_log_p_2, log_Odds_2)

bubble_table_2 <- data.frame(rep("This is the analysis- SHORT"), bubble_table_2)

bubble_table_2[, 2] <- factor(bubble_table_2[, 2], levels = bubble_table_2[, 2]) 

colnames(bubble_table_2) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")





P_3 <- as.numeric(Non_Nested_ATAC_pval$Test_3)

log_Odds_3 <- log2(as.numeric(Non_Nested_ATAC_OR$Test_3))

neg_log_p_3 <- -log2(P_3)

bubble_table_3 <- cbind(rownames(Non_Nested_ATAC_OR), Non_Nested_ATAC_pval$Test_3, neg_log_p_3, log_Odds_3)

bubble_table_3 <- data.frame(rep("This is the analysis- LONG"), bubble_table_3)

bubble_table_3[, 2] <- factor(bubble_table_3[, 2], levels = bubble_table_3[, 2])

colnames(bubble_table_3) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3)




find_size <- c()


for(i in 1:nrow(c_tests)){
  

  size_point <- 14
  
  if(as.numeric(c_tests[i, 3]) < 0.1){

    size_point <- 50
    
  }

  find_size <- c(find_size, size_point)
}


ATAC_c_combined_tests <- cbind(c_tests, find_size)




guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")




ATAC_c_combined_tests$Analysis <- factor(ATAC_c_combined_tests$Analysis, levels = unique(ATAC_c_combined_tests$Analysis))

# Figure 4 B

ATAC_lengths_domain <-
  
  ggplot(ATAC_c_combined_tests, aes(y=Description, x=as.numeric(log_odds)))+
  geom_point(aes(color=Analysis, size=find_size))+
  expand_limits(x=c(-3, 3))+
  scale_x_continuous(breaks = seq(-3, 3, by = 1))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue", "aquamarine2", "darkcyan"))+
  
 
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  

  scale_alpha(name = paste("Transparency of bubble represents \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value < 0.05", "p value > 0.05"))+
  
  
  ggtitle(paste("Figure 4 B"))



ATAC_lengths_domain  <- 
  
  ATAC_lengths_domain  + labs(y="Chromatin Domains\n", x=paste('log2(Odds Ratio) for enrichment'))+
  theme(axis.title=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))


ATAC_lengths_domain



#---------------------------------------
# small RNA


# When we are considering the loci with chromatin domain annotations we can only consider the loci which map to Ahringer as this is where we get the domain annotations from. 

# The lists of length categorised epimutations must be restricted to genes which have Ahringer annotations 


small_RNA_background <- small_RNA_annotation_table[small_RNA_annotation_table$maps_to_Ahr==1, ]


all_Ahr_small_RNA_genes <- unique(small_RNA_background[, 2])

Ahr_small_RNA_not_inherit <- intersect(small_RNA_Muts_Not_Inherited, all_Ahr_small_RNA_genes)
Ahr_small_RNA_short <- intersect(small_RNA_Muts_Inherited_not_long, all_Ahr_small_RNA_genes)
Ahr_small_RNA_long <- intersect(small_RNA_Top_genes, all_Ahr_small_RNA_genes)


# small_RNA_background are genes which map to Ahringer

Active_targets <- unique(small_RNA_background[small_RNA_background$Gene_domain == "A", 2])

Regulated_targets <- unique(small_RNA_background[small_RNA_background$Gene_domain == "R", 2])

Chr_X_targets <- unique(small_RNA_background[small_RNA_background$Gene_domain == ".", 2])

piRNA_targets <- unique(small_RNA_background[small_RNA_background$piRNA_target==1, 2])

# the background data set therefore is length(unique(small_RNA_background$Target_gene))


small_RNA_test_list <- list(Ahr_small_RNA_not_inherit, Ahr_small_RNA_short, Ahr_small_RNA_long)

small_RNA_Features <- list(Regulated_targets, Active_targets, piRNA_targets, Chr_X_targets)




Non_Nested_small_RNA_OR <- c()

Non_Nested_small_RNA_pval <- c()


for(d in 1:length(small_RNA_test_list)){
  
  test <- small_RNA_test_list[[d]]
  
  temp_pval <- matrix(0, ncol = 1, nrow = 4)
  temp_oR <- matrix(0, ncol = 1, nrow = 4)
  
  for (i in 1:length(small_RNA_Features)) {
    common_in <- length(intersect(small_RNA_Features[[i]], test))
    Epimutation_total_in <- length(test)
    Feature_total_in <- length(small_RNA_Features[[i]])
    
    output <- epiIntersect(
      common = common_in,
      Epimutation_total = Epimutation_total_in,
      Feature_total = Feature_total_in,
      overall_total = 6535  # Total number of small RNA gene targets mapping to Ahringer
    )
    temp_pval[i, ] <- output[1]
    temp_oR[i, ] <- output[2]
  }
  
  
  
  temp_pval <- as.matrix(p.adjust(temp_pval, "bonferroni", 4))
  
  
  colnames(temp_pval) <-
    c("pval")
  rownames(temp_pval) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  
  colnames(temp_oR) <-
    c("oR")
  rownames(temp_oR) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  Non_Nested_small_RNA_OR <- cbind(Non_Nested_small_RNA_OR, temp_oR)
  
  Non_Nested_small_RNA_pval <- cbind(Non_Nested_small_RNA_pval, temp_pval)
  
}


colnames(Non_Nested_small_RNA_OR) <- c("Not_inherited", "Short", "Long")
colnames(Non_Nested_small_RNA_pval) <- c("Not_inherited", "Short", "Long")



#-----------------

# 22G-RNA epimutations distributions in distinct chromatin domains

ordered_sequence <- c("piRNA", "ChrX", "Regulated", "Active")

colnames(Non_Nested_small_RNA_OR) <- c("Test_1", "Test_2", "Test_3")

colnames(Non_Nested_small_RNA_pval) <- c("Test_1", "Test_2", "Test_3")


Non_Nested_small_RNA_OR <- as.data.frame(Non_Nested_small_RNA_OR[ordered_sequence, ])

Non_Nested_small_RNA_pval <- as.data.frame(Non_Nested_small_RNA_pval[ordered_sequence, ])



P_1 <- as.numeric(Non_Nested_small_RNA_pval$Test_1)

log_Odds_1 <- log2(as.numeric(Non_Nested_small_RNA_OR$Test_1))

neg_log_p_1 <- -log2(P_1)

bubble_table_1 <- cbind(rownames(Non_Nested_small_RNA_OR), Non_Nested_small_RNA_pval$Test_1, neg_log_p_1, log_Odds_1)

bubble_table_1 <- data.frame(rep("This is the analysis- NON INHERIT"), bubble_table_1)

bubble_table_1[, 2] <- factor(bubble_table_1[, 2], levels = bubble_table_1[, 2])

colnames(bubble_table_1) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")





P_2 <- as.numeric(Non_Nested_small_RNA_pval$Test_2)

log_Odds_2 <- log2(as.numeric(Non_Nested_small_RNA_OR$Test_2))

neg_log_p_2 <- -log2(P_2)

bubble_table_2 <- cbind(rownames(Non_Nested_small_RNA_OR), Non_Nested_small_RNA_pval$Test_2, neg_log_p_2, log_Odds_2)

bubble_table_2 <- data.frame(rep("This is the analysis- SHORT"), bubble_table_2)

bubble_table_2[, 2] <- factor(bubble_table_2[, 2], levels = bubble_table_2[, 2])

colnames(bubble_table_2) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")





P_3 <- as.numeric(Non_Nested_small_RNA_pval$Test_3)

log_Odds_3 <- log2(as.numeric(Non_Nested_small_RNA_OR$Test_3))

neg_log_p_3 <- -log2(P_3)

bubble_table_3 <- cbind(rownames(Non_Nested_small_RNA_OR), Non_Nested_small_RNA_pval$Test_3, neg_log_p_3, log_Odds_3)

bubble_table_3 <- data.frame(rep("This is the analysis- LONG"), bubble_table_3)

bubble_table_3[, 2] <- factor(bubble_table_3[, 2], levels = bubble_table_3[, 2])

colnames(bubble_table_3) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3)



find_size <- c()

for(i in 1:nrow(c_tests)){
  

  size_point <- 14
  
  if(as.numeric(c_tests[i, 3]) < 0.1){

    size_point <- 50
    
  }
  
  find_size <- c(find_size, size_point)
}


smallRNA_c_combined_tests <- cbind(c_tests, find_size)




guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")


smallRNA_c_combined_tests$Analysis <- factor(smallRNA_c_combined_tests$Analysis, levels = unique(smallRNA_c_combined_tests$Analysis))



# Figure 4 C 

# Bubble plot to show distribution of small RNAs in chromatin domains according to length


small_RNA_lengths_domain <-
  
  ggplot(smallRNA_c_combined_tests, aes(y=Description, x=as.numeric(log_odds)))+
  geom_point(aes(color=Analysis, size= find_size))+
  expand_limits(x=c(-3, 3))+
  scale_x_continuous(breaks = seq(-3, 3, by = 1))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue", "aquamarine2", "darkcyan"))+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  
  
  scale_alpha(name = paste("Transparency of bubble represents \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value < 0.1", "p value > 0.1"))+
  
  
  ggtitle(paste("Figure 4 C"))



small_RNA_lengths_domain <- 
  
  small_RNA_lengths_domain + labs(y="Chromatin Domains\n", x=paste('log2(Odds Ratio) for enrichment'))+
  theme(axis.title=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

small_RNA_lengths_domain

#------------------------------------------------------------------------------------------
# Supplementary Figure 5

# look at whether the UP/DOWN direction drives enrichment or depletions in different chromatin domains

# Approach: derive long-lived, short-lived and non-inherited categories again but from UP and DOWN epimutations separately 



library(lme4)


# RNA 

RNA_ratio <- (as.numeric(RNA_minus_suspect$length)/as.numeric(RNA_minus_suspect$number_transitions))

RNA_ratio[RNA_ratio %in% NaN] <- 0

RNA_tab <- data.frame(RNA_minus_suspect$gene, RNA_ratio, RNA_minus_suspect$is_up, RNA_minus_suspect$is_down)

colnames(RNA_tab) <- c("ID","meanLength", "up", "down")


# For UP

#filter so that only UP epimutations are considered.

UP_RNA_InputTab<-RNA_tab[which(RNA_tab[,3]>0),]


RNAModel<-lmer(meanLength~(1|ID),data=UP_RNA_InputTab)

lengthEst <-ranef(RNAModel)[[1]]


# K means clustering

Clusters <- kmeans(lengthEst, 2)

RNA_KMeans_Cluster_1 <- names(Clusters$cluster[Clusters$cluster == 1])

RNA_KMeans_Cluster_2 <- names(Clusters$cluster[Clusters$cluster == 2])



if(length(RNA_KMeans_Cluster_1) < length(RNA_KMeans_Cluster_2)){
  
  UP_RNA_longest_genes <- RNA_KMeans_Cluster_1
}else{UP_RNA_longest_genes <- RNA_KMeans_Cluster_2}




UP_RNA_Top_genes <- unique(All_lins_integrated_table[All_lins_integrated_table$gene %in% UP_RNA_longest_genes, 3])





# Get the different gene categories with RNA expression changes


# Now find the UP epimutations which have a length of 1

UP_epimutation_length_one <- unique(RNA_minus_suspect[RNA_minus_suspect$length == 1 & RNA_minus_suspect$is_up==1, 1])

# Now find the epimutations which have a length > 1

UP_epimutation_length_two_plus <- unique(RNA_minus_suspect[RNA_minus_suspect$length > 1 & RNA_minus_suspect$is_up==1, 1])

# now set difference between inherited and long lasting

UP_RNA_Muts_Inherited_not_long <- setdiff(UP_epimutation_length_two_plus, UP_RNA_Top_genes)

# Now set difference between non inherited and inherited not long

UP_RNA_Muts_Not_Inherited <- setdiff(UP_epimutation_length_one, UP_RNA_Muts_Inherited_not_long)

# Now set difference between non inherited and long lasting

UP_RNA_Muts_Not_Inherited <- setdiff(UP_RNA_Muts_Not_Inherited, UP_RNA_Top_genes)






# Get the gene names so we can intersect them with gene names known to be in different chromatin domains

UP_RNA_Muts_Inherited_not_long_genes <- c()

for(i in 1:length(UP_RNA_Muts_Inherited_not_long)){
  
  name <-  strsplit(UP_RNA_Muts_Inherited_not_long[[i]], ":")[[1]][4]
  UP_RNA_Muts_Inherited_not_long_genes <- c(UP_RNA_Muts_Inherited_not_long_genes, name)
  
}



UP_RNA_Muts_Not_Inherited_genes <- c()

for(i in 1:length(UP_RNA_Muts_Not_Inherited)){
  
  name <-  strsplit(UP_RNA_Muts_Not_Inherited[[i]], ":")[[1]][4]
  UP_RNA_Muts_Not_Inherited_genes <- c(UP_RNA_Muts_Not_Inherited_genes, name)
  
}




UP_RNA_longest_genes <- c()

for(i in 1:length(UP_RNA_Top_genes)){
  
  name <-  strsplit(UP_RNA_Top_genes[[i]], ":")[[1]][4]
  UP_RNA_longest_genes <- c(UP_RNA_longest_genes, name)
  
}




# For DOWN

DOWN_RNA_InputTab<-RNA_tab[which(RNA_tab[,4]>0),]

RNAModel<-lmer(meanLength~(1|ID),data=DOWN_RNA_InputTab)

lengthEst <-ranef(RNAModel)[[1]]


# K means clustering

Clusters <- kmeans(lengthEst, 2)

RNA_KMeans_Cluster_1 <- names(Clusters$cluster[Clusters$cluster == 1])

RNA_KMeans_Cluster_2 <- names(Clusters$cluster[Clusters$cluster == 2])



if(length(RNA_KMeans_Cluster_1) < length(RNA_KMeans_Cluster_2)){
  
  DOWN_RNA_longest_genes <- RNA_KMeans_Cluster_1
}else{DOWN_RNA_longest_genes <- RNA_KMeans_Cluster_2}



DOWN_RNA_Top_genes <- unique(All_lins_integrated_table[All_lins_integrated_table$gene %in% DOWN_RNA_longest_genes, 3])



# Get the different gene categories with RNA expression changes


# Now find the UP epimutations which have a length of 1

DOWN_epimutation_length_one <- unique(RNA_minus_suspect[RNA_minus_suspect$length == 1 & RNA_minus_suspect$is_down==1, 1])

# Now find the epimutations which have a length > 1

DOWN_epimutation_length_two_plus <- unique(RNA_minus_suspect[RNA_minus_suspect$length > 1 & RNA_minus_suspect$is_down==1, 1])

# now set difference between inherited and long lasting

DOWN_RNA_Muts_Inherited_not_long <- setdiff(DOWN_epimutation_length_two_plus, DOWN_RNA_Top_genes)

# Now set difference between non inherited and inherited not long

DOWN_RNA_Muts_Not_Inherited <- setdiff(DOWN_epimutation_length_one, DOWN_RNA_Muts_Inherited_not_long)

# Now set difference between non inherited and long lasting

DOWN_RNA_Muts_Not_Inherited <- setdiff(DOWN_RNA_Muts_Not_Inherited, DOWN_RNA_Top_genes)






# Get the gene names

DOWN_RNA_Muts_Inherited_not_long_genes <- c()

for(i in 1:length(DOWN_RNA_Muts_Inherited_not_long)){
  
  name <-  strsplit(DOWN_RNA_Muts_Inherited_not_long[[i]], ":")[[1]][4]
  DOWN_RNA_Muts_Inherited_not_long_genes <- c(DOWN_RNA_Muts_Inherited_not_long_genes, name)
  
}



DOWN_RNA_Muts_Not_Inherited_genes <- c()

for(i in 1:length(DOWN_RNA_Muts_Not_Inherited)){
  
  name <-  strsplit(DOWN_RNA_Muts_Not_Inherited[[i]], ":")[[1]][4]
  DOWN_RNA_Muts_Not_Inherited_genes <- c(DOWN_RNA_Muts_Not_Inherited_genes, name)
  
}




DOWN_RNA_longest_genes <- c()

for(i in 1:length(DOWN_RNA_Top_genes)){
  
  name <-  strsplit(DOWN_RNA_Top_genes[[i]], ":")[[1]][4]
  DOWN_RNA_longest_genes <- c(DOWN_RNA_longest_genes, name)
  
}





# When we are considering the genes with chromatin domain annotations we can only consider the genes which map to Ahringer as this is where we get the domain annotations from. 

all_Ahr_RNA_genes <- unique(RNA_centric_All_table[, 2])


# The lists of length categorised expression changes must be restricted to genes which have Ahringer annotations 


UP_Ahr_RNA_not_inherit <- intersect(UP_RNA_Muts_Not_Inherited_genes, all_Ahr_RNA_genes)
UP_Ahr_RNA_short <- intersect(UP_RNA_Muts_Inherited_not_long_genes, all_Ahr_RNA_genes)
UP_Ahr_RNA_long <- intersect(UP_RNA_longest_genes, all_Ahr_RNA_genes)



DOWN_Ahr_RNA_not_inherit <- intersect(DOWN_RNA_Muts_Not_Inherited_genes, all_Ahr_RNA_genes)
DOWN_Ahr_RNA_short <- intersect(DOWN_RNA_Muts_Inherited_not_long_genes, all_Ahr_RNA_genes)
DOWN_Ahr_RNA_long <- intersect(DOWN_RNA_longest_genes, all_Ahr_RNA_genes)




RNA_test_list <- list(UP_Ahr_RNA_not_inherit, UP_Ahr_RNA_short, UP_Ahr_RNA_long, 
                      DOWN_Ahr_RNA_not_inherit, DOWN_Ahr_RNA_short, DOWN_Ahr_RNA_long)

RNA_Features <- list(Regulated_genes, Active_genes, piRNA_cluster_genes, Chr_X_genes)






# Define the function

epiIntersect <-
  function(common,
           Epimutation_total,
           Feature_total,
           overall_total) {
    contingency_table <-
      rbind(c(common, Epimutation_total - common),
            c(
              Feature_total - common,
              overall_total - (Epimutation_total + Feature_total - common)
            ))
    FT_out <- fisher.test(contingency_table)
    pval <- FT_out$p.value
    OddsRatio <- FT_out$estimate
    return(c(pval, OddsRatio))
  }





UP_DOWN_Non_Nested_RNA_OR <- c()

UP_DOWN_Non_Nested_RNA_pval <- c()


for(d in 1:length(RNA_test_list)){
  
  test <- RNA_test_list[[d]]
  
  temp_pval <- matrix(0, ncol = 1, nrow = 4)
  temp_oR <- matrix(0, ncol = 1, nrow = 4)
  
  for (i in 1:length(RNA_Features)) {
    common_in <- length(intersect(RNA_Features[[i]], test))
    Epimutation_total_in <- length(test)
    Feature_total_in <- length(RNA_Features[[i]])
    
    
    if(common_in == 0){
      
      common_in <- common_in +1
      Epimutation_total_in <- Epimutation_total_in +1
      Feature_total_in <- Feature_total_in +1
    } 
    
    output <- epiIntersect(
      common = common_in,
      Epimutation_total = Epimutation_total_in,
      Feature_total = Feature_total_in,
      overall_total = 12617  # Total number of RNA genes mapping to Ahringer
    )
    temp_pval[i, ] <- output[1]
    temp_oR[i, ] <- output[2]
  }
  
  
  
  temp_pval <- as.matrix(p.adjust(temp_pval, "bonferroni", 4))
  
  
  colnames(temp_pval) <-
    c("pval")
  rownames(temp_pval) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  
  colnames(temp_oR) <-
    c("oR")
  rownames(temp_oR) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  UP_DOWN_Non_Nested_RNA_OR <- cbind(UP_DOWN_Non_Nested_RNA_OR, temp_oR)
  
  UP_DOWN_Non_Nested_RNA_pval <- cbind(UP_DOWN_Non_Nested_RNA_pval, temp_pval)
  
}


colnames(UP_DOWN_Non_Nested_RNA_OR) <- c("UP_Not_inherited", "UP_Short", "UP_Long", 
                                         "DOWN_Not_inherited", "DOWN_Short", "DOWN_Long")
colnames(UP_DOWN_Non_Nested_RNA_pval) <- c("UP_Not_inherited", "UP_Short", "UP_Long", 
                                           "DOWN_Not_inherited", "DOWN_Short", "DOWN_Long")



#---------------------------------------------------------------------------------------------------------------------
# Now for ATAC

# ATAC length categorisation 

ATAC_ABC_ep_ratio <- (as.numeric(ATAC_ep_minus_suspect$length)/as.numeric(ATAC_ep_minus_suspect$number_transitions))

ATAC_ABC_ep_ratio[ATAC_ABC_ep_ratio  %in% NaN] <- 0

ATAC_ABC_ep_tab <- data.frame(ATAC_ep_minus_suspect$names, ATAC_ABC_ep_ratio, ATAC_ep_minus_suspect$is_up, ATAC_ep_minus_suspect$is_down)

colnames(ATAC_ABC_ep_tab) <- c("ID","meanLength","up", "down")




# For UP

# Just select the coordinates which have UP epimutations

UP_ATAC_ABC_ep_tab <- ATAC_ABC_ep_tab[which(ATAC_ABC_ep_tab[,3]>0),]

# the intercept (shown by '1' can vary by ID)

UP_ATACModel<-lmer(meanLength~(1|ID),data=UP_ATAC_ABC_ep_tab) 

#this makes the mean length depend on the particular gene as a random factor.  

# IE does the identity of the gene influence how long it lasts for on average?

#for each gene an estimate of the influence that it makes on the mean length is generated.  This information is what we want.
#to extract this:

lengthEst <-ranef(UP_ATACModel)[[1]]



# K means clustering

Clusters<-kmeans(lengthEst,2)

ATAC_KMeans_Cluster_1 <- names(Clusters$cluster[Clusters$cluster == 1])

ATAC_KMeans_Cluster_2 <- names(Clusters$cluster[Clusters$cluster == 2])


if(length(ATAC_KMeans_Cluster_1) < length(ATAC_KMeans_Cluster_2)){
  
  UP_ATAC_ep_ABC_Top <- ATAC_KMeans_Cluster_1
}else{UP_ATAC_ep_ABC_Top <- ATAC_KMeans_Cluster_2}






# Get the categories of genes with chromatin epimutations

# find never epimutated genes

never_epimutated <- unique(ATAC_ep_All[ATAC_ep_All$length ==0, 1])


# Now find the epimutations which have a length of 1

UP_epimutation_length_one <- unique(ATAC_ep_minus_suspect[ATAC_ep_minus_suspect$length == 1 & ATAC_ep_minus_suspect$is_up==1, 1])


# Now find the epimutations which have a length > 1

UP_epimutation_length_two_plus <- unique(ATAC_ep_minus_suspect[ATAC_ep_minus_suspect$length > 1 & ATAC_ep_minus_suspect$is_up==1, 1])


# now set difference between inherited and long lasting

UP_ATAC_ep_Muts_Inherited_not_long <- setdiff(UP_epimutation_length_two_plus, UP_ATAC_ep_ABC_Top)


# Now set difference between non inherited and inherited not long & long

UP_ATAC_ep_Muts_Not_Inherited <- setdiff(UP_epimutation_length_one, c(UP_ATAC_ep_Muts_Inherited_not_long, UP_ATAC_ep_ABC_Top))





# For DOWN

# Just select the coordinates which have DOWN epimutations

DOWN_ATAC_ABC_ep_tab <- ATAC_ABC_ep_tab[which(ATAC_ABC_ep_tab[,4]>0),]

# the intercept (shown by '1' can vary by ID)

DOWN_ATACModel<-lmer(meanLength~(1|ID),data=DOWN_ATAC_ABC_ep_tab) 

#this makes the mean length depend on the particular gene as a random factor.  

# IE does the identity of the gene influence how long it lasts for on average?

#for each gene an estimate of the influence that it makes on the mean length is generated.  This information is what we want.
#to extract this:

lengthEst <-ranef(DOWN_ATACModel)[[1]]




# K means clustering

Clusters<-kmeans(lengthEst,2)

ATAC_KMeans_Cluster_1 <- names(Clusters$cluster[Clusters$cluster == 1])

ATAC_KMeans_Cluster_2 <- names(Clusters$cluster[Clusters$cluster == 2])


if(length(ATAC_KMeans_Cluster_1) < length(ATAC_KMeans_Cluster_2)){
  
  DOWN_ATAC_ep_ABC_Top <- ATAC_KMeans_Cluster_1
}else{DOWN_ATAC_ep_ABC_Top <- ATAC_KMeans_Cluster_2}






# Get the categories of genes with chromatin epimutations


# Now find the epimutations which have a length of 1

DOWN_epimutation_length_one <- unique(ATAC_ep_minus_suspect[ATAC_ep_minus_suspect$length == 1 & ATAC_ep_minus_suspect$is_down==1, 1])


# Now find the epimutations which have a length > 1

DOWN_epimutation_length_two_plus <- unique(ATAC_ep_minus_suspect[ATAC_ep_minus_suspect$length > 1 & ATAC_ep_minus_suspect$is_down==1, 1])


# now set difference between inherited and long lasting

DOWN_ATAC_ep_Muts_Inherited_not_long <- setdiff(DOWN_epimutation_length_two_plus, DOWN_ATAC_ep_ABC_Top)


# Now set difference between non inherited and inherited not long & long

DOWN_ATAC_ep_Muts_Not_Inherited <- setdiff(DOWN_epimutation_length_one, c(DOWN_ATAC_ep_Muts_Inherited_not_long, DOWN_ATAC_ep_ABC_Top))






ATAC_test_list <- list(UP_ATAC_ep_Muts_Not_Inherited, UP_ATAC_ep_Muts_Inherited_not_long, UP_ATAC_ep_ABC_Top, 
                       DOWN_ATAC_ep_Muts_Not_Inherited, DOWN_ATAC_ep_Muts_Inherited_not_long, DOWN_ATAC_ep_ABC_Top)

ep_Features <- list(ep_Active_domains, ep_Regulated_domains, ep_piRNA_domains, ep_X_domains)



UP_DOWN_Non_Nested_ATAC_OR <- c()

UP_DOWN_Non_Nested_ATAC_pval <- c()


# For the background use the number of loci in the binarised ATAC table as these all have ahringer annotations 

# nrow(ATAC_binarised_table_A_2.25_ehancer_promoter)



for(d in 1:length(ATAC_test_list)){
  
  test <- ATAC_test_list[[d]]
  
  temp_pval <- matrix(0, ncol = 1, nrow = 4)
  temp_oR <- matrix(0, ncol = 1, nrow = 4)
  
  for (i in 1:length(ep_Features)) {
    common_in <- length(intersect(ep_Features[[i]], test))
    Epimutation_total_in <- length(test)
    Feature_total_in <- length(ep_Features[[i]])
    
    if(common_in == 0){
      
      common_in <- common_in +1
      Epimutation_total_in <- Epimutation_total_in +1
      Feature_total_in <- Feature_total_in +1
    } 
    
    output <- epiIntersect(
      common = common_in,
      Epimutation_total = Epimutation_total_in,
      Feature_total = Feature_total_in,
      overall_total = 31618
    )
    temp_pval[i, ] <- output[1]
    temp_oR[i, ] <- output[2]
  }
  
  # Bonferroni adjust:
  
  temp_pval <- as.matrix(p.adjust(temp_pval, "bonferroni", 4))
  
  
  colnames(temp_pval) <-
    c("pval")
  rownames(temp_pval) <-
    c("Active", "Regulated","piRNA", "ChrX")
  
  
  colnames(temp_oR) <-
    c("oR")
  rownames(temp_oR) <-
    c("Active", "Regulated", "piRNA", "ChrX")
  
  UP_DOWN_Non_Nested_ATAC_OR <- cbind(UP_DOWN_Non_Nested_ATAC_OR, temp_oR)
  
  UP_DOWN_Non_Nested_ATAC_pval <- cbind(UP_DOWN_Non_Nested_ATAC_pval, temp_pval)
  
}


colnames(UP_DOWN_Non_Nested_ATAC_OR) <- c("UP_Not_inherited", "UP_Short", "UP_Long", 
                                          "DOWN_Not_inherited", "DOWN_Short", "DOWN_Long")

colnames(UP_DOWN_Non_Nested_ATAC_pval) <- c("UP_Not_inherited", "UP_Short", "UP_Long", 
                                            "DOWN_Not_inherited", "DOWN_Short", "DOWN_Long")





# For 22G RNAs


# For UP 


small_RNA_ratio <- (as.numeric(small_RNA_minus_suspect_muts$length)/as.numeric(small_RNA_minus_suspect_muts$number_transitions))

small_RNA_tab <- data.frame(small_RNA_minus_suspect_muts$gene_name, small_RNA_ratio, small_RNA_minus_suspect_muts$is_up, small_RNA_minus_suspect_muts$is_down)

colnames(small_RNA_tab) <- c("ID","meanLength", "up", "down")



#filter so that only UP epimutations are considered.
UP_small_RNA_InputTab<-small_RNA_tab[which(small_RNA_tab[,3]>0),]

smallRNAModel<-lmer(meanLength~(1|ID),data=UP_small_RNA_InputTab)



lengthEst <-ranef(smallRNAModel)[[1]]


# K means clustering

Clusters<-kmeans(lengthEst,2)

smallRNA_KMeans_Cluster_1 <- names(Clusters$cluster[Clusters$cluster == 1])

smallRNA_KMeans_Cluster_2 <- names(Clusters$cluster[Clusters$cluster == 2])


if(length(smallRNA_KMeans_Cluster_1) < length(smallRNA_KMeans_Cluster_2)){
  
  UP_small_RNA_Top_genes <- smallRNA_KMeans_Cluster_1
}else{UP_small_RNA_Top_genes <- smallRNA_KMeans_Cluster_2}


# Now find the epimutations which have a length of 1
UP_epimutation_length_one <- unique(small_RNA_minus_suspect[small_RNA_minus_suspect$length == 1 & small_RNA_minus_suspect$is_up==1, 1])

# Now find the epimutations which have a length > 1
UP_epimutation_length_two_plus <- unique(small_RNA_minus_suspect[small_RNA_minus_suspect$length > 1 & small_RNA_minus_suspect$is_up==1, 1])

# now set difference between inherited and long lasting
UP_small_RNA_Muts_Inherited_not_long <- setdiff(UP_epimutation_length_two_plus, UP_small_RNA_Top_genes)

# Now set difference between non inherited and inherited not long
UP_small_RNA_Muts_Not_Inherited <- setdiff(UP_epimutation_length_one, UP_small_RNA_Muts_Inherited_not_long)

# Now set difference between non inherited and long lasting
UP_small_RNA_Muts_Not_Inherited <- setdiff(UP_small_RNA_Muts_Not_Inherited, UP_small_RNA_Top_genes)


UP_Ahr_small_RNA_not_inherit <- intersect(UP_small_RNA_Muts_Not_Inherited, all_Ahr_small_RNA_genes)
UP_Ahr_small_RNA_short <- intersect(UP_small_RNA_Muts_Inherited_not_long, all_Ahr_small_RNA_genes)
UP_Ahr_small_RNA_long <- intersect(UP_small_RNA_Top_genes, all_Ahr_small_RNA_genes)




# For DOWN

#filter so that only DOWN epimutations are considered.

DOWN_small_RNA_InputTab<-small_RNA_tab[which(small_RNA_tab[,4]>0),]

smallRNAModel<-lmer(meanLength~(1|ID),data=DOWN_small_RNA_InputTab)



lengthEst <-ranef(smallRNAModel)[[1]]


# K means clustering

Clusters<-kmeans(lengthEst,2)

smallRNA_KMeans_Cluster_1 <- names(Clusters$cluster[Clusters$cluster == 1])

smallRNA_KMeans_Cluster_2 <- names(Clusters$cluster[Clusters$cluster == 2])


if(length(smallRNA_KMeans_Cluster_1) < length(smallRNA_KMeans_Cluster_2)){
  
  DOWN_small_RNA_Top_genes <- smallRNA_KMeans_Cluster_1
}else{DOWN_small_RNA_Top_genes <- smallRNA_KMeans_Cluster_2}


# Now find the epimutations which have a length of 1
DOWN_epimutation_length_one <- unique(small_RNA_minus_suspect[small_RNA_minus_suspect$length == 1 & small_RNA_minus_suspect$is_down==1, 1])

# Now find the epimutations which have a length > 1
DOWN_epimutation_length_two_plus <- unique(small_RNA_minus_suspect[small_RNA_minus_suspect$length > 1 & small_RNA_minus_suspect$is_down==1, 1])

# now set difference between inherited and long lasting
DOWN_small_RNA_Muts_Inherited_not_long <- setdiff(DOWN_epimutation_length_two_plus, DOWN_small_RNA_Top_genes)

# Now set difference between non inherited and inherited not long
DOWN_small_RNA_Muts_Not_Inherited <- setdiff(DOWN_epimutation_length_one, DOWN_small_RNA_Muts_Inherited_not_long)

# Now set difference between non inherited and long lasting
DOWN_small_RNA_Muts_Not_Inherited <- setdiff(DOWN_small_RNA_Muts_Not_Inherited, DOWN_small_RNA_Top_genes)


DOWN_Ahr_small_RNA_not_inherit <- intersect(DOWN_small_RNA_Muts_Not_Inherited, all_Ahr_small_RNA_genes)
DOWN_Ahr_small_RNA_short <- intersect(DOWN_small_RNA_Muts_Inherited_not_long, all_Ahr_small_RNA_genes)
DOWN_Ahr_small_RNA_long <- intersect(DOWN_small_RNA_Top_genes, all_Ahr_small_RNA_genes)



small_RNA_test_list <- list(UP_Ahr_small_RNA_not_inherit, UP_Ahr_small_RNA_short, UP_Ahr_small_RNA_long, 
                            DOWN_Ahr_small_RNA_not_inherit, DOWN_Ahr_small_RNA_short, DOWN_Ahr_small_RNA_long)

small_RNA_Features <- list(Regulated_targets, Active_targets, piRNA_targets, Chr_X_targets)




UP_DOWN_Non_Nested_small_RNA_OR <- c()

UP_DOWN_Non_Nested_small_RNA_pval <- c()


for(d in 1:length(small_RNA_test_list)){
  
  test <- small_RNA_test_list[[d]]
  
  temp_pval <- matrix(0, ncol = 1, nrow = 4)
  temp_oR <- matrix(0, ncol = 1, nrow = 4)
  
  for (i in 1:length(small_RNA_Features)) {
    common_in <- length(intersect(small_RNA_Features[[i]], test))
    Epimutation_total_in <- length(test)
    Feature_total_in <- length(small_RNA_Features[[i]])
    
    if(common_in == 0){
      
      common_in <- common_in +1
      Epimutation_total_in <- Epimutation_total_in +1
      Feature_total_in <- Feature_total_in +1
    } 
    
    output <- epiIntersect(
      common = common_in,
      Epimutation_total = Epimutation_total_in,
      Feature_total = Feature_total_in,
      overall_total = 6535  # Total number of small RNA gene targets mapping to Ahringer
    )
    temp_pval[i, ] <- output[1]
    temp_oR[i, ] <- output[2]
  }
  
  
  
  temp_pval <- as.matrix(p.adjust(temp_pval, "bonferroni", 4))
  
  
  colnames(temp_pval) <-
    c("pval")
  rownames(temp_pval) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  
  colnames(temp_oR) <-
    c("oR")
  rownames(temp_oR) <-
    c("Regulated", "Active", "piRNA", "ChrX")
  
  UP_DOWN_Non_Nested_small_RNA_OR <- cbind(UP_DOWN_Non_Nested_small_RNA_OR, temp_oR)
  
  UP_DOWN_Non_Nested_small_RNA_pval <- cbind(UP_DOWN_Non_Nested_small_RNA_pval, temp_pval)
  
}


colnames(UP_DOWN_Non_Nested_small_RNA_OR) <- c("UP_Not_inherited", "UP_Short", "UP_Long", 
                                               "DOWN_Not_inherited", "DOWN_Short", "DOWN_Long")

colnames(UP_DOWN_Non_Nested_small_RNA_pval) <- c("UP_Not_inherited", "UP_Short", "UP_Long", 
                                                 "DOWN_Not_inherited", "DOWN_Short", "DOWN_Long")






# Make a bubble plot for RNA 


ordered_sequence <- c("piRNA", "ChrX", "Regulated", "Active")

colnames(UP_DOWN_Non_Nested_RNA_OR) <- c("UP_non", "UP_short", "UP_long", "DOWN_non", "DOWN_short", "DOWN_long")

colnames(UP_DOWN_Non_Nested_RNA_pval) <- c("UP_non", "UP_short", "UP_long", "DOWN_non", "DOWN_short", "DOWN_long")


UP_DOWN_Non_Nested_RNA_OR <- as.data.frame(UP_DOWN_Non_Nested_RNA_OR[ordered_sequence, ])

UP_DOWN_Non_Nested_RNA_pval <- as.data.frame(UP_DOWN_Non_Nested_RNA_pval[ordered_sequence, ])



P_1 <- as.numeric(UP_DOWN_Non_Nested_RNA_pval$UP_non)

log_Odds_1 <- log2(as.numeric(UP_DOWN_Non_Nested_RNA_OR$UP_non))

neg_log_p_1 <- -log2(P_1)

bubble_table_1 <- cbind(rownames(UP_DOWN_Non_Nested_RNA_OR), UP_DOWN_Non_Nested_RNA_pval$UP_non, neg_log_p_1, log_Odds_1)

bubble_table_1 <- data.frame(rep("This is the analysis- UP NON INHERIT"), bubble_table_1)

bubble_table_1[, 2] <- factor(bubble_table_1[, 2], levels = bubble_table_1[, 2])

colnames(bubble_table_1) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_2 <- as.numeric(UP_DOWN_Non_Nested_RNA_pval$UP_short)

log_Odds_2 <- log2(as.numeric(UP_DOWN_Non_Nested_RNA_OR$UP_short))

neg_log_p_2 <- -log2(P_2)

bubble_table_2 <- cbind(rownames(UP_DOWN_Non_Nested_RNA_OR), UP_DOWN_Non_Nested_RNA_pval$UP_short, neg_log_p_2, log_Odds_2)

bubble_table_2 <- data.frame(rep("This is the analysis- UP SHORT"), bubble_table_2)

bubble_table_2[, 2] <- factor(bubble_table_2[, 2], levels = bubble_table_2[, 2])

colnames(bubble_table_2) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_3 <- as.numeric(UP_DOWN_Non_Nested_RNA_pval$UP_long)

log_Odds_3 <- log2(as.numeric(UP_DOWN_Non_Nested_RNA_OR$UP_long))

neg_log_p_3 <- -log2(P_3)

bubble_table_3 <- cbind(rownames(UP_DOWN_Non_Nested_RNA_OR), UP_DOWN_Non_Nested_RNA_pval$UP_long, neg_log_p_3, log_Odds_3)

bubble_table_3 <- data.frame(rep("This is the analysis- UP LONG"), bubble_table_3)

bubble_table_3[, 2] <- factor(bubble_table_3[, 2], levels = bubble_table_3[, 2])

colnames(bubble_table_3) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_4 <- as.numeric(UP_DOWN_Non_Nested_RNA_pval$DOWN_non)

log_Odds_4 <- log2(as.numeric(UP_DOWN_Non_Nested_RNA_OR$DOWN_non))

neg_log_p_4 <- -log2(P_4)

bubble_table_4 <- cbind(rownames(UP_DOWN_Non_Nested_RNA_OR), UP_DOWN_Non_Nested_RNA_pval$DOWN_non, neg_log_p_4, log_Odds_4)

bubble_table_4 <- data.frame(rep("This is the analysis- DOWN NON INHERIT"), bubble_table_4)

bubble_table_4[, 2] <- factor(bubble_table_4[, 2], levels = bubble_table_4[, 2])

colnames(bubble_table_4) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_5 <- as.numeric(UP_DOWN_Non_Nested_RNA_pval$DOWN_short)

log_Odds_5 <- log2(as.numeric(UP_DOWN_Non_Nested_RNA_OR$DOWN_short))

neg_log_p_5 <- -log2(P_5)

bubble_table_5 <- cbind(rownames(UP_DOWN_Non_Nested_RNA_OR), UP_DOWN_Non_Nested_RNA_pval$DOWN_short, neg_log_p_5, log_Odds_5)

bubble_table_5 <- data.frame(rep("This is the analysis- DOWN SHORT"), bubble_table_5)

bubble_table_5[, 2] <- factor(bubble_table_5[, 2], levels = bubble_table_5[, 2])

colnames(bubble_table_5) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_6 <- as.numeric(UP_DOWN_Non_Nested_RNA_pval$DOWN_long)

log_Odds_6 <- log2(as.numeric(UP_DOWN_Non_Nested_RNA_OR$DOWN_long))

neg_log_p_6 <- -log2(P_6)

bubble_table_6 <- cbind(rownames(UP_DOWN_Non_Nested_RNA_OR), UP_DOWN_Non_Nested_RNA_pval$DOWN_long, neg_log_p_6, log_Odds_6)

bubble_table_6 <- data.frame(rep("This is the analysis- DOWN LONG"), bubble_table_6)

bubble_table_6[, 2] <- factor(bubble_table_6[, 2], levels = bubble_table_6[, 2])

colnames(bubble_table_6) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")




c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3, 
                 bubble_table_4, bubble_table_5, bubble_table_6)



find_size <- c()

for(i in 1:nrow(c_tests)){
  
  
  size_point <- 14
  
  if(as.numeric(c_tests[i, 3]) < 0.1){
    
    size_point <- 50
    
  }
  
  find_size <- c(find_size, size_point)
}


UP_DOWN_RNA_c_combined_tests <- cbind(c_tests, find_size)




guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")


UP_DOWN_RNA_c_combined_tests$Analysis <- factor(UP_DOWN_RNA_c_combined_tests$Analysis, levels = unique(UP_DOWN_RNA_c_combined_tests$Analysis))



# Supplemental figure A

# Bubble plot to show distribution of RNA based expression changes in chromatin domains according to length & direction


UP_DOWN_RNA_lengths_domain <-
  
  ggplot(UP_DOWN_RNA_c_combined_tests, aes(y=Description, x=as.numeric(log_odds)))+
  geom_point(aes(color=Analysis, size= find_size))+
  expand_limits(x=c(-4, 4))+
  scale_x_continuous(breaks = seq(-4, 4, by = 1))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue", "forestgreen", "magenta4", "darkgoldenrod2", "maroon1", "navy" ))+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  
  
  scale_alpha(name = paste("Transparency of bubble represents \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value < 0.1", "p value > 0.1"))+
  
  
  ggtitle(paste("RNA UP/DOWN/Domain"))



UP_DOWN_RNA_lengths_domain  <- 
  
  UP_DOWN_RNA_lengths_domain + labs(y="Chromatin Domains\n", x=paste('log2(Odds Ratio) for enrichment'))+
  theme(axis.title=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

UP_DOWN_RNA_lengths_domain 






# Make a bubble plot for ATAC


ordered_sequence <- c("piRNA", "ChrX", "Regulated", "Active")

colnames(UP_DOWN_Non_Nested_ATAC_OR) <- c("UP_non", "UP_short", "UP_long", "DOWN_non", "DOWN_short", "DOWN_long")

colnames(UP_DOWN_Non_Nested_ATAC_pval) <- c("UP_non", "UP_short", "UP_long", "DOWN_non", "DOWN_short", "DOWN_long")


UP_DOWN_Non_Nested_ATAC_OR <- as.data.frame(UP_DOWN_Non_Nested_ATAC_OR[ordered_sequence, ])

UP_DOWN_Non_Nested_ATAC_pval <- as.data.frame(UP_DOWN_Non_Nested_ATAC_pval[ordered_sequence, ])



P_1 <- as.numeric(UP_DOWN_Non_Nested_ATAC_pval$UP_non)

log_Odds_1 <- log2(as.numeric(UP_DOWN_Non_Nested_ATAC_OR$UP_non))

neg_log_p_1 <- -log2(P_1)

bubble_table_1 <- cbind(rownames(UP_DOWN_Non_Nested_ATAC_OR), UP_DOWN_Non_Nested_ATAC_pval$UP_non, neg_log_p_1, log_Odds_1)

bubble_table_1 <- data.frame(rep("This is the analysis- UP NON INHERIT"), bubble_table_1)

bubble_table_1[, 2] <- factor(bubble_table_1[, 2], levels = bubble_table_1[, 2])

colnames(bubble_table_1) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_2 <- as.numeric(UP_DOWN_Non_Nested_ATAC_pval$UP_short)

log_Odds_2 <- log2(as.numeric(UP_DOWN_Non_Nested_ATAC_OR$UP_short))

neg_log_p_2 <- -log2(P_2)

bubble_table_2 <- cbind(rownames(UP_DOWN_Non_Nested_ATAC_OR), UP_DOWN_Non_Nested_ATAC_pval$UP_short, neg_log_p_2, log_Odds_2)

bubble_table_2 <- data.frame(rep("This is the analysis- UP SHORT"), bubble_table_2)

bubble_table_2[, 2] <- factor(bubble_table_2[, 2], levels = bubble_table_2[, 2])

colnames(bubble_table_2) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_3 <- as.numeric(UP_DOWN_Non_Nested_ATAC_pval$UP_long)

log_Odds_3 <- log2(as.numeric(UP_DOWN_Non_Nested_ATAC_OR$UP_long))

neg_log_p_3 <- -log2(P_3)

bubble_table_3 <- cbind(rownames(UP_DOWN_Non_Nested_ATAC_OR), UP_DOWN_Non_Nested_ATAC_pval$UP_long, neg_log_p_3, log_Odds_3)

bubble_table_3 <- data.frame(rep("This is the analysis- UP LONG"), bubble_table_3)

bubble_table_3[, 2] <- factor(bubble_table_3[, 2], levels = bubble_table_3[, 2])

colnames(bubble_table_3) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_4 <- as.numeric(UP_DOWN_Non_Nested_ATAC_pval$DOWN_non)

log_Odds_4 <- log2(as.numeric(UP_DOWN_Non_Nested_ATAC_OR$DOWN_non))

neg_log_p_4 <- -log2(P_4)

bubble_table_4 <- cbind(rownames(UP_DOWN_Non_Nested_ATAC_OR), UP_DOWN_Non_Nested_ATAC_pval$DOWN_non, neg_log_p_4, log_Odds_4)

bubble_table_4 <- data.frame(rep("This is the analysis- DOWN NON INHERIT"), bubble_table_4)

bubble_table_4[, 2] <- factor(bubble_table_4[, 2], levels = bubble_table_4[, 2])

colnames(bubble_table_4) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_5 <- as.numeric(UP_DOWN_Non_Nested_ATAC_pval$DOWN_short)

log_Odds_5 <- log2(as.numeric(UP_DOWN_Non_Nested_ATAC_OR$DOWN_short))

neg_log_p_5 <- -log2(P_5)

bubble_table_5 <- cbind(rownames(UP_DOWN_Non_Nested_ATAC_OR), UP_DOWN_Non_Nested_ATAC_pval$DOWN_short, neg_log_p_5, log_Odds_5)

bubble_table_5 <- data.frame(rep("This is the analysis- DOWN SHORT"), bubble_table_5)

bubble_table_5[, 2] <- factor(bubble_table_5[, 2], levels = bubble_table_5[, 2])

colnames(bubble_table_5) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_6 <- as.numeric(UP_DOWN_Non_Nested_ATAC_pval$DOWN_long)

log_Odds_6 <- log2(as.numeric(UP_DOWN_Non_Nested_ATAC_OR$DOWN_long))

neg_log_p_6 <- -log2(P_6)

bubble_table_6 <- cbind(rownames(UP_DOWN_Non_Nested_ATAC_OR), UP_DOWN_Non_Nested_ATAC_pval$DOWN_long, neg_log_p_6, log_Odds_6)

bubble_table_6 <- data.frame(rep("This is the analysis- DOWN LONG"), bubble_table_6)

bubble_table_6[, 2] <- factor(bubble_table_6[, 2], levels = bubble_table_6[, 2])

colnames(bubble_table_6) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")




c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3, 
                 bubble_table_4, bubble_table_5, bubble_table_6)



find_size <- c()

for(i in 1:nrow(c_tests)){
  
  
  size_point <- 14
  
  if(as.numeric(c_tests[i, 3]) < 0.1){
    
    size_point <- 50
    
  }
  
  find_size <- c(find_size, size_point)
}


UP_DOWN_ATAC_c_combined_tests <- cbind(c_tests, find_size)




guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")


UP_DOWN_ATAC_c_combined_tests$Analysis <- factor(UP_DOWN_ATAC_c_combined_tests$Analysis, levels = unique(UP_DOWN_ATAC_c_combined_tests$Analysis))



# Supplemental figure A

# Bubble plot to show distribution of RNA based expression changes in chromatin domains according to length & direction


UP_DOWN_ATAC_lengths_domain <-
  
  ggplot(UP_DOWN_ATAC_c_combined_tests, aes(y=Description, x=as.numeric(log_odds)))+
  geom_point(aes(color=Analysis, size= find_size))+
  expand_limits(x=c(-4, 4))+
  scale_x_continuous(breaks = seq(-4, 4, by = 1))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue", "forestgreen", "magenta4", "darkgoldenrod2", "maroon1", "navy" ))+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  
  
  scale_alpha(name = paste("Transparency of bubble represents \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value < 0.1", "p value > 0.1"))+
  
  
  ggtitle(paste("ATAC UP/DOWN/Domain"))



UP_DOWN_ATAC_lengths_domain  <- 
  
  UP_DOWN_ATAC_lengths_domain + labs(y="Chromatin Domains\n", x=paste('log2(Odds Ratio) for enrichment'))+
  theme(axis.title=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

UP_DOWN_ATAC_lengths_domain 






# Make a bubble plot for small RNA 


ordered_sequence <- c("piRNA", "ChrX", "Regulated", "Active")

colnames(UP_DOWN_Non_Nested_small_RNA_OR) <- c("UP_non", "UP_short", "UP_long", "DOWN_non", "DOWN_short", "DOWN_long")

colnames(UP_DOWN_Non_Nested_small_RNA_pval) <- c("UP_non", "UP_short", "UP_long", "DOWN_non", "DOWN_short", "DOWN_long")


UP_DOWN_Non_Nested_small_RNA_OR <- as.data.frame(UP_DOWN_Non_Nested_small_RNA_OR[ordered_sequence, ])

UP_DOWN_Non_Nested_small_RNA_pval <- as.data.frame(UP_DOWN_Non_Nested_small_RNA_pval[ordered_sequence, ])



P_1 <- as.numeric(UP_DOWN_Non_Nested_small_RNA_pval$UP_non)

log_Odds_1 <- log2(as.numeric(UP_DOWN_Non_Nested_small_RNA_OR$UP_non))

neg_log_p_1 <- -log2(P_1)

bubble_table_1 <- cbind(rownames(UP_DOWN_Non_Nested_small_RNA_OR), UP_DOWN_Non_Nested_small_RNA_pval$UP_non, neg_log_p_1, log_Odds_1)

bubble_table_1 <- data.frame(rep("This is the analysis- UP NON INHERIT"), bubble_table_1)

bubble_table_1[, 2] <- factor(bubble_table_1[, 2], levels = bubble_table_1[, 2])

colnames(bubble_table_1) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_2 <- as.numeric(UP_DOWN_Non_Nested_small_RNA_pval$UP_short)

log_Odds_2 <- log2(as.numeric(UP_DOWN_Non_Nested_small_RNA_OR$UP_short))

neg_log_p_2 <- -log2(P_2)

bubble_table_2 <- cbind(rownames(UP_DOWN_Non_Nested_small_RNA_OR), UP_DOWN_Non_Nested_small_RNA_pval$UP_short, neg_log_p_2, log_Odds_2)

bubble_table_2 <- data.frame(rep("This is the analysis- UP SHORT"), bubble_table_2)

bubble_table_2[, 2] <- factor(bubble_table_2[, 2], levels = bubble_table_2[, 2])

colnames(bubble_table_2) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_3 <- as.numeric(UP_DOWN_Non_Nested_small_RNA_pval$UP_long)

log_Odds_3 <- log2(as.numeric(UP_DOWN_Non_Nested_small_RNA_OR$UP_long))

neg_log_p_3 <- -log2(P_3)

bubble_table_3 <- cbind(rownames(UP_DOWN_Non_Nested_small_RNA_OR), UP_DOWN_Non_Nested_small_RNA_pval$UP_long, neg_log_p_3, log_Odds_3)

bubble_table_3 <- data.frame(rep("This is the analysis- UP LONG"), bubble_table_3)

bubble_table_3[, 2] <- factor(bubble_table_3[, 2], levels = bubble_table_3[, 2])

colnames(bubble_table_3) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_4 <- as.numeric(UP_DOWN_Non_Nested_small_RNA_pval$DOWN_non)

log_Odds_4 <- log2(as.numeric(UP_DOWN_Non_Nested_small_RNA_OR$DOWN_non))

neg_log_p_4 <- -log2(P_4)

bubble_table_4 <- cbind(rownames(UP_DOWN_Non_Nested_small_RNA_OR), UP_DOWN_Non_Nested_small_RNA_pval$DOWN_non, neg_log_p_4, log_Odds_4)

bubble_table_4 <- data.frame(rep("This is the analysis- DOWN NON INHERIT"), bubble_table_4)

bubble_table_4[, 2] <- factor(bubble_table_4[, 2], levels = bubble_table_4[, 2])

colnames(bubble_table_4) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_5 <- as.numeric(UP_DOWN_Non_Nested_small_RNA_pval$DOWN_short)

log_Odds_5 <- log2(as.numeric(UP_DOWN_Non_Nested_small_RNA_OR$DOWN_short))

neg_log_p_5 <- -log2(P_5)

bubble_table_5 <- cbind(rownames(UP_DOWN_Non_Nested_small_RNA_OR), UP_DOWN_Non_Nested_small_RNA_pval$DOWN_short, neg_log_p_5, log_Odds_5)

bubble_table_5 <- data.frame(rep("This is the analysis- DOWN SHORT"), bubble_table_5)

bubble_table_5[, 2] <- factor(bubble_table_5[, 2], levels = bubble_table_5[, 2])

colnames(bubble_table_5) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")



P_6 <- as.numeric(UP_DOWN_Non_Nested_small_RNA_pval$DOWN_long)

log_Odds_6 <- log2(as.numeric(UP_DOWN_Non_Nested_small_RNA_OR$DOWN_long))

neg_log_p_6 <- -log2(P_6)

bubble_table_6 <- cbind(rownames(UP_DOWN_Non_Nested_small_RNA_OR), UP_DOWN_Non_Nested_small_RNA_pval$DOWN_long, neg_log_p_6, log_Odds_6)

bubble_table_6 <- data.frame(rep("This is the analysis- DOWN LONG"), bubble_table_6)

bubble_table_6[, 2] <- factor(bubble_table_6[, 2], levels = bubble_table_6[, 2])

colnames(bubble_table_6) <- c("Analysis",  "Description", "p_val", "neg_log_P", "log_odds")




c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3, 
                 bubble_table_4, bubble_table_5, bubble_table_6)



find_size <- c()

for(i in 1:nrow(c_tests)){
  
  
  size_point <- 14
  
  if(as.numeric(c_tests[i, 3]) < 0.1){
    
    size_point <- 50
    
  }
  
  find_size <- c(find_size, size_point)
}


UP_DOWN_small_RNA_c_combined_tests <- cbind(c_tests, find_size)




guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")


UP_DOWN_small_RNA_c_combined_tests$Analysis <- factor(UP_DOWN_small_RNA_c_combined_tests$Analysis, levels = unique(UP_DOWN_small_RNA_c_combined_tests$Analysis))



# Supplemental figure A

# Bubble plot to show distribution of RNA based expression changes in chromatin domains according to length & direction


UP_DOWN_small_RNA_lengths_domain <-
  
  ggplot(UP_DOWN_small_RNA_c_combined_tests, aes(y=Description, x=as.numeric(log_odds)))+
  geom_point(aes(color=Analysis, size= find_size))+
  expand_limits(x=c(-4, 4))+
  scale_x_continuous(breaks = seq(-4, 4, by = 1))+
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values = c("dodgerblue", "forestgreen", "magenta4", "darkgoldenrod2", "maroon1", "navy"  ))+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  
  
  scale_alpha(name = paste("Transparency of bubble represents \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value < 0.1", "p value > 0.1"))+
  
  
  ggtitle(paste("Small RNA UP/DOWN/Domain"))



UP_DOWN_small_RNA_lengths_domain  <- 
  
  UP_DOWN_small_RNA_lengths_domain + labs(y="Chromatin Domains\n", x=paste('log2(Odds Ratio) for enrichment'))+
  theme(axis.title=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

UP_DOWN_small_RNA_lengths_domain 


#-------------------------------------------------------------------------------------------
# FIGURE 5


# Gene ontology analyses to be performed here:


# Figure 5 A: RNA seq


# Figure 5 B: ATAC seq


# Figure 5 C: 22G-RNA


# Supplementary Figure 7: RNA + simultaneous ATAC


#-------------------------
# Load the required packages

# here list packages that will be needed
packages = c("stringr",
             "biomaRt",
             "enrichR")

# Now load all (and install if necessary)
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# tell the enrichR package to use worms
setEnrichrSite("WormEnrichr")

# find gene set databases available
dbs <- listEnrichrDbs()

# We need to tell enrichR which databases (from the selection in dbs) we would like to query.
# We can start with KEGG and the three 2018 GO databases

chosendbs <- c("KEGG_2019",
               "GO_Cellular_Component_2018",
               "GO_Molecular_Function_2018",
               "GO_Biological_Process_2018", 
               "InterPro_Domains_2019")


#--------------------------

# Figure 5 A. 

# RNA-seq Gene Set Enrichment Analysis

 RNA_epimutated_genes_list <- list(RNA_longest_genes, RNA_Muts_Inherited_not_long_genes, RNA_Muts_Not_Inherited_genes)



RNA_epimutated_enriched_list <- lapply(RNA_epimutated_genes_list, function(x){
  
  enriched <- enrichr(unlist(x), chosendbs)
  
})



# Sample the Never epimutated list 50 times at a sample size of 8000 genes on each occasion 


RNA_never_epimutated_enriched_list <- list()

for(i in 1:50){
  
  never_mutated_sample <- sample(RNA_No_Muts_genes, 8000, replace = F)
  
  
  never_mutated_sample <- as.data.frame(never_mutated_sample)
  
  
  never_epimutated_enriched  <- enrichr(never_mutated_sample$never_mutated_sample, chosendbs)
  
  
  RNA_never_epimutated_enriched_list[[i]] <- never_epimutated_enriched
  
  
}




combined_sample_result <- lapply(RNA_never_epimutated_enriched_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(RNA_never_epimutated_enriched_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})




# then we combine them and put them in order of significance
combined_NEVER_epimutated_enriched <- lapply(combined_sample_result, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})




# Now we have 2 lists, epimutated_enriched_list has 3 results, and never_epimutated_enriched_list has 50 results 
# We want to compare the three epimutated results with each of the 100 never epimutated sample results 


# and let's combine the results from the 5 libraries 

# we annotate them according to the library


combined_result <- lapply(RNA_epimutated_enriched_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(RNA_epimutated_enriched_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})




# then we combine them and put them in order of significance
combined_epimutated_enriched <- lapply(combined_result, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})




sample_list <- combined_NEVER_epimutated_enriched


split_value_list <- list()

for(i in 1:length(sample_list)){
  
  split <- strsplit(sample_list[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  sample_list[[i]]$Overlap <- split_value
  
}




# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# mean(sample_term) is mean overlap
# sample_not is 8000 minus mean overlap

# then add 1 to all if either test_term OR sample_term = 0

# then only pursue if test_term &/OR sample_term > 4


# let the terms of interest be the terms for the enrichr result for the long_lasting genes

terms_of_interest <- combined_epimutated_enriched[[1]]$Term

# and if the term is not present in the other gene lists then assume there is no enrichment for that term in those lists

long_results_terms_of_interest <- combined_epimutated_enriched[[1]]
short_results_terms_of_interest <- combined_epimutated_enriched[[2]]
non_inherited_results_terms_of_interest <- combined_epimutated_enriched[[3]]

test_list <- list(long_results_terms_of_interest, short_results_terms_of_interest, non_inherited_results_terms_of_interest)


split_value_list <- list()

for(i in 1:length(test_list)){
  
  split <- strsplit(test_list[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
    
    split_value_list[[i]] <- split_value
    
    
  }}


test_list[[1]]$Overlap <- split_value_list[[1]]

test_list[[2]]$Overlap <- split_value_list[[2]]

test_list[[3]]$Overlap <- split_value_list[[3]]




names(test_list) <- c("RNA_longest_results", "RNA_short_results", "RNA_not_inherited_results")


totals_list <- list(length(RNA_longest_genes), length(RNA_Muts_Inherited_not_long_genes), length(RNA_Muts_Not_Inherited_genes))



table_comp <- list(
  LL_table <- c(),
  SL_table <- c(),
  non_inherit_table <- c()
)



for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  
  
  mean_sample_term <- c()
  
  
  
  for(f in 1:length(sample_list)){
    
    sample_term <- 0
    
    if(Term_select %in% sample_list[[f]]$Term == T){
      
      sample_term <- sample_list[[f]][which(sample_list[[f]]$Term %in% Term_select), 2]}
    
    mean_sample_term <- c(mean_sample_term, sample_term)
    
  }
  
  mean_Sample_term <- sum(mean_sample_term)/length(sample_list)
  
  # don't round the mean sample term
  
  
  
  for(p in 1:length(test_list)){
    
    test_term <- 0
    
    if(Term_select %in% test_list[[p]]$Term == T){
      
      test_term <- test_list[[p]][which(test_list[[p]]$Term %in% Term_select), 2]}
    
    test_not <- totals_list[[p]] - test_term
    
    sample_not <- 8000 - mean_Sample_term
    
    
    # save the components in a row in a new table
    
    component_row <- data.frame(test_term, test_not, mean_Sample_term, sample_not)
    
    rownames(component_row) <- Term_select
    
    table_comp[[p]] <- rbind(table_comp[[p]], component_row)
    colnames(table_comp[[p]]) <- c("test_term", "test_not", "mean_Sample_term", "sample_not")
    
  }
  
  
}


for(y in 1:length(table_comp)){
  
  
  # Now modify the table to remove any rows where test term or mean sample term are < 5
  # Now apply the function only if at least one category has 5 or more genes
  # It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category
  
  
  table_comp[[y]] <- table_comp[[y]][(table_comp[[y]]$test_term <5 & table_comp[[y]]$mean_Sample_term < 5)==F, ]
  
  # Now we will adjust  all values in the whole table by adding 1 if there are any zeros
  
  
  zeros <- colSums(table_comp[[y]] == 0)
  
  if(sum(zeros) > 0){
    
    table_comp[[y]]  <-  table_comp[[y]] + 1
  }}


# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list


terms_in_2 <- intersect(rownames(table_comp[[1]]), rownames(table_comp[[2]]))
terms_in_3 <- intersect(terms_in_2, rownames(table_comp[[3]]))

Pval_Results_frame <- c()
X_Results_frame <- c()

for(q in 1:length(table_comp)){
  
  Pval_col <- c()
  X_col <- c()
  
  test_table <- table_comp[[q]][terms_in_3, ]
  
  for(w in 1:nrow(test_table)){
    
    test_term <- test_table[w, ][, 1] 
    test_not <- test_table[w, ][, 2]
    mean_Sample_term <- test_table[w, ][, 3]
    sample_not <- test_table[w, ][, 4]
    
    contingency_table <-
      rbind(c(test_term, test_not),
            c(mean_Sample_term, sample_not))
    X_out <- chisq.test(contingency_table, correct = F)
    pval <- X_out$p.value
    X_squared <- X_out$statistic[[1]]
    
    Pval_col <- c(Pval_col, pval)
    X_col<- c(X_col, X_squared)
    
  }    
  
  Pval_Results_frame <- cbind(Pval_Results_frame, Pval_col)
  X_Results_frame  <- cbind(X_Results_frame, X_col)
  
}


colnames(Pval_Results_frame) <- names(test_list)
rownames(Pval_Results_frame) <- terms_in_3

colnames(X_Results_frame) <- names(test_list)
rownames(X_Results_frame) <- terms_in_3



RNA_LL_enriched_results_pval <- as.data.frame(Pval_Results_frame)
RNA_LL_enriched_results_X <- as.data.frame(X_Results_frame)



# bonferroni correction

for(y in 1:ncol(RNA_LL_enriched_results_pval)){
  RNA_LL_enriched_results_pval[, y] <- p.adjust(RNA_LL_enriched_results_pval[, y], method="bonferroni", n = length(terms_of_interest)) 
  
}



# keep only data for which at least 1 result is significant

test <- RNA_LL_enriched_results_pval < 0.1

keep_all <- c()

for(i in 1:nrow(test)){
  if(TRUE %in% test[i, ]){
    keep_all <- c(keep_all, rownames(RNA_LL_enriched_results_pval[i, ]))
  }
}



RNA_LL_enriched_results_pval <- RNA_LL_enriched_results_pval[keep_all, ] 

RNA_LL_enriched_results_X <- RNA_LL_enriched_results_X[keep_all, ]




# Make plots for this


# i) In each plot order them by the largest significant enrichment out of all of the comparisons otherwise it looks very jumbled and disorganised and is difficult to interpret.


RNA_LL_Sig_X_ordered <- RNA_LL_enriched_results_X[order(RNA_LL_enriched_results_X[, 1], decreasing = T), ]

RNA_LL_Sig_Pval_ordered <- RNA_LL_enriched_results_pval[rownames(RNA_LL_Sig_X_ordered), ]



# Make a bubble plot 

P_1 <- as.numeric(RNA_LL_Sig_Pval_ordered[, 3])

X_1 <- log10(as.numeric(RNA_LL_Sig_X_ordered[, 3]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- cbind(rownames(RNA_LL_Sig_X_ordered), X_1,  P_1, neg_log_p_1)

bubble_table_1 <- data.frame(rep("NON INHERIT"), bubble_table_1)

colnames(bubble_table_1) <- c("Category", "Term", "X_stat", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))



P_2 <- as.numeric(RNA_LL_Sig_Pval_ordered[, 2])

X_2 <- log10(as.numeric(RNA_LL_Sig_X_ordered[, 2]))

neg_log_p_2 <- -log(P_2)

bubble_table_2 <- cbind(rownames(RNA_LL_Sig_X_ordered), X_2,  P_2, neg_log_p_2)

bubble_table_2 <- data.frame(rep("SHORT"), bubble_table_2)

colnames(bubble_table_2) <- c("Category", "Term", "X_stat", "p_val", "neg_log_P")

bubble_table_2$Term <- factor(bubble_table_2$Term, levels = rev(bubble_table_2$Term))



P_3 <- as.numeric(RNA_LL_Sig_Pval_ordered[, 1])

X_3 <- log10(as.numeric(RNA_LL_Sig_X_ordered[, 1]))

neg_log_p_3 <- -log(P_3)

bubble_table_3 <- cbind(rownames(RNA_LL_Sig_X_ordered), X_3,  P_3, neg_log_p_3)

bubble_table_3 <- data.frame(rep("LONG"), bubble_table_3)

colnames(bubble_table_3) <- c("Category", "Term", "X_stat", "p_val", "neg_log_P")

bubble_table_3$Term <- factor(bubble_table_3$Term, levels = rev(bubble_table_3$Term))




c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3)


# Identify which are significant 

find_size <- c()
find_alpha <- c()

for(i in 1:nrow(c_tests)){
  
  alpha <- 0.3
  size_point <- 14
  
  if(as.numeric(c_tests[i, 4]) < 0.1){
    alpha <- 1
    size_point <- 50
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}


RNA_c_combined_tests <- cbind(c_tests, find_alpha, find_size)



# We want to reduce the plot size so only retain the rows for the top 10 enrichment results 


Top_ten <- RNA_c_combined_tests[1:10, 2]

RNA_c_combined_tests <- RNA_c_combined_tests[RNA_c_combined_tests$Term %in% Top_ten, ]



# Trunc terms manually entered as Y axis labels in Adobe Illustrator


 trunc_terms <- c(
   "response to ionizing radiation",                                                                                                                                                             
   "nucleus",
   "defense response to Gram-positive bacterium",
   "detection of chemical stimulus involved in sensory perception",  
   "cytosol",
   "sensory perception of chemical stimulus", 
   "nematode larval development",
   "mitochondrion", 
   "regulation of gene expression", 
   "cellular protein modification process"
 )
 



guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")



RNA_c_combined_tests$Category <- factor(RNA_c_combined_tests$Category, levels = c( "NON INHERIT","LONG", "SHORT"))




# Figure 5 A


bubble_enrichr_non_nest_RNA <-
  
  ggplot(RNA_c_combined_tests, aes(y=Term, x=as.numeric(X_stat)))+
  geom_point(aes(color=Category, size=as.numeric(find_size), alpha=find_alpha))+
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4)) +
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values=c("dodgerblue4", "turquoise4", "mediumaquamarine"))+
  

  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+
  
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.1", "p value not significant > 0.1"))+
  
  
  # ggtitle(paste("Enrichment of ontology terms in genes with altered expression levels \ncompared to genes with no alteration in expression level"))
  
  ggtitle(paste("Figure 5 A"))

bubble_enrichr_non_nest_RNA <- bubble_enrichr_non_nest_RNA + labs(y="Gene Ontology Terms\n", x = "log10(Chi-square statistic for enrichment)")+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(trunc_terms))+
  # scale_y_discrete(labels = " ")+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))


bubble_enrichr_non_nest_RNA

#------------------------------------
# Figure 5 B 

# ATAC Gene Set Enrichment Analysis

# Run the ATAC gene lists through enrichr


ATAC_ep_epimutated_genes_list <- list(ATAC_ep_ABC_Top_genes, ATAC_ep_Muts_Inherited_not_long_genes, ATAC_ep_Muts_Not_Inherited_genes)


ATAC_ep_epimutated_enriched_list <- lapply(ATAC_ep_epimutated_genes_list, function(x){
  
  enriched <- enrichr(unlist(x), chosendbs)
  
})


# and let's combine the results from the 5 libraries 

# we annotate them according to the library


combined_result <- lapply(ATAC_ep_epimutated_enriched_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(ATAC_ep_epimutated_enriched_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})




# then we combine them and put them in order of significance
combined_epimutated_enriched <- lapply(combined_result, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})



# Sample the Never epimutated list 50 times at a sample size of 8000 genes on each occasion 
# and run through enrichr


ATAC_ep_never_epimutated_enriched_list <- list()

for(i in 1:50){
  
  never_mutated_sample <- sample(ATAC_ep_No_Muts_genes, 8000, replace = F)
  
  
  never_mutated_sample <- as.data.frame(never_mutated_sample)
  
  
  never_epimutated_enriched  <- enrichr(never_mutated_sample$never_mutated_sample, chosendbs)
  
  
  ATAC_ep_never_epimutated_enriched_list[[i]] <- never_epimutated_enriched
  
  
}


# and let's combine the results from the 5 libraries 

# we annotate them according to the library


combined_NEVER_result <- lapply(ATAC_ep_never_epimutated_enriched_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(ATAC_ep_never_epimutated_enriched_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})



# then we combine them and put them in order of significance
combined_NEVER_epimutated_enriched <- lapply(combined_NEVER_result, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})




# Now we have 2 lists, combined_epimutated_enriched has 3 results, and combined_never_epimutated_enriched has 50 results 
# We want to compare the three epimutated results with each of the 50 never epimutated sample results 

# Find the overlaps for each category

sample_list <- combined_NEVER_epimutated_enriched


split_value_list <- list()

for(i in 1:length(sample_list)){
  
  split <- strsplit(sample_list[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  sample_list[[i]]$Overlap <- split_value
  
}




# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is the number of genes in the test set which overlap with the term of interest
# test_not is length of gene set minus overlap

# mean(sample_term) is mean number of genes which overlap with term of interest from sample set
# sample_not is 8000 minus mean overlap

# then add 1 to all if either test_term OR sample_term = 0

# then only pursue if test_term &/OR sample_term > 4 (there has to be a minimum of 5 genes in either category)


# let the terms of interest be the terms for the enrichr result for the long_lasting genes

terms_of_interest <- combined_epimutated_enriched[[1]]$Term

# and if the term is not present in the other gene lists then assume there is no enrichment for that term in those lists

long_results_terms_of_interest <- combined_epimutated_enriched[[1]]
short_results_terms_of_interest <- combined_epimutated_enriched[[2]]
non_inherited_results_terms_of_interest <- combined_epimutated_enriched[[3]]

test_list <- list(long_results_terms_of_interest, short_results_terms_of_interest, non_inherited_results_terms_of_interest)


split_value_list <- list()

for(i in 1:length(test_list)){
  
  split <- strsplit(test_list[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
    
    split_value_list[[i]] <- split_value
    
    
  }}


test_list[[1]]$Overlap <- split_value_list[[1]]

test_list[[2]]$Overlap <- split_value_list[[2]]

test_list[[3]]$Overlap <- split_value_list[[3]]




names(test_list) <- c("ATAC_ep_longest_results", "ATAC_ep_short_results", "ATAC_ep_not_inherited_results")


totals_list <- list(length(ATAC_ep_ABC_Top_genes), 
                    length(ATAC_ep_Muts_Inherited_not_long_genes), 
                    length(ATAC_ep_Muts_Not_Inherited_genes))


# We are using the chisq.test because we do not know the exact frequencies in the background population as we are taking the mean from samples
# Do  not round the mean, use the exact mean value 


table_comp <- list(
  LL_table <- c(),
  SL_table <- c(),
  non_inherit_table <- c()
)



for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  
  
  mean_sample_term <- c()
  
  
  
  for(f in 1:length(sample_list)){
    
    sample_term <- 0
    
    if(Term_select %in% sample_list[[f]]$Term == T){
      
      sample_term <- sample_list[[f]][which(sample_list[[f]]$Term %in% Term_select), 2]}
    
    mean_sample_term <- c(mean_sample_term, sample_term)
    
  }
  
  mean_Sample_term <- sum(mean_sample_term)/length(sample_list)
  
  # don't round the mean sample term
  
  
  for(p in 1:length(test_list)){
    
    test_term <- 0
    
    if(Term_select %in% test_list[[p]]$Term == T){
      
      test_term <- test_list[[p]][which(test_list[[p]]$Term %in% Term_select), 2]}
    
    test_not <- totals_list[[p]] - test_term
    
    sample_not <- 8000 - mean_Sample_term
    
    
    # save the components in a row in a new table
    
    component_row <- data.frame(test_term, test_not, mean_Sample_term, sample_not)
    
    rownames(component_row) <- Term_select
    
    table_comp[[p]] <- rbind(table_comp[[p]], component_row)
    colnames(table_comp[[p]]) <- c("test_term", "test_not", "mean_Sample_term", "sample_not")
    
  }
  
  
}


for(y in 1:length(table_comp)){
  
  
  # Now modify the table to remove any rows where test term or mean sample term are < 5
  # Now apply the function only if at least one category has 5 or more genes
  # It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category
  
  
  table_comp[[y]] <- table_comp[[y]][(table_comp[[y]]$test_term <5 & table_comp[[y]]$mean_Sample_term < 5)==F, ]
  
  # Now we will adjust  all values by adding 1 if there are any zeros
  
  
  zeros <- colSums(table_comp[[y]] == 0)
  
  if(sum(zeros) > 0){
    
    table_comp[[y]]  <-  table_comp[[y]] + 1
  }}


# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list


terms_in_2 <- intersect(rownames(table_comp[[1]]), rownames(table_comp[[2]]))
terms_in_3 <- intersect(terms_in_2, rownames(table_comp[[3]]))

Pval_Results_frame <- c()
X_Results_frame <- c()

for(q in 1:length(table_comp)){
  
  Pval_col <- c()
  X_col <- c()
  
  test_table <- table_comp[[q]][terms_in_3, ]
  
  for(w in 1:nrow(test_table)){
    
    test_term <- test_table[w, ][, 1] 
    test_not <- test_table[w, ][, 2]
    mean_Sample_term <- test_table[w, ][, 3]
    sample_not <- test_table[w, ][, 4]
    
    contingency_table <-
      rbind(c(test_term, test_not),
            c(mean_Sample_term, sample_not))
    X_out <- chisq.test(contingency_table, correct = F)
    pval <- X_out$p.value
    X_squared <- X_out$statistic[[1]]
    
    Pval_col <- c(Pval_col, pval)
    X_col<- c(X_col, X_squared)
    
  }    
  
  Pval_Results_frame <- cbind(Pval_Results_frame, Pval_col)
  X_Results_frame  <- cbind(X_Results_frame, X_col)
  
}


colnames(Pval_Results_frame) <- names(test_list)
rownames(Pval_Results_frame) <- terms_in_3

colnames(X_Results_frame) <- names(test_list)
rownames(X_Results_frame) <- terms_in_3


ATAC_ep_LL_enriched_results_pval <- as.data.frame(Pval_Results_frame)
ATAC_ep_LL_enriched_results_X <- as.data.frame(X_Results_frame)



# Bonferroni correction

for(y in 1:ncol(ATAC_ep_LL_enriched_results_pval)){
  ATAC_ep_LL_enriched_results_pval[, y] <- p.adjust(ATAC_ep_LL_enriched_results_pval[, y], method="bonferroni", n = length(terms_in_3)) 
  
}



# keep only data for which at least 1 result is significant

test <- ATAC_ep_LL_enriched_results_pval < 0.1

keep_all <- c()

for(i in 1:nrow(test)){
  if(TRUE %in% test[i, ]){
    keep_all <- c(keep_all, rownames(ATAC_ep_LL_enriched_results_pval[i, ]))
  }
}



ATAC_ep_LL_enriched_results_pval <- ATAC_ep_LL_enriched_results_pval[keep_all, ] 

ATAC_ep_LL_enriched_results_X <- ATAC_ep_LL_enriched_results_X[keep_all, ]



ATAC_LL_ep_Sig_X_ordered <- ATAC_ep_LL_enriched_results_X[order(ATAC_ep_LL_enriched_results_X[, 1], decreasing = T), ]

ATAC_LL_ep_Sig_Pval_ordered <- ATAC_ep_LL_enriched_results_pval[rownames(ATAC_LL_ep_Sig_X_ordered), ]



# Make a bubble plot 

P_1 <- as.numeric(ATAC_LL_ep_Sig_Pval_ordered[, 3])

X_1 <- log10(as.numeric(ATAC_LL_ep_Sig_X_ordered[, 3]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- cbind(rownames(ATAC_LL_ep_Sig_X_ordered), X_1,  P_1, neg_log_p_1)

bubble_table_1 <- data.frame(rep("NON INHERIT"), bubble_table_1)

colnames(bubble_table_1) <- c("Category", "Term", "X_stat", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))



P_2 <- as.numeric(ATAC_LL_ep_Sig_Pval_ordered[, 2])

X_2 <- log10(as.numeric(ATAC_LL_ep_Sig_X_ordered[, 2]))

neg_log_p_2 <- -log(P_2)

bubble_table_2 <- cbind(rownames(ATAC_LL_ep_Sig_X_ordered), X_2,  P_2, neg_log_p_2)

bubble_table_2 <- data.frame(rep("SHORT"), bubble_table_2)

colnames(bubble_table_2) <- c("Category", "Term", "X_stat", "p_val", "neg_log_P")

bubble_table_2$Term <- factor(bubble_table_2$Term, levels = rev(bubble_table_2$Term))



P_3 <- as.numeric(ATAC_LL_ep_Sig_Pval_ordered[, 1])

X_3 <- log10(as.numeric(ATAC_LL_ep_Sig_X_ordered[, 1]))

neg_log_p_3 <- -log(P_3)

bubble_table_3 <- cbind(rownames(ATAC_LL_ep_Sig_X_ordered), X_3,  P_3, neg_log_p_3)

bubble_table_3 <- data.frame(rep("LONG"), bubble_table_3)

colnames(bubble_table_3) <- c("Category", "Term", "X_stat", "p_val", "neg_log_P")

bubble_table_3$Term <- factor(bubble_table_3$Term, levels = rev(bubble_table_3$Term))



c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3)



find_size <- c()
find_alpha <- c()

for(i in 1:nrow(c_tests)){
  
  alpha <- 0.3
  size_point <- 14
  
  if(as.numeric(c_tests[i, 4]) < 0.1){
    alpha <- 1
    size_point <- 50
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}


ATAC_ep_c_combined_tests <- cbind(c_tests, find_alpha, find_size)



Top_ten <- ATAC_ep_c_combined_tests[1:10, 2]

ATAC_ep_c_combined_tests <- ATAC_ep_c_combined_tests[ATAC_ep_c_combined_tests$Term %in% Top_ten, ]


# Trunc terms manually entered as Y axis labels in Adobe Illustrator

 
 trunc_terms <- c(
   "regulation of defense response to bacterium" ,            
   "SET domain",                                                            
   "translation initiation factor activity",                
   "translation factor activity, RNA binding",                
   "establishment of protein localization to membrane",       
   "cyclin-dependent protein kinase activity",               
   "cyclin-dependent protein serine/threonine kinase activity",
   "regulation of dendrite development",                      
   "RNA binding",                                            
   "nuclear lumen")



guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")




# FIGURE 5 B

ATAC_ep_c_combined_tests$Category <- factor(ATAC_ep_c_combined_tests$Category, levels = c( "NON INHERIT","LONG", "SHORT"))



bubble_enrichr_non_nest_ATAC_ep <-
  
  ggplot(ATAC_ep_c_combined_tests, aes(y=Term, x=as.numeric(X_stat)))+
  geom_point(aes(color=Category, size=as.numeric(find_size), alpha=find_alpha))+
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4)) +
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values=c("dodgerblue4", "turquoise4", "mediumaquamarine"))+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+  
  
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.1", "p value not significant > 0.1"))+
  
  
  # ggtitle(paste("Enrichment of ontology terms in genes with altered chromatin states \ncompared to genes with no alteration in chromatin state"))
  ggtitle(paste("Figure 5 B"))

bubble_enrichr_non_nest_ATAC_ep <- bubble_enrichr_non_nest_ATAC_ep + labs(y="Gene Ontology Terms\n", x="log10(Chi-square statistic for enrichment)")+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(trunc_terms))+
  # scale_y_discrete(labels = " ")+
  theme(axis.text.x = element_text(color="#000000", size = 10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))


bubble_enrichr_non_nest_ATAC_ep 


#----------------------

# Figure 5 C 
# Small RNA Gene Set Enrichment Analysis 

small_RNA_epimutated_bg_genes_list <- list(small_RNA_Top_genes, small_RNA_Muts_Inherited_not_long, small_RNA_Muts_Not_Inherited, small_RNA_No_Muts)




# Send the test lists to enrichr

small_RNA_epimutated_bg_enriched_list <- lapply(small_RNA_epimutated_bg_genes_list, function(x){
  
  enriched <- enrichr(unlist(x), chosendbs)
  
})


combined_sample_bg_result <- lapply(small_RNA_epimutated_bg_enriched_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(small_RNA_epimutated_bg_enriched_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})




# then we combine them and put them in order of significance
combined_sample_bg_RESULT <- lapply(combined_sample_bg_result, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})



# We want to compare the three epimutated results with the never epimutated sample results 


test_list <- list(combined_sample_bg_RESULT[[1]], 
                  combined_sample_bg_RESULT[[2]],
                  combined_sample_bg_RESULT[[3]],
                  combined_sample_bg_RESULT[[4]])



split_value_list <- list()

for(i in 1:length(test_list)){
  
  split <- strsplit(test_list[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
    
    split_value_list[[i]] <- split_value
    
    
  }}


test_list[[1]]$Overlap <- split_value_list[[1]]

test_list[[2]]$Overlap <- split_value_list[[2]]

test_list[[3]]$Overlap <- split_value_list[[3]]

test_list[[4]]$Overlap <- split_value_list[[4]]


names(test_list) <- c("smallRNA_mut_longest_results", "smallRNA_mut_short_results", "smallRNA_mut_not_inherited_results", "smallRNA_mut_not_mutated_results")


totals_list <- list(length(small_RNA_Top_genes), 
                    length(small_RNA_Muts_Inherited_not_long), 
                    length(small_RNA_Muts_Not_Inherited))


BG_list <- test_list[[4]]


# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap

# then add 1 to all if either test_term OR sample_term = 0

# then only pursue if test_term &/OR sample_term > 4


# let the terms of interest be the terms for the enrichr result for the list of genes with long lasting small RNA changes

terms_of_interest <- combined_sample_bg_RESULT[[1]]$Term

# and if the term is not present in the other gene lists then assume there is no enrichment for that term in those lists

combined_test_RESULT   <- test_list[1:3]

BG_List <- test_list[[4]]



table_comp <- list(
  LL_table <- c(),
  SL_table <- c(),
  non_inherit_table <- c()
)




for(j in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[j]]
  
  
  
  for(i in 1:length(combined_test_RESULT)){
    
    Test_List <- combined_test_RESULT[[i]]
    
    
    
    
    BG_term <- 0
    
    if(Term_select %in% BG_List$Term == T){
      
      BG_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
    
    BG_not <- length(RNA_No_Muts) - BG_term
    
    
    test_term <- 0
    
    if(Term_select %in% Test_List$Term == T){
      
      test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]}
    test_not <- totals_list[[i]] - test_term
    
    
    # save the components in a row in a new table
    
    component_row <- data.frame(test_term, test_not, BG_term, BG_not)
    
    rownames(component_row) <- Term_select
    
    table_comp[[i]] <- rbind(table_comp[[i]], component_row)
    colnames(table_comp[[i]]) <- c("test_term", "test_not", "BG_term", "BG_not")
  }
  
}



for(y in 1:length(table_comp)){
  
  # Now modify the table to remove any rows where test term or mean sample term are < 5
  # Now apply the function only if at least one category has 5 or more genes
  # It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category
  
  
  table_comp[[y]] <- table_comp[[y]][(table_comp[[y]]$test_term <5 & table_comp[[y]]$BG_term < 5)==F, ]
  
  # Now we will adjust  all values by adding 1 if there are any zeros
  
  
  zeros <- colSums(table_comp[[y]] == 0)
  
  if(sum(zeros) > 0){
    
    table_comp[[y]]  <-  table_comp[[y]] + 1
  }}


# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list


terms_in_2 <- intersect(rownames(table_comp[[1]]), rownames(table_comp[[2]]))
terms_in_3 <- intersect(terms_in_2, rownames(table_comp[[3]]))

Pval_Results_frame <- c()
OR_Results_frame <- c()

for(q in 1:length(table_comp)){
  
  Pval_col <- c()
  OR_col <- c()
  
  test_table <- table_comp[[q]][terms_in_3, ]
  
  for(w in 1:nrow(test_table)){
    
    test_term <- test_table[w, ][, 1] 
    test_not <- test_table[w, ][, 2]
    BG_term <- test_table[w, ][, 3]
    BG_not <- test_table[w, ][, 4]
    
    contingency_table <-
      rbind(c(test_term, test_not),
            c(BG_term, BG_not))
    FT_out <- fisher.test(contingency_table)
    pval <- FT_out$p.value
    OddsRatio <- FT_out$estimate
    
    Pval_col <- c(Pval_col, pval)
    OR_col<- c(OR_col, OddsRatio)
    
  }    
  
  Pval_Results_frame <- cbind(Pval_Results_frame, Pval_col)
  OR_Results_frame  <- cbind(OR_Results_frame, OR_col)
  
}


colnames(Pval_Results_frame) <- c("Long_lasting", "Short_lasting", "Not_inherited")
rownames(Pval_Results_frame) <- terms_in_3

colnames(OR_Results_frame) <- c("Long_lasting", "Short_lasting", "Not_inherited")
rownames(OR_Results_frame) <- terms_in_3



smallRNA_LL_enriched_results_pval <- as.data.frame(Pval_Results_frame)
smallRNA_LL_enriched_results_OR <- as.data.frame(OR_Results_frame)



# bonferroni correction

for(y in 1:ncol(smallRNA_LL_enriched_results_pval)){
  smallRNA_LL_enriched_results_pval[, y] <- p.adjust(smallRNA_LL_enriched_results_pval[, y], method="bonferroni", n = length(terms_in_3)) 
  
}



# keep only data for which at least 1 result is significant

test <- smallRNA_LL_enriched_results_pval < 0.1

keep_all <- c()

for(i in 1:nrow(test)){
  if(TRUE %in% test[i, ]){
    keep_all <- c(keep_all, rownames(smallRNA_LL_enriched_results_pval[i, ]))
  }
}



smallRNA_LL_enriched_results_pval <- smallRNA_LL_enriched_results_pval[keep_all, ] 

smallRNA_LL_enriched_results_OR <- smallRNA_LL_enriched_results_OR[keep_all, ]



# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 

smallRNA_LL_ep_OR_ordered <- smallRNA_LL_enriched_results_OR[order(smallRNA_LL_enriched_results_OR[, 1], decreasing = T), ]

smallRNA_LL_ep_Sig_Pval_ordered <- smallRNA_LL_enriched_results_pval[rownames(smallRNA_LL_ep_OR_ordered), ]



# Make a bubble plot 

P_1 <- as.numeric(smallRNA_LL_ep_Sig_Pval_ordered[, 3])

X_1 <- log10(as.numeric(smallRNA_LL_ep_OR_ordered[, 3]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- cbind(rownames(smallRNA_LL_ep_OR_ordered), X_1,  P_1, neg_log_p_1)

bubble_table_1 <- data.frame(rep("NON INHERIT"), bubble_table_1)

colnames(bubble_table_1) <- c("Category", "Term", "OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))



P_2 <- as.numeric(smallRNA_LL_ep_Sig_Pval_ordered[, 2])

X_2 <- log10(as.numeric(smallRNA_LL_ep_OR_ordered[, 2]))

neg_log_p_2 <- -log(P_2)

bubble_table_2 <- cbind(rownames(smallRNA_LL_ep_OR_ordered), X_2,  P_2, neg_log_p_2)

bubble_table_2 <- data.frame(rep("SHORT"), bubble_table_2)

colnames(bubble_table_2) <- c("Category", "Term", "OR", "p_val", "neg_log_P")

bubble_table_2$Term <- factor(bubble_table_2$Term, levels = rev(bubble_table_2$Term))



P_3 <- as.numeric(smallRNA_LL_ep_Sig_Pval_ordered[, 1])

X_3 <- log10(as.numeric(smallRNA_LL_ep_OR_ordered[, 1]))

neg_log_p_3 <- -log(P_3)

bubble_table_3 <- cbind(rownames(smallRNA_LL_ep_OR_ordered), X_3,  P_3, neg_log_p_3)

bubble_table_3 <- data.frame(rep("LONG"), bubble_table_3)

colnames(bubble_table_3) <- c("Category", "Term", "OR", "p_val", "neg_log_P")

bubble_table_3$Term <- factor(bubble_table_3$Term, levels = rev(bubble_table_3$Term))


c_tests <- rbind(bubble_table_1, bubble_table_2, bubble_table_3)




find_size <- c()
find_alpha <- c()

for(i in 1:nrow(c_tests)){
  
  alpha <- 0.3
  size_point <- 14
  
  if(as.numeric(c_tests[i, 4]) < 0.1){
    alpha <- 1
    size_point <- 50
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}


smallRNA_ep_c_combined_tests <- cbind(c_tests, find_alpha, find_size)




Top_ten <- smallRNA_ep_c_combined_tests[1:10, 2]



smallRNA_ep_c_combined_tests <- smallRNA_ep_c_combined_tests[smallRNA_ep_c_combined_tests$Term %in% Top_ten, ]


# Trunc terms manually entered as Y axis labels in Adobe Illustrator

 trunc_terms <- c(
   "modulation of chemical synaptic transmission",
   "ion transport",                         
   "regulation of eating behavior",              
   "protein targeting",                       
   "regulation of pharyngeal pumping",         
   "inductive cell migration",                    
   "microbody",                                 
   "peroxisome",                                  
   "chloride transmembrane transport",           
   "chloride transport" 
 )





guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")



# Figure 5 C

smallRNA_ep_c_combined_tests$Category <- factor(smallRNA_ep_c_combined_tests$Category, levels = c( "NON INHERIT","LONG", "SHORT"))


bubble_enrichr_non_nest_smallRNA_ep <-
  
  ggplot(smallRNA_ep_c_combined_tests, aes(y=Term, x=as.numeric(OR)))+
  geom_point(aes(color=Category, size=as.numeric(find_size), alpha=find_alpha))+
  scale_x_continuous(breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), limits = c(-4, 4)) +
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  scale_color_manual(values=c("dodgerblue4", "turquoise4", "mediumaquamarine"))+
  
  
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+  
  
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.1", "p value not significant > 0.1"))+
  
  
  # ggtitle(paste("Enrichment of ontology terms in genes targeted by small RNA epimutations \ncompared to genes without small RNA epimutations"))
  
  
  ggtitle(paste("Figure 5 C"))


bubble_enrichr_non_nest_smallRNA_ep <- bubble_enrichr_non_nest_smallRNA_ep + labs(y="Gene Ontology Terms\n", x="log10(Odds Ratio for enrichment)")+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(trunc_terms))+
  # scale_y_discrete(labels= " ")+
  theme(axis.text.x = element_text(color="#000000", size = 10))+
  theme(axis.text.y = element_text(color="#000000", size = 8))+
  guides(size = guide_legend(order = 3), colour = guide_legend(order = 1), alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

bubble_enrichr_non_nest_smallRNA_ep 



#----------------------------

# Supplementary Figure 7

# GSEA for genes with inherited RNA expression changes AND simultaneous ATAC

# Do a Gene Ontology analysis of genes which have simultaneous RNAseq and Chromatin epimutations 
# for which RNA seq changes are inherited
# For the background use genes with non inherited RNAseq changes (of any kind) and non simultaneous chromatin epimutations 

Sim_test_list <- unique(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1 & All_lins_integrated_table$is_RNA_inherited==1 & All_lins_integrated_table$is_ep_epimut==1 & All_lins_integrated_table$has_time_matched==1,  2])

Sim_bg_list <- unique(All_lins_integrated_table[All_lins_integrated_table$gene_maps_to_ep==1&All_lins_integrated_table$is_RNA_exp_change==1 & All_lins_integrated_table$has_time_matched==0,  2])

# Remove any genes from the bg list that feature in the test list

intersect <- intersect(Sim_bg_list, Sim_test_list)
Sim_bg_list <- Sim_bg_list[Sim_bg_list %in% intersect == F]

Sim_epimutated_enriched_list  <- enrichr(Sim_test_list, chosendbs)
Sim_epimutated_bg_list  <- enrichr(Sim_bg_list, chosendbs)

Sim_test_bg_list <- list(Sim_epimutated_enriched_list, Sim_epimutated_bg_list)

# Annotate with name of libraries

Sim_test_bg_results <- lapply(Sim_test_bg_list, function(x){
  
  lapply(1:5, function(j){
    
    if(nrow(x[[j]]) > 0){
      x[[j]][1:nrow(x[[j]]), "library"] <- names(Sim_test_bg_list[[1]])[j]
      return(x[[j]])
    } 
    
  })
  
})




# then we combine them and put them in order of significance
Sim_test_bg_RESULTS <- lapply(Sim_test_bg_results, function(x){
  
  tempdf <- do.call(rbind, x)
  tempdf[order(tempdf$Old.Adjusted.P.value), ]
  
})





split_value_list <- list()

for(i in 1:length(Sim_test_bg_RESULTS)){
  
  split <- strsplit(Sim_test_bg_RESULTS[[i]][, 2], "/")
  
  split_value <- c()
  
  for(j in 1:length(split)){
    
    save <- split[[j]][1]
    
    split_value <- as.numeric(c(split_value, save))
  }
  
  Sim_test_bg_RESULTS[[i]]$Overlap <- split_value
  
}




# get the test_term, the test_not, the mean(sample_term) and the sample_not

# test_term is overlap
# test_not is length of gene set minus overlap


# then add 1 to all if either test_term OR sample_term = 0
# then only pursue if test_term &/OR sample_term > 4
# let the terms of interest be the terms for the enrichr result for the test list, i.e. with simultaneously epimutated genes

terms_of_interest <- Sim_test_bg_RESULTS[[1]]$Term

# and if the term is not present in the other gene lists then assume there is no enrichment for that term in those lists



Test_List <- Sim_test_bg_RESULTS[[1]]
BG_List <- Sim_test_bg_RESULTS[[2]]




table_comp <- c()



for(i in 1:length(terms_of_interest)){
  
  Term_select <- terms_of_interest[[i]]
  
  
  sample_term <- 0
  
  if(Term_select %in% BG_List$Term == T){
    
    sample_term <- BG_List[which(BG_List$Term %in% Term_select), 2]}
  sample_not <- length(Sim_bg_list) - sample_term
  
  test_term <- Test_List[which(Test_List$Term %in% Term_select), 2]
  test_not <- length(Sim_test_list) - test_term
  
  
  # save the components in a row in a new table
  
  component_row <- data.frame(test_term, test_not, sample_term, sample_not)
  
  rownames(component_row) <- Term_select
  
  table_comp <- rbind(table_comp, component_row)
  colnames(table_comp) <- c("test_term", "test_not", "sample_term", "sample_not")
  
}



# Now we will adjust  all values by adding 1 if there are any zeros


zeros <- colSums(table_comp == 0)

if(sum(zeros) > 0){
  
  table_comp <-  table_comp + 1
}


# Now modify the table to remove any rows where test term or mean sample term are < 5
# Now apply the function only if at least one category has 5 or more genes
# It is only reasonable to assess for relative enrichment/depletion if either the test list or sample list has > 5 genes in that ontology category


table_comp <- table_comp[(table_comp$test_term < 5 & table_comp$sample_term < 5)==F, ]

# Now we have found terms where there are at least 5 genes in 1 test or sample category
# we have adjusted the data for each gene set for any 0 values relative to the specific test list
# Now we can do the analysis but only on terms present in each table_comp list


Pval_col <- c()
OR_col <- c()

for(q in 1:nrow(table_comp)){
  
  test_term <- table_comp[q, 1] 
  test_not <- table_comp[q, 2]
  mean_Sample_term <- table_comp[q, 3]
  sample_not <- table_comp[q, 4]
  
  contingency_table <-
    rbind(c(test_term, test_not),
          c(mean_Sample_term, sample_not))
  FT_out <- fisher.test(contingency_table)
  pval <- FT_out$p.value
  OR <- FT_out$estimate
  
  Pval_col <- c(Pval_col, pval)
  OR_col<- c(OR_col, OR)
  
}    


OR_pval_frame <- data.frame(Pval_col, OR_col)

colnames(OR_pval_frame) <- c("Pval", "OR")
rownames(OR_pval_frame) <- rownames(table_comp)


# Bonferroni correction

OR_pval_frame[, 1] <- p.adjust(OR_pval_frame[, 1], method="bonferroni", n = nrow(OR_pval_frame)) 


# Make plots for this

# i) In each plot order them by the largest significant enrichment out of all of the comparisons 


OR_pval_frame_ordered <- OR_pval_frame[order(OR_pval_frame[, 2], decreasing = T), ]



# Make a bubble plot 

P_1 <- as.numeric(OR_pval_frame_ordered[, 1])

OR_1 <- log10(as.numeric(OR_pval_frame_ordered[, 2]))

neg_log_p_1 <- -log(P_1)

bubble_table_1 <- data.frame(rownames(OR_pval_frame_ordered), OR_1,  P_1, neg_log_p_1)

colnames(bubble_table_1) <- c("Term", "log_OR", "p_val", "neg_log_P")

bubble_table_1$Term <- factor(bubble_table_1$Term, levels = rev(bubble_table_1$Term))





find_size <- c()
find_alpha <- c()

for(i in 1:nrow(bubble_table_1)){
  
  alpha <- 0.3
  size_point <- 14
  
  if(as.numeric(bubble_table_1[i, 3]) < 0.1){
    alpha <- 1
    size_point <- 50
    
  }
  find_alpha <- c(find_alpha, alpha)  
  find_size <- c(find_size, size_point)
}


bubble_table <- cbind(bubble_table_1, find_alpha, find_size)




# just plot first 10

bubble_table <- bubble_table[1:10, ]
trunc_terms <- as.character(bubble_table$Term)


# Trunc terms manually entered as Y axis labels in Adobe Illustrator

trunc_terms <- c(

 "defense response to Gram-positive bacterium (GO:0050830)",
 "Autophagy",                                               
 "nematode larval development (GO:0002119)",                
 "MAPK signaling pathway",                                  
 "ATPase activity (GO:0016887)",                            
 "defense response to bacterium (GO:0042742)",              
 "EGF-like domain",                                       
 "cellular response to organic substance (GO:0071310)",     
 "defense response to fungus (GO:0050832)",                 
 "aspartic-type peptidase activity (GO:0070001)"
)





guides_merge <- function(gdefs) {
  gdefs <- lapply(gdefs, function(g) { g$hash <- paste(g$order, g$hash, sep = "z"); g})
  tapply(gdefs, sapply(gdefs, function(g)g$hash), function(gs)Reduce(guide_merge, gs))
}
environment(guides_merge) <- environment(ggplot)
assignInNamespace("guides_merge", guides_merge, pos = "package:ggplot2")


# Supplementary Figure 7

bubble_Sim_RNA_chrom <-
  
  ggplot(bubble_table, aes(y=Term, x=as.numeric(log_OR)))+
  geom_point(aes(color=Term, size=as.numeric(find_size), alpha=find_alpha))+
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2), limits = c(-2, 2)) +
  geom_vline(xintercept = 0, colour = "grey")+
  theme_bw()+
  
 
  scale_size_continuous(range = c(2, 15), breaks = c(14, 50), 
                        limits = c(12, 51))+  
  
  
  scale_alpha(name = paste("Transparency indicates \nsignificance of enrichment"), range = c(0.3, 1),
              breaks = seq(1, 0.3, length = 2),
              limits = c(0.29, 1.05), 
              labels = c("p value significant < 0.1", "p value not significant > 0.1"))+
  
  
  ggtitle(paste("Supplementary Figure 7"))

bubble_Sim_RNA_chrom <- bubble_Sim_RNA_chrom + labs(y="Gene Ontology Terms\n", x = "log10(Odds Ratio for enrichment)")+
  theme(axis.title.x=element_text(face = "bold", size=14))+
  theme(axis.title.y=element_text(face = "bold", size=14))+
  theme(plot.title = element_text(face = "bold", size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_discrete(labels= rev(trunc_terms))+
  # scale_y_discrete(labels= "")+
  theme(axis.text.x = element_text(color="#000000", size=10))+
  theme(axis.text.y = element_text(color="#000000", size = 14))+
  guides(size = guide_legend(order = 3), colour = "none", alpha = guide_legend(order = 2))+
  theme(legend.text=element_text(color="#000000", size=12))

bubble_Sim_RNA_chrom



#-------------------------------
# GSEA for genes with RNA expression changes AND simultaneous small RNA epimutations 
# No significant results found
#-------------------------------


# Figure 6 and Supplementary Figure 8

# PGP genes

# RNA pgp

pgp_rows <- rownames(RNA_Z_score_table_A)[grep("pgp-", rownames(RNA_Z_score_table_A))]

RNA_A_pgp <- RNA_Z_score_table_A[pgp_rows,]


RNA_A_pgp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 1]),
    cbind(rep(4, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 2]),
    cbind(rep(6, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 3]),
    cbind(rep(8, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 4]),
    cbind(rep(10, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 5]),
    cbind(rep(12, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 6]),
    cbind(rep(14, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 7]),
    cbind(rep(16, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 8]),
    cbind(rep(18, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 9]),
    cbind(rep(20, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 10]))


RNA_pgp_genes <- strsplit(pgp_rows, split = ":")

pgp_genes <- c()

for(i in 1:length(RNA_pgp_genes)){
  
  pgp_genes <- c(pgp_genes, RNA_pgp_genes[[i]][4])
  
}


RNA_A_pgp_Z_scores <- cbind(rownames(RNA_A_pgp_Z_scores), rep(pgp_genes, 10), data.frame(RNA_A_pgp_Z_scores, row.names=NULL)) 

colnames(RNA_A_pgp_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")

pgp_trunc_names <- c("pgp-5", "pgp-6", "pgp-7", "pgp-8")

RNA_A_pgp_Z_scores_trunc <- RNA_A_pgp_Z_scores[RNA_A_pgp_Z_scores$gene_name %in% pgp_trunc_names, ]


# -------------
# ATAC pgp

ATAC_ep_pgp_rows <- unique(All_lins_ATAC_table[grep("pgp-", All_lins_ATAC_table$Gene), 3])

ATAC_A_pgp <- ATAC_Z_score_table_A_enhancer_promoter[ATAC_ep_pgp_rows,]


Pgp_ep_loci <- rownames(ATAC_A_pgp)

Pgp_ep_genes <- All_lins_ATAC_table[All_lins_ATAC_table$Locus %in% Pgp_ep_loci & All_lins_ATAC_table$Lineage == "A", 2]
Pgp_ep_genes <- Pgp_ep_genes[grep("pgp-", Pgp_ep_genes)]


ATAC_A_pgp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 1]),
    cbind(rep(4, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 2]),
    cbind(rep(6, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 3]),
    cbind(rep(8, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 4]),
    cbind(rep(10, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 5]),
    cbind(rep(12, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 6]),
    cbind(rep(14, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 7]),
    cbind(rep(16, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 8]),
    cbind(rep(18, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 9]),
    cbind(rep(20, length = nrow(ATAC_A_pgp)), ATAC_A_pgp[, 10]))

pgp_ep_loci <- paste(rownames(ATAC_A_pgp_Z_scores), Pgp_ep_genes, sep = ":")

ATAC_A_pgp_Z_scores  <- data.frame(rep(Pgp_ep_genes, 10), ATAC_A_pgp_Z_scores, pgp_ep_loci)

colnames(ATAC_A_pgp_Z_scores) <- c("gene_name", "Generation", "Expression", "ep_locus")


ATAC_A_pgp_Z_scores_trunc <- ATAC_A_pgp_Z_scores[ATAC_A_pgp_Z_scores$gene_name %in% pgp_trunc_names, ]

#------------------------------------------------------

# small RNA pgp

pgp_genes <- rownames(smallRNA_Z_score_table_A)[grep("pgp-", rownames(smallRNA_Z_score_table_A))]

smallRNA_A_pgp <- smallRNA_Z_score_table_A[pgp_genes,]


smallRNA_A_pgp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 1]),
    cbind(rep(4, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 2]),
    cbind(rep(6, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 3]),
    cbind(rep(8, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 4]),
    cbind(rep(10, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 5]),
    cbind(rep(12, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 6]),
    cbind(rep(14, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 7]),
    cbind(rep(16, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 8]),
    cbind(rep(18, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 9]),
    cbind(rep(20, length = nrow(smallRNA_A_pgp)), smallRNA_A_pgp[, 10]))



smallRNA_A_pgp_Z_scores <- cbind(paste("small_RNA_targeted", rownames(smallRNA_A_pgp_Z_scores), sep = ":"), rep(pgp_genes, 10), data.frame(smallRNA_A_pgp_Z_scores, row.names=NULL)) 

colnames(smallRNA_A_pgp_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")

pgp_trunc_names <- c("pgp-5", "pgp-6", "pgp-7", "pgp-8")

smallRNA_A_pgp_Z_scores_trunc <- smallRNA_A_pgp_Z_scores[smallRNA_A_pgp_Z_scores$gene_name %in% pgp_trunc_names, ]




# -------------------------------

# Figure 6 


# Figure 6 A is this plot but with tracks for promoter chromatin state aside from pgp-6 removed


ggplot(RNA_A_pgp_Z_scores_trunc, aes(x=Generation, y=Expression, group=gene_name, color=gene_name)) +
  geom_line()+
  # scale_color_manual(values = mycolors)+
  geom_line(data = ATAC_A_pgp_Z_scores_trunc, aes(x=Generation, y=Expression, group=ep_locus, color=ep_locus))+
  ggtitle("Figure 6 A") +
  labs(x = "Generations)",
       y = "Z score for expression/chromatin",
       color = "pgp genes")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=c(-2.25, 2.25), col = "black")+
  scale_x_continuous(breaks=seq(2, 20, 2), limits = c(2, 20))+
  scale_y_continuous(breaks = c(-5, -4, -3, -2.25, -2, -1, 0, 1, 2, 2.25, 3, 4, 5), limits = c(-5, 5),
                     labels = c("-5", "-4", "-3", "Epimutations below cut off", "-2", "-1", "0", "1", "2", "Epimutations above cut off", "3", "4", "5"))+
  # geom_rect(aes(xmin=13.5, xmax=14.5, ymin=-Inf, ymax=Inf), colour = NA, fill="gray", alpha=0.01)+
  theme(plot.margin = margin(1,1,1,1,"cm"))





# Fig 6 B 
# Lineage A pgp RNA and 22G-RNA 

ggplot(RNA_A_pgp_Z_scores_trunc, aes(x=Generation, y=Expression, group=gene_name, color=gene_name)) +
  geom_line()+
  # scale_color_manual(values = mycolors)+
  geom_line(data = smallRNA_A_pgp_Z_scores_trunc, aes(x=Generation, y=Expression, group=locus_name, color=locus_name))+
  ggtitle("Figure 6 B") +
  labs(x = "Generations)",
       y = "Z score for expression/chromatin",
       color = "pgp genes")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=c(-2.25, 2.25), col = "black")+
  scale_x_continuous(breaks=seq(2, 20, 2), limits = c(2, 20))+
  scale_y_continuous(breaks = c(-5, -4, -3, -2.25, -2, -1, 0, 1, 2, 2.25, 3, 4, 5), limits = c(-5, 5),
                     labels = c("-5", "-4", "-3", "Epimutations below cut off", "-2", "-1", "0", "1", "2", "Epimutations above cut off", "3", "4", "5"))+
  # geom_rect(aes(xmin=13.5, xmax=14.5, ymin=-Inf, ymax=Inf), colour = NA, fill="gray", alpha=0.01)+
  theme(plot.margin = margin(1,1,1,1,"cm"))





# We want to see if apparent minor chromatin changes on the X chromosome are significant
# Generate Z scores just for the subset of genes on the X chromosome and then plot these 

# Start with the ATACseq changes

# Produce the binarised tables through using the linear model to identify epimutations accorrding to z score cut off of 2.25
# But restrict the distribution to X chromosome genes

# Lineage A

# Find Z scores for ATAC in the context of the X chromosome only

ChrX_rows <- grep("chrX", rownames(normalized_table_Lin_A))

ChrX_normalized_table_Lin_A <- normalized_table_Lin_A[ChrX_rows, ]

ChrX_normalized_table_A2_to_A20 <- ChrX_normalized_table_Lin_A[, 2:11]

ChrX_ATAC_Z_score_table_A_enhancer_promoter <- matrix(0, ncol = 10, nrow = nrow(ChrX_normalized_table_Lin_A))

for (i in 1:ncol(ChrX_normalized_table_A2_to_A20)) {
  LM_temp <-
    lm(log2(ChrX_normalized_table_A2_to_A20[, i]) ~ log2(ChrX_normalized_table_Lin_A[, 1]))
  Residuals <- residuals(LM_temp)
  ChrX_ATAC_Z_score_table_A_enhancer_promoter[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}


rownames(ChrX_ATAC_Z_score_table_A_enhancer_promoter) <- rownames(ChrX_normalized_table_A2_to_A20)
colnames(ChrX_ATAC_Z_score_table_A_enhancer_promoter) <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)


ATAC_ep_pgp_rows <- unique(All_lins_ATAC_table[grep("pgp-", All_lins_ATAC_table$Gene), 3])

chrX_ATAC_A_pgp_rows <- ATAC_ep_pgp_rows[grep("chrX:", ATAC_ep_pgp_rows)]

ChrX_ATAC_A_pgp_z_scores <- ChrX_ATAC_Z_score_table_A_enhancer_promoter[chrX_ATAC_A_pgp_rows, ]


genes <- All_lins_ATAC_table[All_lins_ATAC_table$Locus %in% rownames(ChrX_ATAC_A_pgp_z_scores) & All_lins_ATAC_table$Lineage == "A", ]

pgp_gene_frame <- genes[grep("pgp-", genes$Gene), ]


ChromX_ATAC_A_pgp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 1]),
    cbind(rep(4, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 2]),
    cbind(rep(6, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 3]),
    cbind(rep(8, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 4]),
    cbind(rep(10, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 5]),
    cbind(rep(12, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 6]),
    cbind(rep(14, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 7]),
    cbind(rep(16, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 8]),
    cbind(rep(18, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 9]),
    cbind(rep(20, length = nrow(ChrX_ATAC_A_pgp_z_scores)), ChrX_ATAC_A_pgp_z_scores[, 10]))


pgp_ep_loci <- paste(rownames(ChromX_ATAC_A_pgp_Z_scores), pgp_gene_frame$Gene, sep = ":")

ChromX_ATAC_A_pgp_Z_scores  <- data.frame(rep(pgp_gene_frame$Gene, 10), ChromX_ATAC_A_pgp_Z_scores, pgp_ep_loci)

colnames(ChromX_ATAC_A_pgp_Z_scores) <- c("gene_name", "Generation", "Expression", "ep_locus")


pgp_trunc_names <- c("pgp-6", "pgp-7", "pgp-8")

ChromX_ATAC_A_pgp_Z_scores_trunc <- ChromX_ATAC_A_pgp_Z_scores[ChromX_ATAC_A_pgp_Z_scores$gene_name %in% pgp_trunc_names, ]

#-----

# Normalise RNA in context of X chromosome only

ChrX_Lineage_A_RNA_norm <-  
  
  Lineage_A_RNA_norm[grep("chrX:", rownames(Lineage_A_RNA_norm)), ]

ChrX_RNA_Z_score_table_A <- matrix(0, ncol = 10, nrow = nrow(ChrX_Lineage_A_RNA_norm))

ChrX_normalized_table_A2_to_A20 <- ChrX_Lineage_A_RNA_norm[, 2:11]

for (i in 1:ncol(ChrX_normalized_table_A2_to_A20)) {
  LM_temp <-
    lm(log2(ChrX_normalized_table_A2_to_A20[, i]) ~ log2(ChrX_Lineage_A_RNA_norm[, 1]))
  Residuals <- residuals(LM_temp)
  ChrX_RNA_Z_score_table_A[, i] <- (Residuals - mean(Residuals)) / sd(Residuals)
}


rownames(ChrX_RNA_Z_score_table_A) <- rownames(ChrX_normalized_table_A2_to_A20)
colnames(ChrX_RNA_Z_score_table_A) <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)


ChrX_RNA_A_Z_score_pgp_rows <-
  ChrX_RNA_Z_score_table_A[grep("pgp-", rownames(ChrX_RNA_Z_score_table_A)), ]


ChrX_pgp_genes <- c("pgp-10", "pgp-12", "pgp-13", "pgp-14", "pgp-3", "pgp-4", "pgp-5", "pgp-6", "pgp-7", "pgp-8")



ChromX_RNA_A_pgp_Z_scores_rows <-
  
  rbind(
    cbind(rep(2, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 1]),
    cbind(rep(4, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 2]),
    cbind(rep(6, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 3]),
    cbind(rep(8, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 4]),
    cbind(rep(10, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 5]),
    cbind(rep(12, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 6]),
    cbind(rep(14, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 7]),
    cbind(rep(16, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 8]),
    cbind(rep(18, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 9]),
    cbind(rep(20, length = nrow(ChrX_RNA_A_Z_score_pgp_rows)), ChrX_RNA_A_Z_score_pgp_rows[, 10]))




ChrX_RNA_A_pgp_Z_scores <- cbind(rownames(ChromX_RNA_A_pgp_Z_scores_rows), rep(ChrX_pgp_genes, 10), data.frame(ChromX_RNA_A_pgp_Z_scores_rows, row.names=NULL)) 

colnames(ChrX_RNA_A_pgp_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")

pgp_trunc_names <- c("pgp-5", "pgp-6", "pgp-7", "pgp-8")

ChrX_RNA_A_pgp_Z_scores_trunc <- ChrX_RNA_A_pgp_Z_scores[ChrX_RNA_A_pgp_Z_scores$gene_name %in% pgp_trunc_names, ]





# Supplementary Figure 8

# Plot for X chromosome restricted Z scores ATAC and RNA Lineage A

# library(ggplot2)

# test <-

ggplot(ChrX_RNA_A_pgp_Z_scores_trunc, aes(x=Generation, y=Expression, group=gene_name, color=gene_name)) +
  geom_line()+
  # scale_color_manual(values = mycolors)+
  geom_line(data = ChromX_ATAC_A_pgp_Z_scores_trunc, aes(x=Generation, y=Expression, group=ep_locus, color=ep_locus))+
  ggtitle("Supplementary Figure 8") +
  labs(x = "Generations)",
       y = "Z score for expression/chromatin",
       color = "pgp genes")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=c(-2.25, 2.25), col = "black")+
  scale_x_continuous(breaks=seq(2, 20, 2), limits = c(2, 20))+
  scale_y_continuous(breaks = c(-5, -4, -3, -2.25, -2, -1, 0, 1, 2, 2.25, 3, 4, 5, 6), limits = c(-5, 6),
                     labels = c("-5", "-4", "-3", "Epimutations below cut off", "-2", "-1", "0", "1", "2", "Epimutations above cut off", "3", "4", "5", "6"))+
  # geom_rect(aes(xmin=13.5, xmax=14.5, ymin=-Inf, ymax=Inf), colour = NA, fill="gray", alpha=0.01)+
  theme(plot.margin = margin(1,1,1,1,"cm"))



#---------------------------------------------------------------------------------------------------------

# Supplementary Figure 9

# Investigating fluctuations in Xenobiotic defence gene families and comparing to housekeeping genes


#-----------------------------------

# NHR Gene Family

nhr_genes <- rownames(RNA_Z_score_table_A)[grep("nhr-", rownames(RNA_Z_score_table_A))]

RNA_A_nhr <- RNA_Z_score_table_A[nhr_genes,]

RNA_A_nhr_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_A_nhr)), RNA_A_nhr[, 1]),
    cbind(rep(4, length = nrow(RNA_A_nhr)), RNA_A_nhr[, 2]),
    cbind(rep(6, length = nrow(RNA_A_nhr)), RNA_A_nhr[, 3]),
    cbind(rep(8, length = nrow(RNA_A_nhr)), RNA_A_nhr[, 4]),
    cbind(rep(10, length = nrow(RNA_A_nhr)), RNA_A_nhr[, 5]),
    cbind(rep(12, length = nrow(RNA_A_nhr)), RNA_A_nhr[, 6]),
    cbind(rep(16, length = nrow(RNA_A_nhr)), RNA_A_nhr[, 7]),
    cbind(rep(18, length = nrow(RNA_A_nhr)), RNA_A_nhr[, 8]),
    cbind(rep(20, length = nrow(RNA_A_nhr)), RNA_A_nhr[, 9]))

RNA_A_nhr_genes <- strsplit(nhr_genes, split = ":")

nhr_genes <- c()

for(i in 1:length(RNA_A_nhr_genes)){
  
  nhr_genes <- c(nhr_genes, RNA_A_nhr_genes[[i]][4])
  
}



RNA_A_nhr_Z_scores <- cbind(paste("RNA", rownames(RNA_A_nhr_Z_scores), sep = ":"), rep(nhr_genes, 9), data.frame(RNA_A_nhr_Z_scores, row.names=NULL)) 

colnames(RNA_A_nhr_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")







nhr_genes <- rownames(RNA_Z_score_table_B)[grep("nhr-", rownames(RNA_Z_score_table_B))]


RNA_B_nhr <- RNA_Z_score_table_B[nhr_genes,]


RNA_B_nhr_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_B_nhr)), RNA_B_nhr[, 1]),
    cbind(rep(6, length = nrow(RNA_B_nhr)), RNA_B_nhr[, 2]),
    cbind(rep(8, length = nrow(RNA_B_nhr)), RNA_B_nhr[, 3]),
    cbind(rep(10, length = nrow(RNA_B_nhr)), RNA_B_nhr[, 4]),
    cbind(rep(12, length = nrow(RNA_B_nhr)), RNA_B_nhr[, 5]),
    cbind(rep(14, length = nrow(RNA_B_nhr)), RNA_B_nhr[, 6]),
    cbind(rep(16, length = nrow(RNA_B_nhr)), RNA_B_nhr[, 7]),
    cbind(rep(18, length = nrow(RNA_B_nhr)), RNA_B_nhr[, 8]),
    cbind(rep(20, length = nrow(RNA_B_nhr)), RNA_B_nhr[, 9]))


RNA_B_nhr_genes <- strsplit(nhr_genes, split = ":")

nhr_genes <- c()

for(i in 1:length(RNA_B_nhr_genes)){
  
  nhr_genes <- c(nhr_genes, RNA_B_nhr_genes[[i]][4])
  
}


RNA_B_nhr_Z_scores <- cbind(rownames(RNA_B_nhr_Z_scores), rep(nhr_genes, 9), data.frame(RNA_B_nhr_Z_scores, row.names=NULL)) 

colnames(RNA_B_nhr_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")







nhr_genes <- rownames(RNA_Z_score_table_C)[grep("nhr-", rownames(RNA_Z_score_table_C))]

RNA_C_nhr <- RNA_Z_score_table_C[nhr_genes,]


RNA_C_nhr_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_C_nhr)), RNA_C_nhr[, 1]),
    cbind(rep(4, length = nrow(RNA_C_nhr)), RNA_C_nhr[, 2]),
    cbind(rep(6, length = nrow(RNA_C_nhr)), RNA_C_nhr[, 3]),
    cbind(rep(8, length = nrow(RNA_C_nhr)), RNA_C_nhr[, 4]),
    cbind(rep(10, length = nrow(RNA_C_nhr)), RNA_C_nhr[, 5]),
    cbind(rep(12, length = nrow(RNA_C_nhr)), RNA_C_nhr[, 6]),
    cbind(rep(16, length = nrow(RNA_C_nhr)), RNA_C_nhr[, 7]),
    cbind(rep(18, length = nrow(RNA_C_nhr)), RNA_C_nhr[, 8]),
    cbind(rep(20, length = nrow(RNA_C_nhr)), RNA_C_nhr[, 9]))


RNA_C_nhr_genes <- strsplit(nhr_genes, split = ":")

nhr_genes <- c()

for(i in 1:length(RNA_C_nhr_genes)){
  
  nhr_genes <- c(nhr_genes, RNA_C_nhr_genes[[i]][4])
  
}


RNA_C_nhr_Z_scores <- cbind(rownames(RNA_C_nhr_Z_scores), rep(nhr_genes, 9), data.frame(RNA_C_nhr_Z_scores, row.names=NULL)) 

colnames(RNA_C_nhr_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")





tmp_A <- data.frame(rep("A"), RNA_A_nhr_Z_scores)

colnames(tmp_A) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_B <- data.frame(rep("B"), RNA_B_nhr_Z_scores)

colnames(tmp_B) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_C <- data.frame(rep("C"), RNA_C_nhr_Z_scores)

colnames(tmp_C) <- c("Lin", "locus", "gene", "Generation", "Expression")


nhr_all_lins_tmp <- rbind(tmp_A, tmp_B, tmp_C)


get_direction <- c()

for(i in 1:nrow(nhr_all_lins_tmp)){
  
  direction <- "baseline"
  
  if(as.numeric(nhr_all_lins_tmp[i, 5]) > 2.25){
    
    direction <- "Up"
  }
  
  if(as.numeric(nhr_all_lins_tmp[i, 5]) < -2.25){
    direction <- "Down"
  }
  
  get_direction <- c(get_direction, direction)
  
}

nhr_all_lins_tmp <- cbind(nhr_all_lins_tmp, get_direction)


nhr_all_lins_tmp <- data.frame(paste(nhr_all_lins_tmp$Lin, nhr_all_lins_tmp$get_direction, sep = "_"), nhr_all_lins_tmp)

colnames(nhr_all_lins_tmp) <- c("Lin_direction", "Lin", "locus", "gene", "Generation", "Expression", "direction")




nhr_all_lins_tmp$Lin_direction <- factor(nhr_all_lins_tmp$Lin_direction, levels = c("A_Down", "A_baseline", "A_Up",
                                                                                    "B_Down", "B_baseline", "B_Up",
                                                                                    "C_Down", "C_baseline", "C_Up"))


#-----------------------------------

# CYP Gene Family


cyp_genes <- rownames(RNA_Z_score_table_A)[grep("cyp-", rownames(RNA_Z_score_table_A))]

RNA_A_cyp <- RNA_Z_score_table_A[cyp_genes,]

RNA_A_cyp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_A_cyp)), RNA_A_cyp[, 1]),
    cbind(rep(4, length = nrow(RNA_A_cyp)), RNA_A_cyp[, 2]),
    cbind(rep(6, length = nrow(RNA_A_cyp)), RNA_A_cyp[, 3]),
    cbind(rep(8, length = nrow(RNA_A_cyp)), RNA_A_cyp[, 4]),
    cbind(rep(10, length = nrow(RNA_A_cyp)), RNA_A_cyp[, 5]),
    cbind(rep(12, length = nrow(RNA_A_cyp)), RNA_A_cyp[, 6]),
    cbind(rep(16, length = nrow(RNA_A_cyp)), RNA_A_cyp[, 7]),
    cbind(rep(18, length = nrow(RNA_A_cyp)), RNA_A_cyp[, 8]),
    cbind(rep(20, length = nrow(RNA_A_cyp)), RNA_A_cyp[, 9]))

RNA_A_cyp_genes <- strsplit(cyp_genes, split = ":")

cyp_genes <- c()

for(i in 1:length(RNA_A_cyp_genes)){
  
  cyp_genes <- c(cyp_genes, RNA_A_cyp_genes[[i]][4])
  
}



RNA_A_cyp_Z_scores <- cbind(paste("RNA", rownames(RNA_A_cyp_Z_scores), sep = ":"), rep(cyp_genes, 9), data.frame(RNA_A_cyp_Z_scores, row.names=NULL)) 

colnames(RNA_A_cyp_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")





cyp_genes <- rownames(RNA_Z_score_table_B)[grep("cyp-", rownames(RNA_Z_score_table_B))]


RNA_B_cyp <- RNA_Z_score_table_B[cyp_genes,]


RNA_B_cyp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_B_cyp)), RNA_B_cyp[, 1]),
    cbind(rep(6, length = nrow(RNA_B_cyp)), RNA_B_cyp[, 2]),
    cbind(rep(8, length = nrow(RNA_B_cyp)), RNA_B_cyp[, 3]),
    cbind(rep(10, length = nrow(RNA_B_cyp)), RNA_B_cyp[, 4]),
    cbind(rep(12, length = nrow(RNA_B_cyp)), RNA_B_cyp[, 5]),
    cbind(rep(14, length = nrow(RNA_B_cyp)), RNA_B_cyp[, 6]),
    cbind(rep(16, length = nrow(RNA_B_cyp)), RNA_B_cyp[, 7]),
    cbind(rep(18, length = nrow(RNA_B_cyp)), RNA_B_cyp[, 8]),
    cbind(rep(20, length = nrow(RNA_B_cyp)), RNA_B_cyp[, 9]))


RNA_B_cyp_genes <- strsplit(cyp_genes, split = ":")

cyp_genes <- c()

for(i in 1:length(RNA_B_cyp_genes)){
  
  cyp_genes <- c(cyp_genes, RNA_B_cyp_genes[[i]][4])
  
}


RNA_B_cyp_Z_scores <- cbind(rownames(RNA_B_cyp_Z_scores), rep(cyp_genes, 9), data.frame(RNA_B_cyp_Z_scores, row.names=NULL)) 

colnames(RNA_B_cyp_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")




cyp_genes <- rownames(RNA_Z_score_table_C)[grep("cyp-", rownames(RNA_Z_score_table_C))]

RNA_C_cyp <- RNA_Z_score_table_C[cyp_genes,]


RNA_C_cyp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_C_cyp)), RNA_C_cyp[, 1]),
    cbind(rep(4, length = nrow(RNA_C_cyp)), RNA_C_cyp[, 2]),
    cbind(rep(6, length = nrow(RNA_C_cyp)), RNA_C_cyp[, 3]),
    cbind(rep(8, length = nrow(RNA_C_cyp)), RNA_C_cyp[, 4]),
    cbind(rep(10, length = nrow(RNA_C_cyp)), RNA_C_cyp[, 5]),
    cbind(rep(12, length = nrow(RNA_C_cyp)), RNA_C_cyp[, 6]),
    cbind(rep(16, length = nrow(RNA_C_cyp)), RNA_C_cyp[, 7]),
    cbind(rep(18, length = nrow(RNA_C_cyp)), RNA_C_cyp[, 8]),
    cbind(rep(20, length = nrow(RNA_C_cyp)), RNA_C_cyp[, 9]))


RNA_C_cyp_genes <- strsplit(cyp_genes, split = ":")

cyp_genes <- c()

for(i in 1:length(RNA_C_cyp_genes)){
  
  cyp_genes <- c(cyp_genes, RNA_C_cyp_genes[[i]][4])
  
}


RNA_C_cyp_Z_scores <- cbind(rownames(RNA_C_cyp_Z_scores), rep(cyp_genes, 9), data.frame(RNA_C_cyp_Z_scores, row.names=NULL)) 

colnames(RNA_C_cyp_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")





tmp_A <- data.frame(rep("A"), RNA_A_cyp_Z_scores)

colnames(tmp_A) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_B <- data.frame(rep("B"), RNA_B_cyp_Z_scores)

colnames(tmp_B) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_C <- data.frame(rep("C"), RNA_C_cyp_Z_scores)

colnames(tmp_C) <- c("Lin", "locus", "gene", "Generation", "Expression")


cyp_all_lins_tmp <- rbind(tmp_A, tmp_B, tmp_C)


get_direction <- c()

for(i in 1:nrow(cyp_all_lins_tmp)){
  
  direction <- "baseline"
  
  if(as.numeric(cyp_all_lins_tmp[i, 5]) > 2.25){
    
    direction <- "Up"
  }
  
  if(as.numeric(cyp_all_lins_tmp[i, 5]) < -2.25){
    direction <- "Down"
  }
  
  get_direction <- c(get_direction, direction)
  
}

cyp_all_lins_tmp <- cbind(cyp_all_lins_tmp, get_direction)


cyp_all_lins_tmp <- data.frame(paste(cyp_all_lins_tmp$Lin, cyp_all_lins_tmp$get_direction, sep = "_"), cyp_all_lins_tmp)

colnames(cyp_all_lins_tmp) <- c("Lin_direction", "Lin", "locus", "gene", "Generation", "Expression", "direction")


#-----------------------------------

# UGT Gene Family

ugt_genes <- rownames(RNA_Z_score_table_A)[grep("ugt-", rownames(RNA_Z_score_table_A))]

RNA_A_ugt <- RNA_Z_score_table_A[ugt_genes,]

RNA_A_ugt_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_A_ugt)), RNA_A_ugt[, 1]),
    cbind(rep(4, length = nrow(RNA_A_ugt)), RNA_A_ugt[, 2]),
    cbind(rep(6, length = nrow(RNA_A_ugt)), RNA_A_ugt[, 3]),
    cbind(rep(8, length = nrow(RNA_A_ugt)), RNA_A_ugt[, 4]),
    cbind(rep(10, length = nrow(RNA_A_ugt)), RNA_A_ugt[, 5]),
    cbind(rep(12, length = nrow(RNA_A_ugt)), RNA_A_ugt[, 6]),
    cbind(rep(16, length = nrow(RNA_A_ugt)), RNA_A_ugt[, 7]),
    cbind(rep(18, length = nrow(RNA_A_ugt)), RNA_A_ugt[, 8]),
    cbind(rep(20, length = nrow(RNA_A_ugt)), RNA_A_ugt[, 9]))

RNA_A_ugt_genes <- strsplit(ugt_genes, split = ":")

ugt_genes <- c()

for(i in 1:length(RNA_A_ugt_genes)){
  
  ugt_genes <- c(ugt_genes, RNA_A_ugt_genes[[i]][4])
  
}



RNA_A_ugt_Z_scores <- cbind(paste("RNA", rownames(RNA_A_ugt_Z_scores), sep = ":"), rep(ugt_genes, 9), data.frame(RNA_A_ugt_Z_scores, row.names=NULL)) 

colnames(RNA_A_ugt_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")



ugt_genes <- rownames(RNA_Z_score_table_B)[grep("ugt-", rownames(RNA_Z_score_table_B))]


RNA_B_ugt <- RNA_Z_score_table_B[ugt_genes,]


RNA_B_ugt_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_B_ugt)), RNA_B_ugt[, 1]),
    cbind(rep(6, length = nrow(RNA_B_ugt)), RNA_B_ugt[, 2]),
    cbind(rep(8, length = nrow(RNA_B_ugt)), RNA_B_ugt[, 3]),
    cbind(rep(10, length = nrow(RNA_B_ugt)), RNA_B_ugt[, 4]),
    cbind(rep(12, length = nrow(RNA_B_ugt)), RNA_B_ugt[, 5]),
    cbind(rep(14, length = nrow(RNA_B_ugt)), RNA_B_ugt[, 6]),
    cbind(rep(16, length = nrow(RNA_B_ugt)), RNA_B_ugt[, 7]),
    cbind(rep(18, length = nrow(RNA_B_ugt)), RNA_B_ugt[, 8]),
    cbind(rep(20, length = nrow(RNA_B_ugt)), RNA_B_ugt[, 9]))


RNA_B_ugt_genes <- strsplit(ugt_genes, split = ":")

ugt_genes <- c()

for(i in 1:length(RNA_B_ugt_genes)){
  
  ugt_genes <- c(ugt_genes, RNA_B_ugt_genes[[i]][4])
  
}


RNA_B_ugt_Z_scores <- cbind(rownames(RNA_B_ugt_Z_scores), rep(ugt_genes, 9), data.frame(RNA_B_ugt_Z_scores, row.names=NULL)) 

colnames(RNA_B_ugt_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")







ugt_genes <- rownames(RNA_Z_score_table_C)[grep("ugt-", rownames(RNA_Z_score_table_C))]

RNA_C_ugt <- RNA_Z_score_table_C[ugt_genes,]


RNA_C_ugt_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_C_ugt)), RNA_C_ugt[, 1]),
    cbind(rep(4, length = nrow(RNA_C_ugt)), RNA_C_ugt[, 2]),
    cbind(rep(6, length = nrow(RNA_C_ugt)), RNA_C_ugt[, 3]),
    cbind(rep(8, length = nrow(RNA_C_ugt)), RNA_C_ugt[, 4]),
    cbind(rep(10, length = nrow(RNA_C_ugt)), RNA_C_ugt[, 5]),
    cbind(rep(12, length = nrow(RNA_C_ugt)), RNA_C_ugt[, 6]),
    cbind(rep(16, length = nrow(RNA_C_ugt)), RNA_C_ugt[, 7]),
    cbind(rep(18, length = nrow(RNA_C_ugt)), RNA_C_ugt[, 8]),
    cbind(rep(20, length = nrow(RNA_C_ugt)), RNA_C_ugt[, 9]))


RNA_C_ugt_genes <- strsplit(ugt_genes, split = ":")

ugt_genes <- c()

for(i in 1:length(RNA_C_ugt_genes)){
  
  ugt_genes <- c(ugt_genes, RNA_C_ugt_genes[[i]][4])
  
}


RNA_C_ugt_Z_scores <- cbind(rownames(RNA_C_ugt_Z_scores), rep(ugt_genes, 9), data.frame(RNA_C_ugt_Z_scores, row.names=NULL)) 

colnames(RNA_C_ugt_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")






tmp_A <- data.frame(rep("A"), RNA_A_ugt_Z_scores)

colnames(tmp_A) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_B <- data.frame(rep("B"), RNA_B_ugt_Z_scores)

colnames(tmp_B) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_C <- data.frame(rep("C"), RNA_C_ugt_Z_scores)

colnames(tmp_C) <- c("Lin", "locus", "gene", "Generation", "Expression")


ugt_all_lins_tmp <- rbind(tmp_A, tmp_B, tmp_C)


get_direction <- c()

for(i in 1:nrow(ugt_all_lins_tmp)){
  
  direction <- "baseline"
  
  if(as.numeric(ugt_all_lins_tmp[i, 5]) > 2.25){
    
    direction <- "Up"
  }
  
  if(as.numeric(ugt_all_lins_tmp[i, 5]) < -2.25){
    direction <- "Down"
  }
  
  get_direction <- c(get_direction, direction)
  
}

ugt_all_lins_tmp <- cbind(ugt_all_lins_tmp, get_direction)


ugt_all_lins_tmp <- data.frame(paste(ugt_all_lins_tmp$Lin, ugt_all_lins_tmp$get_direction, sep = "_"), ugt_all_lins_tmp)

colnames(ugt_all_lins_tmp) <- c("Lin_direction", "Lin", "locus", "gene", "Generation", "Expression", "direction")


#-----------------------------------

# PGP Gene Family


pgp_genes <- rownames(RNA_Z_score_table_A)[grep("pgp-", rownames(RNA_Z_score_table_A))]

RNA_A_pgp <- RNA_Z_score_table_A[pgp_genes,]

RNA_A_pgp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 1]),
    cbind(rep(4, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 2]),
    cbind(rep(6, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 3]),
    cbind(rep(8, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 4]),
    cbind(rep(10, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 5]),
    cbind(rep(12, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 6]),
    cbind(rep(16, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 7]),
    cbind(rep(18, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 8]),
    cbind(rep(20, length = nrow(RNA_A_pgp)), RNA_A_pgp[, 9]))

RNA_A_pgp_genes <- strsplit(pgp_genes, split = ":")

pgp_genes <- c()

for(i in 1:length(RNA_A_pgp_genes)){
  
  pgp_genes <- c(pgp_genes, RNA_A_pgp_genes[[i]][4])
  
}



RNA_A_pgp_Z_scores <- cbind(paste("RNA", rownames(RNA_A_pgp_Z_scores), sep = ":"), rep(pgp_genes, 9), data.frame(RNA_A_pgp_Z_scores, row.names=NULL)) 

colnames(RNA_A_pgp_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")







pgp_genes <- rownames(RNA_Z_score_table_B)[grep("pgp-", rownames(RNA_Z_score_table_B))]


RNA_B_pgp <- RNA_Z_score_table_B[pgp_genes,]


RNA_B_pgp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_B_pgp)), RNA_B_pgp[, 1]),
    cbind(rep(6, length = nrow(RNA_B_pgp)), RNA_B_pgp[, 2]),
    cbind(rep(8, length = nrow(RNA_B_pgp)), RNA_B_pgp[, 3]),
    cbind(rep(10, length = nrow(RNA_B_pgp)), RNA_B_pgp[, 4]),
    cbind(rep(12, length = nrow(RNA_B_pgp)), RNA_B_pgp[, 5]),
    cbind(rep(14, length = nrow(RNA_B_pgp)), RNA_B_pgp[, 6]),
    cbind(rep(16, length = nrow(RNA_B_pgp)), RNA_B_pgp[, 7]),
    cbind(rep(18, length = nrow(RNA_B_pgp)), RNA_B_pgp[, 8]),
    cbind(rep(20, length = nrow(RNA_B_pgp)), RNA_B_pgp[, 9]))


RNA_B_pgp_genes <- strsplit(pgp_genes, split = ":")

pgp_genes <- c()

for(i in 1:length(RNA_B_pgp_genes)){
  
  pgp_genes <- c(pgp_genes, RNA_B_pgp_genes[[i]][4])
  
}


RNA_B_pgp_Z_scores <- cbind(rownames(RNA_B_pgp_Z_scores), rep(pgp_genes, 9), data.frame(RNA_B_pgp_Z_scores, row.names=NULL)) 

colnames(RNA_B_pgp_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")







pgp_genes <- rownames(RNA_Z_score_table_C)[grep("pgp-", rownames(RNA_Z_score_table_C))]

RNA_C_pgp <- RNA_Z_score_table_C[pgp_genes,]


RNA_C_pgp_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_C_pgp)), RNA_C_pgp[, 1]),
    cbind(rep(4, length = nrow(RNA_C_pgp)), RNA_C_pgp[, 2]),
    cbind(rep(6, length = nrow(RNA_C_pgp)), RNA_C_pgp[, 3]),
    cbind(rep(8, length = nrow(RNA_C_pgp)), RNA_C_pgp[, 4]),
    cbind(rep(10, length = nrow(RNA_C_pgp)), RNA_C_pgp[, 5]),
    cbind(rep(12, length = nrow(RNA_C_pgp)), RNA_C_pgp[, 6]),
    cbind(rep(16, length = nrow(RNA_C_pgp)), RNA_C_pgp[, 7]),
    cbind(rep(18, length = nrow(RNA_C_pgp)), RNA_C_pgp[, 8]),
    cbind(rep(20, length = nrow(RNA_C_pgp)), RNA_C_pgp[, 9]))


RNA_C_pgp_genes <- strsplit(pgp_genes, split = ":")

pgp_genes <- c()

for(i in 1:length(RNA_C_pgp_genes)){
  
  pgp_genes <- c(pgp_genes, RNA_C_pgp_genes[[i]][4])
  
}


RNA_C_pgp_Z_scores <- cbind(rownames(RNA_C_pgp_Z_scores), rep(pgp_genes, 9), data.frame(RNA_C_pgp_Z_scores, row.names=NULL)) 

colnames(RNA_C_pgp_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")



tmp_A <- data.frame(rep("A"), RNA_A_pgp_Z_scores)

colnames(tmp_A) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_B <- data.frame(rep("B"), RNA_B_pgp_Z_scores)

colnames(tmp_B) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_C <- data.frame(rep("C"), RNA_C_pgp_Z_scores)

colnames(tmp_C) <- c("Lin", "locus", "gene", "Generation", "Expression")


pgp_all_lins_tmp <- rbind(tmp_A, tmp_B, tmp_C)


get_direction <- c()

for(i in 1:nrow(pgp_all_lins_tmp)){
  
  direction <- "baseline"
  
  if(as.numeric(pgp_all_lins_tmp[i, 5]) > 2.25){
    
    direction <- "Up"
  }
  
  if(as.numeric(pgp_all_lins_tmp[i, 5]) < -2.25){
    direction <- "Down"
  }
  
  get_direction <- c(get_direction, direction)
  
}

pgp_all_lins_tmp <- cbind(pgp_all_lins_tmp, get_direction)


pgp_all_lins_tmp <- data.frame(paste(pgp_all_lins_tmp$Lin, pgp_all_lins_tmp$get_direction, sep = "_"), pgp_all_lins_tmp)

colnames(pgp_all_lins_tmp) <- c("Lin_direction", "Lin", "locus", "gene", "Generation", "Expression", "direction")


#------------------------------------------------------------------------

# Housekeeping Genes 


# 13 Validated housekeeping genes obtained from Tao et al. 2020


# Lineage A

rpl_genes <- c("chrIV:12390229:12391053:rps-23", 
               "chrI:14759929:14760626:rps-26", 
               "chrV:103393:104064:rps-27", 
               "chrV:15000017:15000594:rps-16", 
               "chrIV:7925297:7926391:rps-2", 
               "chrIV:7083686:7084682:rps-4",
               "chrI:6220092:6220766:rps-17", 
               "chrI:4585114:4586184:rpl-24.1",
               "chrI:1834877:1835438:rpl-27", 
               "chrII:7105562:7106316:rpl-33",
               "chrIII:7180242:7180647:rpl-36",
               "chrIII:7855133:7855646:rpl-35",
               "chrIV:653427:654576:rpl-15" 
)




RNA_A_rpl <- RNA_Z_score_table_A[rpl_genes,]

RNA_A_rpl_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_A_rpl)), RNA_A_rpl[, 1]),
    cbind(rep(4, length = nrow(RNA_A_rpl)), RNA_A_rpl[, 2]),
    cbind(rep(6, length = nrow(RNA_A_rpl)), RNA_A_rpl[, 3]),
    cbind(rep(8, length = nrow(RNA_A_rpl)), RNA_A_rpl[, 4]),
    cbind(rep(10, length = nrow(RNA_A_rpl)), RNA_A_rpl[, 5]),
    cbind(rep(12, length = nrow(RNA_A_rpl)), RNA_A_rpl[, 6]),
    cbind(rep(16, length = nrow(RNA_A_rpl)), RNA_A_rpl[, 7]),
    cbind(rep(18, length = nrow(RNA_A_rpl)), RNA_A_rpl[, 8]),
    cbind(rep(20, length = nrow(RNA_A_rpl)), RNA_A_rpl[, 9]))


RNA_A_rpl_genes <- strsplit(rpl_genes, split = ":")

rpl_genes <- c()

for(i in 1:length(RNA_A_rpl_genes)){
  
  rpl_genes <- c(rpl_genes, RNA_A_rpl_genes[[i]][4])
  
}




RNA_A_rpl_Z_scores <- cbind(paste("RNA", rownames(RNA_A_rpl_Z_scores), sep = ":"), rep(rpl_genes, 9), data.frame(RNA_A_rpl_Z_scores, row.names=NULL)) 

colnames(RNA_A_rpl_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")





# Lineage B

rpl_genes <- c("chrIV:12390229:12391053:rps-23", 
               "chrI:14759929:14760626:rps-26", 
               "chrV:103393:104064:rps-27", 
               "chrV:15000017:15000594:rps-16", 
               "chrIV:7925297:7926391:rps-2", 
               "chrIV:7083686:7084682:rps-4",
               "chrI:6220092:6220766:rps-17", 
               "chrI:4585114:4586184:rpl-24.1",
               "chrI:1834877:1835438:rpl-27", 
               "chrII:7105562:7106316:rpl-33",
               "chrIII:7180242:7180647:rpl-36",
               "chrIII:7855133:7855646:rpl-35",
               "chrIV:653427:654576:rpl-15" 
)


RNA_B_rpl <- RNA_Z_score_table_B[rpl_genes,]


RNA_B_rpl_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_B_rpl)), RNA_B_rpl[, 1]),
    cbind(rep(6, length = nrow(RNA_B_rpl)), RNA_B_rpl[, 2]),
    cbind(rep(8, length = nrow(RNA_B_rpl)), RNA_B_rpl[, 3]),
    cbind(rep(10, length = nrow(RNA_B_rpl)), RNA_B_rpl[, 4]),
    cbind(rep(12, length = nrow(RNA_B_rpl)), RNA_B_rpl[, 5]),
    cbind(rep(14, length = nrow(RNA_B_rpl)), RNA_B_rpl[, 6]),
    cbind(rep(16, length = nrow(RNA_B_rpl)), RNA_B_rpl[, 7]),
    cbind(rep(18, length = nrow(RNA_B_rpl)), RNA_B_rpl[, 8]),
    cbind(rep(20, length = nrow(RNA_B_rpl)), RNA_B_rpl[, 9]))


RNA_B_rpl_genes <- strsplit(rpl_genes, split = ":")

rpl_genes <- c()

for(i in 1:length(RNA_B_rpl_genes)){
  
  rpl_genes <- c(rpl_genes, RNA_B_rpl_genes[[i]][4])
  
}


RNA_B_rpl_Z_scores <- cbind(rownames(RNA_B_rpl_Z_scores), rep(rpl_genes, 9), data.frame(RNA_B_rpl_Z_scores, row.names=NULL)) 

colnames(RNA_B_rpl_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")




# Lineage C

rpl_genes <- c("chrIV:12390229:12391053:rps-23", 
               "chrI:14759929:14760626:rps-26", 
               "chrV:103393:104064:rps-27", 
               "chrV:15000017:15000594:rps-16", 
               "chrIV:7925297:7926391:rps-2", 
               "chrIV:7083686:7084682:rps-4",
               "chrI:6220092:6220766:rps-17", 
               "chrI:4585114:4586184:rpl-24.1",
               "chrI:1834877:1835438:rpl-27", 
               "chrII:7105562:7106316:rpl-33",
               "chrIII:7180242:7180647:rpl-36",
               "chrIII:7855133:7855646:rpl-35",
               "chrIV:653427:654576:rpl-15" 
)

RNA_C_rpl <- RNA_Z_score_table_C[rpl_genes,]


RNA_C_rpl_Z_scores <-
  
  rbind(
    cbind(rep(2, length = nrow(RNA_C_rpl)), RNA_C_rpl[, 1]),
    cbind(rep(4, length = nrow(RNA_C_rpl)), RNA_C_rpl[, 2]),
    cbind(rep(6, length = nrow(RNA_C_rpl)), RNA_C_rpl[, 3]),
    cbind(rep(8, length = nrow(RNA_C_rpl)), RNA_C_rpl[, 4]),
    cbind(rep(10, length = nrow(RNA_C_rpl)), RNA_C_rpl[, 5]),
    cbind(rep(12, length = nrow(RNA_C_rpl)), RNA_C_rpl[, 6]),
    cbind(rep(16, length = nrow(RNA_C_rpl)), RNA_C_rpl[, 7]),
    cbind(rep(18, length = nrow(RNA_C_rpl)), RNA_C_rpl[, 8]),
    cbind(rep(20, length = nrow(RNA_C_rpl)), RNA_C_rpl[, 9]))


RNA_rpl_genes <- strsplit(rpl_genes, split = ":")

rpl_genes <- c()

for(i in 1:length(RNA_rpl_genes)){
  
  rpl_genes <- c(rpl_genes, RNA_rpl_genes[[i]][4])
  
}


RNA_C_rpl_Z_scores <- cbind(rownames(RNA_C_rpl_Z_scores), rep(rpl_genes, 9), data.frame(RNA_C_rpl_Z_scores, row.names=NULL)) 

colnames(RNA_C_rpl_Z_scores) <- c("locus_name", "gene_name", "Generation", "Expression")




tmp_A <- data.frame(rep("A"), RNA_A_rpl_Z_scores)

colnames(tmp_A) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_B <- data.frame(rep("B"), RNA_B_rpl_Z_scores)

colnames(tmp_B) <- c("Lin", "locus", "gene", "Generation", "Expression")

tmp_C <- data.frame(rep("C"), RNA_C_rpl_Z_scores)

colnames(tmp_C) <- c("Lin", "locus", "gene", "Generation", "Expression")


rpl_all_lins_tmp <- rbind(tmp_A, tmp_B, tmp_C)


get_direction <- c()

for(i in 1:nrow(rpl_all_lins_tmp)){
  
  direction <- "baseline"
  
  if(as.numeric(rpl_all_lins_tmp[i, 5]) > 2.25){
    
    direction <- "Up"
  }
  
  if(as.numeric(rpl_all_lins_tmp[i, 5]) < -2.25){
    direction <- "Down"
  }
  
  get_direction <- c(get_direction, direction)
  
}

rpl_all_lins_tmp <- cbind(rpl_all_lins_tmp, get_direction)


rpl_all_lins_tmp <- data.frame(paste(rpl_all_lins_tmp$Lin, rpl_all_lins_tmp$get_direction, sep = "_"), rpl_all_lins_tmp)

colnames(rpl_all_lins_tmp) <- c("Lin_direction", "Lin", "locus", "gene", "Generation", "Expression", "direction")


rpl_all_lins_tmp$Lin_direction <- factor(rpl_all_lins_tmp$Lin_direction, levels = c("A_Down", "A_baseline", "A_Up",
                                                                                    "B_Down", "B_baseline", "B_Up",
                                                                                    "C_Down", "C_baseline", "C_Up"))


#----------------


# Create data object to produce figure with

Housekeeping_frame <- data.frame(rep("Housekeeping"), rpl_all_lins_tmp)

colnames(Housekeeping_frame)[[1]] <- "Category"

NHR_frame <- data.frame(rep("NHR"), nhr_all_lins_tmp)

colnames(NHR_frame)[[1]] <- "Category"

PGP_frame <- data.frame(rep("PGP"), pgp_all_lins_tmp)

colnames(PGP_frame)[[1]] <- "Category"

UGT_frame <- data.frame(rep("UGT"), ugt_all_lins_tmp)

colnames(UGT_frame)[[1]] <- "Category"

CYP_frame <- data.frame(rep("CYP"), nhr_all_lins_tmp)

colnames(CYP_frame)[[1]] <- "Category"


Xeno_HK_frame <- rbind(Housekeeping_frame, 
                       NHR_frame, 
                       PGP_frame, 
                       UGT_frame, 
                       CYP_frame)





Xeno_HK_frame$Category <- factor(Xeno_HK_frame$Category, levels = c("Housekeeping","PGP","NHR","UGT","CYP"))


# Supplementary Figure 9 A


Xeno_all_Lineage_compare <- ggplot(Xeno_HK_frame, aes(x=Category, y=Expression, fill = Category)) + 
  # geom_boxplot(fatten = 1, lwd = 1)+
  geom_violin()+
  #scale_color_manual(values=c("dodgerblue4", "turquoise4", "mediumaquamarine"))+
  ylim(-5, 5)+
  labs(y = " ")+
  geom_boxplot(width=0.05, fill = "white") + theme_minimal()+
  # geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center', dotsize=0.1)+
  theme_bw()+
  theme_classic()+
  theme(legend.position="none")+
  geom_hline(yintercept=c(-2.25, 2.25), col = "black")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  ggtitle(paste("Supplementary Figure 9 A"))

Xeno_all_Lineage_compare





pgp_trunc_names <- c("pgp-5", "pgp-6", "pgp-7", "pgp-8")

PGP_non_op_frame <- PGP_frame[!PGP_frame$gene %in% pgp_trunc_names, ]

Xeno_non_op_HK_frame <- rbind(Housekeeping_frame, 
                              NHR_frame, 
                              PGP_non_op_frame, 
                              UGT_frame, 
                              CYP_frame)


Xeno_non_op_HK_frame$Category <- factor(Xeno_non_op_HK_frame$Category, levels = c("Housekeeping","PGP","NHR","UGT","CYP"))


# Supplementary Figure 9 B


Xeno_non_op_all_Lineage_compare <- ggplot(Xeno_non_op_HK_frame, aes(x=Category, y=Expression, fill = Category)) + 
  # geom_boxplot(fatten = 1, lwd = 1)+
  geom_violin()+
  #scale_color_manual(values=c("dodgerblue4", "turquoise4", "mediumaquamarine"))+
  ylim(-5, 5)+
  labs(y = " ")+
  geom_boxplot(width=0.05, fill = "white") + theme_minimal()+
  # geom_dotplot(binaxis='y', binwidth = 0.25, stackdir='center', dotsize=0.1)+
  theme_bw()+
  theme_classic()+
  theme(legend.position="none")+
  geom_hline(yintercept=c(-2.25, 2.25), col = "black")+
  theme(axis.text.x = element_text(color = "black", size = 20, angle = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, face = "bold"),
        axis.title.y = element_text(color = "black", size = 16, angle = 90, face = "bold"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))+  
  ggtitle(paste("Supplementary Figure 9 B"))

Xeno_non_op_all_Lineage_compare






#----------------------------------------------------
# Calculations for data reported at end of paper
#----------------------------------------------------


# "Changes in 22G-RNAs may account for around 26 % of the heritable RNA seq changes."  


# "Although significantly greater than expected by chance, changes in chromatin accessibility explain only around 12 % of the remaining heritable changes."



# We can only calculate this from the set of genes with inherited RNAseq changes which have small RNA data and ATAC data available

RNA_inherit_mut_genes <- smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$RNA_mut==1 & smallRNA_All_lins_integrated_table$RNA_inherited==1 & smallRNA_All_lins_integrated_table$maps_to_smallRNA==1 & smallRNA_All_lins_integrated_table$maps_to_Ahr==1, 2]


# Considering all 3 lineages combined so don't correct for duplicates 


# Of which- proportion 'accounted' for by time matched small RNA changes (could also have ATAC?)

mut_genes_smallRNA <- smallRNA_All_lins_integrated_table[smallRNA_All_lins_integrated_table$RNA_mut==1 & smallRNA_All_lins_integrated_table$RNA_inherited==1 & smallRNA_All_lins_integrated_table$maps_to_smallRNA==1 & smallRNA_All_lins_integrated_table$maps_to_Ahr==1 & smallRNA_All_lins_integrated_table$smallRNA_mut==1 & smallRNA_All_lins_integrated_table$time_matched_to_RNA==1, 2]

library(vecsets)

small_RNA_accounts <- vintersect(mut_genes_smallRNA, RNA_inherit_mut_genes)

proportion_smallRNA <- length(small_RNA_accounts)/length(RNA_inherit_mut_genes)*100


# What is the remainder? i.e genes with inherited RNAseq changes with ATAC and small RNA data which don't have time matched small RNA epimutations



test_small_RNA_accounts <- small_RNA_accounts

test_RNA_inherit <- RNA_inherit_mut_genes


for(i in 1:length(test_small_RNA_accounts)){
  
  name <- test_small_RNA_accounts[[i]]
  
  if(name %in% test_RNA_inherit){
    
    positions <- which(test_RNA_inherit == name)
    
    test_RNA_inherit <- test_RNA_inherit[-positions[1]]
    
  }
  
}

Remainder <- test_RNA_inherit

# Of which- proportion of the remainder is 'accounted' for by time matched ATAC changes?

# To get this we have to use All_lins_integrated table

temp_1 <- All_lins_integrated_table[All_lins_integrated_table$gene %in% Remainder & All_lins_integrated_table$is_RNA_exp_change==1 & All_lins_integrated_table$is_RNA_inherited==1 & All_lins_integrated_table$has_time_matched==1, 2]


proportion_ATAC <- length(temp_1)/length(RNA_inherit_mut_genes)*100


