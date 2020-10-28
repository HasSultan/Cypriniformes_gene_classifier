#Correct Project connected to git
#### library load/ include section for downloads (for sally etc)----
#install.packages("rentrez")
library(rentrez)
#install.packages("seqinr")
library(seqinr)
#BiocManager::install("Biostrings")
library(Biostrings)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("randomForest")
library(randomForest)

#### get data ----
entrez_dbs()
entrez_db_summary("nuccore")
entrez_db_searchable("nuccore")

?entrez_search
Cyp_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag1[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])", retmax=150)

maxHits <- Cyp_search$count
maxHits

Cyp_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag1[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])", retmax = 150)

Cyp_search
Cyp_search$ids
length(Cyp_search$ids)

Rag1_summary <- entrez_summary(db = "nuccore", id = Cyp_search$ids)
Rag1_summary

class(Rag1_summary)
class(Rag1_summary[[1]])

Rag1_summary$`1910926681`$title
Rag1_summary$`1910926681`$organism
Rag1_summary$`1910926681`$taxid
Rag1_summary$`1910926681`$genome

extract_from_esummary(Rag1_summary, "organism")
?entrez_fetch

# Let's specify that we want the return type ("rettype") in FASTA format.
Rag1_fetch <- entrez_fetch(db = "nuccore", id = Cyp_search$ids, rettype = "fasta")

# What class is it?
class(Rag1_fetch)

# Its an extremely long character vector.
head(Rag1_fetch)

write(Rag1_fetch, "Rag1_fetch.fasta", sep = "\n")
#check in text editor and looks correct... talk about more...

# Read it back in as DNA StringSet using the readDNAStringSet() function.
Rag1_stringSet <- readDNAStringSet("Rag1_fetch.fasta")

Rag1df <- data.frame(Rag1_Title = names(Rag1_stringSet), Rag1_Sequence = paste(Rag1_stringSet))

Rag1df$Species_Name <- word(Rag1df$Rag1_Title, 2L, 3L)
Rag1df <- Rag1df[, c("Rag1_Title", "Species_Name", "Rag1_Sequence")]
view(Rag1df)


##### Rag1 ----
Rag2_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag2[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])", retmax = 150, retmax = 100)
Rag2_fetch <- entrez_fetch(db = "nuccore", id = Rag2_search$ids, rettype = "fasta")
write(Rag2_fetch, "Rag2_fetch.fasta", sep = "\n") 
head(Rag2_fetch)


Rag2stringSet <- readDNAStringSet("Rag2_fetch.fasta")

Rag2df <- data.frame(Rag2_Title=names(Rag2stringSet), Rag2_Sequence = paste(Rag2stringSet))
# Clean up species names.

Rag2df$Species_Name <- word(Rag2df$Rag2_Title, 2L, 3L)
# Rearrange the columns.

Rag2df <- Rag2df[, c("Rag2_Title", "Species_Name", "Rag2_Sequence")]
# Let's look at the dataframe.
View(Rag2df)


#how many unique species
length(unique(Rag1df$Species_Name))

#df with only unique species
Rag1df_Subset <- Rag1df %>% 
  group_by(Species_Name) %>%  ## Group by Name.
  sample_n(1) 

#check if correct
all.equal(length(unique(Rag1df$Species_Name)), nrow(Rag1df_Subset))
#TRUE

#how many unique species
length(unique(Rag2df$Species_Name))

#df with only unique species
Rag2df_Subset <- Rag2df %>% 
  group_by(Species_Name) %>%  ## Group by Name.
  sample_n(1) 

#check if correct
all.equal(length(unique(Rag2df$Species_Name)), nrow(Rag2df_Subset))
#TRUE

Ragoverlap_df <- merge(Rag1df_Subset, Rag2df_Subset, by = "Species_Name", all = F)

hist(nchar(Ragoverlap_df$Rag1_Sequence), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of Rag1 Sequence Lengths")
summary(nchar(Ragoverlap_df$Rag1_Sequence))
