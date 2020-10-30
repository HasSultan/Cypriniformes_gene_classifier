##### - Introduction ----

#### - Loading Relevant Packages ----
#install.packages("tidyverse")
library(tidyverse)
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
#install.packages("ggvis")
library(ggvis)
#install.packages("DECIPHER")
library(DECIPHER)
#install.packages("shiny")
library(shiny)

#### - Obtaining Gene Data from NCBI and Creating Dataframe ----
# We will be using the nuccore database to obtain sequence information for the genes Rag1 and Rag2 from the taxanomic group Cyprinidae. Cyprinidae is a family of fish with over 2000 species.

# We can get some summary information about the nuccore database to see when it was last updated and the name of the database
entrez_db_summary("nuccore")

# We can also look at the various search fields within the nuccore database
entrez_db_searchable("nuccore")

# A search of the Rag1 gene in UniProt, a database for protein sequences, tells us that the gene is just over 1500bp in some species. Similarly, UniProt tells us that Rag2 is just over 500pb. We can restrict our search to between 400 and 2000bp to ensure we are getting the correct data and not overloading R. 

# Lets get an overview of how to use entrez_search
?entrez_search

# We can search NCBI and the nuccore database for our gene of interest for our taxanomic group. We have included restrictions on the length of sequences to be between 400 to 2000bp
Cyp_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag1[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])")

# We can check the total number of hits
Cyp_search$count

# The maximum number of hits is restricted to 150 to make it easier to work with and easily retrieved through R.
Cyp_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag1[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])", retmax=150)

Cyp_search

# We can check to see the correct number of records are retrieved
length(Cyp_search$ids)

# Using the entrez_summary function, we can get summary information about the sequence records to ensure we have the correct gene and taxanomic group 
Rag1_summary <- entrez_summary(db = "nuccore", id = Cyp_search$ids)
Rag1_summary

# Now we can take a closer look at some information about the first record with the id 1910926681
Rag1_summary$`1910926681`$title
Rag1_summary$`1910926681`$organism
Rag1_summary$`1910926681`$taxid
Rag1_summary$`1910926681`$genome

# This data contains sequence information on the recombination activating protein 1 (Rag1) gene of the organism within the olive barb species 

# Now we can obtain sequence information using entrez_fetch
?entrez_fetch

# Here we are specifying that we want to obtain the Rag1 gene information from Cyp_search in FASTA format 
Rag1_fetch <- entrez_fetch(db = "nuccore", id = Cyp_search$ids, rettype = "fasta")

# Check what class it is
class(Rag1_fetch)

# Preview of the data
head(Rag1_fetch)

# Let's write this file to disk so we can look at it in a text editor
write(Rag1_fetch, "Rag1_fetch.fasta", sep = "\n")

# When looking at the Rag1 FASTA file in a text editor, I notice that it looks as it should with no exessively long sequences. Furthermore, it looks very clean and free of Ns.

# Now we can read it back in as a stringset 
Rag1_stringSet <- readDNAStringSet("Rag1_fetch.fasta")

##FIXLets create a dataframe and clean up the titles
Rag1df <- data.frame(Rag1_Title = names(Rag1_stringSet), Rag1_Sequence = paste(Rag1_stringSet))
Rag1df$Species_Name <- word(Rag1df$Rag1_Title, 2L, 3L)
Rag1df <- Rag1df[, c("Rag1_Title", "Species_Name", "Rag1_Sequence")]

# Check the dataframe again to make sure it looks as it should
view(Rag1df)

# We can also make another column specifying that this database contains the Rag1 genes. This will be helpful later when building our random forest model
Rag1df$Gene <- 'Rag1'

# Check to make sure it looks correct and the new column has been added
view(Rag1df)

# Lets get rid of some of the data that we dont need in our environment
rm(Cyp_search, Rag1_summary)

# We can retrieve our second gene of interest, Rag2, in the same way by searching and then fetching the data in FASTA format. We are keeping our sequence length and number of hits consistent with the first gene, Rag1.
Rag2_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag2[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])", retmax = 150, retmax = 100)
Rag2_fetch <- entrez_fetch(db = "nuccore", id = Rag2_search$ids, rettype = "fasta")

# Similarly, we can write the file to disk and take a look at it in a text editor
write(Rag2_fetch, "Rag2_fetch.fasta", sep = "\n") 

# The file looks correctly with the no exessively long sequences. The file is also clean without any Ns

# Now we can read it back as a stringset
Rag2_stringSet <- readDNAStringSet("Rag2_fetch.fasta")

###FIXX
Rag2df <- data.frame(Rag2_Title=names(Rag2_stringSet), Rag2_Sequence = paste(Rag2_stringSet))
Rag2df$Species_Name <- word(Rag2df$Rag2_Title, 2L, 3L)
Rag2df <- Rag2df[, c("Rag2_Title", "Species_Name", "Rag2_Sequence")]

# Let's look at the dataframe.
View(Rag2df)

# Similarly, we can add another column specifying the gene 
Rag2df$Gene <- "Rag2"

# Let's check to see if it has been done correctly
view(Rag2df)

# In order to merge our dataframe, the dataframe titles need to be the same. The column names are being changed so that they are the same for Rag1df_name and Rag2_name
Rag1df_name <- Rag1df %>%
  rename(Title = Rag1_Title, Sequence = Rag1_Sequence)

Rag2df_name <- Rag2df %>%
  rename(Title = Rag2_Title, Sequence = Rag2_Sequence)

# Let's take a look at both dataframes to make sure the names have been changed
view(Rag1df_name)
view(Rag2df_name)

# Now we can merge the two dataframes so that we have one dataframe called Rag_genes
Rag_genes <- rbind(Rag1df_name, Rag2df_name)

# Lets remove the data that we dont need in our environment
rm(Rag1df, Rag2df, Rag1df_name, Rag2df_name, Rag2_search)

#### - Preliminary Exploration ----

# Lets ensure we have the correct number of observations. We have 150 Rag1 obervations and 150 Rag2 observations so we should have a total of 300 observations
dim(Rag_genes)

### We are building a classification model ur sequence features may be affected by sequence lengthWe can create histograms to look at the sequence length distribution for Rag1 and Rag2. This is important because we are looking sequence features when we build our classification model. This could be affected by sequence lengths 

# We are layering our distributions on the same graph to get a better comparison of the sequence length distributions
# We can plot the distribution of the Rag1 sequence lengths first specifying the range of x and y values. We have a 150 sequences for each gene with lengths ranging from 400 to 2000. 
hist(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag1"]), col=rgb(1,0,0,0.5), xlim=(c(400, 2000)), ylim=(c(0, 150)), xlab="Sequence Length", 
     ylab="Frequency", main="Distribution of Rag1 and Rag2 Sequence Lengths" )

# Now we can plot the distribution for Rag2 sequence lengths specifying that we want to add this plot on top of the previous
hist(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag2"]), col=rgb(0,0,1,0.5), add=T)

# Lets add a legend so we can distinguish between the two distributions
legend("topright", legend=c("Rag1","Rag2"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15, title = "Genes" )

###talk about plot

# Lets look at some summary information regarding the min, max, median, first and third quartiles for the sequence lengths
summary(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag1"]))
summary(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag2"]))

##talk about results         

#### - Calculating Sequence Features ----

# Lets convert our nucleotides to a DNAStringSet so that we can use the Biostrings package for calculating the sequences 
Rag_genes$Sequence <- DNAStringSet(Rag_genes$Sequence)

# We can calculate our letter frquency and append it to our dataframe
Rag_genes <- cbind(Rag_genes, as.data.frame(letterFrequency(Rag_genes$Sequence, letters = c("A", "C","G", "T"))))

# Lets view our dataframe to make sure it has been done correctly
view(Rag_genes)

# We can also calculate and add our A, T, G proportions to our dataframe
Rag_genes$Aprop <- (Rag_genes$A) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)

Rag_genes$Gprop <- (Rag_genes$G) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)

Rag_genes$Tprop <- (Rag_genes$T) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)


# We can look at our dataframe again to make sure it has been done correctly
view(Rag_genes)

# Here we are adding our dinucleotide frequency proportions 
Rag_genes <- cbind(Rag_genes, as.data.frame(dinucleotideFrequency(Rag_genes$Sequence, as.prob = TRUE)))

# Lets also add our trinucleotide frequency proportions 
Rag_genes <- cbind(Rag_genes, as.data.frame(trinucleotideFrequency(Rag_genes$Sequence, as.prob = TRUE)))

# Lets look at our dataframe again to confirm the dinucleotide and trinucleotide frequence have been added
view(Rag_genes)

# Before we start building our model, it will be interesting to see if any variables (Nucleotides) are correlated
# Lets create a scatterplot to see the correlation between the frequency of A and G between the two genes of interest
Rag_genes %>% 
  ggvis(~A, ~G, fill = ~Gene, shape = ~Gene, size := input_slider(10, 100)) %>% 
  layer_points() %>%
  add_axis("x", title = "Frequency of A") %>%
  add_axis("y", title = "Frequency of G") %>%
  add_axis("x", orient = "top", title = "Correlation Between Frequency of A and G for Rag1 and Rag2 Genes")

#### - Building Classification Model ----

# We can convert our data back to character format to use tidyverse functions
Rag_genes$Sequence <- as.character(Rag_genes$Sequence)

# We are using randomization in the generation of the trees when we build our random forest model so it is important to set our seed to ensure reproducibility. Using the same seed should give us the same result.
set.seed(213)

# We have a sample size of 150 for each gene so here we will set aside 35/150 (23%) for both genes to serve as the validation data set
Rag_validation <- Rag_genes %>%
  group_by(Gene) %>%
  sample_n(35)

# We have our validation data set with 35 Rag1 gene sequences and 35 Rag2 gene sequences
table(Rag_validation$Gene)

#
set.seed(23)

# Now we are creating a training dataset that excludes the samples used for the validation dataset. 
Rag_training <- Rag_genes %>%
  filter(!Title %in% Rag_validation$Title) %>%
  group_by(Gene) %>%
  sample_n(115)

# Lets make sure we have equal records of each gene
table(Rag_training$Gene)

# Lets take a look at variable names to see which variables to use for our classifier
names(Rag_training)

# We can start off by using our A, G and T proportions as predictors for classifying the genes
Rag_gene_classifier_simple <- randomForest(x = Rag_training[, 9:11], y = as.factor(Rag_training$Gene), ntree = 1000, importance = TRUE)

# Lets take a look at our results
Rag_gene_classifier_simple

# This worked pretty well! Only one of the Rag2 genes was misclassified. 

#We can also dive into the specifics, such as by looking at the following. See the documentation for the randomForest() function for more information.
Rag_gene_classifier$importance
Rag_gene_classifier$oob.times

# Lets use our dinucleotide sequences as well to see how it affects our performance
Rag_gene_classifier2 <- randomForest(x = Rag_training[, 9:27], y = as.factor(Rag_training$Gene), ntree = 1000, importance = TRUE)

Rag_gene_classifier2
# Perfect performance

# Lets look at just A and G
Rag_gene_classifierAG <- randomForest(x = Rag_training[, 9:10], y = as.factor(Rag_training$Gene), ntree = 1000, importance = TRUE)

Rag_gene_classifierAG

# We can take a look to see how the classifier for the dinucleotide sequences performs on unseen data 
predic_valid <- predict(Rag_gene_classifier2, Rag_validation[, c(4, 9:27)])

#Let's have a look at the data format here. There are the predictions.
predic_valid
class(predic_valid)
length(predic_valid)

#Create our own confusion matrix.
Rag_validation$Gene <- as.factor(Rag_validation$Gene)

table(observed = Rag_validation$Gene, predicted = predic_valid)

class(predic_valid)
class()

#Yes, this classifier also works on the unseen data.
varImpPlot(Rag_gene_classifier2,  main = "Variable Importance for Rag Gene Classifier", bg = "red", color= "black")

#### - KNN ----


##### - Conclusion ----