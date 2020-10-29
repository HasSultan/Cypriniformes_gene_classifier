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

#### - Obtaining Gene Data from NCBI and Creating Dataframe ----
#We will be using the nuccore database to obtain sequence information for the genes Rag1 and Rag2 from the taxanomic group Cyprinidae. Cyprinidae is a family of fish with over 2000 species.

#We can get some summary information about the nuccore database to see when it was last updated and the name of the database
entrez_db_summary("nuccore")

#We can also look at the various search fields within the nuccore database
entrez_db_searchable("nuccore")

#A search of the Rag1 gene in UniProt, a database for protein sequences, tells us that the gene is just over 1500bp in some species. Similarly, UniProt tells us that Rag2 is just over 500pb. We can restrict our search to between 400 and 2000bp to ensure we are getting the correct data and not overloading R. 

#Lets get an overview of how to use entrez_search
?entrez_search

#We can search NCBI and the nuccore database for our gene of interest for our taxanomic group. We have included restrictions on the length of sequences to be between 400 to 2000bp
Cyp_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag1[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])")

#We can check the total number of hits
Cyp_search$count

#The maximum number of hits is restricted to 150 to make it easier to work with and easily retrieved through R.
Cyp_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag1[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])", retmax=150)

#We can check to see the correct number of records are retrieved
length(Cyp_search$ids)

#Using the entrez_summary function, we can get summary information about the sequence records to ensure we have the correct gene and taxanomic group 
Rag1_summary <- entrez_summary(db = "nuccore", id = Cyp_search$ids)
Rag1_summary

#Now we can take a closer look at some information about the first record with the id 1910926681
Rag1_summary$`1910926681`$title
Rag1_summary$`1910926681`$organism
Rag1_summary$`1910926681`$taxid
Rag1_summary$`1910926681`$genome

#This data contains sequence information on the recombination activating protein 1 (Rag1) gene of the organism within the olive barb species 

#Now we can obtain sequence information using entrez_fetch
?entrez_fetch

#Here we are specifying that we want to obtain the Rag1 gene information from Cyp_search in FASTA format 
Rag1_fetch <- entrez_fetch(db = "nuccore", id = Cyp_search$ids, rettype = "fasta")

#Check what class it is
class(Rag1_fetch)

#Preview of the dat
head(Rag1_fetch)

#Let's write this file to disk so we can look at it in a text editor
write(Rag1_fetch, "Rag1_fetch.fasta", sep = "\n")

#When looking at the Rag1 FASTA file in a text editor, I notice that it looks as it should with no exessively long sequences. Furthermore, it looks very clean and free of Ns.

#Now we can read it back in as a stringset 
Rag1_stringSet <- readDNAStringSet("Rag1_fetch.fasta")

##FIXLets create a dataframe and clean up the titles
Rag1df <- data.frame(Rag1_Title = names(Rag1_stringSet), Rag1_Sequence = paste(Rag1_stringSet))
Rag1df$Species_Name <- word(Rag1df$Rag1_Title, 2L, 3L)
Rag1df <- Rag1df[, c("Rag1_Title", "Species_Name", "Rag1_Sequence")]

#Check the dataframe again to make sure it looks as it should
view(Rag1df)

#We can also make another column specifying that this database contains the Rag1 genes. This will be helpful later when building our random forest model
Rag1df$Gene <- 'Rag1'

#Check to make sure it looks correct and the new column has been added
view(Rag1df)

#We can retrieve our second gene of interest, Rag2, in the same way by searching and then fetching the data in FASTA format
Rag2_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag2[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])", retmax = 150, retmax = 100)
Rag2_fetch <- entrez_fetch(db = "nuccore", id = Rag2_search$ids, rettype = "fasta")

#Similarly, we can write the file to disk and take a look at it in a text editor
write(Rag2_fetch, "Rag2_fetch.fasta", sep = "\n") 
#The file looks correctly with the length and is clean without any Ns

#Now we can read it back as a stringset
Rag2stringSet <- readDNAStringSet("Rag2_fetch.fasta")


Rag2df <- data.frame(Rag2_Title=names(Rag2stringSet), Rag2_Sequence = paste(Rag2stringSet))
Rag2df$Species_Name <- word(Rag2df$Rag2_Title, 2L, 3L)
Rag2df <- Rag2df[, c("Rag2_Title", "Species_Name", "Rag2_Sequence")]

# Let's look at the dataframe.
View(Rag2df)

#Similarly, we can add another column specifying the gene 
Rag2df$Gene <- "Rag2"

#Let's check to see if it has been done correctly
view(Rag2df)

#In order to merge our dataframe, the dataframe titles need to be the same. The column names are being changed so that they are the same for Rag1df_name and Rag2_name
Rag1df_name <- Rag1df %>%
  rename(Title = Rag1_Title, Sequence = Rag1_Sequence)

Rag2df_name <- Rag2df %>%
  rename(Title = Rag2_Title, Sequence = Rag2_Sequence)

#Let's take a look at both dataframes to make sure the names have been changed
view(Rag1df_name)
view(Rag2df_name)

#Now we can merge the two dataframes so that we have one dataframe called Rag_genes
Rag_genes <- rbind(Rag1df_name, Rag2df_name)

#### - Preliminary Exploration ----
### to use histograms or not?
hist(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag1"]), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of Rag1 Sequence Lengths")
summary(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag1"]))

hist(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag2"]), xlab = "Sequence Length", ylab = "
     Frequency", main = "Frequency Histogram of Rag2 Sequence Lengths")
summary(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag2"]))

view(Rag_genes)
library(ggvis)
Rag_genes %>% ggvis(~A, ~C, fill = ~Gene) %>% layer_points()

library(ggvis)
Rag_genes %>% ggvis(~G, ~T, fill = ~Gene) %>% layer_points()


dim(Rag_genes)
sum(is.na(Rag_genes$Rag1_Sequence))
sum(is.na(Rag_genes$Rag2_Sequence))
summary(str_count(Rag_genes$Rag1_Sequence))
summary(str_count(Rag_genes$Rag2_Sequence))
         

#### - Calculating Sequence Features ----
Rag_genes <- as.data.frame(Rag_genes)
Rag_genes$Sequence <- DNAStringSet(Rag_genes$Sequence)

BrowseSeqs(Rag_genes$Sequence[Rag_genes$Gene == "Rag1"])

BrowseSeqs(Rag_genes$Sequence[Rag_genes$Gene == "Rag2"])

summary(nchar(dfBOLD$nucleotides[dfBOLD$markercode == "28S"]))
Rag_genes <- cbind(Rag_genes, as.data.frame(letterFrequency(Rag_genes$Sequence, letters = c("A", "C","G", "T"))))
view(Rag_genes)

#
Rag_genes$Aprop <- (Rag_genes$A) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)

Rag_genes$Tprop <- (Rag_genes$T) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)

Rag_genes$Gprop <- (Rag_genes$G) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)

view(Rag_genes)

Rag_genes <- cbind(Rag_genes, as.data.frame(dinucleotideFrequency(Rag_genes$Sequence, as.prob = TRUE)))
view(Rag_genes)

Rag_genes <- cbind(Rag_genes, as.data.frame(trinucleotideFrequency(Rag_genes$Sequence, as.prob = TRUE)))
view(Rag_genes)

Rag_genes$Sequence <- as.character(Rag_genes$Sequence)

table(Rag_genes$Gene)

#### - Building Classification Model ----
#As our maximum sample size for these two genes is 82 for the 28S gene in dfCotesia, we will below sample 20 individuals (about 25% of the total for the 28S gene) from each gene to serve as the validation data set. These samples will be entirely set aside during the model training process. These samples should NOT be used until validating the model later. sample_n() is a very useful function for sampling a specific number of rows from each group. Here, we are sampling 20 rows from each group (i.e. each marker code). See the resulting object to confirm.

#setting seed, as we are using randomization in next step, and we want results to be reproducible.

set.seed(213)

Rag_validation <- Rag_genes %>%
  group_by(Gene) %>%
  sample_n(35)

#Checking object. make prediction about what we should see:
table(Rag_validation$Gene)

#Now, we are creating a training dataset that does NOT overlap with the validation set. To do this, we will first remove processids that were among the samples randomly selected for the validation dataset. We can do this by asking for processid's that are not (!) in the validation set using %in%. Second, we will pick 62 individuals with each gene to serve in the training set. Note that this is a small sample size; this script is an EXAMPLE. The sample sizes don't necessarily need to be exactly equal, but in general we need to aim for "class balance" (i.e. similar representation between classes, or groups). If there is a large imbalance in sample size, as would be the case here if we used all COI sequences, the classifier could achieve good performance simply by saying that everything is COI! That is NOT what we want... we want to train the classifier to recognize COI vs. 28S sequences.
set.seed(23)
Rag_training <- Rag_genes %>%
  filter(!Title %in% Rag_validation$Title) %>%
  group_by(Gene) %>%
  sample_n(115)

#Checking we should have 62 records of each gene
table(Rag_training$Gene)

#Checking our variable names and the column numbers.
names(Rag_training)

#Next, we are building a classifier to separate the COI and 28S genes in these datasets, using first the A, T, and G proportion (columns 86-88) as predictors. Then, if needed, we can see if it is helpful to add more complex features (e.g. dinucleotide frequencies, 3-mers, etc.). The response variable is markercode; we are trying to predict which gene a sequence belongs to, on the basis of simple sequence features alone.
Rag_gene_classifier <- randomForest(x = Rag_training[, 9:11], y = as.factor(Rag_training$Gene), ntree = 1000, importance = TRUE)

#Let's look at results.
Rag_gene_classifier

plot(Rag_gene_classifier)

Rag_gene_classifier <- randomForest(x = Rag_training[, 9:27], y = as.factor(Rag_training$Gene), ntree = 1000, importance = TRUE)
Rag_gene_classifier







#### - KNN ----
##Generate a random number that is 90% of the total number of rows in dataset.
rand_num <- sample(1:nrow(Rag_genes), 0.9 * nrow(Rag_genes)) 

##the normalization function is created
normalization_func <-function(x) { (x -min(x))/(max(x)-min(x))   }

##Run nomalization on first 4 coulumns of dataset because they are the predictors
Rag_norm <- as.data.frame(lapply(Rag_genes[,c(9,10,11)], normalization_func))


##extract training set
R_train <- Rag_norm[ran,] 
##extract testing set
R_test <- Rag_norm[-ran,] 
##extract 5th column of train dataset because it will be used as 'cl' argument in knn function.
R_target_category <- iris[ran,5]
##extract 5th column if test dataset to measure the accuracy
R_test_category <- iris[-ran,5]
##load the package class
library(class)
##run knn function
pr <- knn(iris_train,iris_test,cl=iris_target_category,k=13)

##create confusion matrix
tab <- table(pr,iris_test_category)

##this function divides the correct predictions by total number of predictions that tell us how accurate teh model is.

accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
accuracy(tab)
## [1] 80


##### - Conclusion ----