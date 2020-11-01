##### - Introduction ----
#Genes are an integral part of living organisms that contain information for creating proteins. These proteins then perform various tasks essential for functioning organisms. For the most part, genes encode instructions for well-being, however mutations disrupting a gene’s function can happen, resulting in genetic disorders. Ultimately, understanding how genes function is essential in targeting disorders for prevention and treatment. In this classification project, the recombination-activating genes (RAGs), RAG1 and RAG2 will be explored. These protein coding genes initiate V(D)J recombination, a process that generates unique antigen receptors for the recognition of various molecules (Alberts, 2002). This is an important part of the immune system and the development of T and B cells. Furthermore, it has been reported that a disruption in either gene can result in severe combined immunodeficiency (SCID) (Delmonte, 2018). The closely linked recombination-activating genes are involved in the same process; however, their distinguishing features are less known (Delmonte, 2018). In this project, the supervised machine learning algorithm, random forest is used to build a classifier using nucleotide proportions and k-mer frequencies as features. It is hypothesized that due to their similarity and relatedness in function, the RAG genes would also share similar sequence properties, making them difficult to classify. 

#In this classification project, the taxonomic group Cyprinidae, consisting of the species Danio rerio, commonly Zebrafish, is used. In recent years, Zebrafish are increasingly being used as a model organism for studying human diseases (Lieschke, 2007). Zebrafish share a high genetic similarity with humans making them an excellent model for studying gene function and the mechanisms contributing to the expression of those genes (Lieschke, 2007). The RAG genes are orthologs in humans and zebrafish, making zebrafish a potential model for studying RAG associated diseases. Furthermore, the ability to classify these genes, in the long term can help in studying their unique functions, determining the cause of mutations resulting in disorders, and opportunity for treatment. 

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
#install.packages("data.table")
library(data.table)

#### - Obtaining Gene Data from NCBI and Creating Dataframe ----

# We will be using the nuccore database to obtain sequence information for the genes RAG1 and RAG2 from the taxonomic group Cyprinidae. Cyprinidae is a family of fish with over 2000 species. The specie Danio rerio, commonly known as zebrafish, are a part of this taxonomic group. They are commonly used as model organisms for studying human diseases as 70% of human genes are also found in zebrafish, including the RAG genes (Howe, 2013). 
# In order to ensure sufficient data is obtained, our search is broadened to the family Cyprinidae, encompassing many species of freshwater fish (Cyprinid, n.d).

# We can get some summary information about the nuccore database to see when it was last updated and the name of the database
entrez_db_summary("nuccore")

# We can also look at the various search fields within the nuccore database
entrez_db_searchable("nuccore")

# A search of the Rag1 gene in UniProt, a database for protein sequences, tells us that the gene is just over 1500bp in some species. Similarly, UniProt tells us that Rag2 is just over 500bp. We can restrict our search to between 400 and 2000bp to ensure we are getting the correct data and not overloading R. 

# Let's get an overview of how to use entrez_search
?entrez_search

# We can search NCBI and the nuccore database for our gene of interest for our taxonomic group. We have included restrictions on the length of sequences to be between 400 to 2000bp.
Cyp_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag1[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])")

# We can check the total number of hits
Cyp_search$count

# The sample size required for building a machine learning model is largely dependent on how different the sequences are between the genes and how many features are used. We can start off by using a 150 RAG1 genes and 150 RAG2 genes so that the data is easily retrieved through R. If this poses to be a problem in our classification, we can consider increasing our sample size.
Cyp_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag1[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])", retmax=150)

Cyp_search

# We can check to see if the correct number of records are retrieved
length(Cyp_search$ids)

# Using the entrez_summary function, we can get summary information about the sequence records to ensure we have the correct gene and taxonomic group 
Rag1_summary <- entrez_summary(db = "nuccore", id = Cyp_search$ids)
Rag1_summary

# Now we can take a closer look at some information about the first record with the id 1910926681
Rag1_summary$`1910926681`$title
Rag1_summary$`1910926681`$organism
Rag1_summary$`1910926681`$taxid
Rag1_summary$`1910926681`$genome

# This data contains sequence information on the recombination activating protein 1 (Rag1) gene of the organism within the olive barb species. Here we can confirm that the gene and taxonomic group is correct.

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

# When looking at the Rag1 FASTA file in a text editor, I notice that it looks as it should with no excessively long sequences. Furthermore, it looks clean and free missing characters denoted by N.

# Now we can read it back in as a stringset 
Rag1_stringSet <- readDNAStringSet("Rag1_fetch.fasta")

# We can then create a data frame and organize our data frame more clearly by creating a column for our species name. We should now have a column with the title of the sequence, species name and sequence.
Rag1df <- data.frame(Rag1_Title = names(Rag1_stringSet), Rag1_Sequence = paste(Rag1_stringSet))
Rag1df$Species_Name <- word(Rag1df$Rag1_Title, 2L, 3L)
Rag1df <- Rag1df[, c("Rag1_Title", "Species_Name", "Rag1_Sequence")]

# Check the data frame again to make sure it looks as it should
view(Rag1df)

# We can also make another column specifying that this database contains the RAG1 genes. This will be helpful later when building our random forest model
Rag1df$Gene <- 'Rag1'

# Check to make sure it looks correct and the new column has been added
view(Rag1df)

# Let's get rid of some of the data that we don't need in our environment
rm(Cyp_search, Rag1_summary)

# We can retrieve our second gene of interest, RAG2, in the same way by searching and then fetching the data in FASTA format. We are keeping our sequence length and number of hits consistent with the gene RAG1.
Rag2_search <- entrez_search(db = "nuccore", term = "(Cyprinidae[ORGN] AND Rag2[Gene] AND 400:2000[SLEN]) NOT (genome[TITL])", retmax = 150, retmax = 100)
Rag2_fetch <- entrez_fetch(db = "nuccore", id = Rag2_search$ids, rettype = "fasta")

# Similarly, we can write the file to disk and take a look at it in a text editor
write(Rag2_fetch, "Rag2_fetch.fasta", sep = "\n") 

# The file looks correctly with the no excessively long sequences. The file is also clean without any missing characters denoted by N.

# Now we can read it back as a stringset
Rag2_stringSet <- readDNAStringSet("Rag2_fetch.fasta")

# Similarly, we can create a data frame and organize our data frame more clearly by creating a column for our species name. We should now have a column with the title of the sequence, species name and sequence.
Rag2df <- data.frame(Rag2_Title=names(Rag2_stringSet), Rag2_Sequence = paste(Rag2_stringSet))
Rag2df$Species_Name <- word(Rag2df$Rag2_Title, 2L, 3L)
Rag2df <- Rag2df[, c("Rag2_Title", "Species_Name", "Rag2_Sequence")]

# Let's look at the data frame.
View(Rag2df)

# We can add another column specifying the gene 
Rag2df$Gene <- "Rag2"

# Let's check to see if it has been done correctly
view(Rag2df)

# In order to merge our data frames, the data frame titles need to be the same. The column names are changed for both data frames so that they are the same and generalizable to both genes.
setnames(Rag1df, old = c('Rag1_Title','Rag1_Sequence'), new = c('Title','Sequence'))

setnames(Rag2df, old = c('Rag2_Title','Rag2_Sequence'), new = c('Title','Sequence'))

# Let's take a look at both data frames to make sure the names have been changed
view(Rag1df)
view(Rag2df)

# Now we can merge the two data frames so that we have one data frame called Rag_genes
Rag_genes <- rbind(Rag1df, Rag2df)

# Let's remove the data that we dont need in our environment
rm(Rag1df, Rag2df, Rag2_search)

#### - Preliminary Exploration ----

# Let's ensure we have the correct number of observations. We have 150 Rag1 observations and 150 Rag2 observations so we should have a total of 300 observations
dim(Rag_genes)

# The length of DNA sequences can influence our k-mer frequency. Our search was previously limited to 400 to 2000bp, but we can also create a histogram to compare the distribution of sequence lengths between the two genes.

# We are layering our distributions on the same graph to get a better comparison of the sequence length distributions.
# We can plot the distribution of the Rag1 sequence lengths first specifying the range of x and y values. We have 150 sequences for each gene with lengths ranging from 400 to 2000. 
hist(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag1"]), col=rgb(1,0,0,0.5), xlim=(c(400, 2000)), ylim=(c(0, 150)), xlab="Sequence Length", 
     ylab="Frequency", main="Distribution of Rag1 and Rag2 Sequence Lengths" )

# Now we can plot the distribution for Rag2 sequence lengths specifying that we want to add this plot on top of the previous
hist(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag2"]), col=rgb(0,0,1,0.5), add=T)

# Let's add a legend so we can distinguish between the two distributions
legend("topright", legend=c("Rag1","Rag2"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15, title = "Genes" )


# Here we notice that a majority of the RAG1 and Rag2 genes are  between 1000 to 1500bp. There is a bit of variability in the RAG1 genes with sequence lengths present between 500 to 1000bp, while none of the RAG2 genes are under 1000bp. This is something to keep in mind when we conduct our downstream analysis as the variability in sequence lengths of the RAG1 gene could influence our classification. Furthermore, there seems to be an outlier in the RAG2 gene sequence lengths. While outliers can be a result of improper data collection or storage, they can also be informative and indicate variability. For now we can keep this range and explore how if it affects our data.

# Let's look at some summary information regarding the min, max, median, first and third quartiles for the sequence lengths
summary(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag1"]))
summary(nchar(Rag_genes$Sequence[Rag_genes$Gene == "Rag2"]))

# We can consider creating a more stringent filter, however, a majority of the sequences among RAG1 and RAG2 fall between a comparable range. For now, we can keep this variability and observe how it affects our classification. We can also consider applying our filter of 1000 to 1500bp if we have trouble with our classification later on. 

#### - Calculating Sequence Features ----

# We can now generate our features will be used for our random forest model. The features we will be exploring are nucleotide frequencies and k-mer lengths of 2 and 3. We can also consider adding more complex features of k-mer lengths of 4 and 5 to increase our accuracy if needed.

# Let's convert our nucleotides to a DNAStringSet so that we can use the Biostrings package for calculating the sequences 
Rag_genes$Sequence <- DNAStringSet(Rag_genes$Sequence)

# We can calculate our letter frquency and append it to our dataframe
Rag_genes <- cbind(Rag_genes, as.data.frame(letterFrequency(Rag_genes$Sequence, letters = c("A", "C","G", "T"))))

# Let's view our dataframe to make sure it has been done correctly
view(Rag_genes)

# We can also calculate and add our A, T, G proportions to our data frame. Proportions can be helpful to use as they consider the variability in the length of the sequence rather than the raw count.
Rag_genes$Aprop <- (Rag_genes$A) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)

Rag_genes$Gprop <- (Rag_genes$G) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)

Rag_genes$Tprop <- (Rag_genes$T) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)

Rag_genes$Cprop <- (Rag_genes$C) / (Rag_genes$A + Rag_genes$T + Rag_genes$C + Rag_genes$G)

# We can look at our data frame again to make sure it has been done correctly
view(Rag_genes)

# Here we are adding our dinucleotide frequency proportions 
Rag_genes <- cbind(Rag_genes, as.data.frame(dinucleotideFrequency(Rag_genes$Sequence, as.prob = TRUE)))

# Let's also add our trinucleotide frequency proportions 
Rag_genes <- cbind(Rag_genes, as.data.frame(trinucleotideFrequency(Rag_genes$Sequence, as.prob = TRUE)))

# Let's look at our data frame again to confirm the dinucleotide and trinucleotide frequence have been added
view(Rag_genes)

# Before we start building our model, it will be interesting to see if any variables (nucleotide proportions) are correlated.
# Let's create a scatterplot to see the correlation between the frequency of A and G between the two genes of interest.
Rag_genes %>% 
  ggvis(~Aprop, ~Gprop, fill = ~Gene, shape = ~Gene) %>% 
  layer_points() %>%
  add_axis("x", title = "A Proportions") %>%
  add_axis("y", title = "G Proportions", properties=axis_props(title=list(dy=-30))) %>%
  add_axis("x", orient = "top", title = "Correlation Between A and G Proportions for RAG1 and RAG2 Genes")

# Interestingly, there does seem to be a high correlation between proportion of A and G of the RAG1 and RAG2 genes. The clustering of the two groups gives us an indication that there is a defining feature among the two genes. 

#### - Building the Classification Model ----

# We can convert our data back to character format to use tidyverse functions
Rag_genes$Sequence <- as.character(Rag_genes$Sequence)

# We are using randomization in the generation of the trees when we build our random forest model, so it is important to set our seed to ensure reproducibility. Using the same seed should give us the same result.
set.seed(213)

# We have a sample size of 150 for each gene so here we will set aside 35/150 (23%) for both genes to serve as the validation data set.
Rag_validation <- Rag_genes %>%
  group_by(Gene) %>%
  sample_n(35)

# We have our validation data set with 35 Rag1 gene sequences and 35 Rag2 gene sequences
table(Rag_validation$Gene)

# We can set our seed again for the for reproducibility
set.seed(23)

# Now we are creating a training dataset that excludes the samples used for the validation dataset. 
Rag_training <- Rag_genes %>%
  filter(!Title %in% Rag_validation$Title) %>%
  group_by(Gene) %>%
  sample_n(115)

# Let's make sure we have equal records of each gene
table(Rag_training$Gene)

# Let's take a look at variable names to see which variables to use for our classifier
names(Rag_training)

# We can start off by using our A, G, T and C proportions as predictors for classifying the genes
Rag_gene_classifier_simple <- randomForest(x = Rag_training[, 9:12], y = as.factor(Rag_training$Gene), ntree = 1000, importance = TRUE)

# Let's take a look at our results
Rag_gene_classifier_simple

# This worked pretty well! Only one of the RAG2 genes was misclassified. 

# Let's use our dinucleotide sequences as well to see how it affects our performance
Rag_gene_classifier2 <- randomForest(x = Rag_training[, 9:28], y = as.factor(Rag_training$Gene), ntree = 1000, importance = TRUE)

Rag_gene_classifier2
# Perfect performance

# We can take a look to see how the classifier for the dinucleotide sequences performs on unseen data 
predic_valid <- predict(Rag_gene_classifier2, Rag_validation[, c(4, 9:28)])

# Let's take a look at the class of our predic_valid object
predic_valid
class(predic_valid)
length(predic_valid)

# We can convert our Gene column in the validation dataset to a factor so that we can reference it when building our table for our validation results.
Rag_validation$Gene <- as.factor(Rag_validation$Gene)

table(observed = Rag_validation$Gene, predicted = predic_valid)
# In our validation table, we see that the performance of our model including nucleotide proportions and dinucleotide frequency as features had a 100% accuracy in classifying unseen data. From this we can conclude there are distinguishing features between the RAG gene sequences allowing us to accurately distinguish and classify them.

# Our classification model was able to accurately classify the genes using nucleotide proportions and dinucleotide frequencies. Let's look a little deeper into which features (variables) are the most important using the variable importance plots. 
varImpPlot(Rag_gene_classifier2, main = "Feature Importance for Rag Gene Classifier", bg = "red", color= "black")

# The VarImpPlot function gives us the MeanDecreaseAccuracy and MeanDecreaseGini plots with features on the y-axis and measures on the x-axis. On the left we have our MeanDecreaseAccuracy plot which tells us the decrease in the accuracy of our model is most affected by the removal of the feature GA. Similarly, the MeanDecreaseGini is another measure of variable importance for estimating a target variable.A higher MeanDecreaseGini indicates higher feature importance. Both plots indicate that GA and C proportion are the most important features in classification of the RAG genes. 

# Let's look at just the GA and Cproportion features
Rag_gene_classifierGA <- randomForest(x = Rag_training[, 21], y = as.factor(Rag_training$Gene), ntree = 1000, importance = TRUE)

Rag_gene_classifierGA
# This model with only GA as a feature had 0% error rate!

Rag_gene_classifierC <- randomForest(x = Rag_training[, 12], y = as.factor(Rag_training$Gene), ntree = 1000, importance = TRUE)

Rag_gene_classifierC
# Similarly, the model with only C proportions had an error rate of 0.43%

##### - Conclusion ----

#The immune system is vital to any organism’s survival. Its main tasks are to fight pathogens and disease-causing substance in the body and recognize and counteract harmful substances outside of the body (NCBI, 2020). Genes including the recombination-activating genes, play a role in the function of a healthy immune system, however mutations in these genes can also arise, resulting in immune disorders. The RAG genes are known to be responsible for the development of T and B cells, however their unique function is still being studied. It was hypothesized that their similarity in function would result in similar sequence features. In practice, it was observed that the RAG1 and RAG2 genes are quite dissimilar in their sequences. A variety of features were explored, including nucleotide proportions, and k-mer proportions of 2 and 3. Additionally, it was revealed through the variable importance plot that the most important feature for classifying the gene sequences was the 2-mer, GA. A model was then built using this feature alone, resulting in 100% accuracy. Furthermore, the discovery of this feature can bring insight into its functional significance.

#Currently, little is known about the distinguishing functions of the RAG1 and RAG2 genes in the immune system. The ability to distinguish between related genes is an important step for ultimately being able to study their individual functions. The ability to look at important sequence features for genes can bring insight into their consequent biological function. In future studies, these distinguishing features can be considered further to analyse their potential role in disorders and the research for treatment.

##### - Acknowledgements ----

#I spoke to Abi about our projects and potential plots that we could include to add to the exploration of our classification. We spoke about the use of ROC curves and feature importance plots. I explored the use of ROC curves, however decided to opt for a feature importance plot in my assignment instead. 
##### - References ----

#UniProt. (2018): a worldwide hub of protein knowledge. Nucleic Acids Research, 47(D1), D506–D515. https://doi.org/10.1093/nar/gky1049

#Alberts, B., Lewis A. J., et al. (2002). The Generation of Antibody Diversity. Molecular Biology of the Cell 4th edition. Retreived from https://www.ncbi.nlm.nih.gov/books/NBK26860/
  
#Delmonte, O. M., Schuetz, C., & Notarangelo, L. D. (2018). RAG Deficiency: Two Genes, Many Diseases. Journal of Clinical Immunology, 38(6), 646–655. https://doi.org/10.1007/s10875-018-0537-4

#Howe, K., Clark, M. D., Torroja, C. F., Torrance, J., et al. (2013). The zebrafish reference genome sequence and its relationship to the human genome. Nature, 496(7446), 498–503. https://doi.org/10.1038/nature12111

#Lieschke, G. J., & Currie, P. D. (2007). Animal models of human disease: zebrafish swim into view. Nature Reviews Genetics, 8(5), 353–367. https://doi.org/10.1038/nrg2091

#NCBI. (2020, April 23). How does the immune system work? Retrieved from https://www.ncbi.nlm.nih.gov/books/NBK279364/
  
#New World Encyclopedia. (n.d.). Cyprinid. Retrieved from https://www.newworldencyclopedia.org/entry/Cyprinid

#Willett, C. E., Cherry, J. J., & Steiner, L. A. (1997). Characterization and expression of the recombination activating genes (rag1 and rag2) of zebrafish. Immunogenetics, 45(6), 394–404. https://doi.org/10.1007/s002510050221



