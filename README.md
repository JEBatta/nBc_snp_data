# nBc_snp_data
The code use a database of 347 single nucleotide polymorphisms (SNPs) expression (number of mutant alleles) in Mexican adults to perform a binomial test (to score the contribution of a single predictor on class membership) and a naive Bayes classification of individuals on the following classes: Obesity, Abdominal Obesity and High TGB (1:have the condition, 0:does not have the condition). The programm consists in a 10x10 fold validation of the same data to show predictability of models with AUCs in ROC analysis.

We include a list of variants (SNPs) on metabolic pathways to evaluate groups of SNPs as class predictors. These are included in the database, with the code employed to score pathway lists on individuals. 

In this repository you can find only the codes. An database example is included to verify the format of the data.
