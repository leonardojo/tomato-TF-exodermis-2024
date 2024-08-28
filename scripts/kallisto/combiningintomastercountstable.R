library(tidyverse)

## each count table will be stored in a individual folder (i.e. myb2KO-1)
## output of kallisto is an pseudocount table named abundance.tsv
## this R script uses a for loop to open all the abundance.tsv files and create a master counts table
## the final output is a table in which the first column in the feature(gene name) and adjacent columns are each


## Listing the directory names
folder.names <- list.dirs(full.names = F, recursive = F,path="KALLISTO_COUNTS")

## Creating a final.df object to store all the count files
final.df <- NULL
k <- folder.names[1]
final.df <- read.table(paste0("KALLISTO_COUNTS/",k,"/abundance.tsv"),header=T) %>%
  select(target_id,est_counts)
names(final.df)[2] <- k
for (x in folder.names[-1]) {
  df <- read.table(paste0("KALLISTO_COUNTS/",x,"/abundance.tsv"),header=T) %>%
    select(target_id,est_counts)
  names(df)[2] <- x
  final.df <- full_join(final.df,df,by="target_id")
}

## saving the master table as a csv
## Removing the Undetermined sample
## Arranging the genes alphabetically
final.df %>%
  arrange(target_id) %>%
  select(-Undetermined) %>%
  write.csv("mastercounts.csv",row.names=F)

## saving the master table as a tsv
## Removing the Undetermined sample
## Arranging the genes alphabetically
final.df %>%
  arrange(target_id) %>%
  select(-Undetermined) %>%
  write_tsv("mastercounts.tsv")
