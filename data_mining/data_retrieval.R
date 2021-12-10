source("required_functions.R")

# Install packages if not already installed
install.packages("easyPubMed")
install.packages("devtools")
devtools::install_github("dami82/easyPubMed")
devtools::install_github("dami82/businessPubMed")

#load library packages R
library(devtools)
library(easyPubMed)
library(rentrez)
library(kableExtra)
library(dplyr)
library(businessPubMed)
library(githubinstall)

# Design the query for NCBI
my.query <- paste('Kenya[Affiliation]')
my.query <- paste(my.query, 'AND ("1980"[EDAT]:"2019"[EDAT])')

# view the query
cat(my.query)

# define the regex filter for countries other than Kenya 
data("countries")
countries <- countries[countries != "Kenya"]
countries.filter <- gsub("[[:punct:]]", "[[:punct:]]", toupper(countries))
countries.filter <- gsub(" ", "[[:space:]]", countries.filter)
countries.filter <- paste("(",countries.filter,")", sep = "", collapse = "|")
countries.filter <- paste("(UK)|", countries.filter, sep = "")
cat(countries.filter)

# Show filter string
cat(paste(substr(countries.filter, 1, 85), "...", sep = ""))


#intial time 
t_init <- Sys.time()

#extract data from PubMed
my.data3 <- extract_pubMed_data(pubMed_query = my.query, 
                                batch_size = 1000,
                                affi_regex_exclude = countries.filter)
# Check time t_final
t_final <- Sys.time()

#few data set acquired 
head(my.data3$data, 10)

# exporting data from the R console to the desktop
write.table(my.data3$data, "Kenyapapers.csv", row.names=F, sep ="\t")

#Bash script for extracting pmid
x <- "bash ../script/extract_pmid.sh"

#importing the pmid id in R
Kenya_pmid.csv <- read.csv("Kenyapmid.csv")

# number of id available
count(Kenya_pmid.csv)


#identifying the pmids 
head(Kenya_pmid.csv)

#Extract data from pmid in xml format 
kenya_pmid.xml <- fetch_PMID_data (my.data3$data$pmid)

# Retrival of the pmcid 
Kenyapmcid2.csv <- extract_article_ids(kenya_pmid.xml)

# # exporting data from the R console to the desktop
write.table(Kenyapmcid2.csv, "Kenyapmcid2.csv", row.names=F, sep ="\t")

#Bash script for designing the Kenyan article table
y <- "bash ../script/data_creation.sh"

#importing the data in R
PMID_PMC_Journal_Year_Kenya2.csv<- read.csv("PMID_PMC_Journal_Year_Kenya2.csv")

head(PMID_PMC_Journal_Year_Kenya2.csv)

#Bash script for data manipulation
z <- "bash ../script/script.sh"

#importing the data in R
PMID_PMC_Journal_Year_Kenya3.csv<- read.csv("PMID_PMC_Journal_Year_Kenya3.csv")

head(PMID_PMC_Journal_Year_Kenya3.csv)

