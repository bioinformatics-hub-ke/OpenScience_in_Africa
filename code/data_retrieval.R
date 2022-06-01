source("required_functions.R")

# Install packages if not already installed
# install.packages("easyPubMed")
# install.packages("devtools")
# devtools::install_github("dami82/easyPubMed")
# devtools::install_github("dami82/businessPubMed")


#load library packages R
library(devtools)
library(easyPubMed)
library(rentrez)
library(kableExtra)
library(dplyr)
library(businessPubMed)
library(githubinstall)
library(ggplot2)

#-------------------------------------------------------------
# Edit these parameters
# Country to query
country <- "Cameroon"

# Design the query for NCBI -- Change the year as necessary
my.query <- paste0(country,"[Affiliation]")
my.query <- paste(my.query, 'AND ("2011"[EDAT]:"2021"[EDAT])')

# for plotting/visualization purposes:
# set the minimum publications per journal to filter out those with few publications
minimum_publications <- 10
#-------------------------------------------------------------

# view the query
cat(my.query)

# define the regex filter for countries other than Kenya 
data("countries")
countries <- countries[countries != country]
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
pubmed_data <- extract_pubMed_data(pubMed_query = my.query, 
                                batch_size = 1000,
                                affi_regex_exclude = countries.filter)
# Check time t_final
t_final <- Sys.time()

# check how much time it took for the query
t_final - t_init

# save the pubmed data
saveRDS(pubmed_data, here::here("data", "processed", paste0(country,"_pubmed_data.RDS")))

pubmed_data <- readRDS(here::here("data", "processed", "Cameroon_pubmed_data.RDS"))
# pubmed_data <- readRDS(here::here("data", "processed", "Ghana_pubmed_data.RDS"))


# get the pubmed data from the query result
pubmed_data <- pubmed_data$data

#view the data set acquired 
View(pubmed_data)

# get the pubmed data for affiliations of the queried country.
country_pubmed_data <- dplyr::filter(pubmed_data, grepl(paste(country),address))

# subset the data to remove duplicated pubmed ids
country_pubmed_data_unique <- country_pubmed_data[!duplicated(country_pubmed_data$pmid),]

# number of pubmed id available
count(country_pubmed_data_unique)

#Extract data from PMID in XML format 
country_PMID <- fetch_PMID_data(country_pubmed_data_unique$pmid)

# Save the PMID data for the country used in the query
saveRDS(country_PMID, here::here("data", "processed", paste0(country,"_PMID_data.RDS")))

country_PMID <- readRDS(here::here("data", "processed","Cameroon_PMID_data.RDS"))
# country_PMID <- readRDS(here::here("data", "processed", "Ghana_PMID_data.RDS"))

# Retrival of the PMCID form the XML format
country_PMCID <- extract_article_ids(country_PMID)

View(country_PMCID)

# remove dulicates from the PMCID data
country_PMCID_filtered <- country_PMCID[!duplicated(country_PMCID$PMID),]

# filter to get only PMID in the original PMID query
country_PMCID_filtered <- country_PMCID_filtered[(country_PMCID_filtered$PMID 
                                                  %in% country_pubmed_data_unique$pmid),]

# merge the PMCID data with the queried PMID data
country_PMID_PMCID_data <- merge(country_pubmed_data_unique, 
                                 country_PMCID_filtered, by.x = "pmid", by.y = "PMID")

# View the data 
View(country_PMID_PMCID_data)

# add the publication status based on the presence or absence of PMCID
publications_open_access_status <- country_PMID_PMCID_data %>% 
  mutate(Status = ifelse(is.na(country_PMID_PMCID_data$PMCID), "closed", "open"))

# view the data
View(publications_open_access_status)

#######################################################
# Plots
#######################################################

# summarize to the number of publications by year
publications_by_year <- publications_open_access_status %>%
  group_by(year) %>%
  summarize(n())

# rename the columns
colnames(publications_by_year) <- c("year", "publications_number")

View(publications_by_year)

# Plot to visualize the number of papers by year
ggplot(publications_by_year, aes(x=year, y=publications_number, fill=year)) + 
  geom_bar(stat="identity") +
  ggtitle(paste("Number of Papers Per Year by", country, "Authors")) +
  xlab("Year") +
  ylab("Number of Publications") + 
  theme(legend.position = "none")
ggsave(here::here("data", "processed", paste0("papers_per_year_",country,".png")))


#----------------------------------------------------------------------
# summarize to get the number of papers per journal
journals_by_year <- publications_open_access_status %>%
  group_by(journal) %>%
  summarize(n())

# rename the columns
colnames(journals_by_year) <- c("journal", "publications_number")

# get journals with more than 50 published papers
journals_by_year_filtered <- journals_by_year[journals_by_year$publications_number > minimum_publications,]

tidy_name <- function(name, n_char) {
  ifelse(nchar(name) > (n_char - 2), 
         {substr(name, 1, n_char) %>% paste0(., "..")},
         name)
}

# abbreviate the long journals names.
journals_by_year_filtered <- journals_by_year_filtered %>% 
  mutate(journal_abbreviation = journals_by_year_filtered$journal %>% tidy_name(20))

# sort the data based on the number of publications
journals_by_year_filtered <- journals_by_year_filtered[order(-journals_by_year_filtered$publications_number),]

View(journals_by_year_filtered)

# Plot to visualize the number of paper per journal
ggplot(journals_by_year_filtered, aes(x=journal_abbreviation, y=publications_number)) + 
  geom_bar(stat="identity") +
  ggtitle("Number of Papers Per Journal") +
  xlab("Journal") +
  ylab("Number of Publications") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(here::here("data", "processed", paste0("papers_per_journal_",country,".png")))

#-----------------------------------------------------------------------------
# Summarize open access status of publication per year 
open_access_status <- as.data.frame(table(publications_open_access_status$year, publications_open_access_status$Status))

# rename columns
colnames(open_access_status) <- c("year", "status", "publications_number")

# view the data
View(open_access_status)

ggplot(data=open_access_status, aes(x=year, y=publications_number, fill=status)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  xlab("Year") +
  ylab("Number of Papers") +
  ggtitle("Relationship between years and no. of papers by openess")
ggsave(here::here("data", "processed", paste0("open_access_status_",country,".png")))