# R Workshop, March 21-22 2024
# Ty Wagner, USGS

# Exercise 2: Introduction to R - plotting
# Description: To gain experience using ggplot to visualize data
# Dataset: Mauch Chunk Lake (Big Bass reg). Nighttime Boat Electrofishing


# Load R packages
library(dplyr) # data manipulation
library(tidyverse) # data manipulation
library(lubridate) # work with dates
library(kableExtra) # make tables
# New libraries for Ex 2
library(ggplot2) # plotting data and spatial (simple feature; sf) objects
library(sf) # map creation (simple features)
library(spData) # provides access to polygons of US states

# Exercise 1 --------------------------------------------------------------



# Read in data
# The "././" syntax is backing out 2 directories from our current file
# so we can navigate into the 02_Data folder where our data are located
dat <- read_csv("././02_Data/Mauch_lake/Mauch_lake_surveys_NBE.csv")

str(dat)


# Data notes:
# Note: in 1986 surveys 35946, 35947 are two seprate entries in the ARDB but have 
# exact catch/effort - likely duplicated. Thus only use one (i.e., 35946) and 
# not include the other 						
# Note: in 2000 surveys 35970 and 35971 are likely duplicated in the ARDB, 
# thus on retain one for analysis (i.e., 35970) and not inculde the other.   						


# Clean data
dat_clean <- dat %>% 
  select(-c(`Water Name`, `Unique Water Site ID`,
            `Water Section ID`, `Group Size Fish Length`,
            `Survey Purpose`, `Water Site Comment`,
            `Water Site Survey Comment`, Month, Description)) %>% # remove unwanted columns
  rename(Water_site_survey_ID  = `Water Site Survey ID`,
         Lat = `Survey Site Lat DD`,
         Long = `Survey Site Lon DD`,
         Species = `Fish Species Name`,
         Number_caught = `Number Caught`,
         Effort = `Effort Hours`) %>% # rename columns
  filter(Water_site_survey_ID != 35947 & Water_site_survey_ID != 35971)  %>% # remove duplicate surveys
  mutate(Date = mdy_hm(`Site Date`)) %>% # Convert Site Date to date format
  select(-c(`Site Date`)) %>% 
  # filter(Species == "Largemouth Bass") %>% # Select Largemouth bass
  rename_with(tolower) # make all column names lower case
  
  
str(dat_clean)

head(dat_clean, 20)

# Calculate total catch for each year, survey, and species
dat_tot_catch <- dat_clean %>% 
  group_by(year, water_site_survey_id, species) %>% # group data 
  mutate(total_catch = sum(number_caught)) %>% # sum over catch for each variable in group_by
  ungroup() %>% # ungroup data
  distinct(year, water_site_survey_id, species, .keep_all=TRUE) %>% # retain distinct combos since we don't need size-specific numbers here
  select(-c(number_caught)) # remove number caught, no longer needed


# Species that were not caught in a given survey and year are not recorded as zero catch,
# but we would like to do so. Thus, we have to input missing species names into surveys
# and years where they were not recorded and give them a zero catch values

dat_tot_catch2 <- dat_tot_catch %>% 
  select(water_site_survey_id, year, species, total_catch, effort) %>% # Select columns of interest for summarizing
  complete(nesting(water_site_survey_id, year), species,  fill = list(total_catch=0)) %>% # Input missing species for surveys and years
  arrange(year, species) # sort by year and species

# Look at first few rows of tot_catch2
head(dat_tot_catch2, 30)


# Lets replace the NA values for imputed species efforts to the actual effort of the survey
dat_tot_catch2 <- dat_tot_catch2 %>% 
  group_by(water_site_survey_id) %>% 
  mutate(effort = replace_na(mean(effort, na.rm=T))) %>% 
  ungroup()

tail(dat_tot_catch2, 30)

# Summarize mean catch (across surveys for each year and species)
# sample size, and mean effort.
table1 <- dat_tot_catch2 %>% 
  group_by(year, species) %>% 
 summarize(n = n(), total_catch = sum(total_catch, na.rm=T),
           total_effort = sum(effort, na.rm=T)) %>% 
  arrange(year, species)

# Create a table of the catch summary
table1 %>% 
  kbl(caption = "Catch and effort (hrs) summary table.", digits=2,
      col.names = c("Year",
                    "Species",
                    "n (surveys)",
                    "Total catch",
                    "Total effort (hrs)")) %>%
  row_spec(0,bold=TRUE) %>% 
  kable_classic(full_width = F, html_font = "Cambria") 


# Exercise 2 --------------------------------------------------------------


