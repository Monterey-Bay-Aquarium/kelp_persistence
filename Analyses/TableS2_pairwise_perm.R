#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

librarian::shelf(tidyverse, here, vegan, ggplot2, cluster, ggforce, 
                 pairwiseAdonis, broom)




################################################################################
#set directories and load data
basedir <- here::here("output")

mv_dat <- load(file.path(basedir, "monitoring_data/processed/multivariate_data.Rdata")) 

tab_dir <- here::here("tables")



################################################################################
#prepare data 

#----------------process standardized data--------------------------------------

#define group vars
stan_group_vars <- stan_group_vars %>% 
                    mutate(site_period = paste(site, MHW),
                           outbreak_period = ifelse(year <2014, "before","after"),
                           site_outbreak = paste(site, ifelse(year <2014, "before","after")))


################################################################################
#build pairwise PERMANOVA for regional comparison

#pariwise permanova for regional
set.seed(1985)

pair_perm <- pairwise.adonis2(stan_max_distmat ~ outbreak_period, 
                                    data = stan_group_vars, 
                                    permutations = 999,
                                    num_cores = 20)

##### EXTRACT OUTPUT
# Remove the different list item
pair_perm$parent_call <- NULL

# make it a data frame
regional_result_tab <- pair_perm %>% map_dfr(~tidy(.x), .id="name")%>%
  filter(term == "outbreak_period")


################################################################################
#build pairwise PERMANOVA for site-level comparison

#pariwise permanova for regional
set.seed(1985)

pair_perm <- pairwise.adonis2(stan_max_distmat ~ site_outbreak, 
                              data = stan_group_vars, 
                              permutations = 999,
                              num_cores = 20)

##### EXTRACT OUTPUT
# Remove the different list item
pair_perm$parent_call <- NULL

# make it a data frame
site_result_tab <- pair_perm %>% map_dfr(~tidy(.x), .id="name")%>%
  filter(term == "site_outbreak") %>%
  mutate(group_1 = str_extract(name, ".*(?=_vs)"), #extract everything before "_vs"
         group_2 = str_extract(name, "(?<=vs_).*"), #extract everything after "vs_"
         site_1 = str_extract(group_1, "^[^ ]+"), #extract everything before " "
         site_2 = str_extract(group_2, "^[^ ]+"), #extract everything before " "
         period_1 = str_extract(group_1, "\\w+$"),
         period_2 = str_extract(group_2, "\\w+$")) %>%
  #filter matrix
  filter(site_1 == site_2,
         period_1 == "before")%>%
  dplyr::select(Site = site_1,
                DF = df, `Sum of Sq.` = SumOfSqs, `R-squared` = R2, `F` = statistic, `P(perm)` = p.value) %>%
  mutate(Site = gsub("_", " ", Site),
         Site = tolower(Site),
         Site = tools::toTitleCase(Site),
         Site = sub("(.*)( Dc| Uc)$", "\\1\\U\\2", Site, perl = TRUE))



write.csv(site_result_tab, file.path(tab_dir, "TableS2_pairwise_permanova.csv"), row.names = FALSE)


