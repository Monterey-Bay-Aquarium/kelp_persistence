#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, ggplot2, mvabund)


################################################################################
#set directories and load data
basedir <- here::here("output")
figdir <- here::here("figures")
tabdir <- here::here("tables")

#load raw dat
swath_raw <- read.csv(file.path(basedir, "monitoring_data/processed/kelp_swath_counts_CC.csv"))

upc_raw <- read.csv(file.path(basedir, "monitoring_data/processed/kelp_upc_cov_CC.csv")) 

fish_raw <- read.csv(file.path(basedir, "monitoring_data/processed/kelp_fish_counts_CC.csv"))

#load species attribute table
spp_attribute <- read.csv(file.path(tabdir,"TableS1_spp_table.csv")) %>% janitor::clean_names() %>%
                    mutate(taxa = tolower(gsub(" ", "_", taxa)),
                          #fix names
                          taxa = ifelse(taxa == "loxorhynchus_crispatus__or_scyra_acutifrons",
                                        "loxorhynchus_crispatus_scyra_acutifrons",
                                        taxa))


################################################################################

#drop species that were never encountered
swath_build1 <- swath_raw %>% dplyr::select(where(~ any(. != 0)))
upc_build1 <- upc_raw %>% dplyr::select(where(~ any(. != 0)))
fish_build1 <- fish_raw %>% dplyr::select(where(~ any(. != 0)))

################################################################################
#reshape data and join traits
swath_long <- swath_raw %>% 
              pivot_longer(cols = 12:70, names_to = "taxa", values_to = "density") %>%
              mutate(survey_method = "Swath") %>%
              #join traits
              left_join(spp_attribute, by = c("survey_method","taxa")) %>%
              #fix taxa
              mutate(trophic_ecology = case_when(
                taxa == "cancridae" ~ "Macroinvertivore",
                taxa == "linckia_columbiae" ~ "Detritivore (algal)",
                taxa == "solaster_dawsoni" ~ "Macroinvertivore",
                TRUE ~ trophic_ecology
              )) 

swath_L <- swath_long %>%
          dplyr::select(!(c(common_name, taxonomic_level, trophic_ecology)))%>%
          pivot_wider(names_from = "taxa",values_from = "density")%>%
          dplyr::select(13:71)%>%
          as.data.frame() %>%
          mutate(across(everything(), ~coalesce(., 0))) %>% sample_n(size = 25, replace = FALSE)


swath_R <- swath_long %>%
           dplyr::select(!(c(common_name, taxonomic_level, trophic_ecology)))%>%
           pivot_wider(names_from = "taxa",values_from = "density")%>%
           #select envr vars, note mvabund does not take function formula
            dplyr::select(MHW)%>%
            as.data.frame() %>% sample_n(size = 25, replace = FALSE)


swath_Q <- swath_long %>%
          dplyr::select(taxa, trophic_ecology) %>%
          distinct(taxa, .keep_all = TRUE)%>%
          #set spp to row name
          column_to_rownames(var = "taxa")%>%
          mutate(trophic_ecology = factor(trophic_ecology))%>%
          as.data.frame() 

dput(swath_L)

################################################################################
# trait based model

#============================CCFRP==============================================

swath_mod<- traitglm(swath_L, swath_R, swath_Q,
                         method="glm1path", family = "negative.binomial")

summary(swath_mod)

#test plot
swath_mod$fourth

#check residuals
plot(swath_mod)

library(lattice)

a        = max( abs(swath_mod$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th.anom = levelplot(t(as.matrix(swath_mod$fourth.corner)), xlab="Environmental Variables",
                          ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                          scales = list( x= list(rot = 45)))
print(plot.4th.anom)



################################################################################
#Tidy coef for anom

coef_swath <- swath_mod$fourth %>%
  as.data.frame()%>%
  rownames_to_column(var="Trophic ecology")%>%
  rename("After" = MHWafter,
         "Before" = MHWbefore,
         "During" = MHWduring)%>%
  mutate(`Trophic ecology` = recode(`Trophic ecology`,
                                     "trophic_ecologyAutotroph" = "Autotroph",
                                     "trophic_ecologyDetritivore (algal)" = "Detritivore (algal)",
                                    "trophic_ecologyHerbivore" = "Herbivore",
                                    "trophic_ecologyMacroinvertivore" = "Macroinvertivore",
                                    "trophic_ecologyMicroinvertivore" = "Microinvertivore",
                                    "trophic_ecologyPlanktivore" = "Planktivore"
                                    ))%>%
  mutate(method = "Swath")







dput(swath_L)








# Calculate the total density for each site at each year
site_year_total_density <- swath_long %>%
  mutate(transition_site = ifelse(site == "HOPKINS_UC" | site == "CANNERY_UC" |
                                    site == "SIREN" | site == "CANNERY_DC","no","yes"))%>%
  group_by(year, site, trophic_ecology, transition_site) %>%
  summarize(total_density = sum(density, na.rm=TRUE))

# Calculate the mean total density for each MHW level
mean_total_density <- site_year_total_density %>%
  group_by(year, trophic_ecology, transition_site) %>%
  summarize(mean_density = mean(total_density, na.rm=TRUE))

# Calculate the proportion using the mean total density
proportion_data <- mean_total_density %>%
  group_by(year, transition_site) %>%
  mutate(proportion = mean_density / sum(mean_density))

# Create the bar plot
ggplot(proportion_data, aes(x = year, y = proportion, fill = trophic_ecology)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Proportion of Trophic Ecology Levels by MHW",
       x = "MHW",
       y = "Proportion") +
  facet_wrap(~transition_site)+
  scale_fill_mba(exhibit = "mba2", n_colors = 6, type = "discrete") + # You can choose a different color palette
  theme_minimal()


