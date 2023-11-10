#Communtiy change analyses
#Joshua G. Smith; jossmith@mbayaq.org

rm(list=ls())

librarian::shelf(tidyverse, here, ggplot2, mvabund, lattice)


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

#--------------------------------swath----------------------------------------#

swath_long <- swath_build1 %>% 
              pivot_longer(cols = 12:68, names_to = "taxa", values_to = "density") %>%
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
          dplyr::select(13:69)%>%
          as.data.frame() %>%
          mutate(across(everything(), ~coalesce(., 0))) 


swath_R <- swath_long %>%
           dplyr::select(!(c(common_name, taxonomic_level, trophic_ecology)))%>%
           pivot_wider(names_from = "taxa",values_from = "density")%>%
           #select envr vars, note mvabund does not take function formula
            dplyr::select(MHW)%>%
            as.data.frame()


swath_Q <- swath_long %>%
          dplyr::select(taxa, trophic_ecology) %>%
          distinct(taxa, .keep_all = TRUE)%>%
          #set spp to row name
          column_to_rownames(var = "taxa")%>%
          mutate(trophic_ecology = factor(trophic_ecology))%>%
          as.data.frame() 


#----------------------------------UPC-----------------------------------------#


upc_long <- upc_build1 %>% 
  pivot_longer(cols = 12:50, names_to = "taxa", values_to = "density") %>%
  mutate(survey_method = "UPC") %>%
  #join traits
  left_join(spp_attribute, by = c("survey_method","taxa")) %>% # %>% filter(is.na(trophic_ecology)) %>% distinct(taxa)
  mutate(density = coalesce(density, 0))

upc_L <- upc_long %>%
  dplyr::select(!(c(common_name, taxonomic_level, trophic_ecology)))%>%
  pivot_wider(names_from = "taxa",values_from = "density")%>%
  dplyr::select(13:51)%>%
  replace(is.na(.),0) %>%
  as.data.frame() %>% ceiling(.)


upc_R <- upc_long %>%
  dplyr::select(!(c(common_name, taxonomic_level, trophic_ecology)))%>%
  pivot_wider(names_from = "taxa",values_from = "density")%>%
  #select envr vars, note mvabund does not take function formula
  dplyr::select(MHW)%>%
  as.data.frame()


upc_Q <- upc_long %>%
  dplyr::select(taxa, trophic_ecology) %>%
  distinct(taxa, .keep_all = TRUE)%>%
  #set spp to row name
  column_to_rownames(var = "taxa")%>%
  mutate(trophic_ecology = factor(trophic_ecology))%>%
  as.data.frame() 


#----------------------------------FISH----------------------------------------#


fish_long <- fish_build1 %>% 
  pivot_longer(cols = 12:64, names_to = "taxa", values_to = "density") %>%
  mutate(survey_method = "Fish") %>%
  #join traits
  left_join(spp_attribute, by = c("survey_method","taxa")) # %>% filter(is.na(trophic_ecology)) %>% distinct(taxa)


fish_L <- fish_long %>%
  dplyr::select(!(c(common_name, taxonomic_level, trophic_ecology)))%>%
  pivot_wider(names_from = "taxa",values_from = "density")%>%
  dplyr::select(13:65)%>%
  as.data.frame() %>%
  mutate(across(everything(), ~coalesce(., 0))) 


fish_R <- fish_long %>%
  dplyr::select(!(c(common_name, taxonomic_level, trophic_ecology)))%>%
  pivot_wider(names_from = "taxa",values_from = "density")%>%
  #select envr vars, note mvabund does not take function formula
  dplyr::select(MHW)%>%
  as.data.frame()


fish_Q <- fish_long %>%
  dplyr::select(taxa, trophic_ecology) %>%
  distinct(taxa, .keep_all = TRUE)%>%
  #set spp to row name
  column_to_rownames(var = "taxa")%>%
  mutate(trophic_ecology = factor(trophic_ecology))%>%
  as.data.frame() 



################################################################################
# trait based models

#------------------------------swath-------------------------------------------#

####WARNING: Slow
swath_mod<- traitglm(swath_L, swath_R, swath_Q,
                         method="glm1path", family = "negative.binomial")

#glance at results
a        = max( abs(swath_mod$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th.anom = levelplot(t(as.matrix(swath_mod$fourth.corner)), xlab="Environmental Variables",
                          ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                          scales = list( x= list(rot = 45)))
print(plot.4th.anom)

#--------------------------------upc-------------------------------------------#

####WARNING: Slow 
upc_mod<- traitglm(upc_L, upc_R, upc_Q,
                     method="glm1path", family = "negative.binomial")

#glance at results
a        = max( abs(upc_mod$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th.anom = levelplot(t(as.matrix(upc_mod$fourth.corner)), xlab="Environmental Variables",
                          ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                          scales = list( x= list(rot = 45)))
print(plot.4th.anom)


#-------------------------------fish-------------------------------------------#

####WARNING: Slow 
fish_mod<- traitglm(fish_L, fish_R, fish_Q,
                   method="glm1path", family = "negative.binomial")


a        = max( abs(fish_mod$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th.anom = levelplot(t(as.matrix(fish_mod$fourth.corner)), xlab="Environmental Variables",
                          ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                          scales = list( x= list(rot = 45)))
print(plot.4th.anom)


################################################################################
#Tidy coefs

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
                                    "trophic_ecologyPlanktivore" = "Planktivore"),
         method = "swath")
                                 

coef_fish <- fish_mod$fourth %>%
  as.data.frame()%>%
  rownames_to_column(var="Trophic ecology")%>%
  rename("After" = MHWafter,
         "Before" = MHWbefore,
         "During" = MHWduring)%>%
  mutate(`Trophic ecology` = recode(`Trophic ecology`,
                                    "trophic_ecologyHerbivore" = "Herbivore",
                                    "trophic_ecologyMacroinvertivore" = "Macroinvertivore",
                                    "trophic_ecologyMicroinvertivore" = "Microinvertivore",
                                    "trophic_ecologyPiscivore" = "Piscivore",
                                    "trophic_ecologyPlanktivore" = "Planktivore"),
         method = "fish")


coef_upc <- upc_mod$fourth %>%
  as.data.frame()%>%
  rownames_to_column(var="Trophic ecology")%>%
  rename("After" = MHWafter,
         "Before" = MHWbefore,
         "During" = MHWduring)%>%
  mutate(`Trophic ecology` = recode(`Trophic ecology`,
                                    "trophic_ecologyPlanktivore" = "Planktivore"),
         method = "upc")


fourth_dat <- rbind(coef_swath, coef_fish, coef_upc) %>%
              pivot_longer(2:4, names_to = "Heatwave period",values_to = "Coef") %>%
              #center and scale
              mutate(beta_sd=scale(Coef, center=F, scale=T),
                     `Heatwave period` = factor(`Heatwave period`,levels = c("Before","During","After"))) 


################################################################################
#plot


# Plot all four corner results
ggplot(fourth_dat, aes(x = `Heatwave period`, y = `Trophic ecology`, fill = beta_sd)) +
  facet_col(~ method, scales = "free_y", space = "free",
            strip.position = "top") +
  geom_tile() +
  # Plot 0s (no interaction term)
  geom_point(data=fourth_dat %>% filter(Coef==0), shape="x") +
  # Labels
  labs(x="Indicator", y="Thermal affinity", tag="") +
  scale_fill_mba("drifters2",n_colors=50, type = "continuous", name = "Effect") +
  # Legend
  #scale_fill_gradient2(name="Effect",
   #                    midpoint=0,
    #                   breaks=c(-1.5, 0, 2), 
     #                  labels=c("-", "0", "+"),
      #                 low="navy", high="darkred", mid="white") +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black")) +
  # Theme
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(-5,0,0,0),
        axis.title.y = element_text(vjust = 4), 
        # axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1))




















