#Mixed effects model of patch drivers
#Joshua G. Smith
#May 9, 2023

rm(list=ls())
librarian::shelf(tidyverse, here, vegan, cowplot, ggpubr, mgcv, gratia)


################################################################################
#set directories and load data
basedir <- here::here("output")
figdir <- here::here("figures")
tabdir <- here::here("tables")

#load sampling data
swath_raw <- read.csv(file.path(basedir, "monitoring_data/processed/kelp_swath_counts_CC.csv")) %>%
  #convert 60m2 kelp density to m2
  mutate(stipe_den_m2 = macrocystis_pyrifera/60)

#load model predictors
mod_predict <- readRDS(file.path(basedir, "environmental_data/processed/envr_at_pisco_sites.Rds"))

#load published urchin behavior data
#accessed from "https://github.com/joshgsmith/PatchDynamics/blob/main/Data/raw_data/patch_data.csv"
patch_raw <- read.csv(file.path(basedir, "monitoring_data/processed/patch_data.csv"))


################################################################################
#process environmental data

mod_predict_build1 <- mod_predict %>%
  dplyr::filter(year>=2007 & year <= 2020)%>%
  #calculate annual means
  group_by(site, year)%>%
  summarize(across(8:46,mean, na.rm=TRUE)) %>%
  #drop monthly statistics
  dplyr::select(!(c(cuti_month_baseline, cuti_month_sd,
                    beuti_month_baseline, beuti_month_baseline_sd,
                    sst_month_baseline,
                    sst_month_baseline_sd, sst_month_anom, 
  ))) %>% #these are calculated at monthly intervals so irrelevant here. 
  dplyr::mutate_all(~ ifelse(is.nan(.), NA, .))

################################################################################
#calculate baseline kelp density

kelp_baseline <- swath_raw %>% dplyr::select(year, site, zone, transect, macrocystis_pyrifera)%>%
  filter(year <= 2013)%>%
  group_by(site)%>%
  summarize(baseline_kelp = mean(macrocystis_pyrifera, na.rm=TRUE),
            baseline_kelp_cv = (sd(macrocystis_pyrifera, na.rm = TRUE) / mean(macrocystis_pyrifera, na.rm = TRUE)) * 100
  )

################################################################################
#calculate mean urchin densities

urchin_density <- swath_raw %>% dplyr::select(year, site, zone, transect, strongylocentrotus_purpuratus)%>%
  group_by(year, site)%>%
  summarize(urchin_density = mean(strongylocentrotus_purpuratus, na.rm=TRUE))

################################################################################
#simulate behavior

##Note: prop_pur_exp is average across replicate transects, not calculated directly 
#from the density averages. 

#examine data

ggplot(data = patch_raw, aes(x = mac_stipe_den, y = prop_pur_exp)) +
  geom_point() +
  xlab("mac_stipe_den") +
  ylab("prop_pur_exposed") +
  ggtitle("Scatter Plot of prop_pur_exposed vs. mac_stipe_den")

#remove six outliers

exp_build1 <- patch_raw %>% filter(!(prop_pur_exp > 0.25 & mac_stipe_den > 4))

# Remove rows with NA values in mac_stipe_den column
exp_build1 <- exp_build1[!is.na(exp_build1$mac_stipe_den), ]

# Fit modified exponential decay model using nls()
fit <- nls(prop_pur_exp ~ a * (1 - exp(-b * mac_stipe_den)) + c, data = exp_build1, start = list(a = 1, b = 1, c = 0))

# Create line data using the cleaned exp_build1 data frame
line_data <- data.frame(mac_stipe_den = seq(min(exp_build1$mac_stipe_den), max(exp_build1$mac_stipe_den), length.out = 100))
line_data$prop_pur_exposed <- predict(fit, newdata = line_data)

# Plot with fitted line
ggplot(data = exp_build1, aes(x = mac_stipe_den, y = prop_pur_exp)) +
  geom_point() +
  geom_line(data = line_data, aes(x = mac_stipe_den, y = prop_pur_exposed), color = "red") +
  xlab("mac_stipe_den") +
  ylab("prop_pur_exposed") +
  ggtitle("Scatter Plot of prop_pur_exposed vs. mac_stipe_den")

# Calculate R-squared value
residuals <- residuals(fit)
y_mean <- mean(exp_build1$prop_pur_exp, na.rm=TRUE)
ss_total <- sum((exp_build1$prop_pur_exp - y_mean)^2, na.rm=TRUE)
ss_residual <- sum(residuals^2)
r_squared <- 1 - (ss_residual / ss_total)

# Get the formula
formula_text <- as.character(as.formula(fit))


################################################################################
#process response variable and exposed urchin 
# Simulate behavior for 'swath_raw' using the fitted model
swath_raw$predicted_values <- predict(fit, newdata = data.frame(mac_stipe_den = swath_raw$stipe_den_m2)) + sample(residuals(fit), size = 1) 

swath_raw <- swath_raw %>% mutate(line_fit = -0.6487 * (1 - exp(-5.1438 * (stipe_den_m2))) + 0.7169)

response_vars <- swath_raw %>% dplyr::select(year, site, zone, transect,
                                             macrocystis_pyrifera, predicted_values) %>%
  group_by(year, site) %>%
  dplyr::summarize(stipe_mean = mean(macrocystis_pyrifera, na.rm = TRUE),
                   line_fit = -0.6487 * (1 - exp(-5.1438 * (stipe_mean / 60))) + 0.7169,
                   prop_urch_exp = mean(predicted_values,na.rm=TRUE))

################################################################################
#join sampling data with environmental

mod_dat <- left_join(response_vars, mod_predict_build1, by=c("year","site")) %>%
  #create lagged variables
  mutate(#npp_lag1 = lag(npp,1),
    #npp_lag2 = lag(npp,2),
    #designate resistant site
    resistance = ifelse(site == "HOPKINS_UC" |
                          site == "CANNERY_UC" |
                          site == "SIREN" |
                          site == "CANNERY_DC",
                        #site == "BUTTERFLY_DC",
                        "resistant","transitioned")
  ) %>%
  left_join(kelp_baseline, by="site") %>% 
  left_join(urchin_density, by = c("year","site"))

################################################################################
#build GAM

#  Build full meta-GAM
full_gam <- mgcv::gam(stipe_mean ~ 
                              s(prop_urch_exp) +
                              s(vrm_sum) +
                              s(bat_mean) +
                              s(beuti_month_obs) +
                              s(npp_ann_mean) +
                              s(wave_hs_max)+
                              s(orb_vmax) +
                              s(slope_mean) +
                              s(sst_month_obs) +
                              s(baseline_kelp) +
                            s(baseline_kelp_cv) +
                            s(urchin_density) +
                            s(year, bs = "cc"),
                            data = mod_dat %>% filter(urchin_density < 2500),
                            method = "REML",
                            #use shrinkage selection (Marra and Wood 2011) 
                            #with a double penalty approach
                            select = TRUE) 
            

#drop shrinkage terms
red_gam <- mgcv::gam(stipe_mean ~ 
                       s(prop_urch_exp) +
                       #s(vrm_sum) +
                       #s(bat_mean) +
                       #s(beuti_month_obs) +
                       s(npp_ann_mean) +
                       #s(wave_hs_max)+
                       #s(orb_vmax) +
                       #s(slope_mean) +
                       s(sst_month_obs) +
                       s(baseline_kelp),
                      # s(baseline_kelp_cv) +
                      # s(urchin_density), 
                      # s(year, bs = "cc"),
                     data = mod_dat,
                     method = "REML",
                     #use shrinkage selection (Marra and Wood 2011) 
                     #with a double penalty approach
                     select = TRUE) 


summary(red_gam)
plot(red_gam)

# Get the summary of the GAM model
summary_info <- summary(red_gam)

# Extract terms and p-values
gam_terms <- data.frame(summary_info$s.table) %>%
  tibble::rownames_to_column(var = "smooth") %>%
  rename("EDF" = "edf")

################################################################################
#export model table
library(flextable)
set_flextable_defaults(background.color = "white")

A <- as_flextable(full_gam) 
A <- add_header_lines(A,values = "Full model")

save_as_image(A, path = file.path(figdir, "archive/gam_a.png"), dpi=600)


B <- flextable::as_flextable(red_gam)
B <-flextable::add_header_lines(B,values = "Reduced model")

save_as_image(B, path = file.path(figdir, "archive/gam_b.png"), dpi = 600)


# Load PNG images
img_a <- png::readPNG(file.path(figdir, "archive/gam_a.png"), native = TRUE)
img_b <- png::readPNG(file.path(figdir, "archive/gam_b.png"), native = TRUE)

# Create a blank ggplot
plot <- ggplot() +
  theme_void()

plot <- plot +
  annotation_custom(grid::rasterGrob(img_a), xmin = 0, xmax = 1, ymin = 0.5, ymax = 1) +
  annotation_custom(grid::rasterGrob(img_b), xmin = 0, xmax = 1, ymin = 0, ymax = 0.5)

#print(plot)

ggsave(file.path(tabdir, "TableS4_GAM_results.png"), plot, dpi = 600,
       bg = "white", width = 4, height = 6, units = "in")


################################################################################
#structure data for plotting

sm_dat <- gratia::smooth_estimates(red_gam) %>%
  add_confint()%>%
  pivot_longer(cols = 6:9, names_to = "var", values_to = "var_val") %>%
  left_join(gam_terms, by = "smooth") %>%
  mutate(smooth = str_replace(smooth, "s\\((.+)\\)", "\\1"),
         smooth = str_replace_all(smooth, "_", " "),
         smooth = str_to_sentence(smooth)) %>%
  mutate(smooth = factor(smooth, levels = c("Prop urch exp","Baseline kelp","Sst month obs","Npp ann mean")))

#add residuals to data
resid_dat <- mod_dat %>% 
  ungroup()%>%
  dplyr::select(stipe_mean, prop_urch_exp, npp_ann_mean, sst_month_obs, baseline_kelp) %>%
  na.omit() %>% add_partial_residuals(red_gam) %>%
  #make smooth vars longer
  pivot_longer(cols = 6:ncol(.), names_to = "smooth", values_to = "res_val") %>%
  #make predictors longer
  pivot_longer(cols = c("stipe_mean","prop_urch_exp","npp_ann_mean",
                        "sst_month_obs","baseline_kelp"), 
               names_to = "pred_var", values_to = "pred_val")%>%
  dplyr::select(smooth, res_val, pred_var, pred_val) %>%
  mutate(match_var = paste0("s(",pred_var,")"))%>%
  filter(smooth == match_var) %>%
  dplyr::select(smooth, res_val, pred_val)%>%
  mutate(smooth = str_replace(smooth, "s\\((.+)\\)", "\\1"),
         smooth = str_replace_all(smooth, "_", " "),
         smooth = str_to_sentence(smooth))

################################################################################
#plot GAM

my_theme <-  theme(axis.text=element_text(size=8, color = "black"),
                   axis.text.y = element_text(color = "black",size=7),
                   axis.text.x = element_text(color = "black",size=7),
                   axis.title=element_text(size=9, color = "black"),
                   plot.tag=element_text(size= 7, color = "black"), #element_text(size=8),
                   plot.title =element_text(size=9, face="bold", color = "black"),
                   # Gridlines 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   # Legend
                   legend.key = element_blank(),
                   legend.background = element_rect(fill=alpha('blue', 0)),
                   legend.key.size  = unit(1, "lines"), 
                   legend.text = element_text(size = 10, color = "black"),
                   legend.title = element_text(size = 10, color = "black"),
                   #legend.spacing.y = unit(0.75, "cm"),
                   #facets
                   strip.background = element_blank(),
                   strip.text = element_text(size = 8 , face="plain", color = "black", vjust=-0.5),
                   strip.placement = "outside"
)

pval <- sm_dat %>% dplyr::select(smooth, p.value, EDF) %>% distinct() %>%
  mutate(color = ifelse(p.value <= 0.051, "red","black"),
         label = paste0(ifelse(p.value < 0.001, "p < 0.001", paste("p = ", round(p.value, 3))),
                        ", EDF: ",round(EDF,2))
  )

# Assuming 'smooth' is the variable you want to facet by
sm_dat$smooth <- factor(sm_dat$smooth, levels = c("Baseline kelp", "Npp ann mean", "Prop urch exp", "Sst month obs"))


g <- ggplot(sm_dat) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = var_val), color = "black", linetype = "dashed", fill = NA) +
  #add residuals
  geom_point(aes(x = pred_val, y = res_val),
            data = resid_dat, colour = "steelblue3", alpha = 0.2) +
  geom_line(aes(x = var_val, y = est), lwd = 0.8) +
  facet_wrap(~factor(smooth, levels = c("Prop urch exp", "Baseline kelp","Sst month obs", "Npp ann mean")), scales = "free", nrow = 1,
             labeller = as_labeller(
               c("Prop urch exp" = "Proportion of purple \nsea urchins exposed", "Baseline kelp" = "Baseline kelp density \n(stipes per 60 m²)", "Npp ann mean" = "Net primary productivity \n(annual mean)",
                 "Sst month obs" = "Sea surface \ntemperature (°C)")), 
             strip.position = "bottom"
  )+
  #add summary stats
  geom_text(aes(x = -Inf, y = Inf, 
                label = pval$label),
            hjust = -0.1, vjust = 1.3, size = 2.5, data = pval, color = "black") +
  labs(y = "Partial effect", x  = "", fill = "Scaled residual \ndensity", tag = "A",
       title = "Predictors of kelp density")+
  theme_bw() + my_theme + theme(plot.margin = margin(5,5,0,30),
                                plot.tag.position = c(-.045, 1))


envr_plot <- mod_dat %>% 
             mutate(urchin_m2 = urchin_density/60)%>%
             dplyr::select(year, site, resistance,
                           baseline_kelp, orb_vmax, sst_month_obs, slope_mean,
                           vrm_mean, beuti_month_obs, baseline_kelp_cv, bat_mean,
                           npp_ann_mean, urchin_m2, wave_hs_max, prop_urch_exp)%>%
             pivot_longer(4:15, names_to = "variable",values_to = "value") %>%
             group_by(site, resistance, variable)  %>% summarize(val_mean = mean(value, na.rm=TRUE))


# Perform t-test for each variable
t_test_data <- envr_plot %>%
  dplyr::select(variable, resistance, val_mean, site)

t_test_results <- t_test_data %>%
  group_by(variable) %>%
  do(broom::tidy(t.test(val_mean ~ resistance, data = .)))%>%
  mutate(color = ifelse(p.value <= 0.051, "red","black"),
         label = paste0(ifelse(p.value < 0.001, "p < 0.001", paste("p = ", round(p.value, 3)))))


g2 <- ggplot(data = envr_plot %>% 
               mutate(resistance = ifelse(resistance == "transitioned","Transitioned","Persistent")),
             aes(x = resistance, y = val_mean)) +
  geom_boxplot(color = "black", aes(fill = resistance), show.legend = TRUE) +
  geom_jitter(width = 0.1, height = 0.3, alpha = 0.2, size=1) +
  geom_text(aes(x = -Inf, y = Inf, 
                label = t_test_results$label),
            hjust = -0.1, vjust = 1.3, size = 2.5, data = t_test_results, color = t_test_results$color) +
  facet_wrap(~variable, scales = "free", 
             labeller = as_labeller(
               c(
                 "baseline_kelp" = "Baseline kelp density \n(stipes per 60 m²)",
                 "baseline_kelp_cv" = "Baseline kelp \ncoefficient of variation",
                 "bat_mean" = "Average \ndepth (m)",
                 "beuti_month_obs" = "Biologically effective \nupwelling index",
                 "npp_ann_mean" = "Net primary \nproductivity",
                 "orb_vmax" = "Orbital velocity",
                 "prop_urch_exp" = "Proportion of purple \nsea urchins exposed",
                 "slope_mean" = "Reef \nslope",
                 "sst_month_obs" = "Sea surface \ntemperature (°C)",
                 "urchin_m2" = "Sea urchin density \n(no. per m²)",
                 "vrm_mean" = "Rugosity \n ",
                 "wave_hs_max" = "Wave \nheight (m)"
               )
             ),
             strip.position = "left")+
  #ylim(1,4)+
  labs(x = "", y="", tag = "B", title = "Features of persistent and transitioned forests")+
  MBAcolors::scale_fill_mba("mba2", type = "discrete", rev=TRUE, name = "Site type")+
  theme_classic()+
  my_theme+
  scale_x_discrete(labels = c("Persistent", "Transitioned")) +
  theme(strip.placement = "outside",
        axis.text.x = element_blank(),
        legend.position = "bottom") + 
  #geom_blank(aes(y = y_min)) +
  #geom_blank(aes(y = y_max)) + 
  theme(plot.margin = margin(2,5,5,5))



full_plot <- ggarrange(g,g2, nrow=2,  heights = c(0.4,0.6))
#full_plot




ggsave(full_plot, filename=file.path(figdir, "Fig5_GAM_predictors.png"), 
       width=7.5, height=9, bg="white", units="in", dpi=600,
       device = "png")








