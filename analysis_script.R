library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(DHARMa)
library(RColorBrewer)
library(reshape2)
library(broom)
library(ggsignif)
library(patchwork)
library(cowplot)
library(here)

set.seed(1234)

### Functions

# function to sort data into 5% quantile bins 
twenty_quantiles <- function(x) {quantile(x, seq(0.05, 1, by = 0.05), type = 7, na.rm = TRUE)}

# function to test fit of linear mixed models 
model_fit <- function(x) {simulateResiduals(fittedModel = x, n = 250) %>%
    plot(., asFactor = F)}

# functions to create raincloud plots
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )

theme_anes <- 
  theme_classic() +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),  
        plot.title = element_text(size=14, hjust = 0.5),
        plot.subtitle = element_text(size=12, hjust = 0.5))

# animal treatment group assignments
anesthesia <- c("RJS15", "PS15", "QS15", "UU15", "QH15", "CI15", "UL15", "HT15", "QE15")
control <- c("RAS15", "QU15", "TZ15", "WY15", "TJ15", "OL15", "WJ15", "ME16", "SI15")

# animals by sex
female <- c("TJ15", "SI15", "RAS15", "QU15", "QE15", "QH15", "RJS15", "QS15")
male <- c("WJ15", "WY15", "OL15", "TZ15", "ME16", "CI15", "UL15", "PS15", "HT15", "UU15")

### Functions

# function to add treatment and sex to animal-level data
addgroups <- function(x) {x %>% mutate(group = case_when(
  case %in% control ~ "control",
  case %in% anesthesia ~ "anesthesia"),
  sex = case_when(case %in% female ~ "female",
                  case %in% male ~ "male"))}

## CA1 Synapse Area Analyses

CA1_areas <- readRDS(here("data", "CA1_areas.rds"))

CA1_area_means <- CA1_areas %>% gather(case, area) %>% na.omit() %>% 
  group_by(case) %>% summarize(mean_area = mean(area)) %>% addgroups()

# find anesthesia and control group means
group_mean_area <- CA1_area_means %>% group_by(group) %>% summarize(mean = mean(mean_area), sd = sd(mean_area))
a_mean <- group_mean_area[2][[1]][1]
c_mean <- group_mean_area[2][[1]][2]

# calculate percent decrease in mean synapse size
(c_mean-a_mean)/c_mean*100 

CA1_areas2 <- CA1_areas %>% gather(case, area) %>% na.omit() %>%
  addgroups()

# fit a linear mixed model to the CA1 synapse areas and check group differences with anova 
# accounts for random variance of individual animals
CA1_area_lmer <- lmer(area ~ group*sex + (1|case), data = CA1_areas2)
anova(CA1_area_lmer)

# assign group factor levels to CA1 synapse area means
CA1_area_means$group <- factor(CA1_area_means$group, levels=c("control", "anesthesia"))
group_mean_area$group <- factor(group_mean_area$group, levels=c("control", "anesthesia"))

# sort the CA1 synapse areas into 5% bins, add a column to indicate quantile number, compile into columns
CA1quantiles <- data.frame(apply(CA1_areas, 2, twenty_quantiles))
CA1quantiles$quantile <- seq(1:20)
CA1quantiles <- gather(CA1quantiles, case, area, -quantile)

# add animal treatment & sex to CA1 synapse area means
CA1quantiles <- addgroups(CA1quantiles)

# assign treatments, sex as factor levels to CA1 synapse area quantiles
CA1quantiles$group <- factor(CA1quantiles$group, levels=c("control", "anesthesia"))
CA1quantiles$sex <- factor(CA1quantiles$sex, levels = c("female", "male"))

# analyze treatment & sex diffs across CA1 synapse area quantiles with linear mixed model run through anova
CA1_area_quantile_lmer <- lmer(area ~ group*sex*factor(quantile) + (1|case), data = CA1quantiles)
anova(CA1_area_quantile_lmer)
summary(CA1_area_quantile_lmer)
#CA1quantiles %>% filter(quantile ==19) %>% aov(area ~ group*sex, data = .) %>% summary()

# posthoc pairwise comparisons
CA1_area_quantile_lmer %>% emmeans(pairwise ~ group | quantile)
CA1_area_quantile_lmer %>% emmeans(pairwise ~ group*sex | quantile)
# sig: quantile 19: CF-AF .0047, quantile 20: CF-AF <.0001,  CM-AM .0362

# summarize CA1 synapse areas for quantiles 19 & 20
CA1quantiles %>% filter(quantile == 19) %>% group_by(group, sex) %>% summarize(mean(area), sd(area))
CA1quantiles %>% filter(quantile == 20) %>% group_by(group, sex) %>% summarize(mean(area), sd(area))

## CA1 Synapse Density Analyses

# import CA1 synapse density data & convert to columnar format
CA1density <- readRDS(here("data", "CA1density.rds"))
CA1density <- gather(CA1density, case, density, -pair) 

# add group assignments
CA1density <- addgroups(CA1density)

# analyze treatment & sex diffs for CA1 synapse density with linear mixed model run through anova
CA1_density_lmer <- lmer(density ~ group*sex + (1|case), data = CA1density)
anova(CA1_density_lmer)
CA1density %>% group_by(group, sex) %>% summarize(mean(density))

# create table of CA1 synapse density means per animal with treatment and sex
CA1densmeans <- CA1density %>% group_by(case) %>% summarize(mean_dens = mean(density)) %>%
  addgroups() 

# calculate synapse density means by sex and the percent difference between groups
CA1male <- CA1densmeans %>% filter(sex == "male") %>% summarize(mean(mean_dens), sd(mean_dens))
CA1female <- CA1densmeans %>% filter(sex == "female") %>% summarize(mean(mean_dens), sd(mean_dens))
(CA1male-CA1female)/CA1female *100

# test homogeneity of variances for animal mean CA1 synapse densities by treatment
bartlett.test(mean_dens ~ group, data = CA1densmeans) 

### CA1 Vesicle Docking Analysis 

# import CA1 vesicle counts 
CA1vesicles <- here("data", "CA1_vesicles_labeled.xlsx")

# list all sheets, set column names, read excel sheets into R
CA1vesicle_data_list <- CA1vesicles %>%
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = CA1vesicles)

# specify column names for CA1 vesicle data  
CA1bins <- CA1vesicle_data_list %>% map(~setNames(.x, c("case", "synapse", "vesicle", "bin", "label"))) %>%
  # fill in the case and synapse columns
  map(~fill(.x, case)) %>% map(~fill(.x, synapse)) %>%
  # sort for docking and predocking vesicles 
  map(~filter(.x,label == "0-30 from Pre" | label == "30-60 from Pre")) %>%
  # group by docking and predocking vesicles per synapse per case and count
  map_df(~.x %>% group_by(case,synapse, label) %>% summarize(n=n()) %>% 
           # add columns describing docking and predocking bin labels 
           mutate(bin=case_when(label == "0-30 from Pre" ~ "bin030", label == "30-60 from Pre" ~ "bin3060")) %>%
           # gather synapse, bin, and n columns and redistribute bins per synapse per case into columns with n as variables  
           select(synapse, bin, n) %>% spread(bin, n) %>% 
           # replace NAs with 0 
           replace_na(list(bin030 = 0, bin3060 = 0)) %>%
           # add a ratio column of docking to total vesicles 
           mutate(ratio = 100*(bin030/(bin030+bin3060))) %>%
           separate(case, c("case", NA, NA, NA, NA), "-")) 

# add treatment and sex
CA1bins <- addgroups(CA1bins)

# analyze CA1 vesicle docking ratio distribution by treatment & sex with linear mixed model run through anova
CA1_vesicle_lmer <- lmer(ratio ~ group*sex + (1|case), data = CA1bins)
anova(CA1_vesicle_lmer)

## CA1 Presynaptic Bouton Mitochondria Analyses

# import CA1 mito data
CA1mito <- readRDS(here("data", "CA1mito.rds"))

# analyze CA1 presynaptic bouton mitochondria density by treatment & sex with linear mixed model run through anova
CA1_mito_lmer <- lmer(total_mito_per_um3 ~ group*sex + (1|case), data = CA1mito)
anova(CA1_mito_lmer)

# data for plot: CA1 mean mito density with treatment averages as bars and animals as points 
CA1_mitodensity <- CA1mito %>% group_by(case,group) %>% 
  summarize(mean_mito = mean(total_mito_per_um3)) %>%
  mutate(group = factor(group, levels=c("control", "anesthesia"))) %>% ungroup() 


# analyze CA1 straight mito numbers with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_str_random <- glmer(num_straight ~ (1|case), data = CA1mito, family = poisson())
CA1_str_sex <- glmer(num_straight ~ sex + (1|case), data = CA1mito, family = poisson())
CA1_str_both <- glmer(num_straight ~ group + sex + (1|case), data = CA1mito, family = poisson())
CA1_str_full <- glmer(num_straight ~ group*sex + (1|case), data = CA1mito, family = poisson()) 

#check fit of full model with poisson distribution
CA1_str_full %>% model_fit()

# summarize full model 
summary(CA1_str_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_str_random, CA1_str_sex) # main effect of sex 1 df
anova(CA1_str_sex, CA1_str_both) # main effect of group 1 df (controlling for sex)
anova(CA1_str_full, CA1_str_both) # group x sex interaction 1 df


# analyze CA1 curved mito numbers with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_cur_random <- glmer(num_curved ~ (1|case), data = CA1mito, family = poisson())
CA1_cur_sex <- glmer(num_curved ~ sex + (1|case), data = CA1mito, family = poisson())
CA1_cur_both <- glmer(num_curved ~ group + sex + (1|case), data = CA1mito, family = poisson())
CA1_cur_full <- glmer(num_curved ~ group*sex + (1|case), data = CA1mito, family = poisson()) 

#check fit of full model with poisson distribution
CA1_cur_full %>% model_fit()

# summarize full model 
summary(CA1_cur_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_cur_random, CA1_cur_sex) # main effect of sex 1 df
anova(CA1_cur_sex, CA1_cur_both) # main effect of group 1 df (controlling for sex)
anova(CA1_cur_full, CA1_cur_both) # group x sex interaction 1 df

# summarize significant group numbers for CA1 curved mito 
CA1mito %>% group_by(group) %>% summarize(mean(num_curved), sd(num_curved))


# analyze CA1 curved mito numbers with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_tor_random <- glmer(num_toroid ~ (1 | case), family = "poisson", data = CA1mito, 
                        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
CA1_tor_sex <- glmer(num_toroid ~ sex + (1 | case), family = "poisson", data = CA1mito, 
                     control=glmerControl(optimizer="nlminbwrap", optCtrl=list(maxfun=100000)))
CA1_tor_both <- glmer(num_toroid ~ group + sex + (1|case), data = CA1mito, family = poisson())
CA1_tor_full <- glmer(num_toroid ~ group*sex + (1|case), data = CA1mito, family = poisson()) 

#check fit of full model with poisson distribution
CA1_tor_full %>% model_fit()

# summarize full model 
summary(CA1_tor_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_tor_random, CA1_tor_sex) # main effect of sex 1 df
anova(CA1_tor_sex, CA1_tor_both) # main effect of group 1 df (controlling for sex)
anova(CA1_tor_full, CA1_tor_both) # group x sex interaction 1 df

# find CA1 mito shape count means by treatment group 
CA1_mitoshape <- CA1mito %>% group_by(case, group, sex) %>% 
  summarize(Straight=mean(num_straight),Curved=mean(num_curved),Donut=mean(num_toroid)) %>%
  pivot_longer(Straight:Donut, values_to = "count", names_to = "morphology") %>% 
  mutate(morphology = factor(morphology, levels=c("Straight", "Curved", "Donut")),
         group = factor(group, levels=c("control", "anesthesia"))) 


# calculate and add columns for counts and percent of boutons with 0 mito
CA1mito2 <- CA1mito %>% mutate(mito_0_bouton = total_boutons -mito_1_bouton -mito_2_boutons- mito_3plus_boutons, 
                               pct_0_bouton = mito_0_bouton/total_boutons)

# add bouton mito frequency and other variables as factors
CA1mito_long <- CA1mito2 %>% pivot_longer(c(mito_1_bouton:mito_3plus_boutons, mito_0_bouton), 
                                          values_to = "number_per_bouton", names_to = "bouton_frequency") %>%
  mutate(group = factor(group, levels=c("control", "anesthesia")),
         sex = factor(sex, levels = c("female", "male")), 
         bouton_frequency = factor(bouton_frequency, 
                                   levels = c("mito_0_bouton", 
                                              "mito_1_bouton", 
                                              "mito_2_boutons", 
                                              "mito_3plus_boutons")))

# check Gaussian fit of full model
gaussmod <- lmer(number_per_bouton ~ group*sex*bouton_frequency + (1|case), data = CA1mito_long)
gaussmod %>% model_fit() # it is NOT HAPPY

# analyze CA1 mito frequency per bouton with general linear mixed model with poisson distribution
# analyze group, sex, bouton frequency, and full interaction effects 
CA1_boutons_1 <- glmer(number_per_bouton ~ (1|case), data = CA1mito_long, family = poisson())
CA1_boutons_2 <- glmer(number_per_bouton ~ sex + (1|case), data = CA1mito_long, family = poisson())
CA1_boutons_3 <- glmer(number_per_bouton ~ group + sex + (1|case), data = CA1mito_long, family = poisson())
CA1_boutons_4 <- glmer(number_per_bouton ~ group*sex + (1|case), data = CA1mito_long, family = poisson())
CA1_boutons_5 <- glmer(number_per_bouton ~ bouton_frequency + group*sex + (1|case), data = CA1mito_long, family = poisson())
CA1_boutons_6 <- glmer(number_per_bouton ~ sex:bouton_frequency + bouton_frequency + group*sex + (1|case), data = CA1mito_long, family = poisson())
CA1_boutons_7 <- glmer(number_per_bouton ~ group:bouton_frequency + sex:bouton_frequency + bouton_frequency + group*sex + (1|case), data = CA1mito_long, family = poisson())
CA1_boutons_8 <- glmer(number_per_bouton ~ group*sex*bouton_frequency + (1|case), data = CA1mito_long, family = poisson())

#check fit of full model with poisson distribution
CA1_boutons_8 %>% model_fit()
testOutliers(CA1_boutons_8, type = 'bootstrap')

# summarize full fit model
CA1_boutons_8 %>% summary()

# tabulate group, sex, bouton frequency, and interaction contributions to effects via chi squared tests
anova(CA1_boutons_1, CA1_boutons_2) # sex effects
anova(CA1_boutons_2, CA1_boutons_3) # group effects
anova(CA1_boutons_3, CA1_boutons_4) # group:sex interaction effects
anova(CA1_boutons_4, CA1_boutons_5) # bouton frequency effects
anova(CA1_boutons_5, CA1_boutons_6) # sex:bouton frequency interaction effects
anova(CA1_boutons_6, CA1_boutons_7) # group:bouton frequency interaction effects
anova(CA1_boutons_7, CA1_boutons_8) # group:sex:bouton frequency interaction effects

# summarize CA1 mito / bouton by group, sex, type 
CA1mito_long %>% group_by(bouton_frequency, group, sex) %>% summarize(mean(number_per_bouton), sd(number_per_bouton))


# analyze CA1 numbers of 1 mito/bouton with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_0_random <- glmer(mito_0_bouton ~ (1|case), data = CA1mito2, family = poisson())
CA1_0_sex <- glmer(mito_0_bouton ~ sex + (1|case), data = CA1mito2, family = poisson())
CA1_0_both <- glmer(mito_0_bouton ~ group + sex + (1|case), data = CA1mito2, family = poisson())
CA1_0_full <- glmer(mito_0_bouton ~ group*sex + (1|case), data = CA1mito2, family = poisson()) 

#check fit of full model with poisson distribution
CA1_0_full %>% model_fit()

# summarize full fit model
summary(CA1_0_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_0_random, CA1_0_sex) # main effect of sex 1 df
anova(CA1_0_sex, CA1_0_both) # main effect of group 1 df (controlling for sex)
anova(CA1_0_full, CA1_0_both) # group x sex interaction 1 df

# summarize near significant group numbers for CA1 0/bouton 
CA1mito2 %>% group_by(group,sex) %>% summarize(mean(mito_0_bouton), sd(mito_0_bouton))


# analyze CA1 numbers of 1 mito/bouton with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_1_random <- glmer(mito_1_bouton ~ (1|case), data = CA1mito, family = poisson())
CA1_1_sex <- glmer(mito_1_bouton ~ sex + (1|case), data = CA1mito, family = poisson())
CA1_1_both <- glmer(mito_1_bouton ~ group + sex + (1|case), data = CA1mito, family = poisson())
CA1_1_full <- glmer(mito_1_bouton ~ group*sex + (1|case), data = CA1mito, family = poisson()) 

#check fit of full model with poisson distribution
CA1_1_full %>% model_fit()

# summarize full fit model
summary(CA1_1_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_1_random, CA1_1_sex) # main effect of sex 1 df
anova(CA1_1_sex, CA1_1_both) # main effect of group 1 df (controlling for sex)
anova(CA1_1_full, CA1_1_both) # group x sex interaction 1 df


# analyze CA1 numbers of 2 mito/bouton with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_2_random <- glmer(mito_2_boutons ~ (1|case), data = CA1mito, family = poisson())
CA1_2_sex <- glmer(mito_2_boutons ~ sex + (1|case), data = CA1mito, family = poisson())
CA1_2_both <- glmer(mito_2_boutons ~ group + sex + (1|case), data = CA1mito, family = poisson())
CA1_2_full <- glmer(mito_2_boutons ~ group*sex + (1|case), data = CA1mito, family = poisson()) 

#check fit of full model with poisson distribution
CA1_2_full %>% model_fit()

# summarize full fit model
summary(CA1_2_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_2_random, CA1_2_sex) # main effect of sex 1 df
anova(CA1_2_sex, CA1_2_both) # main effect of group 1 df (controlling for sex)
anova(CA1_2_full, CA1_2_both) # group x sex interaction 1 df


# analyze CA1 numbers of 3+ mito/bouton with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
# singularity appears due to low variance of random effects
CA1_3_random <- glmer(mito_3plus_boutons ~ (1|case), data = CA1mito, family = poisson())
CA1_3_sex <- glmer(mito_3plus_boutons ~ sex + (1|case), data = CA1mito, family = poisson())
CA1_3_both <- glmer(mito_3plus_boutons ~ group + sex + (1|case), data = CA1mito, family = poisson())
CA1_3_full <- glmer(mito_3plus_boutons ~ group*sex + (1|case), data = CA1mito, family = poisson()) 

#check fit of full model with poisson distribution
CA1_3_full %>% model_fit()

# summarize full fit model
summary(CA1_3_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_3_random, CA1_3_sex) # main effect of sex 1 df
anova(CA1_3_sex, CA1_3_both) # main effect of group 1 df (controlling for sex)
anova(CA1_3_full, CA1_3_both) # group x sex interaction 1 df

# summarize near significant group numbers for CA1 3+/bouton 
CA1mito %>% group_by(sex) %>% summarize(mean(mito_3plus_boutons), sd(mito_3plus_boutons))


##### DLPFC analysis

# cases were recoded for DLPFC analyses to maintain blinding to group and sex
# these functions apply case labels that were used for CA1 data analysis
toCA1_case <- function(x) {x %>% 
    mutate(case = case_when(case == "case_1" ~ "TJ15",
                            case == "case_2" ~ "QU15",
                            case == "case_3" ~ "UU15",
                            case == "case_4" ~ "PS15",
                            case == "case_5" ~ "WJ15",
                            case == "case_6" ~ "RAS15",
                            case == "case_7" ~ "QS15",
                            case == "case_8" ~ "RJS15",
                            case == "case_11" ~ "TZ15",
                            case == "case_12" ~ "CI15",
                            case == "case_13" ~ "QH15",
                            case == "case_14" ~ "QE15",
                            case == "case_15" ~ "HT15",
                            case == "case_16" ~ "SI15",
                            case == "case_17" ~ "ME16",
                            case == "case_18" ~ "WY15",
                            case == "case_19" ~ "OL15",
                            case == "case_20" ~ "UL15"))}

toCA1 <- function(x) {x %>% mutate(case = case_when(case == "1" ~ "TJ15",
                                                    case == "2" ~ "QU15",
                                                    case == "3" ~ "UU15",
                                                    case == "4" ~ "PS15",
                                                    case == "5" ~ "WJ15",
                                                    case == "6" ~ "RAS15",
                                                    case == "7" ~ "QS15",
                                                    case == "8" ~ "RJS15",
                                                    case == "11" ~ "TZ15",
                                                    case == "12" ~ "CI15",
                                                    case == "13" ~ "QH15",
                                                    case == "14" ~ "QE15",
                                                    case == "15" ~ "HT15",
                                                    case == "16" ~ "SI15",
                                                    case == "17" ~ "ME16",
                                                    case == "18" ~ "WY15",
                                                    case == "19" ~ "OL15",
                                                    case == "20" ~ "UL15"))}

## dlPFC Synapse Area Analyses

# import dlPFC synapse area data
dl_areas <- readRDS(here("data", "dl_areas.rds"))

# compile dlPFC synapse areas, change cases to CA1 labels, add treatment and sex
dl_areas2 <- dl_areas %>% gather(case, synapse_area) %>% na.omit() %>%
  toCA1_case() %>% addgroups()

# analyze treatment & sex diffs across whole dlPFC synapse area with linear mixed model run through anova
dlpfc_area_lmer <- lmer(synapse_area ~ group*sex + (1|case), data = dl_areas2)
anova(dlpfc_area_lmer)

# sort the dlPFC synapse areas into 5% bins, add a column to indicate quantile number, compile into columns
dl_quantiles <- data.frame(apply(dl_areas, 2, twenty_quantiles))
dl_quantiles$quantile <- seq(1:20)
dl_quantiles <- gather(dl_quantiles, case, area, -quantile)

# change cases to CA1 labels, add treatment and sex
dl_quantiles2 <- dl_quantiles %>% na.omit() %>%
  toCA1_case() %>% addgroups()

# assign treatment groups as factor levels for dlPFC synapse areas
dl_quantiles2$group <- factor(dl_quantiles2$group, levels=c("control", "anesthesia"))

# analyze treatment & sex diffs by dlPFC synapse area quantiles with linear mixed model run through anova
dlpfc_area_quantile_lmer <- lmer(area ~ group*sex*factor(quantile) + (1|case), data = dl_quantiles2)
anova(dlpfc_area_quantile_lmer)
summary(dlpfc_area_quantile_lmer)

# pairwise comparisons
dlpfc_area_quantile_lmer %>% emmeans(pairwise ~ group | quantile)

#only quantile 20: CF-AF (<0.0001) and CM-AM (.0049) sig, CF-CM and AF-AM same

# summarize dlPFC synapse areas for quantile 20
dl_quantiles2 %>% filter(quantile == 20) %>% group_by(group) %>% summarize(mean(area), sd(area))

# calculate quantile 20 synapse area means by treatment and the percent difference between groups
dl_anes <- dl_quantiles2 %>% filter(quantile ==20, group == "anesthesia") %>% summarize(mean(area)) %>% print()
dl_control <- dl_quantiles2 %>% filter(quantile == 20, group == "control") %>% summarize(mean(area)) %>% print()
(dl_control-dl_anes)/dl_control *100

## dlPFC Synapse Density Analyses 

# import dlPFC synapse densities
dl_density <- readRDS(here("data", "dl_density.rds"))

# compile dlPFC synapse densities, change cases to CA1 labels, add treatment and sex, remove suspect cases 
dl_density2 <- gather(dl_density, case, density) %>% na.omit() %>% 
  toCA1() %>% addgroups()

# analyze dlPFC synapse density by treatment & sex with linear mixed model run through anova
dlpfc_density_lmer <- lmer(density ~ group*sex + (1|case), data = dl_density2)
anova(dlpfc_density_lmer)

# create table of dlPFC synapse density means per animal with treatment and sex
dlPFCdensmeans <- dl_density2 %>% group_by(case) %>% summarize(mean_dens = mean(density)) %>%
  addgroups() 
# test homogeneity of variances for animal mean CA1 synapse densities by treatment
bartlett.test(mean_dens ~ group, data = dlPFCdensmeans) 

### dlPFC Vesicle Docking Analysis 

# import dlPFC vesicle counts 
dlpfcvesicles <- here("data", "dlPFC Vesicle Counts.xlsx")

# list all sheets, set column names, read excel sheets into R
dlpfcvesicle_data_list <- dlpfcvesicles %>%
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = dlpfcvesicles)

# specify column names for dlPFC vesicle data  
dlpfcbins <- dlpfcvesicle_data_list %>% map(~setNames(.x, c("case", "synapse", "vesicle", "bin", "label"))) %>%
  # fill in the case and synapse columns
  map(~fill(.x, case)) %>% map(~fill(.x, synapse)) %>%
  # sort for docking and predocking vesicles 
  map(~filter(.x,label == "0-30 from Pre" | label == "30-60 from Pre")) %>%
  # group by docking and predocking vesicles per synapse per case and count
  map_df(~.x %>% group_by(case,synapse, label) %>% summarize(n=n()) %>% 
           # add columns describing docking and predocking bin labels 
           mutate(bin=case_when(label == "0-30 from Pre" ~ "bin030", label == "30-60 from Pre" ~ "bin3060")) %>%
           # gather synapse, bin, and n columns and redistribute bins per synapse per case into columns with n as variables  
           select(synapse, bin, n) %>% spread(bin, n) %>% 
           # replace NAs with 0 
           replace_na(list(bin030 = 0, bin3060 = 0)) %>%
           # add a ratio column of docking to total vesicles 
           mutate(ratio = 100*(bin030/(bin030+bin3060))) %>%
           separate(case, c("case", NA, "stack", NA), "-")) 

# change cases to CA1 labels, add treatment and sex
dlpfcbins <- toCA1(dlpfcbins) %>% addgroups()

# analyze dlPFC vesicle docking ratio distribution by treatment & sex with linear mixed model run through anova
dlpfc_vesicle_lmer <- lmer(ratio ~ group*sex + (1|case), data = dlpfcbins)
anova(dlpfc_vesicle_lmer)

## dlPFC Mitochondria Analyses 

# import dlPFC mitochondria data
dl_mito <- readRDS(here("data", "dl_mito.rds"))

# change cases to CA1 labels, add treatment and sex
dl_mito2 <- toCA1(dl_mito) %>% addgroups()

# analyze dlPFC presynaptic bouton mitochondria density by treatment & sex with linear mixed model run through anova
dlpfc_mito_lmer <- lmer(total_mito_per_um3 ~ group*sex + (1|case), data = dl_mito2)
anova(dlpfc_mito_lmer)

# organize data to plot dlPFC presynaptic bouton mean mitochondrial density by treatment
dl_mitodensity <- dl_mito2 %>% group_by(case,group) %>% summarize(mean_mito = mean(total_mito_per_um3)) %>%
  mutate(group = factor(group, levels=c("control", "anesthesia"))) %>% ungroup()


# analyze dlPFC numbers of straight mito with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_str_random <- glmer(num_straight ~ (1|case), data = dl_mito2, family = poisson())
dl_str_sex <- glmer(num_straight ~ sex + (1|case), data = dl_mito2, family = poisson())
dl_str_both <- glmer(num_straight ~ group + sex + (1|case), data = dl_mito2, family = poisson())
dl_str_full <- glmer(num_straight ~ group*sex + (1|case), data = dl_mito2, family = poisson()) 

#check fit of full model with poisson distribution
dl_str_full %>% model_fit()

# summarize full fit model
summary(dl_str_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_str_random, dl_str_sex) # main effect of sex 1 df
anova(dl_str_sex, dl_str_both) # main effect of group 1 df (controlling for sex)
anova(dl_str_full, dl_str_both) # group x sex interaction 1 df


# analyze dlPFC numbers of curved mito with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_cur_random <- glmer(num_curved ~ (1|case), data = dl_mito2, family = poisson())
dl_cur_sex <- glmer(num_curved ~ sex + (1|case), data = dl_mito2, family = poisson())
dl_cur_both <- glmer(num_curved ~ group + sex + (1|case), data = dl_mito2, family = poisson())
dl_cur_full <- glmer(num_curved ~ group*sex + (1|case), data = dl_mito2, family = poisson()) 

#check fit of full model with poisson distribution
dl_cur_full %>% model_fit()

# summarize full fit model
summary(dl_cur_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_cur_random, dl_cur_sex) # main effect of sex 1 df
anova(dl_cur_sex, dl_cur_both) # main effect of group 1 df (controlling for sex)
anova(dl_cur_full, dl_cur_both) # group x sex interaction 1 df


# analyze dlPFC numbers of toroidal mito with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_tor_random <- glmer(num_toroid ~ (1|case), data = dl_mito2, family = poisson())
dl_tor_sex <- glmer(num_toroid ~ sex + (1|case), data = dl_mito2, family = poisson())
dl_tor_both <- glmer(num_toroid ~ group + sex + (1|case), data = dl_mito2, family = poisson())
dl_tor_full <- glmer(num_toroid ~ group*sex + (1|case), data = dl_mito2, family = poisson()) 

#check fit of full model with poisson distribution
dl_tor_full %>% model_fit()

# summarize full fit model
summary(dl_tor_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_tor_random, dl_tor_sex) # main effect of sex 1 df
anova(dl_tor_sex, dl_tor_both) # main effect of group 1 df (controlling for sex)
anova(dl_tor_full, dl_tor_both) # group x sex interaction 1 df


# summarize CA1 mito by shape per treatment and sex
dl_mito2 %>% group_by(group, sex) %>% 
  summarize(mean_straight = mean(pct_straight), mean_curved = mean(pct_curved), mean_toroid = mean(pct_toroid))

# organize data to plot dlPFC mito shape count means by treatment group 
dl_mitoshape <- dl_mito2 %>% group_by(case, group, sex) %>% 
  summarize(Straight=mean(num_straight),Curved=mean(num_curved),Donut=mean(num_toroid)) %>%
  pivot_longer(Straight:Donut, values_to = "count", names_to = "morphology") %>% 
  mutate(morphology = factor(morphology, levels=c("Straight", "Curved", "Donut")),
         group = factor(group, levels=c("control", "anesthesia"))) 


# calculate and add columns for counts and percent of boutons with 0 mito
dl_mito3 <- dl_mito2 %>% mutate(mito_0_bouton = total_boutons -mito_1_bouton -mito_2_boutons- mito_3plus_boutons, 
                                pct_0_bouton = mito_0_bouton/total_boutons)

# add bouton mito frequency and other variables as factors
dl_mito_long <- dl_mito3 %>% pivot_longer(c(mito_1_bouton:mito_3plus_boutons, mito_0_bouton), 
                                          values_to = "number_per_bouton", names_to = "bouton_frequency") %>%
  mutate(group = factor(group, levels=c("control", "anesthesia")),
         sex = factor(sex, levels = c("female", "male")), 
         bouton_frequency = factor(bouton_frequency, 
                                   levels = c("mito_0_bouton", 
                                              "mito_1_bouton", 
                                              "mito_2_boutons", 
                                              "mito_3plus_boutons")))

# check Gaussian fit of full model
gaussmod <- lmer(number_per_bouton ~ group*sex*bouton_frequency + (1|case), data = dl_mito_long)
gaussmod %>% model_fit() # it is NOT HAPPY

# analyze DLPFC mito frequency per bouton with general linear mixed model with poisson distribution
# analyze group, sex, bouton frequency, and full interaction effects 
dl_boutons_1 <- glmer(number_per_bouton ~ (1|case), data = dl_mito_long, family = poisson())
dl_boutons_2 <- glmer(number_per_bouton ~ sex + (1|case), data = dl_mito_long, family = poisson())
dl_boutons_3 <- glmer(number_per_bouton ~ group + sex + (1|case), data = dl_mito_long, family = poisson())
dl_boutons_4 <- glmer(number_per_bouton ~ group*sex + (1|case), data = dl_mito_long, family = poisson())
dl_boutons_5 <- glmer(number_per_bouton ~ bouton_frequency + group*sex + (1|case), data = dl_mito_long, family = poisson())
dl_boutons_6 <- glmer(number_per_bouton ~ sex:bouton_frequency + bouton_frequency + group*sex + (1|case), data = dl_mito_long, family = poisson())
dl_boutons_7 <- glmer(number_per_bouton ~ group:bouton_frequency + sex:bouton_frequency + bouton_frequency + group*sex + (1|case), data = dl_mito_long, family = poisson())
dl_boutons_8 <- glmer(number_per_bouton ~ group*sex*bouton_frequency + (1|case), data = dl_mito_long, family = poisson())

#check fit of full model with poisson distribution
dl_boutons_8 %>% model_fit()
testOutliers(dl_boutons_8, type = 'bootstrap')

# summarize full fit model
dl_boutons_8 %>% summary()

# tabulate group, sex, bouton frequency, and interaction contributions to effects via chi squared tests
anova(dl_boutons_1, dl_boutons_2) # sex effects
anova(dl_boutons_2, dl_boutons_3) # group effects
anova(dl_boutons_3, dl_boutons_4) # group:sex interaction effects
anova(dl_boutons_4, dl_boutons_5) # bouton frequency effects
anova(dl_boutons_5, dl_boutons_6) # sex:bouton frequency interaction effects
anova(dl_boutons_6, dl_boutons_7) # group:bouton frequency interaction effects
anova(dl_boutons_7, dl_boutons_8) # group:sex:bouton frequency interaction effects

# summarize CA1 mito / bouton by group, sex, type 
dl_mito_long %>% group_by(bouton_frequency, group, sex) %>% summarize(mean(number_per_bouton), sd(number_per_bouton))


# analyze dlPFC numbers of boutons without mito with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_0_random <- glmer(mito_0_bouton ~ (1|case), data = dl_mito3, family = poisson())
dl_0_sex <- glmer(mito_0_bouton ~ sex + (1|case), data = dl_mito3, family = poisson())
dl_0_both <- glmer(mito_0_bouton ~ group + sex + (1|case), data = dl_mito3, family = poisson())
dl_0_full <- glmer(mito_0_bouton ~ group*sex + (1|case), data = dl_mito3, family = poisson()) 

#check fit of full model with poisson distribution
dl_0_full %>% model_fit()

# summarize full fit model
summary(dl_0_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_0_random, dl_0_sex) # main effect of sex 1 df
anova(dl_0_sex, dl_0_both) # main effect of group 1 df (controlling for sex)
anova(dl_0_full, dl_0_both) # group x sex interaction 1 df

# summarize near significant group numbers for dlPFC/bouton 
dl_mito3 %>% group_by(group,sex) %>% summarize(mean(mito_0_bouton), sd(mito_0_bouton))

# pairwise comparisons for 0 mito/bouton
emmeans(dl_0_full, list(pairwise ~ group*sex))


# analyze dlPFC numbers of 1 mito/bouton with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_1_random <- glmer(mito_1_bouton ~ (1|case), data = dl_mito2, family = poisson())
dl_1_sex <- glmer(mito_1_bouton ~ sex + (1|case), data = dl_mito2, family = poisson())
dl_1_both <- glmer(mito_1_bouton ~ group + sex + (1|case), data = dl_mito2, family = poisson())
dl_1_full <- glmer(mito_1_bouton ~ group*sex + (1|case), data = dl_mito2, family = poisson()) 

#check fit of full model with poisson distribution
dl_1_full %>% model_fit()

# summarize full fit model
summary(dl_1_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_1_random, dl_1_sex) # main effect of sex 1 df
anova(dl_1_sex, dl_1_both) # main effect of group 1 df (controlling for sex)
anova(dl_1_full, dl_1_both) # group x sex interaction 1 df

# calculate sex means for 1/bouton counts
dl_mito2 %>% group_by(group,sex) %>% summarize(mean(mito_1_bouton), sd(mito_1_bouton))


# analyze dlPFC numbers of 2 mito/bouton with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_2_random <- glmer(mito_2_boutons ~ (1|case), data = dl_mito2, family = poisson())
dl_2_sex <- glmer(mito_2_boutons ~ sex + (1|case), data = dl_mito2, family = poisson())
dl_2_both <- glmer(mito_2_boutons ~ group + sex + (1|case), data = dl_mito2, family = poisson())
dl_2_full <- glmer(mito_2_boutons ~ group*sex + (1|case), data = dl_mito2, family = poisson()) 

#check fit of full model with poisson distribution
dl_2_full %>% model_fit()

# summarize full fit model
summary(dl_2_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_2_random, dl_2_sex) # main effect of sex 1 df
anova(dl_2_sex, dl_2_both) # main effect of group 1 df (controlling for sex)
anova(dl_2_full, dl_2_both) # group x sex interaction 1 df


# analyze dlPFC numbers of 3+ mito/bouton with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_3_random <- glmer(mito_3plus_boutons ~ (1|case), data = dl_mito2, family = poisson())
dl_3_sex <- glmer(mito_3plus_boutons ~ sex + (1|case), data = dl_mito2, family = poisson())
dl_3_both <- glmer(mito_3plus_boutons ~ group + sex + (1|case), data = dl_mito2, family = poisson())
dl_3_full <- glmer(mito_3plus_boutons ~ group*sex + (1|case), data = dl_mito2, family = poisson()) 

#check fit of full model with poisson distribution
dl_3_full %>% model_fit()

# summarize full fit model
summary(dl_3_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_3_random, dl_3_sex) # main effect of sex 1 df
anova(dl_3_sex, dl_3_both) # main effect of group 1 df (controlling for sex)
anova(dl_3_full, dl_3_both) # group x sex interaction 1 df


# organize data to plot the distribution of dlPFC mito presynaptic bouton frequencies by number
dl_mitofreq <- dl_mito3 %>% group_by(group,sex,case) %>% 
  summarize(Zero = mean(mito_0_bouton), One = mean(mito_1_bouton), Two = mean(mito_2_boutons), "Three or more" = mean(mito_3plus_boutons)) %>% 
  pivot_longer(Zero:"Three or more", values_to = "number", names_to = "counts") %>% 
  mutate(counts = factor(counts, levels=c("Zero", "One", "Two", "Three or more")),
         group = factor(group, levels=c("control", "anesthesia")),
         sex = factor(sex, levels = c("female", "male"))) 



###### 

## CA1 SYNAPSE TYPE COUNTS

## Use Poisson models for everything. 
## Better fit than Gaussian and more appropriate given small counts
## Negative binomial models collapse to Poisson.

# import CA1 synapse type counts
CA1counts <- readRDS(here("data", "CA1counts.rds"))

# add treatment and sex
CA1counts <- CA1counts %>% addgroups()

# calculate case means per synapse type 
CA1countmeans <- CA1counts %>% group_by(case, group, sex) %>% 
  summarize(perf_spin = mean(perf_spin), 
            perf_dend = mean(perf_dend), 
            nonperf_spin = mean(nonperf_spin), 
            nonperf_dend = mean(nonperf_dend))

# add group by sex column 
CA1countmeans <-  CA1countmeans %>% mutate(combined = paste(group,sex)) %>%
  mutate(combined = factor(combined, levels = c("control female",
                                                "control male",
                                                "anesthesia female",
                                                "anesthesia male"))) 


## Is there an interaction of synapse type, group, and sex?

# add synapse type as a factor
CA1counts_long <- CA1counts %>% pivot_longer(perf_spin:nonperf_dend, values_to = "synapse_count", names_to = "synapse_type")

# check Gaussian fit of full model
gaussmod <- lmer(synapse_count ~ group*sex*synapse_type + (1|case), data = CA1counts_long)
gaussmod %>% model_fit() # it is NOT HAPPY

# analyze CA1 synapse counts by type with general linear mixed model with poisson distribution
# analyze group, sex, synapse type, and full interaction effects 
CA1_counts_1 <- glmer(synapse_count ~ (1|case), data = CA1counts_long, family = poisson())
CA1_counts_2 <- glmer(synapse_count ~ sex + (1|case), data = CA1counts_long, family = poisson())
CA1_counts_3 <- glmer(synapse_count ~ group + sex + (1|case), data = CA1counts_long, family = poisson())
CA1_counts_4 <- glmer(synapse_count ~ group*sex + (1|case), data = CA1counts_long, family = poisson())
CA1_counts_5 <- glmer(synapse_count ~ synapse_type + group*sex + (1|case), data = CA1counts_long, family = poisson())
CA1_counts_6 <- glmer(synapse_count ~ sex:synapse_type + synapse_type + group*sex + (1|case), data = CA1counts_long, family = poisson())
CA1_counts_7 <- glmer(synapse_count ~ group:synapse_type + sex:synapse_type + synapse_type + group*sex + (1|case), data = CA1counts_long, family = poisson())
CA1_counts_8 <- glmer(synapse_count ~ group*sex*synapse_type + (1|case), data = CA1counts_long, family = poisson())

#check fit of full model with poisson distribution
CA1_counts_8 %>% model_fit()
testOutliers(CA1_counts_8, type = 'bootstrap')

# summarize full fit model
CA1_counts_8 %>% summary()

# tabulate group, sex, synapse type, and interaction contributions to effects via chi squared tests
anova(CA1_counts_1, CA1_counts_2) # sex effects
anova(CA1_counts_2, CA1_counts_3) # group effects
anova(CA1_counts_3, CA1_counts_4) # group:sex interaction effects
anova(CA1_counts_4, CA1_counts_5) # synapse type effects
anova(CA1_counts_5, CA1_counts_6) # sex:synapse type interaction effects
anova(CA1_counts_6, CA1_counts_7) # group:synapse type interaction effects
anova(CA1_counts_7, CA1_counts_8) # group:sex:synapse type interaction effects

# summarize CA1 synapse counts by group, sex, type 
CA1counts_long %>% group_by(synapse_type, group, sex) %>% summarize(mean(synapse_count), sd(synapse_count))

## What are the effects of anesthesia & sex on synapse types? 

# analyze CA1 perforated spinous synapse counts with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_ps_random <- glmer(perf_spin ~ (1|case), data = CA1counts, family = poisson())
CA1_ps_sex <- glmer(perf_spin ~ sex + (1|case), data = CA1counts, family = poisson())
CA1_ps_both <- glmer(perf_spin ~ group + sex + (1|case), data = CA1counts, family = poisson())
CA1_ps_full <- glmer(perf_spin ~group*sex + (1|case), data = CA1counts, family = poisson()) 

#check fit of full model with poisson distribution
CA1_ps_full %>% model_fit()

# summarize full fit model
summary(CA1_ps_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_ps_random, CA1_ps_sex) # main effect of sex 1 df
anova(CA1_ps_sex, CA1_ps_both) # main effect of group 1 df (controlling for sex)
anova(CA1_ps_full, CA1_ps_both) # group x sex interaction 1 df


# analyze CA1 perforated dendritic synapse counts with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_pd_random <- glmer(perf_dend ~ (1|case), data = CA1counts, family = poisson())
CA1_pd_sex <- glmer(perf_dend ~ sex + (1|case), data = CA1counts, family = poisson())
CA1_pd_both <- glmer(perf_dend ~ group + sex + (1|case), data = CA1counts, family = poisson())
CA1_pd_full <- glmer(perf_dend ~ group*sex + (1|case), data = CA1counts, family = poisson()) 

#check fit of full model with poisson distribution
CA1_pd_full %>% model_fit()

# summarize full fit model
summary(CA1_pd_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_pd_random, CA1_pd_sex) # main effect of sex 1 df
anova(CA1_pd_sex, CA1_pd_both) # main effect of group 1 df (controlling for sex)
anova(CA1_pd_full, CA1_pd_both) # group x sex interaction 1 df

# summarize CA1 perforated dendritic synapse counts by group
CA1counts %>% group_by(group) %>%summarize(mean(perf_dend), sd(perf_dend)) 

# analyze CA1 nonperforated spinous synapse counts with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_ns_random <- glmer(nonperf_spin ~ (1|case), data = CA1counts, family = poisson())
CA1_ns_sex <- glmer(nonperf_spin ~ sex + (1|case), data = CA1counts, family = poisson())
CA1_ns_both <- glmer(nonperf_spin ~ group + sex + (1|case), data = CA1counts, family = poisson())
CA1_ns_full <- glmer(nonperf_spin ~group*sex + (1|case), data = CA1counts, family = poisson())

#check fit of full model with poisson distribution
CA1_ns_full %>% model_fit()

# summarize full fit model
summary(CA1_ns_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_ns_random, CA1_ns_sex) # main effect of sex 1 df
anova(CA1_ns_sex, CA1_ns_both) # main effect of group 1 df (controlling for sex)
anova(CA1_ns_full, CA1_ns_both) # group x sex interaction 1 df


# analyze CA1 nonperforated dendritic synapse counts with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
CA1_nd_random <- glmer(nonperf_dend ~ (1|case), data = CA1counts, family = poisson())
CA1_nd_sex <- glmer(nonperf_dend ~ sex + (1|case), data = CA1counts, family = poisson())
CA1_nd_both <- glmer(nonperf_dend ~ group + sex + (1|case), data = CA1counts, family = poisson())
CA1_nd_full <- glmer(nonperf_dend ~group*sex + (1|case), data = CA1counts, family = poisson()) 

#check fit of full model with poisson distribution
CA1_nd_full %>% model_fit()

# summarize full fit model
summary(CA1_nd_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(CA1_nd_random, CA1_nd_sex) # main effect of sex 1 df
anova(CA1_nd_sex, CA1_nd_both) # main effect of group 1 df (controlling for sex)
anova(CA1_nd_full, CA1_nd_both) # group x sex interaction 1 df

# summarize CA1 nonperforated dendritic synapse counts by group
CA1counts %>% group_by(group) %>%summarize(mean(nonperf_dend), sd(nonperf_dend)) 


### CA1 SYNAPSE TYPE AREAS

# import CA1 synapse type areas file
CA1areas_file <- here("data", "CA1 Synapse Areas.xlsx")

# list all sheets, set column names, read excel sheets into R
CA1areas_data_list <- CA1areas_file %>%
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = CA1areas_file)

# create names column from sheet names
for (i in 1:length(CA1areas_data_list)) {CA1areas_data_list[[i]]$case <- names(CA1areas_data_list)[i]}

# Cyan: nonperforated, dendritic
# Pink: nonperforated, spinous
# X: perforated, spinous
# Yellow: perforated, dendritic 

# set column names, remove NAs from type, combine to a single df, extract case name
CA1synapse_areas <- 
  CA1areas_data_list %>% map(~setNames(.x, c("type", "note", "start", "end", "length", "area", "case"))) %>%
  map(~filter(.x, !is.na(type))) %>%
  bind_rows() %>%
  mutate(case = str_extract(case, "-([^-]+)") %>% str_remove("-")) %>%
  # classify areas as non/perforated and spinous/dendritic
  mutate(perf = case_when(str_detect(type, "cyan") ~ "nonperf",
                          str_detect(type, "pink") ~ "nonperf",
                          str_detect(type, "X") ~ "perf",
                          str_detect(type, "yellow") ~ "perf")) %>%
  mutate(spinous = case_when(str_detect(type, "cyan") ~ "dendritic",
                             str_detect(type, "pink") ~ "spinous",
                             str_detect(type, "X") ~ "spinous",
                             str_detect(type, "yellow") ~ "dendritic")) %>%
  # add treatment and sex, remove suspect cases, select for relevant columns 
  addgroups() %>% select(case,group,sex,perf,spinous,area)

# analyze CA1 synapse areas by treatment & sex & synapse types with linear mixed model run through anova
CA1_synapse_type_area_lmer <- CA1synapse_areas %>% lmer(area ~ group * sex * spinous * perf + (1 | case), data = .)
anova(CA1_synapse_type_area_lmer)

# analyze CA1 synapse areas by type with linear mixed models 
# analyze treatment, sex, and treatment by sex contributions to synapse areas effects by synapse type
overall_area <- lmer(area ~ group * sex * spinous * perf + (1 | case), data = CA1synapse_areas)
reduced_area <- lmer(area ~ group * sex + spinous + perf + spinous:perf + (1 | case), data = CA1synapse_areas)
reduced_area_2 <- lmer(area ~ group * sex + spinous + perf + spinous:perf +
                         sex:spinous + sex:perf + sex:spinous:perf + (1 | case), data = CA1synapse_areas)
reduced_area_3 <- lmer(area ~ group * sex + spinous + perf + spinous:perf +
                         sex:spinous + sex:perf + sex:spinous:perf + 
                         group:spinous + group:perf + group:spinous:perf +(1 | case), data = CA1synapse_areas)

# tabulate treatment, sex, synapse type, and interaction contributions to effects via chi squared tests
anova(reduced_area, reduced_area_2)   # sex effects: sex:spin, sex:perf, and sex:spin:perf interactions
anova(reduced_area_2, reduced_area_3) # group effects: group:spin, group:perf, and group:spin:perf interactions
anova(reduced_area_3, overall_area)   # group by sex effects: group:sex:spin ,group:sex:perf, and group:sex:spin:perf interactions

# analyze synapse areas per type by treatment and sex  
CA1_perf_spin_area_anova <- 
  CA1synapse_areas %>% filter(perf == "perf") %>% filter(spinous == "spinous") %>%
  lmer(area ~ group * sex + (1 | case), data = .) %>% anova()
CA1_perf_dend_area_anova <- 
  CA1synapse_areas %>% filter(perf == "perf") %>% filter(spinous == "dendritic") %>%
  lmer(area ~ group * sex + (1 | case), data = .) %>% anova()
CA1_nonperf_spin_area_anova <- 
  CA1synapse_areas %>% filter(perf == "nonperf") %>% filter(spinous == "spinous") %>%
  lmer(area ~ group * sex + (1 | case), data = .) %>% anova()
CA1_nonperf_dend_area_anova <-
  CA1synapse_areas %>% filter(perf == "nonperf") %>% filter(spinous == "dendritic") %>%
  lmer(area ~ group * sex + (1 | case), data = .) %>% anova()

CA1_perf_spin_area_anova
CA1_perf_dend_area_anova
CA1_nonperf_spin_area_anova
CA1_nonperf_dend_area_anova

# calculate perforated spinous area means by treatment and sex
CA1synapse_areas %>% group_by(group) %>% 
  filter(perf == "perf") %>% filter(spinous == "spinous") %>% 
  summarize(mean_area = mean(area), sd(area))
# calculate nonperforated spinous area means by treatment and sex
CA1synapse_areas %>% group_by(group) %>% 
  filter(perf == "nonperf") %>% filter(spinous == "spinous") %>% 
  summarize(mean_area = mean(area), sd(area))

## 20220809: are synapses abnormally large
## 2 SD criterion

CA1synapse_areas %>% group_by(group, spinous) %>% 
  filter(perf == "nonperf") %>% summarize(mean_area = mean(area), sd(area), n())
(CA1_abnorm <- CA1synapse_areas %>% 
    filter(perf == "nonperf") %>% filter(spinous == "spinous") %>% 
    summarize(mean_area = mean(area), sd_area = sd(area), n()))
CA1_large_synapses <- 
  CA1synapse_areas %>% filter(area > CA1_abnorm$mean_area+2*CA1_abnorm$sd_area) %>% group_by(group) %>% 
  filter(perf == "nonperf") %>% filter(spinous == "spinous")
CA1_large_synapses %>% summarize(mean_area = mean(area), sd(area), n())
# are large nonperforated spinous synapses larger?
t.test(area ~ group, data = CA1_large_synapses)
# are there more abnormally large nonperf spinous in control?
dat_ab <- data.frame(
  "not_large" = c((2646-171), (2652-135)),
  "abnorm_large" = c(171, 135),
  row.names = c("Control", "Anesthesia"),
  stringsAsFactors = FALSE
)
colnames(dat_ab) <- c("Not_abnorm", "Abnorm_large")

dat_ab
fisher.test(dat_ab)

# same thing for dendritic?
(CA1_abnorm_d <- CA1synapse_areas %>% 
    filter(perf == "nonperf") %>% filter(spinous == "dendritic") %>% 
    summarize(mean_area = mean(area), sd_area = sd(area), n()))
CA1_large_synapses_d <- 
  CA1synapse_areas %>% filter(area > CA1_abnorm_d$mean_area+2*CA1_abnorm_d$sd_area) %>% group_by(group) %>% 
  filter(perf == "nonperf") %>% filter(spinous == "dendritic")
CA1_large_synapses_d %>% summarize(mean_area = mean(area), sd(area), n())
# only 8 of each so no

# create graph labels
perf.labs <- c("Perforated", "Nonperforated")
names(perf.labs) <- c("perf", "nonperf")
spinous.labs <- c("Dendritic", "Spinous")
names(spinous.labs) <- c("dendritic", "spinous")
sex.labs <- c("Female", "Male")
names(sex.labs) <- c("female", "male")


### DLPFC SYNAPSE TYPE COUNTS

## Use Poisson models for everything. 
## Better fit than Gaussian and more appropriate given small counts
## Negative binomial models collapse to Poisson.

# import dlPFC synapse type counts
dl_typecounts <- readRDS(here("data", "dl_typecounts.rds"))

# calculate case means per synapse type
dl_typemeans <- dl_typecounts %>% group_by(case, group, sex) %>% 
  summarize(perf_spin = mean(perf_spin), 
            perf_dend = mean(perf_dend), 
            nonperf_spin = mean(nonperf_spin), 
            nonperf_dend = mean(nonperf_dend))

# add group by sex column 
dl_typemeans <-  dl_typemeans %>% mutate(combined = paste(group,sex)) %>%
  mutate(combined = factor(combined, levels = c("control female",
                                                "control male",
                                                "anesthesia female",
                                                "anesthesia male"))) 
## Is there an interaction of synapse type, group, and sex?

# add synapse type as a factor
dl_typecounts_long <- dl_typecounts %>% pivot_longer(perf_spin:nonperf_dend, values_to = "synapse_count", names_to = "synapse_type")

# check Gaussian fit of full model
gaussmod_dl <- lmer(synapse_count ~ group*sex*synapse_type + (1|case), data = dl_typecounts_long)
gaussmod_dl %>% model_fit() # it is NOT HAPPY

# analyze dlPFC synapse counts by type with general linear mixed model with poisson distribution
# analyze group, sex, synapse type, and full interaction effects 
dl_counts_1 <- glmer(synapse_count ~ (1|case), data = dl_typecounts_long, family = poisson())
dl_counts_2 <- glmer(synapse_count ~ sex + (1|case), data = dl_typecounts_long, family = poisson())
dl_counts_3 <- glmer(synapse_count ~ group + sex + (1|case), data = dl_typecounts_long, family = poisson())
dl_counts_4 <- glmer(synapse_count ~ group*sex + (1|case), data = dl_typecounts_long, family = poisson())
dl_counts_5 <- glmer(synapse_count ~ synapse_type + group*sex + (1|case), data = dl_typecounts_long, family = poisson())
dl_counts_6 <- glmer(synapse_count ~ sex:synapse_type + synapse_type + group*sex + (1|case), data = dl_typecounts_long, family = poisson(),
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
dl_counts_7 <- glmer(synapse_count ~ group:synapse_type + sex:synapse_type + synapse_type + group*sex + (1|case), data = dl_typecounts_long, family = poisson(),
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
dl_counts_8 <- glmer(synapse_count ~ group*sex*synapse_type + (1|case), data = dl_typecounts_long, family = poisson(),
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

# check fit of full model with poisson distribution
dl_counts_8 %>% model_fit()
testOutliers(dl_counts_8, type = 'bootstrap')

# summarize full fit model
dl_counts_8 %>% summary()

# tabulate group, sex, synapse type, and interaction contributions to effects via chi squared tests
anova(dl_counts_1, dl_counts_2) # sex effects
anova(dl_counts_2, dl_counts_3) # group effects
anova(dl_counts_3, dl_counts_4) # group:sex interaction effects
anova(dl_counts_4, dl_counts_5) # synapse type effects
anova(dl_counts_5, dl_counts_6) # sex:synapse type interaction effects
anova(dl_counts_6, dl_counts_7) # group:synapse type interaction effects
anova(dl_counts_7, dl_counts_8) # group:sex:synapse type interaction effects

# summarize dlPFC synapse counts by group, sex, type 
dl_typecounts_long %>% group_by(synapse_type, group, sex) %>% summarize(mean(synapse_count), sd(synapse_count))


## What are the effects of anesthesia & sex on synapse types? 

# analyze dlPFC perforated spinous synapse counts with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_ps_random <- glmer(perf_spin ~ (1|case), data = dl_typecounts, family = poisson())
dl_ps_sex <- glmer(perf_spin ~ sex + (1|case), data = dl_typecounts, family = poisson())
dl_ps_both <- glmer(perf_spin ~ group + sex + (1|case), data = dl_typecounts, family = poisson())
dl_ps_full <- glmer(perf_spin ~group*sex + (1|case), data = dl_typecounts, family = poisson()) 

# check fit of full model with poisson distribution
dl_ps_full %>% model_fit()

# summarize full fit model
summary(dl_ps_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_ps_random, dl_ps_sex) # main effect of sex 1 df
anova(dl_ps_sex, dl_ps_both) # main effect of group 1 df (controlling for sex)
anova(dl_ps_full, dl_ps_both) # group x sex interaction 1 df


# analyze dlPFC perforated dendritic synapse counts with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_pd_random <- glmer(perf_dend ~ (1|case), data = dl_typecounts, family = poisson())
dl_pd_sex <- glmer(perf_dend ~ sex + (1|case), data = dl_typecounts, family = poisson())
dl_pd_both <- glmer(perf_dend ~ group + sex + (1|case), data = dl_typecounts, family = poisson())
dl_pd_full <- glmer(perf_dend ~group*sex + (1|case), data = dl_typecounts, family = poisson()) 

# check fit of full model with poisson distribution
dl_pd_full %>% model_fit()

# summarize full fit model
summary(dl_pd_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_pd_random, dl_pd_sex) # main effect of sex 1 df
anova(dl_pd_sex, dl_pd_both) # main effect of group 1 df (controlling for sex)
anova(dl_pd_full, dl_pd_both) # group x sex interaction 1 df


# analyze dlPFC nonperforated spinous synapse counts with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_ns_random <- glmer(nonperf_spin ~ (1|case), data = dl_typecounts, family = poisson())
dl_ns_sex <- glmer(nonperf_spin ~ sex + (1|case), data = dl_typecounts, family = poisson())
dl_ns_both <- glmer(nonperf_spin ~ group + sex + (1|case), data = dl_typecounts, family = poisson())
dl_ns_full <- glmer(nonperf_spin ~group*sex + (1|case), data = dl_typecounts, family = poisson()) 

# check fit of full model with poisson distribution
dl_ns_full %>% model_fit()

# summarize full fit model
summary(dl_ns_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_ns_random, dl_ns_sex) # main effect of sex 1 df
anova(dl_ns_sex, dl_ns_both) # main effect of group 1 df (controlling for sex)
anova(dl_ns_full, dl_ns_both) # group x sex interaction 1 df


# analyze dlPFC nonperforated dendritic synapse counts with general linear mixed model with poisson distribution
# analyze random effects, sex effects, group and sex effects, and full interaction effects
dl_nd_random <- glmer(nonperf_dend ~ (1|case), data = dl_typecounts, family = poisson())
dl_nd_sex <- glmer(nonperf_dend ~ sex + (1|case), data = dl_typecounts, family = poisson())
dl_nd_both <- glmer(nonperf_dend ~ group + sex + (1|case), data = dl_typecounts, family = poisson())
dl_nd_full <- glmer(nonperf_dend ~group*sex + (1|case), data = dl_typecounts, family = poisson()) 

# check fit of full model with poisson distribution
dl_nd_full %>% model_fit()

# summarize full fit model
summary(dl_nd_full)

# tabulate sex, group, and interaction contributions to effects via chi squared tests
anova(dl_nd_random, dl_nd_sex) # main effect of sex 1 df
anova(dl_nd_sex, dl_nd_both) # main effect of group 1 df (controlling for sex)
anova(dl_nd_full, dl_nd_both) # group x sex interaction 1 df


# DLPFC SYNAPSE TYPE AREAS

# import dlPFC synapse type areas file
dl_areas_file <- here("data", "dlPFC Synapse Areas.xlsx")

# list all sheets, set column names, read excel sheets into R
dl_areas_data_list <- dl_areas_file %>%
  excel_sheets() %>% 
  set_names() %>% 
  map(read_excel, path = dl_areas_file)

# create names column from sheet names
for (i in 1:length(dl_areas_data_list)) {dl_areas_data_list[[i]]$case <- names(dl_areas_data_list)[i]}

# Cyan: nonperforated, dendritic
# Pink: nonperforated, spinous
# X: perforated, spinous
# Yellow: perforated, dendritic 

# set column names, remove type NAs, combine into single df, extract case name
dl_synapse_areas <- 
  dl_areas_data_list %>% map(~setNames(.x, c("type", "note", "start", "end", "length", "area", "case"))) %>%
  map(~filter(.x, !is.na(type))) %>%
  bind_rows() %>%
  mutate(case = as.numeric(str_extract(case, "-([^-]+)") %>% str_remove("-"))) %>%
  # classify areas as non/perforated and spinous/dendritic
  mutate(perf = case_when(str_detect(type, "cyan") ~ "nonperf",
                          str_detect(type, "pink") ~ "nonperf",
                          str_detect(type, "X") ~ "perf",
                          str_detect(type, "yellow") ~ "perf")) %>%
  mutate(spinous = case_when(str_detect(type, "cyan") ~ "dendritic",
                             str_detect(type, "pink") ~ "spinous",
                             str_detect(type, "X") ~ "spinous",
                             str_detect(type, "yellow") ~ "dendritic"))

# convert to CA1 case code, add treatment and sex, remove suspect cases, select relevant columns 
dl_synapse_areas <- dl_synapse_areas %>% toCA1() %>% 
  addgroups() %>% select(case,group,sex,perf,spinous,area)

dl_synapse_type_area_lmer <- dl_synapse_areas %>% lmer(area ~ group * sex * spinous * perf + (1 | case), data = .)
anova(dl_synapse_type_area_lmer)

# analyze dlPFC synapse areas by type with linear mixed models 
# analyze treatment, sex, and treatment by sex contributions to synapse areas effects by synapse type
dl_overall_area <- lmer(area ~ group * sex * spinous * perf + (1 | case), data = dl_synapse_areas)
dl_reduced_area <- lmer(area ~ group * sex + spinous + perf + spinous:perf + (1 | case), data = dl_synapse_areas)
dl_reduced_area_2 <- lmer(area ~ group * sex + spinous + perf + spinous:perf +
                            sex:spinous + sex:perf + sex:spinous:perf + (1 | case), data = dl_synapse_areas)
dl_reduced_area_3 <- lmer(area ~ group * sex + spinous + perf + spinous:perf +
                            sex:spinous + sex:perf + sex:spinous:perf + 
                            group:spinous + group:perf + group:spinous:perf +(1 | case), data = dl_synapse_areas)

# tabulate treatment, sex, synapse type, and interaction contributions to effects via chi squared tests
anova(dl_reduced_area, dl_reduced_area_2)   # sex effects: sex:spin, sex:perf, and sex:spin:perf interactions
anova(dl_reduced_area_2, dl_reduced_area_3) # group effects: group:spin, group:perf, and group:spin:perf interactions
anova(dl_reduced_area_3, dl_overall_area)   # group by sex effects: group:sex:spin and group:sex:perf interactions

# analyze synapse areas per type by treatment and sex  
DL_perf_spin_area_anova <- 
  dl_synapse_areas %>% filter(perf == "perf") %>% filter(spinous == "spinous") %>%
  lmer(area ~ group * sex + (1 | case), data = .) %>% anova()
DL_perf_dend_area_anova <- 
  dl_synapse_areas %>% filter(perf == "perf") %>% filter(spinous == "dendritic") %>%
  lmer(area ~ group * sex + (1 | case), data = .) %>% anova() # singular, seems due to low variance of random effects 
DL_nonperf_spin_area_anova <- 
  dl_synapse_areas %>% filter(perf == "nonperf") %>% filter(spinous == "spinous") %>%
  lmer(area ~ group * sex + (1 | case), data = .) %>% anova()
DL_nonperf_dend_area_anova <- 
  dl_synapse_areas %>% filter(perf == "nonperf") %>% filter(spinous == "dendritic") %>%
  lmer(area ~ group * sex + (1 | case), data = .) %>% anova()

DL_perf_spin_area_anova
DL_perf_dend_area_anova
DL_nonperf_spin_area_anova
DL_nonperf_dend_area_anova

# calculate perforated spinous area means by treatment and sex
dl_synapse_areas %>% group_by(group, sex) %>% 
  filter(perf == "perf") %>% filter(spinous == "spinous") %>% 
  summarize(mean_area = mean(area), sd(area))
# calculate nonperforated spinous area means by treatment and sex
dl_synapse_areas %>% group_by(group, sex) %>% 
  filter(perf == "nonperf") %>% filter(spinous == "spinous") %>% 
  summarize(mean_area = mean(area), sd(area))


## 20220809: are synapses abnormally large
## 2 SD criterion

dl_synapse_areas %>% group_by(group, spinous) %>% 
  filter(perf == "nonperf") %>% summarize(mean_area = mean(area), sd(area), n())
(dl_abnorm <- dl_synapse_areas %>% 
    filter(perf == "nonperf") %>% filter(spinous == "spinous") %>% 
    summarize(mean_area = mean(area), sd_area = sd(area), n()))
dl_large_synapses <- 
  dl_synapse_areas %>% filter(area > dl_abnorm$mean_area+2*dl_abnorm$sd_area) %>% group_by(group) %>% 
  filter(perf == "nonperf") %>% filter(spinous == "spinous")
dl_large_synapses %>% summarize(mean_area = mean(area), sd(area), n())
# are large nonperforated spinous synapses larger?
t.test(area ~ group, data = dl_large_synapses)
# are there more abnormally large nonperf spinous in control?
dat_ab_dl <- data.frame(
  "not_large" = c((1701-93), (1902-97)),
  "abnorm_large" = c(93, 97),
  row.names = c("Control", "Anesthesia"),
  stringsAsFactors = FALSE
)
colnames(dat_ab_dl) <- c("Not_abnorm", "Abnorm_large")

dat_ab_dl
fisher.test(dat_ab_dl)

# same thing for dendritic?
(dl_abnorm_d <- dl_synapse_areas %>% 
    filter(perf == "nonperf") %>% filter(spinous == "dendritic") %>% 
    summarize(mean_area = mean(area), sd_area = sd(area), n()))
dl_large_synapses_d <- 
  dl_synapse_areas %>% filter(area > dl_abnorm_d$mean_area+2*dl_abnorm_d$sd_area) %>% group_by(group) %>% 
  filter(perf == "nonperf") %>% filter(spinous == "dendritic")
dl_large_synapses_d %>% summarize(mean_area = mean(area), sd(area), n())
# only ~13-15 of each so no



####### generate statistical summary

hierarch_results <- function(model_object) {
  tibble(Chisq = model_object$Chisq[2], df = model_object$Df[2], pval = model_object$`Pr(>Chisq)`[2]) %>%
    mutate_if(is.numeric, round, 5)
}

### CA1

format_anova_output <- function(anova_call) {
  anova_call %>% tidy() %>% 
    suppressWarnings(classes = "warning") %>% 
    select(term, NumDF, DenDF, statistic, p.value) %>% 
    mutate_if(is.numeric, round, 5) %>%
    pivot_longer(NumDF:p.value, names_to = "parameter", values_to = "number")
}

produce_ANOVA_results_1 <- function(tbl) {
  tbl_mod <- tbl %>% mutate_if(is.numeric, round, 2) %>%
    mutate(group_test = paste0("F(",group_NumDF,", ",group_DenDF,") = ",group_statistic,", p = ",group_p.value),
           sex_test = paste0("F(",sex_NumDF,", ",sex_DenDF,") = ",sex_statistic,", p = ",sex_p.value),
           int_test = paste0("F(",`group:sex_NumDF`,", ",`group:sex_DenDF`,") = ",`group:sex_statistic`,", p = ",`group:sex_p.value`))
  tbl_mod %>% select(type, group_test) %>% print()
  tbl_mod %>% select(type, sex_test) %>% print()
  tbl_mod %>% select(type, int_test) %>% print()
  # for F-statistics (lmer)
}

produce_ANOVA_results_2 <- function(tbl) {
  tbl_mod <- tbl %>% mutate_if(is.numeric, round, 3) %>%
    mutate(group_test = paste0("Chisq(",group_df,") = ",group_Chisq,", p = ",group_pval),
           sex_test = paste0("Chisq(",sex_df,") = ",sex_Chisq,", p = ",sex_pval),
           int_test = paste0("Chisq(",`group x sex_df`,") = ",`group x sex_Chisq`,", p = ",`group x sex_pval`))
  tbl_mod %>% select(type, group_test) %>% print()
  tbl_mod %>% select(type, sex_test) %>% print()
  tbl_mod %>% select(type, int_test) %>% print()
  # for chi-squared (glmer)
}

CA1_area_results <- anova(CA1_area_lmer) %>% format_anova_output() %>% mutate(type = "CA1 area")
CA1_density_results <- anova(CA1_density_lmer) %>% format_anova_output() %>% mutate(type = "CA1 density")
CA1_vesicle_results <- anova(CA1_vesicle_lmer) %>% format_anova_output() %>% mutate(type = "CA1 vesicles")
CA1_mito_results <- anova(CA1_mito_lmer) %>% format_anova_output() %>% mutate(type = "CA1 mitochondria density")

bind_rows(CA1_area_results, CA1_density_results, CA1_vesicle_results, CA1_mito_results) %>%
  pivot_wider(everything(), values_from = number, names_from = c(term, parameter)) %>%
  produce_ANOVA_results_1()

anova(CA1_area_quantile_lmer) %>% format_anova_output() %>% mutate_if(is.numeric, round, 2) %>%
  pivot_wider(everything(), names_from = parameter, values_from = number) %>%
  mutate(analysis = "CA1 area by quantile") %>%
  mutate(result = paste0("F(",NumDF,", ",DenDF,") = ",statistic,", p = ",p.value)) %>%
  select(analysis, term, result) %>% print()

# straight mitochondria
CA1_str <- bind_cols(effect = c("sex", "group", "group x sex"),
                     bind_rows(hierarch_results(anova(CA1_str_random, CA1_str_sex)),
                               hierarch_results(anova(CA1_str_sex, CA1_str_both)),
                               hierarch_results(anova(CA1_str_full, CA1_str_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 straight")
# curved mitochondria
CA1_cur <- bind_cols(effect = c("sex", "group", "group x sex"),
                     bind_rows(hierarch_results(anova(CA1_cur_random, CA1_cur_sex)),
                               hierarch_results(anova(CA1_cur_sex, CA1_cur_both)),
                               hierarch_results(anova(CA1_cur_full, CA1_cur_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 curved")
# toroidal mitochondria
CA1_tor <- bind_cols(effect = c("sex", "group", "group x sex"),
                     bind_rows(hierarch_results(anova(CA1_tor_random, CA1_tor_sex)),
                               hierarch_results(anova(CA1_tor_sex, CA1_tor_both)),
                               hierarch_results(anova(CA1_tor_full, CA1_tor_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 toroid")

bind_rows(CA1_str, CA1_cur, CA1_tor) %>%
  pivot_wider(everything(), values_from = number, names_from = c(effect, parameter)) %>%
  produce_ANOVA_results_2()

# bouton frequency: overall
bind_cols(effect = c("sex", "group", "group x sex", "bouton freq",
                     "sex x bouton freq", "group x bouton freq", 
                     "group x sex x bouton freq"),
          bind_rows(hierarch_results(anova(CA1_boutons_1, CA1_boutons_2)),
                    hierarch_results(anova(CA1_boutons_2, CA1_boutons_3)),
                    hierarch_results(anova(CA1_boutons_3, CA1_boutons_4)),
                    hierarch_results(anova(CA1_boutons_4, CA1_boutons_5)),
                    hierarch_results(anova(CA1_boutons_5, CA1_boutons_6)),
                    hierarch_results(anova(CA1_boutons_6, CA1_boutons_7)),
                    hierarch_results(anova(CA1_boutons_7, CA1_boutons_8)))) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(analysis = "CA1 overall bouton frequency") %>%
  mutate(result = paste0("Chisq(",df,") = ",Chisq,", p = ",pval)) %>%
  select(analysis, effect, result) %>% print()

# 0 mito boutons
CA1_0mito <- bind_cols(effect = c("sex", "group", "group x sex"),
                       bind_rows(hierarch_results(anova(CA1_0_random, CA1_0_sex)),
                                 hierarch_results(anova(CA1_0_sex, CA1_0_both)),
                                 hierarch_results(anova(CA1_0_full, CA1_0_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 0 boutons")
# 1 mito boutons
CA1_1mito <- bind_cols(effect = c("sex", "group", "group x sex"),
                       bind_rows(hierarch_results(anova(CA1_1_random, CA1_1_sex)),
                                 hierarch_results(anova(CA1_1_sex, CA1_1_both)),
                                 hierarch_results(anova(CA1_1_full, CA1_1_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 1 boutons")
# 2 mito boutons
CA1_2mito <- bind_cols(effect = c("sex", "group", "group x sex"),
                       bind_rows(hierarch_results(anova(CA1_2_random, CA1_2_sex)),
                                 hierarch_results(anova(CA1_2_sex, CA1_2_both)),
                                 hierarch_results(anova(CA1_2_full, CA1_2_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 2 boutons")
# 3+ mito boutons
CA1_3mito <- bind_cols(effect = c("sex", "group", "group x sex"),
                       bind_rows(hierarch_results(anova(CA1_3_random, CA1_3_sex)),
                                 hierarch_results(anova(CA1_3_sex, CA1_3_both)),
                                 hierarch_results(anova(CA1_3_full, CA1_3_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 3+ boutons")

bind_rows(CA1_0mito, CA1_1mito, CA1_2mito, CA1_3mito) %>%
  pivot_wider(everything(), values_from = number, names_from = c(effect, parameter)) %>%
  produce_ANOVA_results_2()

### DLPFC

DL_area_results <- anova(dlpfc_area_lmer) %>% format_anova_output() %>% mutate(type = "DLPFC area")
DL_density_results <- anova(dlpfc_density_lmer) %>% format_anova_output() %>% mutate(type = "DLPFC density")
DL_vesicle_results <- anova(dlpfc_vesicle_lmer) %>% format_anova_output() %>% mutate(type = "DLPFC vesicles")
DL_mito_results <- anova(dlpfc_mito_lmer) %>% format_anova_output() %>% mutate(type = "DLPFC mitochondria density")

bind_rows(DL_area_results, DL_density_results, DL_vesicle_results, DL_mito_results) %>%
  pivot_wider(everything(), values_from = number, names_from = c(term, parameter)) %>%
  produce_ANOVA_results_1()

anova(dlpfc_area_quantile_lmer) %>% format_anova_output() %>% mutate_if(is.numeric, round, 2) %>%
  pivot_wider(everything(), names_from = parameter, values_from = number) %>%
  mutate(analysis = "DLPFC area by quantile") %>%
  mutate(result = paste0("F(",NumDF,", ",DenDF,") = ",statistic,", p = ",p.value)) %>%
  select(analysis, term, result) %>% print()

# straight mitochondria
DL_str <- bind_cols(effect = c("sex", "group", "group x sex"),
                    bind_rows(hierarch_results(anova(dl_str_random, dl_str_sex)),
                              hierarch_results(anova(dl_str_sex, dl_str_both)),
                              hierarch_results(anova(dl_str_full, dl_str_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC straight")
# curved mitochondria
DL_cur <- bind_cols(effect = c("sex", "group", "group x sex"),
                    bind_rows(hierarch_results(anova(dl_cur_random, dl_cur_sex)),
                              hierarch_results(anova(dl_cur_sex, dl_cur_both)),
                              hierarch_results(anova(dl_cur_full, dl_cur_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC curved")
# toroidal mitochondria
DL_tor <- bind_cols(effect = c("sex", "group", "group x sex"),
                    bind_rows(hierarch_results(anova(dl_tor_random, dl_tor_sex)),
                              hierarch_results(anova(dl_tor_sex, dl_tor_both)),
                              hierarch_results(anova(dl_tor_full, dl_tor_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC toroid")

bind_rows(DL_str, DL_cur, DL_tor) %>%
  pivot_wider(everything(), values_from = number, names_from = c(effect, parameter)) %>%
  produce_ANOVA_results_2()

# bouton frequency: overall
bind_cols(effect = c("sex", "group", "group x sex", "bouton freq",
                     "sex x bouton freq", "group x bouton freq", 
                     "group x sex x bouton freq"),
          bind_rows(hierarch_results(anova(dl_boutons_1, dl_boutons_2)),
                    hierarch_results(anova(dl_boutons_2, dl_boutons_3)),
                    hierarch_results(anova(dl_boutons_3, dl_boutons_4)),
                    hierarch_results(anova(dl_boutons_4, dl_boutons_5)),
                    hierarch_results(anova(dl_boutons_5, dl_boutons_6)),
                    hierarch_results(anova(dl_boutons_6, dl_boutons_7)),
                    hierarch_results(anova(dl_boutons_7, dl_boutons_8)))) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(analysis = "DLPFC overall bouton frequency") %>%
  mutate(result = paste0("Chisq(",df,") = ",Chisq,", p = ",pval)) %>%
  select(analysis, effect, result) %>% print()

# 0 mito boutons
DL_0mito <- bind_cols(effect = c("sex", "group", "group x sex"),
                      bind_rows(hierarch_results(anova(dl_0_random, dl_0_sex)),
                                hierarch_results(anova(dl_0_sex, dl_0_both)),
                                hierarch_results(anova(dl_0_full, dl_0_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC 0 boutons")
# 1 mito boutons
DL_1mito <- bind_cols(effect = c("sex", "group", "group x sex"),
                      bind_rows(hierarch_results(anova(dl_1_random, dl_1_sex)),
                                hierarch_results(anova(dl_1_sex, dl_1_both)),
                                hierarch_results(anova(dl_1_full, dl_1_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC 1 boutons")
# 2 mito boutons
DL_2mito <- bind_cols(effect = c("sex", "group", "group x sex"),
                      bind_rows(hierarch_results(anova(dl_2_random, dl_2_sex)),
                                hierarch_results(anova(dl_2_sex, dl_2_both)),
                                hierarch_results(anova(dl_2_full, dl_2_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC 2 boutons")
# 3+ mito boutons
DL_3mito <- bind_cols(effect = c("sex", "group", "group x sex"),
                      bind_rows(hierarch_results(anova(dl_3_random, dl_3_sex)),
                                hierarch_results(anova(dl_3_sex, dl_3_both)),
                                hierarch_results(anova(dl_3_full, dl_3_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC 3+ boutons")

bind_rows(DL_0mito, DL_1mito, DL_2mito, DL_3mito) %>%
  pivot_wider(everything(), values_from = number, names_from = c(effect, parameter)) %>%
  produce_ANOVA_results_2()

### CA1 synapse type counts
bind_cols(effect = c("sex", "group", "group x sex", "synapse type",
                     "sex x synapse type", "group x synapse type", 
                     "group x sex x synapse type"),
          bind_rows(hierarch_results(anova(CA1_counts_1, CA1_counts_2)),
                    hierarch_results(anova(CA1_counts_2, CA1_counts_3)),
                    hierarch_results(anova(CA1_counts_3, CA1_counts_4)),
                    hierarch_results(anova(CA1_counts_4, CA1_counts_5)),
                    hierarch_results(anova(CA1_counts_5, CA1_counts_6)),
                    hierarch_results(anova(CA1_counts_6, CA1_counts_7)),
                    hierarch_results(anova(CA1_counts_7, CA1_counts_8)))) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(analysis = "CA1 overall synapse type counts") %>%
  mutate(result = paste0("Chisq(",df,") = ",Chisq,", p = ",pval)) %>%
  select(analysis, effect, result) %>% print()

# CA1 perforated spinous counts
CA1_ps_count <- bind_cols(effect = c("sex", "group", "group x sex"),
                          bind_rows(hierarch_results(anova(CA1_ps_random, CA1_ps_sex)),
                                    hierarch_results(anova(CA1_ps_sex, CA1_ps_both)),
                                    hierarch_results(anova(CA1_ps_full, CA1_ps_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 perforated spinous count")
# CA1 perforated dendritic counts
CA1_pd_count <- bind_cols(effect = c("sex", "group", "group x sex"),
                          bind_rows(hierarch_results(anova(CA1_pd_random, CA1_pd_sex)),
                                    hierarch_results(anova(CA1_pd_sex, CA1_pd_both)),
                                    hierarch_results(anova(CA1_pd_full, CA1_pd_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 perforated dendritic count")
# CA1 nonperforated spinous counts
CA1_ns_count <- bind_cols(effect = c("sex", "group", "group x sex"),
                          bind_rows(hierarch_results(anova(CA1_ns_random, CA1_ns_sex)),
                                    hierarch_results(anova(CA1_ns_sex, CA1_ns_both)),
                                    hierarch_results(anova(CA1_ns_full, CA1_ns_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 nonperforated spinous count")
# CA1 nonperforated dendritic counts
CA1_nd_count <- bind_cols(effect = c("sex", "group", "group x sex"),
                          bind_rows(hierarch_results(anova(CA1_nd_random, CA1_nd_sex)),
                                    hierarch_results(anova(CA1_nd_sex, CA1_nd_both)),
                                    hierarch_results(anova(CA1_nd_full, CA1_nd_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "CA1 nonperforated dendritic count")

bind_rows(CA1_ps_count, CA1_pd_count, CA1_ns_count, CA1_nd_count) %>%
  pivot_wider(everything(), values_from = number, names_from = c(effect, parameter)) %>%
  produce_ANOVA_results_2()

# CA1 synapse type areas
anova(CA1_synapse_type_area_lmer) %>% tidy() %>% 
  suppressWarnings(classes = "warning") %>% 
  mutate_if(is.numeric, round, 2) %>%
  mutate(result = paste0("F(",NumDF,", ",DenDF,") = ",statistic,", p = ",p.value)) %>%
  mutate(analysis = "CA1 synapse type areas") %>%
  select(analysis, term, result) %>% print(n=Inf)

bind_cols(effect = c("sex x type", "group x type", "group x sex x type"),
          bind_rows(hierarch_results(anova(reduced_area, reduced_area_2)),
                    hierarch_results(anova(reduced_area_2, reduced_area_3)),
                    hierarch_results(anova(reduced_area_3, overall_area)))) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(analysis = "CA1 overall synapse type areas") %>%
  mutate(result = paste0("Chisq(",df,") = ",Chisq,", p = ",pval)) %>%
  select(analysis, effect, result) %>% print()

CA1_results_ps <- CA1_perf_spin_area_anova %>% format_anova_output() %>% mutate(type = "CA1 perforated spinous area")
CA1_results_pd <- CA1_perf_dend_area_anova %>% format_anova_output() %>% mutate(type = "CA1 perforated dendritic area")
CA1_results_ns <- CA1_nonperf_spin_area_anova %>% format_anova_output() %>% mutate(type = "CA1 nonperforated spinous area")
CA1_results_nd <- CA1_nonperf_dend_area_anova %>% format_anova_output() %>% mutate(type = "CA1 nonperforated dendritic area")

bind_rows(CA1_results_ps, CA1_results_pd, CA1_results_ns, CA1_results_nd) %>%
  pivot_wider(everything(), values_from = number, names_from = c(term, parameter)) %>%
  produce_ANOVA_results_1()


### DLPFC synapse types
bind_cols(effect = c("sex", "group", "group x sex", "synapse type",
                     "sex x synapse type", "group x synapse type", 
                     "group x sex x synapse type"),
          bind_rows(hierarch_results(anova(dl_counts_1, dl_counts_2)),
                    hierarch_results(anova(dl_counts_2, dl_counts_3)),
                    hierarch_results(anova(dl_counts_3, dl_counts_4)),
                    hierarch_results(anova(dl_counts_4, dl_counts_5)),
                    hierarch_results(anova(dl_counts_5, dl_counts_6)),
                    hierarch_results(anova(dl_counts_6, dl_counts_7)),
                    hierarch_results(anova(dl_counts_7, dl_counts_8)))) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(analysis = "DLPFC overall synapse type counts") %>%
  mutate(result = paste0("Chisq(",df,") = ",Chisq,", p = ",pval)) %>%
  select(analysis, effect, result) %>% print()

# DLPFC perforated spinous counts
DL_ps_count <- bind_cols(effect = c("sex", "group", "group x sex"),
                         bind_rows(hierarch_results(anova(dl_ps_random, dl_ps_sex)),
                                   hierarch_results(anova(dl_ps_sex, dl_ps_both)),
                                   hierarch_results(anova(dl_ps_full, dl_ps_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC perforated spinous count")
# DLPFC perforated dendritic counts
DL_pd_count <- bind_cols(effect = c("sex", "group", "group x sex"),
                         bind_rows(hierarch_results(anova(dl_pd_random, dl_pd_sex)),
                                   hierarch_results(anova(dl_pd_sex, dl_pd_both)),
                                   hierarch_results(anova(dl_pd_full, dl_pd_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC perforated dendritic count")
# DLPFC nonperforated spinous counts
DL_ns_count <- bind_cols(effect = c("sex", "group", "group x sex"),
                         bind_rows(hierarch_results(anova(dl_ns_random, dl_ns_sex)),
                                   hierarch_results(anova(dl_ns_sex, dl_ns_both)),
                                   hierarch_results(anova(dl_ns_full, dl_ns_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC nonperforated spinous count")
# DLPFC nonperforated dendritic counts
DL_nd_count <- bind_cols(effect = c("sex", "group", "group x sex"),
                         bind_rows(hierarch_results(anova(dl_nd_random, dl_nd_sex)),
                                   hierarch_results(anova(dl_nd_sex, dl_nd_both)),
                                   hierarch_results(anova(dl_nd_full, dl_nd_both)))) %>%
  pivot_longer(Chisq:pval, names_to = "parameter", values_to = "number") %>%
  mutate(type = "DLPFC nonperforated dendritic count")

bind_rows(DL_ps_count, DL_pd_count, DL_ns_count, DL_nd_count) %>%
  pivot_wider(everything(), values_from = number, names_from = c(effect, parameter)) %>%
  produce_ANOVA_results_2()

# DLPFC synapse type areas
anova(dl_synapse_type_area_lmer) %>% tidy() %>% 
  suppressWarnings(classes = "warning") %>% 
  mutate_if(is.numeric, round, 2) %>%
  mutate(result = paste0("F(",NumDF,", ",DenDF,") = ",statistic,", p = ",p.value)) %>%
  mutate(analysis = "DLPFC synapse type areas") %>%
  select(analysis, term, result) %>% print(n=Inf)

bind_cols(effect = c("sex x type", "group x type", "group x sex x type"),
          bind_rows(hierarch_results(anova(dl_reduced_area, dl_reduced_area_2)),
                    hierarch_results(anova(dl_reduced_area_2, dl_reduced_area_3)),
                    hierarch_results(anova(dl_reduced_area_3, dl_overall_area)))) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(analysis = "DLPFC overall synapse type areas") %>%
  mutate(result = paste0("Chisq(",df,") = ",Chisq,", p = ",pval)) %>%
  select(analysis, effect, result) %>% print()

DL_results_ps <- DL_perf_spin_area_anova %>% format_anova_output() %>% mutate(type = "DLPFC perforated spinous area")
DL_results_pd <- DL_perf_dend_area_anova %>% format_anova_output() %>% mutate(type = "DLPFC perforated dendritic area")
DL_results_ns <- DL_nonperf_spin_area_anova %>% format_anova_output() %>% mutate(type = "DLPFC nonperforated spinous area")
DL_results_nd <- DL_nonperf_dend_area_anova %>% format_anova_output() %>% mutate(type = "DLPFC nonperforated dendritic area")

bind_rows(DL_results_ps, DL_results_pd, DL_results_ns, DL_results_nd) %>%
  pivot_wider(everything(), values_from = number, names_from = c(term, parameter)) %>%
  produce_ANOVA_results_1()

## correlations with behavior

behav_z <- readRDS(here("data", "behav_z.rds"))

ca1_for_corr <- CA1_area_means %>% left_join(behav_z, by = "case")

dl_mean_areas <- dl_areas2 %>% group_by(case) %>% summarize(mean_dl_area = mean(synapse_area))
ca1_for_corr <- ca1_for_corr %>% left_join(dl_mean_areas, by = "case")

cor.test(ca1_for_corr$mean_area, ca1_for_corr$mean_score_z_con) 
cor.test(ca1_for_corr$mean_dl_area, ca1_for_corr$mean_score_z_con)

ca1_for_corr %>% select(case, group, CA1 = mean_area, DLPFC = mean_dl_area, behav_z = mean_score_z_con) %>% 
  pivot_longer(CA1:DLPFC, names_to = "region", values_to = "mean_area") %>% 
  ggplot(aes(x = behav_z, y = mean_area)) + geom_point(aes(color = group)) + 
  geom_smooth(method = "lm", se = FALSE) + facet_wrap(~ region)

ca1quant19 <- CA1quantiles %>% filter(quantile == 19) %>% select(case, ca1_19 = area)
ca1quant20 <- CA1quantiles %>% filter(quantile == 20) %>% select(case, ca1_20 = area)
DLquant19 <- dl_quantiles2 %>% filter(quantile == 19) %>% select(case, DL_19 = area)
DLquant20 <- dl_quantiles2 %>% filter(quantile == 20) %>% select(case, DL_20 = area)

ca1_for_corr <- ca1_for_corr %>% left_join(ca1quant19, by = "case") %>%
  left_join(ca1quant20, by = "case") %>%
  left_join(DLquant19, by = "case") %>%
  left_join(DLquant20, by = "case")

cor.test(ca1_for_corr$ca1_19, ca1_for_corr$mean_score_z_con)
cor.test(ca1_for_corr$ca1_20, ca1_for_corr$mean_score_z_con)
cor.test(ca1_for_corr$DL_19, ca1_for_corr$mean_score_z_con)
cor.test(ca1_for_corr$DL_20, ca1_for_corr$mean_score_z_con) # ns

ca1_for_corr <- CA1synapse_areas %>% filter(spinous == "spinous") %>% filter(perf == "perf") %>% 
  group_by(case) %>% summarize(CA1_perf_spin_area = mean(area)) %>%
  right_join(ca1_for_corr, by = "case")

ca1_for_corr <- dl_synapse_areas %>% filter(spinous == "spinous") %>% filter(perf == "nonperf") %>% 
  group_by(case) %>% summarize(DL_nonperf_spin_area = mean(area)) %>%
  right_join(ca1_for_corr, by = "case")

cor.test(ca1_for_corr$CA1_perf_spin_area, ca1_for_corr$mean_score_z_con)
cor.test(ca1_for_corr$DL_nonperf_spin_area, ca1_for_corr$mean_score_z_con) # ns



######  FIGURES

# Synapse Area Bar Graphs

CA1areameansgraph <- CA1_area_means %>% 
  ggplot(aes(x=group,y=mean_area)) + 
  geom_jitter(aes(color=group), size=2, 
              position = position_jitter(width=0.1, height = 0,seed = 75)) + 
  geom_bar(data=group_mean_area, aes(x=group,y=mean,fill=group), 
           alpha=0.3, stat="identity") + 
  geom_line(data=tibble(x=c(1,2), y =c(.4,.4)), aes(x=x, y=y)) +
  geom_text(data=tibble(x=c(1.5), y =c(.41)), aes(x=x, y=y, label="**"), size=7) +
  theme_anes +
  #  theme(aspect.ratio=0.5) +
  ggtitle("Mean Synapse Areas", subtitle = "CA1") +
  xlab("") + ylab(bquote('Mean synapse area, '~µm^2)) + 
  guides(fill="none",color="none") + 
  scale_x_discrete(labels=c("Control", "Anesthesia")) +
  scale_y_continuous(limits=c(0,0.45), expand=c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","indianred1")) +
  scale_fill_manual(values=c("dodgerblue3","indianred1"))  

dlareameansgraph <- dl_areas2 %>% group_by(case) %>% 
  summarize(mean = mean(synapse_area), group = unique(group)) %>%
  mutate(group = factor(group, levels=c("control", "anesthesia"))) %>% 
  ungroup() %>%
  ggplot(aes(x=group,y=mean)) + 
  geom_jitter(aes(color=group), size=2, 
              position = position_jitter(width=0.1, height = 0,seed = 60)) + 
  geom_bar(aes(fill=group), position="dodge", 
           alpha=0.3, stat = "summary") + 
  theme_anes +
  #  theme(aspect.ratio=0.5) +
  ggtitle("", subtitle = "dlPFC") +
  xlab("") + ylab(bquote('Mean synapse area, '~µm^2)) + 
  guides(fill="none",color="none") + 
  scale_x_discrete(labels=c("Control", "Anesthesia")) +
  scale_y_continuous(limits=c(0,0.5), expand=c(0,0)) +
  scale_fill_manual(values=c("dodgerblue3","indianred1")) + 
  scale_color_manual(values=c("dodgerblue3","indianred1")) 

# Synapse Area Line Graphs

CA1arealineanes <- CA1quantiles %>% 
  group_by(group,quantile) %>% 
  summarize(mean = mean(area), sd=sd(area)) %>% 
  ggplot(aes(x=quantile,y=mean,color=group)) + 
  geom_point(size=2) + 
  geom_line() + 
  geom_linerange(aes(ymin=mean-sd,ymax=mean+sd)) +
  annotate(geom = "text", x = 20, y = 1.9, label = "**", size = 7) +
  annotate(geom = "text", x = 18.9, y = 1.1, label = "**", size = 7) +
  annotate(geom = "text", x = 18, y = .9, label = "**", size = 7) +
  annotate(geom = "text", x = 17, y = .75, label = "*", size = 7) +
  theme_anes +
  theme(legend.position = c(.1,.8)) +
  ggtitle("Synapse Areas by Quantile", subtitle = "CA1") +
  xlab("Quantile (5% of total synapses)") +  
  ylab(bquote('Mean ± SD synapse area,'~µm^2)) +
  scale_x_continuous(breaks=1:20, labels=1:20) + 
  scale_y_continuous(limit=c(0,2.2)) +
  scale_color_manual(name="", values=c("dodgerblue3","indianred1"),
                     labels=c("Control","Anesthesia"))

ann_text <- data.frame(quantile = 19,wt = 1.25,lab = "**",
                       sex = factor("Female",levels = c("Female", "Male")))

inset <- CA1quantiles %>% 
  group_by(group,sex,quantile) %>% 
  summarize(mean = mean(area), sd=sd(area)) %>% ungroup () %>%
  filter(quantile > 17) %>% 
  mutate(sex = factor(sex, labels = c("Female", "Male"))) %>%
  ggplot(aes(x=quantile,y=mean,
             color=interaction(group,sex),
             shape=interaction(group,sex))) + 
  xlab("Quantile") + ylab(bquote('Synapse area,'~µm^2)) +
  geom_point(size=3) + geom_line() + 
  geom_linerange(aes(ymin=mean-sd,ymax=mean+sd)) +
  #  geom_text(data = ann_text,label = "**") +
  scale_x_continuous(breaks=1:20, labels=1:20) + 
  theme_anes +  
  scale_color_manual(name = "", values=c("dodgerblue3","indianred1", "dodgerblue4", "indianred3"),
                     labels=c("Control Female ","Anesthesia Female","Control Male","Anesthesia Male")) + 
  scale_shape_manual(name = "", values=c(19, 19, 17, 17),
                     labels=c("Control Female ","Anesthesia Female","Control Male","Anesthesia Male")) +
  labs(color = "", shape = "") + guides(color = "none", shape = "none") +
  facet_wrap(sex~.) +
  theme(strip.text.x = element_text(size = 10), 
        strip.placement = "outside", 
        strip.background = element_rect(color = NA),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))

dlarealinegraph <- dl_quantiles2 %>% group_by(group,quantile) %>% 
  summarize(mean = mean(area), sd=sd(area)) %>% 
  ggplot(aes(x=quantile,y=mean,color=group)) + 
  geom_point(size=2) + 
  geom_line() + 
  geom_linerange(aes(ymin=mean-sd,ymax=mean+sd)) +
  annotate(geom = "text", x = 20, y = 1.9, label = "***", size = 7) +
  theme_anes +
  theme(legend.position = c(.15,.8)) +
  ggtitle("", subtitle = "dlPFC") +
  xlab("Quantile (5% of total synapses)") +  
  ylab(bquote('Mean ± SD synapse area,'~µm^2)) +
  scale_x_continuous(breaks=1:20, labels=1:20) + 
  scale_y_continuous(limit=c(0,2.0)) +
  scale_color_manual(name="", values=c("dodgerblue3","indianred1"),
                     labels=c("Control","Anesthesia")) 

meansandline <- 
  CA1areameansgraph + 
  CA1arealineanes + 
  (inset_element(inset, 0.1, 0.4, .5, 1,
                 align_to = 'plot',
                 ignore_tag = T) +
     theme(plot.background = element_rect(linetype = 'solid', colour = 'black'))) + 
  dlareameansgraph + dlarealinegraph +
  plot_layout(widths = c(1,2)) 

# Synapse Density

CA1syndensgraph <- 
  CA1densmeans %>% 
  mutate(combined = paste(group,sex)) %>%
  mutate(combined = factor(combined, levels = c("control female",
                                                "control male",
                                                "anesthesia female",
                                                "anesthesia male"))) %>%
  ggplot(aes(x=combined,y=mean_dens)) + 
  geom_jitter(aes(color=combined, shape=sex), 
              size=2, 
              position = position_jitter(width=0.1, height = 0, seed = 61)) + 
  geom_bar(aes(fill=combined), position="dodge", 
           alpha=0.3, stat = "summary") + 
  geom_line(data=tibble(x=c(1.5,3.5), y =c(2.8,2.8)), aes(x=x, y=y)) +
  geom_line(data=tibble(x=c(1.5,1.5), y =c(2.35,2.8)), aes(x=x, y=y)) +
  geom_line(data=tibble(x=c(3.5,3.5), y =c(2.6,2.8)), aes(x=x, y=y)) +
  geom_line(data=tibble(x=c(1,3), y =c(2.35,2.35)), aes(x=x, y=y)) +
  geom_line(data=tibble(x=c(1,1), y =c(2.25,2.35)), aes(x=x, y=y)) +
  geom_line(data=tibble(x=c(3,3), y =c(2.25,2.35)), aes(x=x, y=y)) +
  geom_line(data=tibble(x=c(2,4), y =c(2.6,2.6)), aes(x=x, y=y)) +
  geom_line(data=tibble(x=c(2,2), y =c(2.5,2.6)), aes(x=x, y=y)) +
  geom_line(data=tibble(x=c(4,4), y =c(2.5,2.6)), aes(x=x, y=y)) +
  geom_text(data=tibble(x=c(2.5), y =c(2.85)), 
            aes(x=x, y=y, label="*"), size=7) +
  theme_anes +
  ggtitle("Synapse Density", subtitle = "CA1") +
  xlab("") + ylab(bquote('Mean synapses per'~µm^3)) + 
  guides(shape="none", fill="none", color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale", 
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,3.1), expand=c(0,0)) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",
                             "indianred1", "indianred3")) + 
  scale_color_manual(name = "", values=c("dodgerblue3","dodgerblue4",
                                         "indianred1", "indianred3"),
                     labels=c("Control Female (CF)","Control Male (CM)",
                              "Anesthesia Female (AF)","Anesthesia Male (AM)")) 

dlsyndensgraph <- dl_density2 %>% 
  group_by(case,group,sex) %>% 
  summarize(mean_dens = mean(density)) %>%
  mutate(combined = paste(group,sex)) %>%
  mutate(combined = factor(combined, 
                           levels = c("control female", 
                                      "control male",
                                      "anesthesia female",
                                      "anesthesia male"))) %>%
  ggplot(aes(x=combined,y=mean_dens)) + 
  geom_jitter(aes(color=combined, shape=sex), 
              size=2, 
              position = position_jitter(width=0.1, height = 0, seed = 61)) + 
  geom_bar(aes(fill=combined), position="dodge", 
           alpha=0.3, stat = "summary") + 
  theme_anes +
  ggtitle("", subtitle = "dlPFC") +
  xlab("") + ylab(bquote('Mean synapses per'~µm^3)) + 
  guides(shape="none",fill="none",color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale", 
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,1.5), expand=c(0,0)) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",
                             "indianred1", "indianred3")) + 
  scale_color_manual(name = "", values=c("dodgerblue3","dodgerblue4",
                                         "indianred1", "indianred3"),
                     labels=c("Control Female (CF)","Control Male (CM)",
                              "Anesthesia Female (AF)","Anesthesia Male (AM)")) 

# Vesicle Docking 

CA1vesiclegraph <- CA1bins %>% 
  ggplot(aes(x=group, y=ratio, fill = group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .5) +
  geom_point(aes(y = ratio, color = group), 
             position = position_jitter(width = .15), size = .5, alpha = 0.6) +
  geom_boxplot(width = .1, guides = "none", outlier.shape = NA, alpha = 0.3) +
  theme_anes + 
  coord_flip() +
  ggtitle("Vesicle Docking", subtitle = "CA1") +
  expand_limits(x = 3) +
  guides(fill="none",color="none") +
  xlab("") + ylab(bquote('Percent of Actively Docking Vesicles at Synapse')) +
  scale_x_discrete(labels=c("Anesthesia", "Control")) +
  scale_fill_manual(values=c("indianred1", "dodgerblue3")) + 
  scale_color_manual(values=c("indianred1", "dodgerblue3")) 

dlvesiclegraph <- dlpfcbins %>% 
  ggplot(aes(x=group, y=ratio, fill = group)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .5) +
  geom_point(aes(y = ratio, color = group), 
             position = position_jitter(width = .15), size = .5, alpha = 0.6) +
  geom_boxplot(width = .1, guides = "none", outlier.shape = NA, alpha = 0.3) +
  theme_anes +
  coord_flip() +
  ggtitle("", subtitle = "dlPFC") +
  expand_limits(x = 3) +
  guides(fill="none",color="none") +
  xlab("") + ylab(bquote('Percent of Actively Docking Vesicles at Synapse')) +
  scale_x_discrete(labels=c("Anesthesia", "Control")) +
  scale_fill_manual(values=c("indianred1", "dodgerblue3")) + 
  scale_color_manual(values=c("indianred1", "dodgerblue3")) 

synandves <- 
  CA1syndensgraph + CA1vesiclegraph + dlsyndensgraph + dlvesiclegraph +
  plot_layout(widths = c(1,2))

syn1 <- meansandline / synandves + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

# Perforated Spinous Synapse Counts

CA1perfspincountgraph <- CA1countmeans %>% 
  ggplot(aes(x=combined,y=perf_spin)) + 
  geom_jitter(aes(color=combined, shape=sex), 
              size=2, 
              position = position_jitter(width=0.1, height = 0,seed = 61)) + 
  geom_bar(aes(fill=combined), position="dodge", 
           alpha=0.3, stat = "summary") + 
  theme_anes + 
  ggtitle("Perforated Spinous\nSynapse Counts", subtitle = "CA1") +
  xlab("") + ylab(bquote('Mean number of synapses')) +
  guides(shape="none", fill="none", color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale",  
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,7), expand=c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","dodgerblue4",  
                              "indianred1","indianred3")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",  
                             "indianred1","indianred3")) 

dlperfspincountgraph <- dl_typemeans %>% 
  ggplot(aes(x=interaction(group,sex),y=perf_spin)) + 
  geom_jitter(aes(color=interaction(group,sex), shape=sex), 
              size=2, 
              position = position_jitter(width=0.2, height = 0,seed = 73)) + 
  geom_bar(aes(fill=interaction(group,sex)), position="dodge", 
           alpha=0.3, stat = "summary") + 
  theme_anes + 
  ggtitle("", subtitle = "dlPFC") +
  xlab("") + ylab(bquote('Mean number of synapses')) +
  guides(shape="none", fill="none", color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale",  
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,8), expand=c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","dodgerblue4",  
                              "indianred1","indianred3")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",  
                             "indianred1","indianred3")) 

# Perforated Dendritic Synapse Counts

signifbars <- read.table(header=TRUE,
                         text =
                           "x xend y yend
                          1.5 3.5 2.5 2.5
                          3.5 3.5 2.3 2.5
                          3.0 4.0 2.3 2.3
                          3.0 3.0 2.2 2.3
                          4.0 4.0 2.2 2.3
                          1.5 1.5 1.9 2.5
                          1.0 2.0 1.9 1.9
                          1.0 1.0 1.8 1.9
                          2.0 2.0 1.8 1.9",
                         stringsAsFactors = FALSE) 

CA1perfdendcountgraph <- 
  CA1countmeans %>% 
  ggplot(aes(x=combined,y=perf_dend)) + 
  geom_jitter(aes(color=combined, shape=sex), 
              size=2, 
              position = position_jitter(width=0.2, height = 0,seed = 69)) + 
  geom_bar(aes(fill=combined), position="dodge", 
           alpha=0.3, stat = "summary") + 
  geom_segment(data=signifbars, 
               aes(x=x, xend=xend, y=y*.29, yend=yend*.29)) +  
  geom_text(data=tibble(x=c(2.5), y =c(.75)), 
            aes(x=x, y=y, label="*"), size=7) +
  theme_anes + 
  ggtitle("Perforated Dendritic\nSynapse Counts", subtitle = "CA1") +
  xlab("") + ylab(bquote('Mean number of synapses')) +
  guides(shape="none", fill="none", color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale",  
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,.8)) +
  scale_color_manual(values=c("dodgerblue3","dodgerblue4",  
                              "indianred1","indianred3")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",  
                             "indianred1","indianred3")) 

dlperfdendcountgraph <- dl_typemeans %>% 
  ggplot(aes(x=combined,y=perf_dend)) + 
  geom_jitter(aes(color=combined, shape=sex), 
              size=2, 
              position = position_jitter(width=0.1, height = 0,seed = 9445)) + 
  geom_bar(aes(fill=combined), position="dodge", 
           alpha=0.3, stat = "summary") + 
  theme_anes + 
  ggtitle("", subtitle = "dlPFC") +
  xlab("") + ylab(bquote('Mean number of synapses')) +
  guides(shape="none", fill="none", color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale",  
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","dodgerblue4",  
                              "indianred1","indianred3")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",  
                             "indianred1","indianred3")) 

# Nonperforated Spinous Synapse Counts

CA1nonperfspincountgraph <- CA1countmeans %>% 
  ggplot(aes(x=combined,y=nonperf_spin)) + 
  geom_jitter(aes(color=combined, shape=sex), 
              size=2, 
              position = position_jitter(width=0.1, height = 0,seed = 61)) + 
  geom_bar(aes(fill=combined), position="dodge", 
           alpha=0.3, stat = "summary") + 
  theme_anes + 
  ggtitle("Nonperforated Spinous\nSynapse Counts", subtitle = "CA1") +
  xlab("") + ylab(bquote('Mean number of synapses')) +
  guides(shape="none", fill="none", color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale",  
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,30), expand=c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","dodgerblue4",  
                              "indianred1","indianred3")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",  
                             "indianred1","indianred3")) 

dlnonperfspincountgraph <- dl_typemeans %>% 
  ggplot(aes(x=combined,y=nonperf_spin)) + 
  geom_jitter(aes(color=combined, shape=sex), 
              size=2, 
              position = position_jitter(width=0.1, height = 0,seed = 65)) + 
  geom_bar(aes(fill=combined), position="dodge", 
           alpha=0.3, stat = "summary") + 
  theme_anes + 
  ggtitle("", subtitle = "dlPFC") +
  xlab("") + ylab(bquote('Mean number of synapses')) +
  guides(shape="none", fill="none", color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale",  
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,20), expand=c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","dodgerblue4",  
                              "indianred1","indianred3")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",  
                             "indianred1","indianred3")) 

# Nonperforated Dendritic Synapse Counts

CA1nonperfdendcountgraph <- 
  CA1countmeans %>% 
  ggplot(aes(x=combined,y=nonperf_dend)) + 
  geom_jitter(aes(color=combined, shape=sex), 
              size=2, 
              position = position_jitter(width=0.1, height = 0,seed = 87)) + 
  geom_bar(aes(fill=combined), position="dodge", 
           alpha=0.3, stat = "summary") + 
  geom_segment(data=signifbars, 
               aes(x=x, xend=xend, y=y*1.32, yend=yend*1.32)) +  
  geom_text(data=tibble(x=c(2.5), y =c(3.4)), 
            aes(x=x, y=y, label="*"), size=7) +
  theme_anes + 
  ggtitle("Nonperforated Dendritic\nSynapse Counts", subtitle = "CA1") +
  xlab("") + ylab(bquote('Mean number of synapses')) +
  guides(shape="none", fill="none", color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale",  
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,3.6), expand=c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","dodgerblue4",  
                              "indianred1","indianred3")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",  
                             "indianred1","indianred3")) 

dlnonperfdendcountgraph <- dl_typemeans %>% 
  ggplot(aes(x=combined,y=nonperf_dend)) + 
  geom_jitter(aes(color=combined, shape=sex), 
              size=2, 
              position = position_jitter(width=0.1, height = 0,seed = 61)) + 
  geom_bar(aes(fill=combined), position="dodge", 
           alpha=0.3, stat = "summary") + 
  theme_anes + 
  ggtitle("", subtitle = "dlPFC") +
  xlab("") + ylab(bquote('Mean number of synapses')) +
  guides(shape="none", fill="none", color="none") + 
  scale_x_discrete(labels=c("Control\nFemale", "Control\nMale",  
                            "Anesthesia\nFemale", "Anesthesia\nMale")) +
  scale_y_continuous(limits=c(0,5), expand=c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","dodgerblue4",  
                              "indianred1","indianred3")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",  
                             "indianred1","indianred3")) 

# Perforated Spinous Synapse Areas

CA1synapse_areas$combined <- paste(CA1synapse_areas$group, 
                                   CA1synapse_areas$sex)

CA1perfspinareagraph <- 
  CA1synapse_areas %>% 
  filter(perf == "perf") %>% filter(spinous == "spinous") %>%
  mutate(combined = factor(combined, levels = c("anesthesia male",
                                                "anesthesia female",
                                                "control male",
                                                "control female"))) %>%
  ggplot(aes(x=combined, y=area)) +
  geom_flat_violin(aes(fill = combined), 
                   position = position_nudge(x = .2, y = 0), alpha = .5) +
  geom_point(aes(y = area, color = combined), 
             position = position_jitter(width = .15), 
             size = .5, alpha = 0.6) +
  geom_segment(data=signifbars, 
               mapping = aes(x=x+.2, xend=xend+.2, y=y*.97, yend=yend*.97)) +
  geom_boxplot(width = .1, guides = "none", 
               outlier.shape = NA, alpha = 0.3) +
  annotate(geom = "text", x = 2.7, y = 2.6, label = "*", size = 7) +
  theme_anes + 
  coord_flip() + 
  ggtitle("Perforated Spinous\nSynapse Areas",subtitle="CA1") +
  xlab("") + ylab(bquote('Synapse Area,'~µm^2)) + 
  guides(shape="none",fill="none",color="none") +
  scale_x_discrete(labels=c("Anesthesia\nMale", "Anesthesia\nFemale", 
                            "Control\nMale", "Control\nFemale")) +
  scale_y_continuous(limits = c(0,2.6), breaks = c(0, .5, 1, 1.5, 2, 2.5)) +
  scale_color_manual(values=c("indianred3","indianred1",
                              "dodgerblue4", "dodgerblue3")) +
  scale_fill_manual(values=c("indianred3","indianred1",
                             "dodgerblue4", "dodgerblue3")) 

dl_synapse_areas$combined <- paste(dl_synapse_areas$group, 
                                   dl_synapse_areas$sex)

dlperfspinareagraph <- dl_synapse_areas %>% 
  filter(perf == "perf") %>% filter(spinous == "spinous") %>%
  mutate(combined = factor(combined, levels = c("anesthesia male",
                                                "anesthesia female",
                                                "control male",
                                                "control female"))) %>%
  ggplot(aes(x=combined, y=area)) +
  geom_flat_violin(aes(fill = combined), 
                   position = position_nudge(x = .2, y = 0), alpha = .5) +
  geom_point(aes(y = area, color = combined), 
             position = position_jitter(width = .15), 
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .1, guides = "none", outlier.shape = NA, alpha = 0.3) +
  theme_anes + 
  coord_flip() + 
  ggtitle("",subtitle="dlPFC") +
  xlab("") + ylab(bquote('Synapse Area,'~µm^2)) + 
  guides(shape="none",fill="none",color="none") +
  scale_x_discrete(labels=c("Anesthesia\nMale", "Anesthesia\nFemale", 
                            "Control\nMale", "Control\nFemale")) +
  scale_y_continuous(limits = c(0,2), breaks = c(0, .5, 1, 1.5, 2)) +
  scale_color_manual(values=c("indianred3","indianred1",
                              "dodgerblue4", "dodgerblue3")) +
  scale_fill_manual(values=c("indianred3","indianred1",
                             "dodgerblue4", "dodgerblue3")) 

# Nonperforated Spinous Synapse Areas

CA1nonperfspinareagraph <- 
  CA1synapse_areas %>% 
  filter(perf == "nonperf") %>% filter(spinous == "spinous") %>%
  mutate(combined = factor(combined, levels = c("anesthesia male",
                                                "anesthesia female",
                                                "control male",
                                                "control female"))) %>%
  ggplot(aes(x=combined, y=area)) +
  geom_flat_violin(aes(fill = combined), 
                   position = position_nudge(x = .2, y = 0), alpha = .5) +
  geom_point(aes(y = area, color = combined), 
             position = position_jitter(width = .15), 
             size = .5, alpha = 0.6) +
  geom_boxplot(width = .1, guides = "none", 
               outlier.shape = NA, alpha = 0.3) +
  theme_anes + 
  coord_flip() + 
  ggtitle("Nonperforated Spinous\nSynapse Areas",subtitle="CA1") +
  xlab("") + ylab(bquote('Synapse Area,'~µm^2)) + 
  guides(shape="none",fill="none",color="none") +
  scale_x_discrete(labels=c("Anesthesia\nMale", "Anesthesia\nFemale", 
                            "Control\nMale", "Control\nFemale")) +
  scale_color_manual(values=c("indianred3","indianred1",
                              "dodgerblue4", "dodgerblue3")) +
  scale_fill_manual(values=c("indianred3","indianred1",
                             "dodgerblue4", "dodgerblue3")) 

signifbars2 <- read.table(header=TRUE,
                          text =
                            "x xend y yend
                          1.5 3.5 2.5 2.5
                          3.5 3.5 2.4 2.5
                          3.0 4.0 2.4 2.4
                          3.0 3.0 2.35 2.4
                          4.0 4.0 2.35 2.4
                          1.5 1.5 2.3 2.5
                          1.0 2.0 2.3 2.3
                          1.0 1.0 2.25 2.3
                          2.0 2.0 2.25 2.3",
                          stringsAsFactors = FALSE) 

dlnonperfspinareagraph <- 
  dl_synapse_areas %>% 
  filter(perf == "nonperf") %>% filter(spinous == "spinous") %>%
  mutate(combined = factor(combined, levels = c("anesthesia male",
                                                "anesthesia female",
                                                "control male",
                                                "control female"))) %>%
  ggplot(aes(x=combined, y=area)) +
  geom_flat_violin(aes(fill = combined), 
                   position = position_nudge(x = .2, y = 0), alpha = .5) +
  geom_point(aes(y = area, color = combined), 
             position = position_jitter(width = .15), 
             size = .5, alpha = 0.6) +
  geom_segment(data=signifbars2, 
               mapping = aes(x=x+.2, xend=xend+.2, y=y*.61, yend=yend*.61)) +
  geom_boxplot(width = .1, guides = "none", 
               outlier.shape = NA, alpha = 0.3) +
  annotate(geom = "text", x = 2.7, y = 1.6, label = "*", size = 7) +
  theme_anes + 
  coord_flip() + 
  ggtitle("",subtitle="dlPFC") +
  xlab("") + ylab(bquote('Synapse Area,'~µm^2)) + 
  guides(shape="none",fill="none",color="none") +
  scale_x_discrete(labels=c("Anesthesia\nMale", "Anesthesia\nFemale", 
                            "Control\nMale", "Control\nFemale")) +
  scale_color_manual(values=c("indianred3","indianred1",
                              "dodgerblue4", "dodgerblue3")) +
  scale_fill_manual(values=c("indianred3","indianred1",
                             "dodgerblue4", "dodgerblue3")) 

col1 <-
  CA1perfdendcountgraph + dlperfdendcountgraph + CA1nonperfdendcountgraph +
  dlnonperfdendcountgraph + 
  plot_layout(nrow=4)

col2 <-
  CA1perfspinareagraph + CA1nonperfspinareagraph +
  dlperfspinareagraph + dlnonperfspinareagraph +
  plot_layout(nrow=2)

syntypesall <- (col1 | plot_spacer() | col2) + plot_layout(widths = c(1.1,.05,2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

# Mitochondria Density 

CA1mitodensgraph <- CA1_mitodensity %>%
  ggplot(aes(x=group,y=mean_mito)) + 
  geom_jitter(aes(color=group),
              position=position_jitter(seed=12,width=0.1, height = 0), 
              size = 2) + 
  geom_bar(aes(fill=group), position = "dodge", alpha=0.3, stat="summary") + 
  theme_anes + 
  ggtitle("Mitochondria Density", subtitle = "CA1") +
  xlab("") + ylab(bquote('Mean mitochondria per'~µm^3)) +
  guides(shape = "none", fill = "none",color = "none") + 
  scale_x_discrete(labels=c("Control", "Anesthesia")) +
  scale_y_continuous(limits=c(0,1), expand = c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","indianred1")) +
  scale_fill_manual(values=c("dodgerblue3","indianred1"))

dlmitodensgraph <- dl_mitodensity %>%
  ggplot(aes(x=group,y=mean_mito)) + 
  geom_jitter(aes(color=group),
              position=position_jitter(seed=12,width=0.1, height = 0), 
              size = 2) + 
  geom_bar(aes(fill=group), position="dodge", alpha=0.3, stat = "summary") + 
  theme_anes + 
  ggtitle("", subtitle = "dlPFC") +
  xlab("") + ylab(bquote('Mean mitochondria per'~µm^3)) + 
  guides(shape="none",fill="none",color="none") + 
  scale_x_discrete(labels=c("Control", "Anesthesia")) +
  scale_y_continuous(limits=c(0,0.61), expand = c(0,0)) +
  scale_color_manual(values=c("dodgerblue3","indianred1")) +
  scale_fill_manual(values=c("dodgerblue3","indianred1"))

# Mitochondria Shape

CA1shapefacets <- read.table(header=TRUE,
                             text = 
                               "morphology ymin ymax
                                Straight 0 40
                                Curved 0 7.7
                                Donut -0.05 1.5", 
                             stringsAsFactors = FALSE)
CA1_shape_facets <- with(CA1shapefacets,
                         data.frame(count=c(ymin,ymax),
                                    morphology=c(morphology, morphology)))

CA1mitoshapegraph <- CA1_mitoshape %>%
  ggplot(aes(x=interaction(group,morphology),y=count)) + 
  geom_jitter(aes(color=group),
              position=position_jitter(seed=582,width=0.2, height = 0), 
              size = 2) + 
  geom_bar(aes(fill=group), position = "dodge", alpha=0.3, stat="summary") + 
  geom_point(data = CA1_shape_facets, x=NA) +
  theme_anes +
  theme(strip.text.x = element_text(size = 12), 
        strip.placement = "outside", 
        strip.background = element_blank()) +
  facet_wrap(~factor(morphology, levels = c("Straight", "Curved", "Donut")), 
             scales = "free") +
  ggtitle("Mitochondria Shape Counts", subtitle = "CA1") +
  xlab("") + ylab(bquote('Mean number of\nmitochondria')) + 
  ggtitle("Mitochondria Shapes") +
  guides(shape="none",fill="none",color="none") + 
  scale_x_discrete(labels=c("Control", "Anesthesia")) +
  scale_fill_manual(values=c("dodgerblue3","indianred1")) + 
  scale_color_manual(values=c("dodgerblue3","indianred1")) +
  geom_signif(data = data.frame(morphology = c("Curved")),
              aes(y_position=c(6.5), xmin = c(1), xmax = c(2), 
                  annotations = c("*")), tip_length = c(0), textsize = c(7),
              manual = T)

dlmitoshapegraph <- dl_mitoshape %>%
  ggplot(aes(x=interaction(group,morphology),y=count)) + 
  geom_jitter(aes(color=group),
              position=position_jitter(seed=594,width=0.2, height = 0), 
              size = 2) + 
  geom_bar(aes(fill=group), position = "dodge", alpha=0.3, stat="summary") + 
  theme_anes +
  theme(strip.text.x = element_text(size = 12), 
        strip.placement = "outside", 
        strip.background = element_blank()) +
  facet_wrap(~morphology, scales = "free") +
  ggtitle("", subtitle = "dlPFC") +
  xlab("") + ylab(bquote('Mean number of\nmitochondria')) + 
  guides(shape="none",fill="none",color="none") + 
  scale_x_discrete(labels=c("Control", "Anesthesia")) +
  #scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=c("dodgerblue3","indianred1")) + 
  scale_color_manual(values=c("dodgerblue3","indianred1")) 

# Mitochondria Per Bouton 

CA1mitofreq <- 
  CA1mito_long %>% select(case,group,sex,bouton_frequency, number_per_bouton) %>% 
  group_by(case, group, sex, bouton_frequency) %>% 
  summarize(number = mean(number_per_bouton)) %>%
  mutate(counts = case_when(bouton_frequency == "mito_0_bouton" ~ "Zero",
                            bouton_frequency == "mito_1_bouton" ~ "One",
                            bouton_frequency == "mito_2_boutons" ~ "Two",
                            bouton_frequency == "mito_3plus_boutons" ~ "Three or more")) %>%
  mutate(counts = factor(counts, levels = c("Zero", "One", "Two", "Three or more"))) %>%
  select(-bouton_frequency) %>% ungroup()

CA1mitoboutongraph <- CA1mitofreq %>% 
  mutate(combined = paste(group,sex)) %>%
  mutate(combined = factor(combined, levels = c("control female",
                                                "control male",
                                                "anesthesia female",
                                                "anesthesia male"))) %>%
  ggplot(aes(x=interaction(combined,counts),y=number)) + 
  geom_jitter(aes(color=combined, 
                  shape=combined),
              position=position_jitter(seed=44456,width=0.2, height = 0), 
              size = 1.5) + 
  geom_bar(aes(fill=combined), 
           position = "dodge", alpha=0.3, stat="summary") + 
  theme_anes + 
  theme(strip.text.x = element_text(size = 12), 
        strip.placement = "outside", 
        strip.background = element_rect(color = NA)) +
  facet_wrap(vars(counts), scales = "free") +
  xlab("") + ylab(bquote('Mean number of mitochondria')) + 
  ggtitle("Mitochondria per Bouton in CA1") +
  guides(fill="none") +
  scale_x_discrete(labels=c("CF", "CM", "AF", "AM")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",
                             "indianred1", "indianred3")) + 
  scale_color_manual(name = "", values=c("dodgerblue3","dodgerblue4",
                                         "indianred1", "indianred3"),
                     labels=c("Control Female (CF)","Control Male (CM)",
                              "Anesthesia Female (AF)","Anesthesia Male (AM)")) +
  scale_shape_manual(name  = "", values=c(19, 17, 19, 17),
                     labels=c("Control Female (CF)", "Control Male (CM)",
                              "Anesthesia Female (AF)", "Anesthesia Male (AM)")) 

dlfreqfacets <- read.table(header=TRUE,
                           text = 
                             "counts ymin ymax
                                Zero 0 37", 
                           stringsAsFactors = FALSE)
dl_freq_facets <- with(dlfreqfacets,
                       data.frame(number=c(ymin,ymax),
                                  counts=c(counts, counts)))

dlmitoboutongraph <- dl_mitofreq %>%
  mutate(combined = paste(group,sex)) %>%
  mutate(combined = factor(combined, levels = c("control female",
                                                "control male",
                                                "anesthesia female",
                                                "anesthesia male"))) %>%
  ggplot(aes(x=interaction(combined,counts),y=number)) + 
  geom_jitter(aes(color=combined, 
                  shape=combined),
              position=position_jitter(seed=44457,width=0.2, height = 0), 
              size = 1.5) + 
  geom_bar(aes(fill=combined), 
           position = "dodge", alpha=0.3, stat="summary") + 
  theme_anes + 
  theme(strip.text.x = element_text(size = 12), 
        strip.placement = "outside", 
        strip.background = element_rect(color = NA)) +
  facet_wrap(~factor(counts, levels = c("Zero", "One", "Two", "Three or more")),
             scales = "free") +
  xlab("") + ylab(bquote('Mean number of mitochondria')) + 
  ggtitle("Mitochondria per Bouton in dlPFC") +
  guides(fill="none") +
  scale_x_discrete(labels=c("CF", "CM", "AF", "AM")) +
  scale_fill_manual(values=c("dodgerblue3","dodgerblue4",
                             "indianred1", "indianred3")) + 
  scale_color_manual(name = "", values=c("dodgerblue3","dodgerblue4",
                                         "indianred1", "indianred3"),
                     labels=c("Control Female (CF)","Control Male (CM)",
                              "Anesthesia Female (AF)","Anesthesia Male (AM)")) +
  scale_shape_manual(name  = "", values=c(19, 17, 19, 17),
                     labels=c("Control Female (CF)", "Control Male (CM)",
                              "Anesthesia Female (AF)", "Anesthesia Male (AM)")) +
  geom_point(data = dl_freq_facets, x=NA) +
  geom_signif(data = data.frame(counts = c("Zero")),
              aes(y_position=c(33), xmin = c(0.8), xmax = c(4.2), 
                  annotations = c("*")), tip_length = c(0), 
              textsize = c(7), vjust=.4, manual = T)

mito1cow <- cowplot::plot_grid(CA1mitodensgraph, CA1mitoshapegraph, 
                               dlmitodensgraph, dlmitoshapegraph, 
                               rel_widths = c(1.5, 3),
                               labels = c('A', 'B', 'C', 'D'))

mito2 <-(CA1mitoboutongraph | plot_spacer() | dlmitoboutongraph) +
  plot_layout(widths = c(1,.05,1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "bottom")

ggsave("Syn1.tiff", plot = syn1, width = 11.5, height = 13.88) # Figure 3
ggsave("SynTypesAll.tiff", syntypesall, width = 11.5, height = 13.88) # Figure 4
ggsave("Mito1.tiff", mito1cow, width = 11.5, height = 6.94) # Figure 5
ggsave("Mito2.tiff", mito2, width = 11.5, height = 6.94) # Figure 6


