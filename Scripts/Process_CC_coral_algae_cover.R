#### Examining composition and nutrient gradient to create and hone MBIO 610 model
####  of species cover through 2 species competition/interaction


##### LOAD LIBRARIES #####
library(tidyverse)
library(here)
library(patchwork)


##### READ IN DATA #####
meta <- read_csv(here("Data", "Full_Metadata.csv"))
comp <- read_csv(here("Data", "Species_Abundances_wide.csv"))
fe <- read_csv(here("Data", "Species_FE.csv"))
traits <- read_csv(here("Data", "Distinct_Taxa.csv"))
chem <- read_csv(here("Data", "Nutrients_Processed_All.csv"))



##### PROCESS & VIEW DATA #####

# which species are most abundant?
comp %>%
  pivot_longer(cols = 2:52, names_to = "Taxa", values_to = "pcover") %>%
  filter(pcover > 10) %>%
  count(Taxa)
# assumption affirmed: porites and turf/turb are the most common top species across reef

# alternative option: coral vs macroalgae: 2 functional group competition rather than 'species'
sub_comp <- comp %>%
  pivot_longer(cols = 2:52, names_to = "Taxa", values_to = "pcover") %>%
  left_join(traits) %>%
  group_by(CowTagID, Taxon_Group, Calcification) %>%
  summarise(mean = mean(pcover),
            sum = sum(pcover))

coral <- sub_comp %>%
  filter(Taxon_Group == "Cnidaria" &
           Calcification == "Hermatypic") %>%
  mutate(functional_group = "Stony Coral")


# check which cnidaria are still high % cover: anemone at 16% at one site and
#  discosoma at 10% another side, but mostly 0's or below 6%
# sub_comp %>%
#   filter(Taxon_Group == "Cnidaria" &
#            Calcification != "Hermatypic")
# View(comp %>%
#   pivot_longer(cols = 2:52, names_to = "Taxa", values_to = "pcover") %>%
#   left_join(traits) %>%
#   filter(Taxon_Group == "Cnidaria" &
#            Calcification != "Hermatypic")
#   )

# check other important groups
sub_comp %>%
  filter(Taxon_Group != "Cnidaria" &
           Taxon_Group != "Cyanobacteria")

## porifera might still be important? so maybe 3 way competition model? coral, macroalgae, sponges?
# eh highest per-site cover is 12% and all else are below 5% so many don't need to include?
sub_comp %>%
  filter(Taxon_Group == "Porifera")


algae <- sub_comp %>%
  filter(Taxon_Group != "Cnidaria" &
           Taxon_Group != "Cyanobacteria" &
           Taxon_Group != "Porifera") %>%
  mutate(functional_group = "ma_turf")


## if I decide to model a three species/functional group competition interaction:
turf <- sub_comp %>%
  filter(Taxon_Group == "Turf") %>%
  mutate(functional_group = "Turf")
ma <- sub_comp %>%
  filter(Taxon_Group != "Cnidaria" &
           Taxon_Group != "Cyanobacteria" &
           Taxon_Group != "Porifera" &
           Taxon_Group != "Turf") %>%
  mutate(functional_group = "Macroalgae")


# combine coral and macroalgae/turf df
# TWO SPECIES
c.ma <- full_join(coral, algae) %>%
  ungroup() %>%
  select(-c(Taxon_Group, Calcification)) %>%
  group_by(CowTagID) %>%
  mutate(total_pcover = sum(sum)) %>%
  arrange(CowTagID) %>%
  ungroup() %>%
  group_by(CowTagID, functional_group) %>%
  mutate(sum_fg = sum(sum)) %>%
  distinct(CowTagID, functional_group, sum_fg, total_pcover) %>%
  ungroup()

# THREE SPECIES
c.ma.t <- full_join(coral, ma) %>%
  full_join(turf) %>%
  ungroup() %>%
  select(-c(Taxon_Group, Calcification)) %>%
  group_by(CowTagID) %>%
  mutate(total_pcover = sum(sum)) %>%
  arrange(CowTagID) %>%
  ungroup() %>%
  group_by(CowTagID, functional_group) %>%
  mutate(sum_fg = sum(sum)) %>%
  distinct(CowTagID, functional_group, sum_fg, total_pcover) %>%
  ungroup()


# Process chem data
cv.chem <- chem %>%
  select(-c(Season, Maximum:CVSeasonal)) %>%
  distinct() %>%
  pivot_wider(names_from = Parameters, values_from = CVAll)


# join comp with chem
# TWO SPECIES
chem.comp2 <- c.ma %>%
  left_join(cv.chem) %>%
  select(-Location)
chem.comp3 <- c.ma.t %>%
  left_join(cv.chem) %>%
  select(-Location)


##### VISUALIZE RELATIONSHIPS #####
# TWO SPECIES
plot1 <- chem.comp2 %>%
  mutate(sum_fg = sum_fg/100) %>% # proportional cover as 0 to 1.0
  filter(CowTagID != "VSEEP") %>%
  select(-c(Ammonia_umolL, HIX,VisibleHumidic_Like, M_C, Tryptophan_Like, Tyrosine_Like)) %>%
  pivot_longer(cols = Salinity:NN_umolL, names_to = "Parameters", values_to = "cv_values") %>%
  ggplot(aes(x = cv_values, y = sum_fg, group = functional_group, color = functional_group)) +
  geom_point(size = 1) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  labs(x = "Coefficient of variation of parameter",
       y = "Proportional benthic cover",
       color = "Functional Group") +
  scale_color_manual(values = c("green4", "orange"),
                     labels = c("Macroalgae+Turf", "Stony coral")) +
  facet_wrap(~Parameters, scales = "free_x", ncol = 3)
plot1
ggsave(here("Output", "pCover_coral_algae_v_SGD_param.png"), plot1, device = "png", width = 8, height = 8)

plot1b <- c.ma %>%
  mutate(sum_fg = sum_fg/100) %>%
  pivot_wider(names_from = functional_group, values_from = sum_fg) %>%
  ggplot(aes(x = ma_turf, y = `Stony Coral`)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Macroalgae+Turf") +
  geom_hline(yintercept = 0, linetype = "dashed")
plot1b
ggsave(here("Output", "pCover_coral_v_algae.png"), plot1b, device = "png", width = 5, height = 5)


# three SPECIES
plot2 <- chem.comp3 %>%
  mutate(sum_fg = sum_fg/100) %>% # proportional cover as 0 to 1.0
  filter(CowTagID != "VSEEP") %>%
  select(-c(M_C, Tryptophan_Like, Tyrosine_Like)) %>%
  pivot_longer(cols = Salinity:VisibleHumidic_Like, names_to = "Parameters", values_to = "cv_values") %>%
  ggplot(aes(x = cv_values, y = sum_fg, group = functional_group, color = functional_group)) +
  geom_point(size = 1) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  labs(x = "Coefficient of variation of parameter",
       y = "Proportional benthic cover",
       color = "Functional Group") +
  scale_color_manual(values = c("green4", "purple3", "yellow3"),
                     labels = c("Macroalgae", "Stony coral", "Turf")) +
  facet_wrap(~Parameters, scales = "free_x")
plot2
ggsave(here("Output", "pCover_coral_algae_turf_v_SGD_param.png"), plot2, device = "png", width = 8, height = 6)

plot2b <- c.ma.t %>%
  mutate(sum_fg = sum_fg/100) %>%
  pivot_wider(names_from = functional_group, values_from = sum_fg) %>%
  ggplot(aes(x = Turf, y = Macroalgae)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")
plot2c <- c.ma.t %>%
  mutate(sum_fg = sum_fg/100) %>%
  pivot_wider(names_from = functional_group, values_from = sum_fg) %>%
  ggplot(aes(x = Macroalgae, y = `Stony Coral`)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")
plot2d <- c.ma.t %>%
  mutate(sum_fg = sum_fg/100) %>%
  pivot_wider(names_from = functional_group, values_from = sum_fg) %>%
  ggplot(aes(x = Turf, y = `Stony Coral`)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")
plot2_patch <- (plot2b + plot2c) / plot2d
ggsave(here("Output", "pCover_coral_v_algae_v_turf.png"), plot2_patch, device = "png", width = 8, height = 8)




# relationship between stony coral and macroalgae+turf
# high correlation!
anova(lm(data = c.ma, sum_fg ~functional_group))
# p ~ 3e-11
summary(lm(data = c.ma, sum_fg ~functional_group))
# r2 ~ 0.69
lm.cma <- c.ma %>% pivot_wider(names_from = functional_group, values_from = sum_fg)
anova(lm(data = lm.cma, `Stony Coral` ~ ma_turf)) # p ~ 2e-9
summary(lm(data = lm.cma, `Stony Coral` ~ ma_turf)) # r2 = 0.87



# relationship between stony coral and macroalgae and turf
# high correlation!
anova(lm(data = c.ma.t, sum_fg ~functional_group))
# p ~ 4.7e4
summary(lm(data = c.ma.t, sum_fg ~functional_group))
# r2 ~ 0.23
lm.cmat <- c.ma.t %>% pivot_wider(names_from = functional_group, values_from = sum_fg)
anova(lm(data = lm.cmat, `Stony Coral` ~ Macroalgae)) # p ~ 0.002
summary(lm(data = lm.cmat, `Stony Coral` ~ Macroalgae)) # r2 = 0.42
#anova(lm(data = lm.cmat, `Stony Coral` ~ Turf)) # p > 0.05
anova(lm(data = lm.cmat, Macroalgae ~ Turf)) # p ~ 1e-6
summary(lm(data = lm.cmat, Macroalgae ~ Turf)) # r2 = 0.72



# save dataframe to use in python?
chem.comp.wide <- chem.comp2 %>%
  filter(CowTagID != "VSEEP") %>%
  mutate(functional_group = if_else(functional_group == "ma_turf", "Macroalgae+Turf", functional_group)) %>%
  select(-c(M_C, Tryptophan_Like, Tyrosine_Like)) %>%
  mutate(sum_fg = sum_fg/100,
         total_pcover = total_pcover/100) # proportional benthic cover from 0 to 1.0
#write_csv(chem.comp.wide, here("Data", "coral_algae_pcover_sgd_param.csv"))


