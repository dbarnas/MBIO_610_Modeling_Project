#### Examining composition and nutrient gradient to create and hone MBIO 610 model
####  of species cover through 2 species competition/interaction


##### LOAD LIBRARIES #####
library(tidyverse)
library(here)


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


# combine coral and macroalgae/turf df
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


# Process chem data
cv.chem <- chem %>%
  select(-c(Season, Maximum:CVSeasonal)) %>%
  distinct() %>%
  pivot_wider(names_from = Parameters, values_from = CVAll)


# join comp with chem
chem.comp <- c.ma %>%
  left_join(cv.chem) %>%
  select(-Location)


##### VISUALIZE RELATIONSHIPS #####
plot1 <- chem.comp %>%
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
  scale_color_manual(values = c("green4", "purple3"),
                     labels = c("Macroalgae+Turf", "Stony coral")) +
  facet_wrap(~Parameters, scales = "free")
plot1

ggsave(here("Output", "pCover_coral_algae_v_SGD_param.png"), plot1, device = "png", width = 8, height = 6)

# relationship between stony coral and macroalgae/turf
# high correlation!
anova(lm(data = c.ma, sum_fg ~functional_group))
# p ~ 3e-11
summary(lm(data = c.ma, sum_fg ~functional_group))
# r2 ~ 0.69

