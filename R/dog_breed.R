library(tidyverse)
library(phytools)
library(here)
library(plotrix)
library(RColorBrewer)
library(ggtree)
library(patchwork)
library(gridExtra)
library(cowplot)

prereg.data <- read.csv(here('data','md1_data_vienna_long.csv')) %>%
  filter(preregistered_data=='after_prereg') #only data that were collected after preregistration

breed.means <- prereg.data %>%
  filter(condition != "odour", condition != "w_2cp_vd") %>% #only test conditions
  group_by(subject_ID, breed, condition) %>%
  summarise(mean.resp = mean(response)) %>%
  ungroup() %>%
  group_by(condition, breed) %>%
  summarise(mean = mean(mean.resp, na.rm = T), sem = std.error(mean.resp, na.rm = T)) %>%
  pivot_wider(id_cols = breed, names_from = condition, values_from = c(mean, sem))

load(here('data','dog_tree.RData')) #cladogram
my.names <- read_csv(here('data','tree_names.csv'))
breed.means$tip <- my.names$tree.name[match(breed.means$breed, my.names$data.name)]
to.drop <- c(setdiff(breed.means$tip, dog.tree$tip.label), setdiff(dog.tree$tip.label, breed.means$tip))
dog.tree <- drop.tip(dog.tree, to.drop)
#dog.tree <- compute.brlen(dog.tree, 1)

theme_set(theme_bw(18))

dog_breed_tree <- ggtree(dog.tree, branch.length='none') +
  theme_tree2() +
  geom_tiplab(align=TRUE, linesize=.5) + 
  ggplot2::xlim(0,9.5) +
  ggplot2::theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())

dog_breed_mean <- prereg.data %>% 
  filter(condition != "w_2cp_vd" & condition != "odour" & breed != "MIX") %>% 
  group_by(subject_ID, breed, condition) %>% summarise(mean = mean(response, na.rm = T)) %>% 
  mutate(breed = factor(breed, levels = c("Australian_Shepherd", "Border_Collie", "Rhodesian_Ridgeback", "Labrador_Retriever", "Golden_Retriever", "White Swiss Shepherd Dog", "JP_Russell_Terrier", "Samoyed", "Siberian_Husky"))) %>% 
  drop_na()

dog_breed_mean$condition[dog_breed_mean$condition == 'ost'] <- 'ostensive'
dog_breed_mean$condition[dog_breed_mean$condition == 'non'] <- 'non-ostensive'
 
cog_plot <- ggplot(dog_breed_mean, aes(x = mean, y = condition, color = condition)) +
  facet_wrap(~breed, ncol = 1, strip.position = "right") +
  geom_jitter(height = 0.1, size = 2, alpha = 0.5) +
  stat_summary(size = 1.05) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = 'red') +
  labs(y = "", x = "proportion correct") +
  scale_color_brewer(palette = 'Set2') +
  theme(axis.text.y = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
 
plot_grid(dog_breed_tree, cog_plot, align = 'hv', axis = 'bt')


ggsave(here("graphs/dog_breed_means.png"), width = 11, height = 8)
