library(ggtree)
load(here('data','dog_tree.RData')) #cladogram
my.names <- read_csv(here('data','tree_names.csv'))
breed.means$tip <- my.names$tree.name[match(breed.means$breed, my.names$data.name)]
to.drop <- c(setdiff(breed.means$tip, dog.tree$tip.label), setdiff(dog.tree$tip.label, breed.means$tip))
dog.tree <- drop.tip(dog.tree, to.drop)
dog.tree <- compute.brlen(dog.tree, 1)

dog_breed_tree <- ggtree(dog.tree, branch.length='none') +
  theme_tree2() +
  geom_tiplab(align=TRUE, linesize=.5)

dog_breed_mean <- prereg.data %>% 
  filter(condition != "w_2cp_vd" & condition != "odour") %>% 
  mutate(breed = factor(breed, levels = c("Australian_Shepherd", "Border_Collie", "Rhodesian_Ridgeback", "Labrador_Retriever", "Golden_Retriever", "White Swiss Shepherd Dog", "JP_Russell_Terrier", "Samoyed", "Siberian_Husky"))) %>% 
  drop_na() %>% 
  ggplot(aes(x = response, y = condition, color = condition)) +
  facet_wrap(~breed, ncol = 1, strip.position = "right") +
  stat_summary()

library(patchwork)
dog_breed_tree + dog_breed_mean
ggsave(here("graphs/dog_breed_means.png"))
