#### Burtoni reproductive opportunity data
## CART 
# https://www.pluralsight.com/guides/explore-r-libraries:-rpart

### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_repro_opportunity/Synchronization/RNAseq/Data/")

#### load libraries ####
library(ggrepel)
library(tidyverse)

# CIT
library(party)
library(partykit)
library(ggparty)

#### load data ####
# load ME z score distance
load("../Figures/WGCNA/ME/MEs.dist.sum.sync.RData")

# load PCA z score dist
load('../Figures/DESEQ2/PCA/pca.dist.sum.sync.RData')

## combine and save
combined.dist.sum.sync = MEs.dist.sum.sync |>
  full_join(pca.dist.sum.sync |> 
              dplyr::rename(SimilarityRank.pca.F = SimilarityRank.F) |> 
              dplyr::rename(SimilarityRank.pca.M = SimilarityRank.M) |> 
              dplyr::select(Observation.id,
                            SimilarityZ.pca.F,
                            SimilarityZ.pca.M,
                            SimilarityRank.pca.F,
                            SimilarityRank.pca.M),
            by = 'Observation.id')  |> 
  relocate(all_of(c('Observation.id',
                    'SimilarityZ.pca.F',
                    'SimilarityZ.pca.M',
                    'SimilarityRank.pca.F',
                    'SimilarityRank.pca.M')))

# save
write.csv(combined.dist.sum.sync,
          '../Figures/combined.dist.sum.sync.csv',
          row.names = F)


#### CIT ME dist ####
### conditional inference tree
## use module eigengene z score distance

#https://stackoverflow.com/questions/29131254/how-to-generate-a-prediction-interval-from-a-regression-tree-rpart-object
#https://www.statmethods.net/advstats/cart.html
#https://martinschweinberger.github.io/TreesUBonn/conditional-inference-trees.html#example-2-prepositions
### Conditional inference trees via party
#"Conditional Inference Trees (CITs) are much better at determining the true effect of a predictor, i.e. the effect of a predictor if all other effects are simultaneously considered"
# dim(train); dim(test)

## rename variables
MEs.dist.sum.sync.rename = MEs.dist.sum.sync %>% 
  dplyr::rename(Time.together = Male.time.female.near.barrier) %>% 
  dplyr::rename(Estradiol = Estradiol_pg.mL.g.female) %>% 
  dplyr::rename(Testosterone = Testosterone_pg.mL.g.male) 

# females
tree_model.female.party = ctree(SimilarityZ.F ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                          data = MEs.dist.sum.sync.rename, 
                          control = ctree_control(minbucket = 5, 
                                                  minsplit = 10,
                                                  testtype = "Teststatistic",
                                                  mincriterion = 0.9))
                          
png('../Figures/WGCNA/ME/CIT/Female party model.png')
plot(tree_model.female.party, main="Female party model")
dev.off()

# plotting
p =
  ggparty(tree_model.female.party) +
  geom_edge() +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", 
                     y = SimilarityZ.F), 
                 coef = Inf, 
                 width = 1, 
                 fill = "#008080"),
    geom_jitter(aes(x = "", 
                    y = SimilarityZ.F),
                width = 0.1,
                height  = 0),
    theme_minimal(),
    xlab(""),
    ylab("Synchronization")
  ),
  shared_axis_labels = TRUE) +
  geom_node_label(mapping = aes(fill = splitvar),
    line_list = list(aes(label = splitvar)),
  line_gpar = list(list(size = 15)), 
  ids = "inner",
  color = 'white') +
  geom_node_label(aes(label = paste0("N = ", nodesize)),
                  ids = "terminal", nudge_y = -0.0, nudge_x = 0.01) +
  geom_edge_label(aes(label = paste(substr(breaks_label, start = 1, stop = 15))),
                  fill = 'white') + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  ggtitle('ME female') +
    scale_fill_manual(
      values = c("Court.total" = '#E69F00',
                 'Testosterone' = '#E69F00',
                 'Estradiol' = "#008080",
                 'Time.together' = 'grey25',
                 'Latency' = "#008080",
                 'GSI.male' = '#E69F00',
                 'GSI.female' = "#008080")) +
    theme(legend.position = "none") 
ggsave('../Figures/WGCNA/ME/CIT/Female party model present.pdf',
       p,
       height = 6.5,
       width = 6.5)


# males
tree_model.male.party = ctree(SimilarityZ.M ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                                data = MEs.dist.sum.sync.rename, 
                                control = ctree_control(minbucket = 5, 
                                                        minsplit = 10,
                                                        testtype = "Teststatistic",
                                                        mincriterion = 0.9)
) 

png('../Figures/WGCNA/ME/CIT/Male party model.png')
plot(tree_model.male.party, main="Male party model")
dev.off()

# plotting
p =
  ggparty(tree_model.male.party) +
  geom_edge() +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", 
                     y = SimilarityZ.M), 
                 coef = Inf, 
                 width = 1, 
                 fill = "#E69F00"),
    geom_jitter(aes(x = "", 
                    y = SimilarityZ.M),
                width = 0.1,
                height  = 0),
    theme_minimal(),
    xlab(""),
    ylab("Synchronization")
  ),
  shared_axis_labels = TRUE) +
  geom_node_label(mapping = aes(fill = splitvar),
                  line_list = list(aes(label = splitvar)),
                  line_gpar = list(list(size = 15)), 
                  ids = "inner",
                  color = 'white') +
  geom_node_label(aes(label = paste0("N = ", nodesize)),
                  ids = "terminal", nudge_y = -0.0, nudge_x = 0.01) +
  geom_edge_label(aes(label = paste(substr(breaks_label, start = 1, stop = 15))),
                  fill = 'white') + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  ggtitle('ME male') +
  scale_fill_manual(
    values = c("Court.total" = '#E69F00',
               'Testosterone' = '#E69F00',
               'Estradiol' = "#008080",
               'Time.together' = 'grey25',
               'Latency' = "#008080",
               'GSI.male' = '#E69F00',
               'GSI.female' = "#008080")) +
  theme(legend.position = "none") 
ggsave('../Figures/WGCNA/ME/CIT/Male party model present.pdf',
       p,
       height = 6.5,
       width = 6.5)


### calculate variable importance
## extract variable importance
variable.importance.ME.df = varimp(tree_model.female.party,
                                   nperm = 1000) %>%
  as.data.frame() %>%
  rownames_to_column('variable') %>%
  mutate(model = 'Female') %>%
  full_join(varimp(tree_model.male.party,
                   nperm = 1000) %>%
              as.data.frame()%>%
              rownames_to_column('variable') %>%
              mutate(model = 'Male')) %>%
  dplyr::rename(variable.importance = '.') %>%
  pivot_wider(names_from = model,
              values_from = variable.importance,
              names_prefix = "variable.importance.") %>%
  replace(is.na(.), 0)

# graph
variable.importance.ME.df %>%
  ggplot(aes(x = variable.importance.Male,
             y = variable.importance.Female,
             label = variable)) +
  geom_point() +
geom_label(hjust = "inward",
             vjust = "inward"
             ) +
  theme_classic() +
  theme(aspect.ratio = 1)+
  ggtitle('ME dist variable importance')
ggsave('../Figures/WGCNA/ME/CIT/female vs male variable importance.png',
       height = 10,
       width = 10)

# paper
variable.importance.ME.df %>%
  ggplot(aes(x = variable.importance.Male,
             y = variable.importance.Female,
             label = variable,
             fill = variable)) +
  geom_point() +
  geom_label(hjust = "inward",
             vjust = "inward",
             color = 'white') +
  theme_classic() +
  theme(aspect.ratio = 1)+
  ggtitle('ME Distance Variable Importance')  +
  scale_fill_manual(
    values = c("Court.total" = '#E69F00',
               'Testosterone' = '#E69F00',
               'Estradiol' = "#008080",
               'Time.together' = 'grey25',
               'Latency' = "#008080",
               'GSI.male' = '#E69F00',
               'GSI.female' = "#008080")) +
  theme(legend.position = "none")  +
  ylab('Female model variable importance') +
  xlab('Male model variable importance') 
ggsave('../Figures/WGCNA/ME/CIT/female vs male variable importance.pdf',
       height = 5,
       width = 5)







#### CIT PCA dist ####
### conditional inference tree
## use weighted PCA z score distance

## rename variables
pca.dist.sum.sync.rename = pca.dist.sum.sync %>% 
  dplyr::rename(Time.together = Male.time.female.near.barrier) %>% 
  dplyr::rename(Estradiol = Estradiol_pg.mL.g.female) %>% 
  dplyr::rename(Testosterone = Testosterone_pg.mL.g.male) 

# females
tree_model.female.party.pca = ctree(SimilarityZ.pca.F ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                                data = pca.dist.sum.sync.rename, 
                                control = ctree_control(minbucket = 5, 
                                                        minsplit = 10,
                                                        testtype = "Teststatistic",
                                                        mincriterion = 0.9))

png('../Figures/DESEQ2/PCA/CIT/Female party model.png')
plot(tree_model.female.party.pca, main="Female party model")
dev.off()

# plotting
p =
  ggparty(tree_model.female.party.pca) +
  geom_edge() +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", 
                     y = SimilarityZ.pca.F), 
                 coef = Inf, 
                 width = 1, 
                 fill = "#008080"),
    geom_jitter(aes(x = "", 
                    y = SimilarityZ.pca.F),
                width = 0.1,
                height  = 0),
    theme_minimal(),
    xlab(""),
    ylab("Synchronization")
  ),
  shared_axis_labels = TRUE) +
  geom_node_label(mapping = aes(fill = splitvar),
                  line_list = list(aes(label = splitvar)),
                  line_gpar = list(list(size = 15)), 
                  ids = "inner",
                  color = 'white') +
  geom_node_label(aes(label = paste0("N = ", nodesize)),
                  ids = "terminal", nudge_y = -0.0, nudge_x = 0.01) +
  geom_edge_label(aes(label = paste(substr(breaks_label, start = 1, stop = 15))),
                  fill = 'white') + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  ggtitle('PCA female') +
  scale_fill_manual(
    values = c("Court.total" = '#E69F00',
               'Testosterone' = '#E69F00',
               'Estradiol' = "#008080",
               'Time.together' = 'grey25',
               'Latency' = "#008080",
               'GSI.male' = '#E69F00',
               'GSI.female' = "#008080")) +
  theme(legend.position = "none") 
ggsave('../Figures/DESEQ2/PCA/CIT/Female party model present.pdf',
       p,
       height = 6.5,
       width = 6.5)



# males
tree_model.male.party.pca = ctree(SimilarityZ.pca.M ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                              data = pca.dist.sum.sync.rename, 
                              control = ctree_control(minbucket = 5, 
                                                      minsplit = 10,
                                                      testtype = "Teststatistic",
                                                      mincriterion = 0.9)
) 

png('../Figures/DESEQ2/PCA/CIT/Male party model.png')
plot(tree_model.male.party.pca, main="Male party model")
dev.off()

# plotting
p =
  ggparty(tree_model.male.party.pca) +
  geom_edge() +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", 
                     y = SimilarityZ.pca.M), 
                 coef = Inf, 
                 width = 1, 
                 fill = "#E69F00"),
    geom_jitter(aes(x = "", 
                    y = SimilarityZ.pca.M),
                width = 0.1,
                height  = 0),
    theme_minimal(),
    xlab(""),
    ylab("Synchronization")
  ),
  shared_axis_labels = TRUE) +
  geom_node_label(mapping = aes(fill = splitvar),
                  line_list = list(aes(label = splitvar)),
                  line_gpar = list(list(size = 15)), 
                  ids = "inner",
                  color = 'white') +
  geom_node_label(aes(label = paste0("N = ", nodesize)),
                  ids = "terminal", nudge_y = -0.0, nudge_x = 0.01) +
  geom_edge_label(aes(label = paste(substr(breaks_label, start = 1, stop = 15))),
                  fill = 'white') + 
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  ggtitle('PCA male') +
  scale_fill_manual(
    values = c("Court.total" = '#E69F00',
               'Testosterone' = '#E69F00',
               'Estradiol' = "#008080",
               'Time.together' = 'grey25',
               'Latency' = "#008080",
               'GSI.male' = '#E69F00',
               'GSI.female' = "#008080")) +
  theme(legend.position = "none") 
ggsave('../Figures/DESEQ2/PCA/CIT/Male party model present.pdf',
       p,
       height = 6.5,
       width = 6.5)


### calculate variable importance
## extract variable importance
variable.importance.pca.df = varimp(tree_model.female.party.pca,
                                    nperm = 1000) %>%
  as.data.frame() %>%
  rownames_to_column('variable') %>%
  mutate(model = 'Female') %>%
  full_join(varimp(tree_model.male.party.pca,
                   nperm = 1000) %>%
              as.data.frame()%>%
              rownames_to_column('variable') %>%
              mutate(model = 'Male')) %>%
  dplyr::rename(variable.importance = '.') %>%
  pivot_wider(names_from = model,
              values_from = variable.importance,
              names_prefix = "variable.importance.") %>%
  replace(is.na(.), 0)

# graph
variable.importance.pca.df %>%
  ggplot(aes(x = variable.importance.Male,
             y = variable.importance.Female,
             label = variable)) +
  geom_point() +
  geom_label(hjust = "inward",
             vjust = "inward"
  ) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggtitle('PCA dist variable importance')  
ggsave('../Figures/DESEQ2/PCA/CIT/female vs male variable importance.png',
       height = 10,
       width = 10)

# paper
variable.importance.pca.df %>%
  ggplot(aes(x = variable.importance.Male,
             y = variable.importance.Female,
             label = variable,
             fill = variable)) +
  geom_point() +
  geom_label(hjust = "inward",
             vjust = "inward",
             color= 'white'
  ) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggtitle('PCA Distance Variable Importance')  +
  scale_fill_manual(
    values = c("Court.total" = '#E69F00',
               'Testosterone' = '#E69F00',
               'Estradiol' = "#008080",
               'Time.together' = 'grey25',
               'Latency' = "#008080",
               'GSI.male' = '#E69F00',
               'GSI.female' = "#008080")) +
  theme(legend.position = "none")  +
  ylab('Female model variable importance') +
  xlab('Male model variable importance') 
ggsave('../Figures/DESEQ2/PCA/CIT/female vs male variable importance.pdf',
       height = 5,
       width = 5)







#### Compare CIT variable between ME and PCA ####
### Combine variable importance across models
variable.importance.all = variable.importance.ME.df |> 
  full_join(variable.importance.pca.df,
            by = 'variable',
            suffix = c('_ME',
                       '_pca'))

## convert NA to 0
variable.importance.all = variable.importance.all%>% 
  mutate(across(everything(), \(x) replace_na(x, 0)))

### graph
## males
variable.importance.all %>%
  ggplot(aes(x = variable.importance.Male_pca,
             y = variable.importance.Male_ME,
             label = variable)) +
  geom_point() +
  geom_label_repel(  xlim = c(0, Inf),
                     ylim = c(0, Inf)) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(x = 'PCA variable',
      y = 'ME variable') +
  ggtitle('Male Distance Variable Importance')
ggsave('../Figures/DESEQ2/PCA/CIT/Male variable importance pca vs ME.png',
       height = 10,
       width = 10)

# paper
variable.importance.all %>%
  ggplot(aes(x = variable.importance.Male_pca,
             y = variable.importance.Male_ME,
             label = variable,
             fill = variable)) +
  geom_point() +
  geom_label_repel(  xlim = c(0, Inf),
                     ylim = c(0, Inf),
                     color = 'white') +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(x = 'PCA variable',
       y = 'ME variable') +
  ggtitle('Male Distance Variable Importance') +
  scale_fill_manual(
    values = c("Court.total" = '#E69F00',
               'Testosterone' = '#E69F00',
               'Estradiol' = "#008080",
               'Time.together' = 'grey25',
               'Latency' = "#008080",
               'GSI.male' = '#E69F00',
               'GSI.female' = "#008080")) +
  theme(legend.position = "none")  +
  labs(x = 'PCA model variable importance',
       y = 'ME model variable importance') 
ggsave('../Figures/DESEQ2/PCA/CIT/Male variable importance pca vs ME.pdf',
       height = 5,
       width = 5)

## females
variable.importance.all %>%
  ggplot(aes(x = variable.importance.Female_pca,
             y = variable.importance.Female_ME,
             label = variable)) +
  geom_point() +
  geom_label_repel(  xlim = c(0, Inf),
                     ylim = c(0, Inf)) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(x = 'PCA variable',
       y = 'ME variable') +
  ggtitle('Female dist variable importance')
ggsave('../Figures/DESEQ2/PCA/CIT/Female variable importance pca vs ME.png',
       height = 10,
       width = 10)

# paper
variable.importance.all %>%
  ggplot(aes(x = variable.importance.Female_pca,
             y = variable.importance.Female_ME,
             label = variable,
             fill = variable)) +
  geom_point() +
  geom_label_repel(  xlim = c(0, Inf),
                     ylim = c(0, Inf),
                     color = 'white') +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(x = 'PCA variable',
       y = 'ME variable') +
  ggtitle('Female Distance Variable Importance') +
  scale_fill_manual(
    values = c("Court.total" = '#E69F00',
               'Testosterone' = '#E69F00',
               'Estradiol' = "#008080",
               'Time.together' = 'grey25',
               'Latency' = "#008080",
               'GSI.male' = '#E69F00',
               'GSI.female' = "#008080")) +
  theme(legend.position = "none")  +
  labs(x = 'PCA model variable importance',
       y = 'ME model variable importance') 
ggsave('../Figures/DESEQ2/PCA/CIT/Female variable importance pca vs ME.pdf',
       height = 5,
       width = 5)



