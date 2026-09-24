#### Burtoni reproductive opportunity data
## transcriptomic synchronization
# > R-4.6.3
## run on ccbbcomp01

### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_repro_opportunity/Synchronization/RNAseq/Data/")


#### load libraries ####
## installation
library(PerformanceAnalytics)
library(tidyverse)
library(distances)

#### load data ####
# gene data
load('vsd.RData')

# create dataframe
vsd.df = assay(vsd)
names(vsd.df) = rownames(vsd.df)

## pca data
# sample PC scores
load("data.behavior.comp.reduce.pca.RData")

# load pca vairance
load("percentVar.df.RData")

# reduce PCs to top 80%
data.behavior.comp.reduce.pca.PCs = data.behavior.comp.reduce.pca |> 
  dplyr::select(c('ID',
                   percentVar.df |> 
                     filter(sum <0.8) |>
                     pull(PC)))

## module eigengene
load('MEs.RData')

# create dataframe
MEs.df = t(MEs)
names(MEs.df) = rownames(MEs.df)

### sample IDs
sample.names = read.csv('../RNA_extraction_QC//DNA-RNA CSV Template_burtoni_repro_IMC.csv')
#replace '_' with '-'
sample.names = sample.names %>% 
  select(Sample.Name) %>% 
  mutate(Sample.Name = str_replace_all(Sample.Name,
                                       '_',
                                       '-'),
         ID = Sample.Name) %>% 
  separate(ID,
           c('Plate.position',
             'Round',
             'Tank',
             'Sex'))

sample.names.select = sample.names %>% 
  mutate(Sample.Name = paste(Plate.position,
                             Round,
                             Tank,
                             Sex,
                             sep = '.'),
         Tank.Round = paste(Tank,
                            Round,
                            sep = '.')) %>% 
  dplyr::select(c(Sample.Name,
                  Sex,
                  Tank.Round))

# get behavior data
# need data.behavior.comp from 'initial_behavior_script.R'
load('../../IMC_Rscripts/data.behavior.comp.RData')
# create reduced dataframe
data.behavior.comp.reduce = data.behavior.comp %>% 
  select(-c(Female.bower.male.not.present,
            Male.time.female.not.near.barrier,
            Male.time.female.not.present.bower,
            Male.time.female.present.bower,
            Male.total.time.in.female.bower,
            lead,
            `lateral display`,
            approach,
            Male.percent.near,
            Female.percent.near,
            Male.total,
            Male.total.time.near.barrier,
            Female.bower.total.bower.time,
            Female.bower.male.present)) %>% 
  dplyr::rename(Estradiol = Estradiol_pg.mL.g.female,
                Testosterone = Testosterone_pg.mL.g.male,
                Time.together = Male.time.female.near.barrier,
                Courtship = Court.total)




#### WGCNA ME distance matrix ####
## create distance matrix
MEs.dist = dist(scale(MEs)) %>% 
  as.matrix()


# add names to matrix
colnames(MEs.dist) = rownames(MEs)
rownames(MEs.dist) = rownames(MEs)

# convert dataframe to long format
MEs.dist = MEs.dist %>% 
  as.data.frame() %>% 
  rownames_to_column('from') %>% 
  pivot_longer(cols = contains('.'),
               names_to = 'to',
               values_to = 'dist') %>% 
  mutate(keep = ifelse(from == to,
                       0,
                       1)) %>% 
  filter(keep == 1) %>% 
  dplyr::select(-c(keep))

# add sample label
MEs.dist = MEs.dist %>%
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.'),
                     from.tank = paste(Round,
                                       Tank,
                                       sep = '.')) %>% 
              dplyr::rename('from' = 'Sample.Name') %>% 
              dplyr::rename('from.sex' = 'Sex') %>% 
              dplyr::select(-c(Round,
                               Tank,
                               Plate.position))) %>% 
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.'),
                     to.tank = paste(Round,
                                     Tank,
                                     sep = '.')) %>% 
              dplyr::rename('to' = 'Sample.Name') %>% 
              dplyr::rename('to.sex' = 'Sex') %>% 
              dplyr::select(-c(Round,
                               Tank,
                               Plate.position))) %>%
  mutate(same = ifelse(from.tank == to.tank,
                       'same',
                       'other'))

#### graph WGCNA ME distance matrix ####
## graph same vs. other
#poster
#females to males
MEs.dist %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>%
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  mutate(dist.same = ifelse(same == 'partner',
                                 dist,
                                 0.01)) %>% 
  group_by(from) %>% 
  mutate(dist.same = max(dist.same)) %>% 
  ungroup() %>% 
  mutate(from = as.numeric(as.factor(from))) %>% 
  ggplot(aes(y = dist,
             x = reorder(from,
                         -dist.same),
             color = same)) +
  geom_boxplot() +
  theme_classic(base_size = 25) +
  # ggtitle('Female rbo score compared to all males') +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Female') +
  ylab('ME Dist') +
  # labs(color = 'Comparison')
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank()) 
ggsave('../Figures/WGCNA/ME/Dist from females to all males paper.pdf',
       height = 5,
       width = 10)

#males to females
MEs.dist %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>%
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  mutate(dist.same = ifelse(same == 'partner',
                                 dist,
                                 0.01)) %>% 
  group_by(from) %>% 
  mutate(dist.same = max(dist.same)) %>% 
  ungroup() %>% 
  mutate(from = as.numeric(as.factor(from))) %>% 
  ggplot(aes(y = dist,
             x = reorder(from,
                         -dist.same),
             color = same)) +
  geom_boxplot() +
  theme_classic(base_size = 25) +
  # ggtitle('Male rbo score compared to all females') +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Male') +
  ylab('ME Dist') +
  # labs(color = 'Comparison')+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank()) 
ggsave('../Figures/WGCNA/ME/Dist from males to all females paper.pdf',
       height = 5,
       width = 10)


#females to males
#poster
MEs.dist %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  ggplot(aes(y = dist,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point( size = 3) +
  theme_classic() +
  scale_color_manual(values = c('black',
                                'red')) +
  ggtitle('Female relative similarity') + 
  theme(text = element_text(size = 20),
        legend.position = 'none') +
  xlab('Pair') +
  ylab('ME Dist')
ggsave('../Figures/WGCNA/ME/Dist from females to all males zscore all poster.pdf',
       height = 5.25,
       width = 5.25)

#males to females
#poster
MEs.dist %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>% 
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  ggplot(aes(y = dist,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point( size = 3) +
  theme_classic() +
  scale_color_manual(values = c('black',
                                'red')) +
  ggtitle('Male relative similarity') + 
  theme(text = element_text(size = 20),
        legend.position = 'none') +
  xlab('Pair') +
  ylab('ME Dist')
ggsave('../Figures/WGCNA/ME/Dist from males to all females zscore all poster.pdf',
       height = 5.25,
       width = 5.25)

### combine with behavior data
## male
MEs.dist.male = MEs.dist %>% 
  filter(same == 'same') %>% 
  filter(from.sex == 'M') %>% 
  dplyr::select(c(from.tank,
                  dist)) %>% 
  dplyr::rename(Tank = from.tank) %>% 
  separate(Tank,
           c('Round',
             'Tank')) %>% 
  mutate(Observation.id = paste(Tank,
                                Round,
                                sep = '_')) %>% 
  dplyr::select(-c(Round, Tank)) %>% 
  full_join(data.behavior.comp.reduce)
  
# check correlation
pdf('../Figures/WGCNA/ME/ME dist and behavior chart correlation male.pdf',
    height = 10,
    width = 10)
chart.Correlation(MEs.dist.male %>% 
                    select(-c(Observation.id)),
                  histogram = TRUE,
                  pch = 19)
dev.off()


## male
MEs.dist.female = MEs.dist %>% 
  filter(same == 'same') %>% 
  filter(from.sex == 'F') %>% 
  dplyr::select(c(from.tank,
                  dist)) %>% 
  dplyr::rename(Tank = from.tank) %>% 
  separate(Tank,
           c('Round',
             'Tank')) %>% 
  mutate(Observation.id = paste(Tank,
                                Round,
                                sep = '_')) %>% 
  dplyr::select(-c(Round, Tank)) %>% 
  full_join(data.behavior.comp.reduce)

# check correlation
pdf('../Figures/WGCNA/ME/ME dist and behavior chart correlation female.pdf',
    height = 10,
    width = 10)
chart.Correlation(MEs.dist.female %>% 
                    select(-c(Observation.id)),
                  histogram = TRUE,
                  pch = 19)
dev.off()






####  ME zscore distance ####
### calculate ME partner zscore and rank
## for z score remove partner from distribution 
MEs.dist.sum = MEs.dist %>%
  group_by(from,
           from.sex,
           to.sex) %>%
  mutate(
    Mean.null = mean(dist[same == "other"]),
    Sd.null   = sd(dist[same == "other"]),
    SimilarityZ = -(dist - Mean.null)/Sd.null,
    SimilarityRank = rank(dist,
                          ties.method = "average")
  ) %>%
  ungroup()

## combine with behavior data
MEs.dist.sum.sync = MEs.dist.sum %>%
  filter(same == "same") %>%
  separate(from.tank,
           into = c("Round","Tank")) %>%
  mutate(
    Observation.id = paste0(Tank,"_",Round)
  ) %>%
  select(
    -c(from,to,
       Round,Tank,
       to.tank,
       to.sex,
       same,
       dist,
       Mean.null,
       Sd.null)
  ) %>%
  pivot_wider(
    names_from = from.sex,
    values_from = c(SimilarityZ,
                    SimilarityRank),
    names_sep = "."
  ) %>%
  full_join(data.behavior.comp)

#### save ME distance data ####
save(MEs.dist.sum.sync,
     file = '../Figures/WGCNA/ME/MEs.dist.sum.sync.RData')

# load('../Figures/WGCNA/Rbo.score/MEs.dist.sum.sync.RData')

#### gene expression PCA distance ####
### calculate weighted PCA distance
### use weighted euclidean distance from "distances::distances" package
## only use top 80% variance PCs weighted
pca.dist = distances::distances(data.behavior.comp.reduce.pca.PCs |> 
                                  dplyr::select(-c(ID)) |> 
                                  as.matrix(),
                             normalize = "none",
                             weights = percentVar.df$percentVar[1:c(ncol(data.behavior.comp.reduce.pca.PCs) - 1)]) |> 
  distance_matrix() |> 
  broom::tidy()

# add all combinations
pca.dist = pca.dist |> 
  rbind(pca.dist |> 
          dplyr::rename(item.tmp = item1) |> 
          transmute(item1 = item2,
                    item2 = item.tmp,
                    distance= distance))

# add sample ID
sample.ID = data.behavior.comp.reduce.pca.PCs |> 
  dplyr::select(ID) |> 
  rownames_to_column('item') 

# rename to sample names
pca.dist = pca.dist |> 
  left_join(sample.ID |> 
              dplyr::rename(from = ID),
            by = c('item1' = 'item'))|> 
  left_join(sample.ID |> 
              dplyr::rename(to = ID),
            by = c('item2' = 'item'))

# add sample information
pca.dist = pca.dist %>%
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.'),
                     from.tank = paste(Round,
                                       Tank,
                                       sep = '.')) %>% 
              dplyr::rename('from' = 'Sample.Name') %>% 
              dplyr::rename('from.sex' = 'Sex') %>% 
              dplyr::select(-c(Round,
                               Tank,
                               Plate.position))) %>% 
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.'),
                     to.tank = paste(Round,
                                     Tank,
                                     sep = '.')) %>% 
              dplyr::rename('to' = 'Sample.Name') %>% 
              dplyr::rename('to.sex' = 'Sex') %>% 
              dplyr::select(-c(Round,
                               Tank,
                               Plate.position))) %>%
  mutate(same = ifelse(from.tank == to.tank,
                       'same',
                       'other')) |> 
  dplyr::rename(dist = distance) |> 
  dplyr::select(-c('item1',
                   'item2'))


#### graph PCA distance matrix ####
## graph same vs. other
#poster
#females to males
pca.dist %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>%
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  mutate(dist.same = ifelse(same == 'partner',
                            dist,
                            0.01)) %>% 
  group_by(from) %>% 
  mutate(dist.same = max(dist.same)) %>% 
  ungroup() %>% 
  mutate(from = as.numeric(as.factor(from))) %>% 
  ggplot(aes(y = dist,
             x = reorder(from,
                         -dist.same),
             color = same)) +
  geom_boxplot() +
  theme_classic(base_size = 25) +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Female') +
  ylab('PCA Dist') +
  # labs(color = 'Comparison')
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank()) 
ggsave('../Figures/DESEQ2/PCA/Dist from females to all males paper.pdf',
       height = 5,
       width = 10)

#males to females
pca.dist %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>%
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  mutate(dist.same = ifelse(same == 'partner',
                            dist,
                            0.01)) %>% 
  group_by(from) %>% 
  mutate(dist.same = max(dist.same)) %>% 
  ungroup() %>% 
  mutate(from = as.numeric(as.factor(from))) %>% 
  ggplot(aes(y = dist,
             x = reorder(from,
                         -dist.same),
             color = same)) +
  geom_boxplot() +
  theme_classic(base_size = 25) +
  # ggtitle('Male rbo score compared to all females') +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Male') +
  ylab('PCA Dist') +
  # labs(color = 'Comparison')+ 
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank()) 
ggsave('../Figures/DESEQ2/PCA/Dist from males to all females paper.pdf',
       height = 5,
       width = 10)


#females to males
#poster
pca.dist %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  ggplot(aes(y = dist,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point( size = 3) +
  theme_classic() +
  scale_color_manual(values = c('black',
                                'red')) +
  ggtitle('Female relative similarity') + 
  theme(text = element_text(size = 20),
        legend.position = 'none') +
  xlab('Pair') +
  ylab('PCA Dist')
ggsave('../Figures/DESEQ2/PCA/Dist from females to all males zscore all poster.pdf',
       height = 5.25,
       width = 5.25)

#males to females
#poster
pca.dist %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>% 
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  ggplot(aes(y = dist,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point( size = 3) +
  theme_classic() +
  scale_color_manual(values = c('black',
                                'red')) +
  ggtitle('Male relative similarity') + 
  theme(text = element_text(size = 20),
        legend.position = 'none') +
  xlab('Pair') +
  ylab('PCA Dist')
ggsave('../Figures/DESEQ2/PCA/Dist from males to all females zscore all poster.pdf',
       height = 5.25,
       width = 5.25)

### combine with behavior data
## male
pca.dist.male = pca.dist %>% 
  filter(same == 'same') %>% 
  filter(from.sex == 'M') %>% 
  dplyr::select(c(from.tank,
                  dist)) %>% 
  dplyr::rename(Tank = from.tank) %>% 
  separate(Tank,
           c('Round',
             'Tank')) %>% 
  mutate(Observation.id = paste(Tank,
                                Round,
                                sep = '_')) %>% 
  dplyr::select(-c(Round, Tank)) %>% 
  full_join(data.behavior.comp.reduce)

# check correlation
pdf('../Figures/DESEQ2/PCA/ME dist and behavior chart correlation male.pdf',
    height = 10,
    width = 10)
chart.Correlation(pca.dist.male %>% 
                    select(-c(Observation.id)),
                  histogram = TRUE,
                  pch = 19)
dev.off()


## female
pca.dist.female = pca.dist %>% 
  filter(same == 'same') %>% 
  filter(from.sex == 'F') %>% 
  dplyr::select(c(from.tank,
                  dist)) %>% 
  dplyr::rename(Tank = from.tank) %>% 
  separate(Tank,
           c('Round',
             'Tank')) %>% 
  mutate(Observation.id = paste(Tank,
                                Round,
                                sep = '_')) %>% 
  dplyr::select(-c(Round, Tank)) %>% 
  full_join(data.behavior.comp.reduce)

# check correlation
pdf('../Figures/DESEQ2/PCA/ME dist and behavior chart correlation female.pdf',
    height = 10,
    width = 10)
chart.Correlation(pca.dist.female %>% 
                    select(-c(Observation.id)),
                  histogram = TRUE,
                  pch = 19)
dev.off()






####  ME zscore distance ####
### calculate ME partner zscore and rank
## for z score remove partner from distribution 
pca.dist.sum = pca.dist %>%
  group_by(from,
           from.sex,
           to.sex) %>%
  mutate(
    Mean.null = mean(dist[same == "other"]),
    Sd.null   = sd(dist[same == "other"]),
    SimilarityZ.pca = -(dist - Mean.null)/Sd.null,
    SimilarityRank = rank(dist,
                          ties.method = "average")
  ) %>%
  ungroup()

## combine with behavior data
pca.dist.sum.sync = pca.dist.sum %>%
  filter(same == "same") %>%
  separate(from.tank,
           into = c("Round","Tank")) %>%
  mutate(
    Observation.id = paste0(Tank,"_",Round)
  ) %>%
  select(
    -c(from,to,
       Round,Tank,
       to.tank,
       to.sex,
       same,
       dist,
       Mean.null,
       Sd.null)
  ) %>%
  pivot_wider(
    names_from = from.sex,
    values_from = c(SimilarityZ.pca,
                    SimilarityRank),
    names_sep = "."
  ) %>% 
  full_join(data.behavior.comp)

## create same sex measure
# combine with behavior data
pca.dist.sum.sex.sync = pca.dist.sum %>%
  filter(from.sex == to.sex) %>%
  separate(from.tank,
           into = c("Round","Tank")) %>%
  separate(to.tank,
           into = c("Round.to","Tank.to")) %>%
  mutate(
    Observation.id = paste0(Tank,"_",Round),
    Observation.id.to = paste0(Tank.to,"_",Round.to)
  ) %>%
  select(
    -c(from,to,
       Round,Tank,
       Tank.to,
       Round.to,
       to.sex,
       same,
       dist,
       Mean.null,
       Sd.null)
  ) %>% head()
  pivot_wider(
    names_from = from.sex,
    values_from = c(SimilarityZ.pca,
                    SimilarityRank),
    names_sep = "."
  ) %>% 
  full_join(data.behavior.comp)

#### save ME distance data ####
save(pca.dist.sum.sync,
     file = '../Figures/DESEQ2/PCA/pca.dist.sum.sync.RData')

#### compare WGCNA ME dist to PCA dist ####
## combine PCA and ME dist

dist.comp.pca.MEs = pca.dist %>% 
  filter(same == 'same') %>% 
  filter(from.sex == 'F') %>% 
  dplyr::select(c(from.tank,
                  dist)) %>% 
  dplyr::rename(Tank = from.tank) %>% 
  separate(Tank,
           c('Round',
             'Tank')) %>% 
  mutate(Observation.id = paste(Tank,
                                Round,
                                sep = '_')) %>% 
  dplyr::select(-c(Round, Tank)) |> 
  dplyr::rename(dist.pca = dist)|> 
  full_join(MEs.dist %>% 
              filter(same == 'same') %>% 
              filter(from.sex == 'F') %>% 
              dplyr::select(c(from.tank,
                              dist)) %>% 
              dplyr::rename(Tank = from.tank) %>% 
              separate(Tank,
                       c('Round',
                         'Tank')) %>% 
              mutate(Observation.id = paste(Tank,
                                            Round,
                                            sep = '_')) %>% 
              dplyr::select(-c(Round, Tank)) |> 
              dplyr::rename(dist.ME = dist)) |> 
  full_join(pca.dist.sum.sync |> 
              dplyr::select(Observation.id,
                            SimilarityZ.pca.M,
                            SimilarityZ.pca.F ))|> 
  full_join(MEs.dist.sum.sync |> 
              dplyr::select(Observation.id,
                            SimilarityZ.M,
                            SimilarityZ.F ))

# check correlation
pdf('../Figures/DESEQ2/PCA/PCA dist and ME dist correlation.pdf',
    height = 10,
    width = 10)
chart.Correlation(dist.comp.pca.MEs %>% 
                    select(c(dist.pca,
                             dist.ME,
                             SimilarityZ.pca.M,
                             SimilarityZ.pca.F,
                             SimilarityZ.M,
                             SimilarityZ.F)),
                  histogram = TRUE,
                  pch = 19)
dev.off()

















