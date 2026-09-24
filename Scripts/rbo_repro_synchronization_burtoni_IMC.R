#### Burtoni reproductive opportunity data
## rank biased overlap
# > R-4.6.3
## run on ccbbcomp01

### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_repro_opportunity/Synchronization/RNAseq/Data/")


#### load libraries ####
## installation
# library(BiocManager)
# BiocManager::install("prada",
#                      version = "3.10")
# install.packages("rlang")
# BiocManager::install("vsn")
# BiocManager::install("cellHTS2")
# # BiocManager::install("cellHTS2",
# #                      version = "3.3")
# BiocManager::install("gespeR")
# 
# library(cellHTS2)
# BiocManager::install("gespeR")

library("gespeR")
library(PerformanceAnalytics)
library(jtools)
library(tidyverse)


#### load data ####
## gene data
# load('vsd.RData')
# 
# # create dataframe
# vsd.df = assay(vsd)
# names(vsd.df) = rownames(vsd.df)

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



#### run rbo across all pairwise Module eigengenes####
#### Module eigengenes 
## create matrix of pairwise rbo scores 
ME.rbo.mat = sapply(1:NCOL(MEs.df), function(i) sapply(1:NCOL(MEs.df), function(j){
  rbo(MEs.df[,i],
      MEs.df[,j],
      p = 0.7,
      k = 5)
}))


# add names to matrix
colnames(ME.rbo.mat) = colnames(MEs.df)
rownames(ME.rbo.mat) = colnames(MEs.df)

# convert dataframe to long format
ME.rbo.df = ME.rbo.mat %>% 
  as.data.frame() %>% 
  rownames_to_column('from') %>% 
  pivot_longer(cols = contains('.'),
               names_to = 'to',
               values_to = 'rbo.score') %>% 
  mutate(keep = ifelse(from == to,
                       0,
                       1)) %>% 
  filter(keep == 1) %>% 
  dplyr::select(-c(keep))
  
# add sample label
ME.rbo.df = ME.rbo.df %>%
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

### create distance statistic
ME.rbo.df.sum = ME.rbo.df %>% 
  group_by(from,
           from.sex,
           to.sex) %>% 
  mutate(Mean = mean(rbo.score),
         Sd = sd(rbo.score)) %>% 
  ungroup() %>% 
  mutate(Zscore = (rbo.score - Mean)/Sd)


#### graph rbo matrix ME ####
## graph same vs. other
#males to females
ME.rbo.df %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = rbo.score,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from males to all females')
ggsave('../Figures/WGCNA/Rbo.score/Rbo from males to all females.png')

#females to males
ME.rbo.df %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = rbo.score,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Rbo from females to all males')
ggsave('../Figures/WGCNA/Rbo.score/Rbo from females to all males.png')

#poster
#females to males
ME.rbo.df %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>%
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  mutate(rbo.score.same = ifelse(same == 'partner',
                                 rbo.score,
                                 0.01)) %>% 
  group_by(from) %>% 
  mutate(rbo.score.same = max(rbo.score.same)) %>% 
  ungroup() %>% 
  mutate(from = as.numeric(as.factor(from))) %>% 
  ggplot(aes(y = rbo.score,
             x = reorder(from,
                         -rbo.score.same),
             color = same)) +
  geom_boxplot() +
  theme_classic(base_size = 30) +
  # ggtitle('Female rbo score compared to all males') +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Female') +
  ylab('rbo score') +
  # labs(color = 'Comparison')
  theme(legend.position = "none")
ggsave('../Figures/WGCNA/Rbo.score/Rbo from females to all males paper.png',
       height = 5,
       width = 10)

#males to females
ME.rbo.df %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>%
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  mutate(rbo.score.same = ifelse(same == 'partner',
                                 rbo.score,
                                 0.01)) %>% 
  group_by(from) %>% 
  mutate(rbo.score.same = max(rbo.score.same)) %>% 
  ungroup() %>% 
  mutate(from = as.numeric(as.factor(from))) %>% 
  ggplot(aes(y = rbo.score,
             x = reorder(from,
                         -rbo.score.same),
             color = same)) +
  geom_boxplot() +
  theme_classic(base_size = 30) +
  # ggtitle('Male rbo score compared to all females') +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Male') +
  ylab('rbo score') +
  # labs(color = 'Comparison')+ 
  theme(legend.position = "none")
ggsave('../Figures/WGCNA/Rbo.score/Rbo from males to all females paper.png',
       height = 5,
       width = 10)

### rbo.score stats 
### z score
## graph same vs. other
#males to females
ME.rbo.df.sum %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Rbo from males to all females zscore')
ggsave('../Figures/WGCNA/Rbo.score/Rbo from males to all females zscore.png')

#females to males
ME.rbo.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Rbo from females to all males zscore')
ggsave('../Figures/WGCNA/Rbo.score/Rbo from females to all males zscore.png')


#males to females
ME.rbo.df.sum %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Rbo from males to all females zscore') 
ggsave('../Figures/WGCNA/Rbo.score/Rbo from males to all females zscore all.png')

#females to males
ME.rbo.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Rbo from females to all males zscore') 
ggsave('../Figures/WGCNA/Rbo.score/Rbo from females to all males zscore all.png')

#females to males
#poster
ME.rbo.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  ggplot(aes(y = Zscore,
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
  ylab('Female rbo zscore')
ggsave('../Figures/WGCNA/Rbo.score/Rbo from females to all males zscore all poster.pdf',
       height = 5.25,
       width = 5.25)

#females to males
#poster
ME.rbo.df.sum %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>% 
  mutate(same = ifelse(same == 'same',
                       'partner',
                       same)) %>% 
  ggplot(aes(y = Zscore,
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
  ylab('Male rbo zscore')
ggsave('../Figures/WGCNA/Rbo.score/Rbo from males to all females zscore all poster.pdf',
       height = 5.25,
       width = 5.25)










#### rbo synchronization ME ####
### add behavior to ME correlation dataframe
ME.rbo.df.sum.sync = ME.rbo.df.sum %>% 
  filter(same == 'same')  %>% 
  separate(from.tank, 
           into = c("Round", "Tank")) %>%
  mutate(Observation.id = paste0(Tank,
                                 "_",
                                 Round)) %>% 
  select(-c(from,
            to,
            Round,
            Tank,
            to.tank,
            to.sex,
            same,
            Mean,
            Sd,
            rbo.score)) %>% 
  pivot_wider(names_from = from.sex,
              values_from = Zscore,
              names_prefix = 'Zscore.') %>% 
  full_join(data.behavior.comp)

### compare zscores
ME.rbo.df.sum.sync %>% 
  ggplot(aes(x = Zscore.M,
             y= Zscore.F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  ggtitle('ME rbo zscore comparison')
ggsave('../Figures/WGCNA/Rbo.score//ME rbo zscore comparison sexes.png',
       height = 10,
       width = 10)
## lm
Zscore.rbo.comp.lm = lm(Zscore.F ~ Zscore.M,
                          data=ME.rbo.df.sum.sync)

summary(Zscore.rbo.comp.lm)
anova(Zscore.rbo.comp.lm)

Zscore.rbo.comp.lm.pvalue = signif(summary(Zscore.rbo.comp.lm)[["coefficients"]]['Zscore.M','Pr(>|t|)'], digits = 1)

Zscore.rbo.comp.lm.rvalue = signif(sqrt(summary(Zscore.rbo.comp.lm)[["r.squared"]]), digits = 2)

# poster
effect_plot(Zscore.rbo.comp.lm,
            pred = Zscore.M,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=ME.rbo.df.sum.sync ,
             aes(x=Zscore.M,
                 y=Zscore.F),
             size = 3)+
  theme_classic()+
  theme(legend.position = 'null')+
  ylab('Zscore.F') +
  xlab('Zscore.M')+
  theme(text = element_text(size = 15)) +
  annotate("text", x=1, y=-1.5, label= paste0('p-value = ',
                                                Zscore.rbo.comp.lm.pvalue,
                                                ', R = ',
                                                Zscore.rbo.comp.lm.rvalue
  ))
ggsave('../Figures/WGCNA/Rbo.score/Comp corr zscore male vs female.png',
       height = 10,
       width = 10)

# check correlation
pdf('../Figures/WGCNA/Rbo.score/ME rbo zscore and behavior chart correlation.pdf',
    height = 10,
    width = 10)
chart.Correlation(ME.rbo.df.sum.sync %>% 
                    select(c(Zscore.M,
                             Zscore.F,
                             Male.time.female.near.barrier,
                             Court.total,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone_pg.mL.g.male,
                             Estradiol_pg.mL.g.female)),
                  histogram = TRUE,
                  pch = 19)
dev.off()

# spearman
pdf('../Figures/WGCNA/Rbo.score/ME rbo zscore and behavior chart correlation spearman.pdf',
    height = 10,
    width = 10)
chart.Correlation(ME.rbo.df.sum.sync %>% 
                    select(c(Zscore.M,
                             Zscore.F,
                             Male.time.female.near.barrier,
                             Court.total,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone_pg.mL.g.male,
                             Estradiol_pg.mL.g.female)),
                  histogram = TRUE,
                  pch = 19,
                  method = 'spearman')
dev.off()

### lm
## time together
Zscore.rbo.F.time.lm = lm(Zscore.F ~ Male.time.female.near.barrier,
                       data=ME.rbo.df.sum.sync)


summary(Zscore.rbo.F.time.lm)
anova(Zscore.rbo.F.time.lm)

Zscore.rbo.F.time.lm.pvalue = signif(summary(Zscore.rbo.F.time.lm)[["coefficients"]]['Male.time.female.near.barrier','Pr(>|t|)'], digits = 1)

Zscore.rbo.F.time.lm.rvalue = signif(sqrt(summary(Zscore.rbo.F.time.lm)[["r.squared"]]), digits = 2)

#graph glm with outlier on graph
# presentation
effect_plot(Zscore.rbo.F.time.lm,
            pred = Male.time.female.near.barrier,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=ME.rbo.df.sum.sync ,
             aes(x=Male.time.female.near.barrier,
                 y=Zscore.F),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  ylab('Zscore female rbo') +
  xlab('Time together')+
  theme(text = element_text(size = 30)) +
  annotate("text", x=450, y=-1.5, label= paste0('p-value = ',
                                                Zscore.rbo.F.time.lm.pvalue
  ))
ggsave('../Figures/WGCNA/Rbo.score/Female corr zscore vs time together.png',
       height = 10,
       width = 10)

# poster
effect_plot(Zscore.rbo.F.time.lm,
            pred = Male.time.female.near.barrier,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=ME.rbo.df.sum.sync ,
             aes(x=Male.time.female.near.barrier,
                 y=Zscore.F),
             size = 3)+
  theme_classic()+
  theme(legend.position = 'null')+
  ylab('Female rbo zscore') +
  xlab('Time together')+
  theme(text = element_text(size = 15)) +
  annotate("text", x=450, y=-1.5, label= paste0('p-value = ',
                                                Zscore.rbo.F.time.lm.pvalue,
                                                ', R = ',
                                                Zscore.rbo.F.time.lm.rvalue
  ))
ggsave('../Figures/WGCNA/Rbo.score/Female corr zscore vs time together.pdf',
       height = 5.25,
       width = 5.25)


## latency
Zscore.rbo.F.latency.lm = lm(Zscore.F ~ Latency,
                          data=ME.rbo.df.sum.sync)


summary(Zscore.rbo.F.latency.lm)
anova(Zscore.rbo.F.latency.lm)

Zscore.rbo.F.latency.lm.pvalue = signif(summary(Zscore.rbo.F.latency.lm)[["coefficients"]]['Latency','Pr(>|t|)'], digits = 1)

Zscore.rbo.F.latency.lm.rvalue = signif(sqrt(summary(Zscore.rbo.F.latency.lm)[["r.squared"]]), digits = 2)

#graph glm with outlier on graph
# presentation
effect_plot(Zscore.rbo.F.latency.lm,
            pred = Latency,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=ME.rbo.df.sum.sync ,
             aes(x=Latency,
                 y=Zscore.F),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  ylab('Zscore female rbo') +
  xlab('Latency')+
  theme(text = element_text(size = 30)) +
  annotate("text", x=450, y=-1.5, label= paste0('p-value = ',
                                                Zscore.rbo.F.latency.lm.pvalue
  ))
ggsave('../Figures/WGCNA/Rbo.score/Female corr zscore vs Latency.png',
       height = 10,
       width = 10)

# poster
effect_plot(Zscore.rbo.F.latency.lm,
            pred = Latency,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=ME.rbo.df.sum.sync ,
             aes(x=Latency,
                 y=Zscore.F),
             size = 3)+
  theme_classic()+
  theme(legend.position = 'null')+
  ylab('Female rbo zscore') +
  xlab('Latency')+
  theme(text = element_text(size = 15)) +
  annotate("text", x=750, y=-1.5, label= paste0('p-value = ',
                                                Zscore.rbo.F.latency.lm.pvalue,
                                                ', R = ',
                                                Zscore.rbo.F.latency.lm.rvalue
  ))
ggsave('../Figures/WGCNA/Rbo.score/Female corr zscore vs Latency.pdf',
       height = 5.25,
       width = 5.25)


## graph residuals
# check correlation
pdf('../Figures/WGCNA/Rbo.score/ME rbo zscore and behavior and female residuals chart correlation.pdf',
    height = 10,
    width = 10)
chart.Correlation(ME.rbo.df.sum.sync %>% 
                    mutate(Time.residual = Zscore.rbo.F.time.lm$residuals,
                           Latency.residual = Zscore.rbo.F.latency.lm$residuals) %>%
                    select(c(Time.residual,
                             Latency.residual,
                             Zscore.M,
                             Zscore.F,
                             Male.time.female.near.barrier,
                             Court.total,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone_pg.mL.g.male,
                             Estradiol_pg.mL.g.female)),
                  histogram = TRUE,
                  pch = 19)
dev.off()




### males lm
## time together
Zscore.rbo.M.time.lm = lm(Zscore.M ~ Male.time.female.near.barrier,
                          data=ME.rbo.df.sum.sync)


summary(Zscore.rbo.M.time.lm)
anova(Zscore.rbo.M.time.lm)

Zscore.rbo.M.time.lm.pvalue = signif(summary(Zscore.rbo.M.time.lm)[["coefficients"]]['Male.time.female.near.barrier','Pr(>|t|)'], digits = 1)

Zscore.rbo.M.time.lm.rvalue = signif(sqrt(summary(Zscore.rbo.M.time.lm)[["r.squared"]]), digits = 2)

#graph glm with outlier on graph
# presentation
effect_plot(Zscore.rbo.M.time.lm,
            pred = Male.time.female.near.barrier,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=ME.rbo.df.sum.sync ,
             aes(x=Male.time.female.near.barrier,
                 y=Zscore.M),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  ylab('Zscore male rbo') +
  xlab('Time together')+
  theme(text = element_text(size = 30)) +
  annotate("text", x=450, y=-1.5, label= paste0('p-value = ',
                                                Zscore.rbo.M.time.lm.pvalue
  ))
ggsave('../Figures/WGCNA/Rbo.score/Male corr zscore vs time together.png',
       height = 10,
       width = 10)

# poster
effect_plot(Zscore.rbo.M.time.lm,
            pred = Male.time.female.near.barrier,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=ME.rbo.df.sum.sync ,
             aes(x=Male.time.female.near.barrier,
                 y=Zscore.M),
             size = 3)+
  theme_classic()+
  theme(legend.position = 'null')+
  ylab('Male rbo zscore') +
  xlab('Time together')+
  theme(text = element_text(size = 15)) +
  annotate("text", x=450, y=-1, label= paste0('p-value = ',
                                                Zscore.rbo.M.time.lm.pvalue,
                                                ', R = ',
                                                Zscore.rbo.M.time.lm.rvalue
  ))
ggsave('../Figures/WGCNA/Rbo.score/Male corr zscore vs time together.pdf',
       height = 5.25,
       width = 5.25)


## court
Zscore.rbo.M.court.lm = lm(Zscore.M ~ Court.total,
                             data=ME.rbo.df.sum.sync)


summary(Zscore.rbo.M.court.lm)
anova(Zscore.rbo.M.court.lm)

Zscore.rbo.M.court.lm.pvalue = signif(summary(Zscore.rbo.M.court.lm)[["coefficients"]]['Court.total','Pr(>|t|)'], digits = 1)

Zscore.rbo.M.court.lm.rvalue = signif(sqrt(summary(Zscore.rbo.M.court.lm)[["r.squared"]]), digits = 2)

#graph glm with outlier on graph
# presentation
effect_plot(Zscore.rbo.M.court.lm,
            pred = Court.total,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=ME.rbo.df.sum.sync ,
             aes(x=Court.total,
                 y=Zscore.M),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  ylab('Zscore male rbo') +
  xlab('Courtship')+
  theme(text = element_text(size = 30)) +
  annotate("text", x=20, y=-1, label= paste0('p-value = ',
                                                Zscore.rbo.M.court.lm.pvalue
  ))
ggsave('../Figures/WGCNA/Rbo.score/Male corr zscore vs Courtship.png',
       height = 10,
       width = 10)

# poster
effect_plot(Zscore.rbo.M.court.lm,
            pred = Court.total,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=ME.rbo.df.sum.sync ,
             aes(x=Court.total,
                 y=Zscore.M),
             size = 3)+
  theme_classic()+
  theme(legend.position = 'null')+
  ylab('Male rbo zscore') +
  xlab('Courtship')+
  theme(text = element_text(size = 15)) +
  annotate("text", x=20, y=-1, label= paste0('p-value = ',
                                                Zscore.rbo.M.court.lm.pvalue,
                                                ', R = ',
                                                Zscore.rbo.M.court.lm.rvalue
  ))
ggsave('../Figures/WGCNA/Rbo.score/Memale corr zscore vs courtship.pdf',
       height = 5.25,
       width = 5.25)

## graph residuals
# check correlation
pdf('../Figures/WGCNA/Rbo.score/ME rbo zscore and behavior and male residuals chart correlation.pdf',
    height = 10,
    width = 10)
chart.Correlation(ME.rbo.df.sum.sync %>% 
                    mutate(Court.residual = Zscore.rbo.M.court.lm$residuals,
                           Time.residual = Zscore.rbo.M.time.lm$residuals) %>% 
                    select(c(Time.residual,
                             Court.residual,
                             Zscore.M,
                             Zscore.F,
                             Male.time.female.near.barrier,
                             Court.total,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone_pg.mL.g.male,
                             Estradiol_pg.mL.g.female)),
                  histogram = TRUE,
                  pch = 19)
dev.off()

#### use rbo score vs MEGs ####
### add behavior to ME correlation dataframe
ME.rbo.df.sync = ME.rbo.df.sum %>% 
  filter(same == 'same')  %>% 
  separate(from.tank, 
           into = c("Round", "Tank")) %>%
  mutate(Observation.id = paste0(Tank,
                                 "_",
                                 Round)) %>% 
  select(-c(from,
            to,
            Round,
            Tank,
            to.tank,
            to.sex,
            same,
            Mean,
            Sd,
            Zscore
            )) %>% 
  pivot_wider(names_from = from.sex,
              values_from = rbo.score,
              names_prefix = 'rbo.score.') %>% 
  full_join(data.behavior.comp)

## compare zscores
ME.rbo.df.sync %>% 
  ggplot(aes(x = rbo.score.M,
             y= rbo.score.F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  ggtitle('ME rbo comparison')
ggsave('../Figures/WGCNA/Rbo.score/ME rbo comparison sexes.png',
       height = 10,
       width = 10)

# check correlation
pdf('../Figures/WGCNA/Rbo.score/ME rbo and behavior chart correlation.pdf',
    height = 10,
    width = 10)
chart.Correlation(ME.rbo.df.sync %>% 
                    select(c(rbo.score.M,
                             rbo.score.F,
                             Male.time.female.near.barrier,
                             Court.total,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone_pg.mL.g.male,
                             Estradiol_pg.mL.g.female)),
                  histogram = TRUE,
                  pch = 19)
dev.off()



#### use rbo score 
### compare to MEGs
ME.rbo.df.sync.MEGs = ME.rbo.df.sum %>% 
  filter(same == 'same')  %>% 
  full_join(MEs %>% 
              rownames_to_column('from')) 

##check correlation
# all
pdf('../Figures/WGCNA/Rbo.score/ME rbo and MEGs chart correlation all.pdf',
    height = 10,
    width = 10)
chart.Correlation(ME.rbo.df.sync.MEGs %>% 
                    filter(from.sex == 'F') %>% 
                    select(-c("from",
                              "to",
                              "rbo.score",
                              "from.sex",
                              "from.tank" ,
                              "to.sex" ,
                              "to.tank",
                              "same",
                              "Mean",
                              "Sd")),
                  histogram = TRUE,
                  pch = 19)
dev.off()

#males
pdf('../Figures/WGCNA/Rbo.score/ME rbo and MEGs chart correlation males.pdf',
    height = 10,
    width = 10)
chart.Correlation(ME.rbo.df.sync.MEGs %>% 
                    filter(from.sex == 'M') %>% 
                    select(-c("from",
                              "to",
                              "rbo.score",
                              "from.sex",
                              "from.tank" ,
                              "to.sex" ,
                              "to.tank",
                              "same",
                              "Mean",
                              "Sd")),
                  histogram = TRUE,
                  pch = 19)
dev.off()

#females
pdf('../Figures/WGCNA/Rbo.score/ME rbo and MEGs chart correlation females.pdf',
    height = 10,
    width = 10)
chart.Correlation(ME.rbo.df.sync.MEGs %>% 
                    filter(from.sex == 'F') %>% 
                    select(-c("from",
                              "to",
                              "rbo.score",
                              "from.sex",
                              "from.tank" ,
                              "to.sex" ,
                              "to.tank",
                              "same",
                              "Mean",
                              "Sd")),
                  histogram = TRUE,
                  pch = 19)
dev.off()


#### create heatmap of comparison ####
library(pheatmap)
library(WGCNA)
##males vs females
#make male dataframe
ME.rbo.df.sync.MEGs.males = ME.rbo.df.sum %>% 
  filter(same == 'same')  %>% 
  full_join(MEs %>% 
              rownames_to_column('from')) %>% 
  filter(from.sex == 'M') %>% 
  select(-c("from",
            "to",
            "rbo.score",
            "from.sex",
            "from.tank" ,
            "to.sex" ,
            "to.tank",
            "same",
            "Mean",
            "Sd")) %>% 
  rename_all(~ paste(., "M", sep = "."))
#make female dataframe
ME.rbo.df.sync.MEGs.females = ME.rbo.df.sum %>% 
  filter(same == 'same')  %>% 
  full_join(MEs %>% 
              rownames_to_column('from')) %>% 
  filter(from.sex == 'F') %>% 
  select(-c("from",
            "to",
            "rbo.score",
            "from.sex",
            "from.tank" ,
            "to.sex" ,
            "to.tank",
            "same",
            "Mean",
            "Sd")) %>% 
  rename_all(~ paste(., "F", sep = "."))

# pairs trait correlation
moduleTraitCor.ME.rbo.MEGs = cor(ME.rbo.df.sync.MEGs.males, 
                                   ME.rbo.df.sync.MEGs.females, 
                                   use = "p")
moduleTraitPvalue.ME.rbo.MEGs = corPvalueStudent(moduleTraitCor.ME.rbo.MEGs, 
                                                   nrow(ME.rbo.df.sync.MEGs.females))

#pairs
#create empty boxes
moduleTraitPvalue.ME.rbo.MEGs.empty = moduleTraitPvalue.ME.rbo.MEGs
moduleTraitPvalue.ME.rbo.MEGs.empty = signif(moduleTraitPvalue.ME.rbo.MEGs.empty, 
                                               1)
moduleTraitPvalue.ME.rbo.MEGs.empty[moduleTraitPvalue.ME.rbo.MEGs.empty > .1] <- ''

#create star 
moduleTraitPvalue.ME.rbo.MEGs.empty.star = moduleTraitPvalue.ME.rbo.MEGs
moduleTraitPvalue.ME.rbo.MEGs.empty.star = signif(moduleTraitPvalue.ME.rbo.MEGs.empty.star, 
                                                    1)
moduleTraitPvalue.ME.rbo.MEGs.empty.star[moduleTraitPvalue.ME.rbo.MEGs.empty.star <= .05] <- '*'
moduleTraitPvalue.ME.rbo.MEGs.empty.star[moduleTraitPvalue.ME.rbo.MEGs.empty.star > .05] <- ''


# poster
pheatmap(moduleTraitCor.ME.rbo.MEGs, 
         # clustering_distance_rows = 'euclidean',
         # clustering_distance_cols = 'euclidean',
         # clustering_method = 'ward.D2',
         cluster_rows = F,
         cluster_cols = F,
         scale = 'none',
         border_color = 'black',
         # cutree_rows = 3,
         # cutree_cols = 3,
         legend = T,
         treeheight_col = 25,
         treeheight_row = 25,
         angle_col = 45,
         display_numbers = moduleTraitPvalue.ME.rbo.MEGs.empty,
         color = blueWhiteRed(50),
         main = 'Comparison between partners',
         fontsize = 15,
         fontsize_number = 20,
         breaks=seq(-0.6, 0.6, length.out=51),
         filename = "../Figures/WGCNA/Rbo.score/Pairs rbo score and MEGs heatmap.pdf",
         width = 10,
         height = 10
)

# ### reorder columns and rows with pvclust
# ## create pvclust
# # males
# male.module.behavior.pvclust = pvclust(moduleTraitCor.tank.behavior,
#                                        method.hclust="ward.D2",
#                                        method.dist="euclidean",
#                                        use.cor="pairwise.complete.obs",
#                                        nboot=1000,
#                                        parallel=T,
#                                        quiet=FALSE)
# 
# #graph
# png(filename = "../Figures/WGCNA/heatmaps/Tank trait and behavior pvclust.png")
# plot(male.module.behavior.pvclust)
# pvrect(male.module.behavior.pvclust, alpha=0.8,  max.only = T)
# dev.off()
# 
# # females
# female.module.behavior.pvclust = pvclust(moduleTraitCor.tank.behavior %>%
#                                            t(),
#                                          method.hclust="ward.D2",
#                                          method.dist="euclidean",
#                                          use.cor="pairwise.complete.obs",
#                                          nboot=1000,
#                                          parallel=T,
#                                          quiet=FALSE)
# 
# #graph
# png(filename = "../Figures/WGCNA/heatmaps/Female trait and behavior pvclust.png")
# plot(female.module.behavior.pvclust)
# pvrect(female.module.behavior.pvclust, alpha=0.8)
# dev.off()


#### save data ####
save(ME.rbo.df.sum.sync,
     file = '../Figures/WGCNA/Rbo.score/ME.rbo.df.sum.sync.RData')

# load('../Figures/WGCNA/Rbo.score/ME.rbo.df.sum.sync.RData')

#### run rbo across all pairwise vsd ####
#### Use 50% most variable genes 
## create matrix of pairwise rbo scores 
# look at the top 1000 genes (gives similar results to default 3336 but runs faster)
vsd.rbo.mat = sapply(1:NCOL(vsd.df), function(i) sapply(1:NCOL(vsd.df), function(j){
  rbo(vsd.df[,i],
      vsd.df[,j],
      p = 1,
      k = 1000
      )
}))

# add names to matrix
colnames(vsd.rbo.mat) = colnames(vsd.df)
rownames(vsd.rbo.mat) = colnames(vsd.df)

# convert dataframe to long format
vsd.rbo.df = vsd.rbo.mat %>% 
  as.data.frame() %>% 
  rownames_to_column('from') %>% 
  pivot_longer(cols = contains('.'),
               names_to = 'to',
               values_to = 'rbo.score') %>% 
  mutate(keep = ifelse(from == to,
                       0,
                       1)) %>% 
  filter(keep == 1) %>% 
  dplyr::select(-c(keep))

# add sample label
vsd.rbo.df = vsd.rbo.df %>%
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

### create distance statistic
vsd.rbo.df.sum = vsd.rbo.df %>% 
  group_by(from,
           from.sex,
           to.sex) %>% 
  mutate(Mean = mean(rbo.score),
         Sd = sd(rbo.score)) %>% 
  ungroup() %>% 
  mutate(Zscore = (rbo.score - Mean)/Sd)


#### graph rbo matrix vsd ####
## graph same vs. other
#males to females
vsd.rbo.df %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = rbo.score,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from males to all females')
ggsave('../Figures/DESEQ2/Rbo/Rbo.score/Rbo from males to all females.png')

#females to males
vsd.rbo.df %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = rbo.score,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Rbo from females to all males')
ggsave('../Figures/DESEQ2/Rbo/Rbo.score/Rbo from females to all males.png')

### rbo.scoreelation stats 
### z score
## graph same vs. other
#males to females
vsd.rbo.df.sum %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Rbo from males to all females zscore')
ggsave('../Figures/DESEQ2/Rbo/Rbo.score/Rbo from males to all females zscore.png')

#females to males
vsd.rbo.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Rbo from females to all males zscore')
ggsave('../Figures/DESEQ2/Rbo/Rbo.score/Rbo from females to all males zscore.png')


#males to females
vsd.rbo.df.sum %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Rbo from males to all females zscore') 
ggsave('../Figures/DESEQ2/Rbo/Rbo.score/Rbo from females to all males zscore all.png')

#females to males
vsd.rbo.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Rbo from females to all males zscore') 
ggsave('../Figures/DESEQ2/Rbo/Rbo.score/Rbo from females to all males zscore all.png')

#females to males
#presentation
vsd.rbo.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point( size = 5) +
  theme_classic() +
  scale_color_manual(values = c('black',
                                'red')) +
  ggtitle('Rbo score females') + 
  theme(text = element_text(size = 30),
        legend.position = 'none') +
  xlab('Pair') +
  ylab('Zscore female rbo.scoreelation')
ggsave('../Figures/DESEQ2/Rbo/Rbo.score/Rbo from females to all males zscore all presentation.png',
       height = 10,
       width = 10)










#### rbo synchronization vsd ####
### add behavior to vsd correlation dataframe
vsd.rbo.df.sum.sync = vsd.rbo.df.sum %>% 
  filter(same == 'same')  %>% 
  separate(from.tank, 
           into = c("Round", "Tank")) %>%
  mutate(Observation.id = paste0(Tank,
                                 "_",
                                 Round)) %>% 
  select(-c(from,
            to,
            Round,
            Tank,
            to.tank,
            to.sex,
            same,
            Mean,
            Sd,
            rbo.score)) %>% 
  pivot_wider(names_from = from.sex,
              values_from = Zscore,
              names_prefix = 'Zscore.') %>% 
  full_join(data.behavior.comp)

## compare zscores
vsd.rbo.df.sum.sync %>% 
  ggplot(aes(x = Zscore.M,
             y= Zscore.F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  ggtitle('vsd rbo zscore comparison')
ggsave('../Figures/DESEQ2/Rbo/Rbo.score//vsd rbo zscore comparison sexes.png',
       height = 10,
       width = 10)

# check correlation
pdf('../Figures/DESEQ2/Rbo/Rbo.score/vsd rbo zscore and behavior chart correlation.pdf',
    height = 10,
    width = 10)
chart.Correlation(vsd.rbo.df.sum.sync %>% 
                    select(c(Zscore.M,
                             Zscore.F,
                             Male.time.female.near.barrier,
                             Court.total,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone_pg.mL.g.male,
                             Estradiol_pg.mL.g.female)),
                  histogram = TRUE,
                  pch = 19)
dev.off()

### lm
## time together
Zscore.rbo.F.time.lm = lm(Zscore.F ~ Male.time.female.near.barrier,
                          data=vsd.rbo.df.sum.sync)


summary(Zscore.rbo.F.time.lm)
anova(Zscore.rbo.F.time.lm)

Zscore.rbo.F.time.lm.pvalue = signif(summary(Zscore.rbo.F.time.lm)[["coefficients"]]['Male.time.female.near.barrier','Pr(>|t|)'], digits = 3)

#graph glm with outlier on graph
# presentation
effect_plot(Zscore.rbo.F.time.lm,
            pred = Male.time.female.near.barrier,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=vsd.rbo.df.sum.sync ,
             aes(x=Male.time.female.near.barrier,
                 y=Zscore.F),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  ylab('Zscore female rbo') +
  xlab('Time together')+
  theme(text = element_text(size = 30)) +
  annotate("text", x=450, y=-1.5, label= paste0('p-value = ',
                                                Zscore.rbo.F.time.lm.pvalue
  ))
ggsave('../Figures/DESEQ2/Rbo/Rbo.score/Female corr zscore vs time together.png',
       height = 10,
       width = 10)


## latency
Zscore.rbo.F.latency.lm = lm(Zscore.F ~ Latency,
                             data=vsd.rbo.df.sum.sync)


summary(Zscore.rbo.F.latency.lm)
anova(Zscore.rbo.F.latency.lm)

Zscore.rbo.F.latency.lm.pvalue = signif(summary(Zscore.rbo.F.latency.lm)[["coefficients"]]['Latency','Pr(>|t|)'], digits = 3)

#graph glm with outlier on graph
# presentation
effect_plot(Zscore.rbo.F.latency.lm,
            pred = Latency,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=vsd.rbo.df.sum.sync ,
             aes(x=Latency,
                 y=Zscore.F),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  ylab('Zscore female rbo') +
  xlab('Latency')+
  theme(text = element_text(size = 30)) +
  annotate("text", x=450, y=-1.5, label= paste0('p-value = ',
                                                Zscore.rbo.F.latency.lm.pvalue
  ))
ggsave('../Figures/DESEQ2/Rbo/Rbo.score/Female corr zscore vs Latency.png',
       height = 10,
       width = 10)

#### use rbo score
### add behavior to ME correlation dataframe
vsd.rbo.df.sync = vsd.rbo.df.sum %>% 
  filter(same == 'same')  %>% 
  separate(from.tank, 
           into = c("Round", "Tank")) %>%
  mutate(Observation.id = paste0(Tank,
                                 "_",
                                 Round)) %>% 
  select(-c(from,
            to,
            Round,
            Tank,
            to.tank,
            to.sex,
            same,
            Mean,
            Sd,
            Zscore
  )) %>% 
  pivot_wider(names_from = from.sex,
              values_from = rbo.score,
              names_prefix = 'rbo.score.') %>% 
  full_join(data.behavior.comp)

# check correlation
pdf('../Figures/DESEQ2/Rbo/Rbo.score.p1.k1000/vsd rbo and behavior chart correlation.pdf',
    height = 10,
    width = 10)
chart.Correlation(vsd.rbo.df.sync %>% 
                    select(c(rbo.score.M,
                             rbo.score.F,
                             Male.time.female.near.barrier,
                             Court.total,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone_pg.mL.g.male,
                             Estradiol_pg.mL.g.female)),
                  histogram = TRUE,
                  pch = 19)
dev.off()

#### test rbo/gespeR matrix ####
#https://rdrr.io/bioc/gespeR/man/rbo.html
## from paper:
# For instance, p = 0.9 means that the first 10 ranks have 86% of the weight of the evaluation; 
# to give the top 50 ranks the same weight involves taking p = 0.98 as the setting.
# raising p value deepens the comparison

##https://ai.plainenglish.io/comparing-top-k-rankings-statistically-9adfc9cfc98b
# weight at a given depth is equal to (1-p)p^(d-1) for a given parameter p and depth d

# so for a dataframe of 10 ME's
test.d = seq(1:10)
test.p = seq(.1:1, 
             by = 0.1)

# set p
test = expand.grid(d = test.d,
                   p = test.p) 

# calculate weight at a given depth
test = test %>% 
  mutate(weight = (1-p)*(p^(d-1))) %>%
  group_by(p) %>%
  mutate(weight.cum = cumsum(weight)) %>% 
  ungroup()

## graph weight results
test %>% 
  ggplot(aes(x = d,
             y = weight,
             group = p,
             color = p,
             label = p)) +
  geom_line() +
  geom_label() +
  theme_classic()


## graph weight results
test %>% 
  ggplot(aes(x = d,
             y = weight.cum,
             group = p,
             color = p,
             label = p)) +
  geom_hline(yintercept = 0.5) +
  geom_line() +
  geom_label() +
  theme_classic()



## larger sample size 
# so for a dataframe of 10 ME's
test.d = seq(from = 1,
             to = 1000,
             by = 10)
test.p = seq(.9:1, 
             by = 0.01)

# set p
test = expand.grid(d = test.d,
                   p = test.p) 

test = test %>% 
  mutate(weight = (1-p)*(p^(d-1))) %>%
  group_by(p) %>%
  mutate(weight.cum = cumsum(weight)) %>% 
  ungroup()


## graph weight results
test %>% 
  filter(d < 250) %>% 
  ggplot(aes(x = d,
             y = weight,
             group = p,
             color = p,
             label = p)) +
  geom_line() +
  geom_point() +
  theme_classic()


## graph weight results
test %>% 
  filter(d < 250) %>% 
  ggplot(aes(x = d,
             y = weight.cum,
             group = p,
             color = p)) +
  geom_hline(yintercept = 0.5) +
  geom_line() +
  geom_point() +
  theme_classic()

###
#test rbo
# genes
rbo(vsd.df[,1],
    vsd.df[,2],
    p = 1,
    k = 1000)

# MEs
rbo(MEs.df[,1],
    MEs.df[,2],
    p = 0.95,
    k= 4)




## test across range of RBO values
rbo.p.values = c(0.5, 0.75, .9,.95,.98,.99,.995,.998,.999,1)
k.value = c(4,8)

tmp = data.frame(V1 = numeric(),
                 V2 = numeric())


tmp = tmp %>% 
  rbind(sapply(k.value, function(j) sapply(rbo.p.values, function(i){
    rbo(MEs.df[,1],
        MEs.df[,2],
        p = i,
        k= j)},
    USE.NAMES = TRUE
  )) %>% 
    as.data.frame() %>% 
    mutate(from = colnames(MEs.df)[1],
           to = colnames(MEs.df)[2],
           partner = 'same',
           rbo.p.values = rbo.p.values))

tmp = tmp %>% 
  rbind(sapply(k.value, function(j) sapply(rbo.p.values, function(i){
    rbo(MEs.df[,1],
        MEs.df[,3],
        p = i,
        k= j)},
    USE.NAMES = TRUE
  )) %>% 
    as.data.frame() %>% 
    mutate(from = colnames(MEs.df)[1],
           to = colnames(MEs.df)[3],
           partner = 'other',
           rbo.p.values = rbo.p.values))

tmp = tmp %>% 
  rbind(sapply(k.value, function(j) sapply(rbo.p.values, function(i){
    rbo(MEs.df[,3],
        MEs.df[,4],
        p = i,
        k= j)},
    USE.NAMES = TRUE
  )) %>% 
    as.data.frame() %>% 
    mutate(from = colnames(MEs.df)[3],
           to = colnames(MEs.df)[4],
           partner = 'same',
           rbo.p.values = rbo.p.values))

tmp = tmp %>% 
  rbind(sapply(k.value, function(j) sapply(rbo.p.values, function(i){
    rbo(MEs.df[,3],
        MEs.df[,2],
        p = i,
        k= j)},
    USE.NAMES = TRUE
  )) %>% 
    as.data.frame() %>% 
    mutate(from = colnames(MEs.df)[3],
           to = colnames(MEs.df)[2],
           partner = 'other',
           rbo.p.values = rbo.p.values))


tmp = tmp %>% 
  dplyr::rename('k.4' = 'V1',
                'k.8' = 'V2')

tmp = tmp %>% 
  pivot_longer(cols = c(k.4,
                        k.8),
               names_to = 'k.value',
               values_to = 'rbo.score')


tmp %>% 
  ggplot(aes(x = rbo.p.values,
             y = rbo.score,
             color = partner,
             shape = k.value)) +
  geom_line(aes(linetype = k.value))+
  geom_point() +
  theme_classic() +
  facet_grid(. ~ from)

####
# rbo.p.values = c(0.5, 0.75, .9,.95,.98,.99,.995,.998,.999,1)
# k.value = c(1000,3336)

rbo.p.values = c(.998,.999,1)
k.value = c(500,1000)

tmp = data.frame(V1 = numeric(),
                 V2 = numeric())
tmp = tmp %>% 
  rbind(sapply(k.value, function(j) sapply(rbo.p.values, function(i){
    rbo(vsd.df[,1],
        vsd.df[,2],
        p = i,
        k= j)},
    USE.NAMES = TRUE
  )) %>% 
    as.data.frame() %>% 
    mutate(from = colnames(vsd.df)[1],
           to = colnames(vsd.df)[2],
           partner = 'same',
           rbo.p.values = rbo.p.values))

tmp = tmp %>% 
  rbind(sapply(k.value, function(j) sapply(rbo.p.values, function(i){
    rbo(vsd.df[,1],
        vsd.df[,3],
        p = i,
        k= j)},
    USE.NAMES = TRUE
  )) %>% 
    as.data.frame() %>% 
    mutate(from = colnames(vsd.df)[1],
           to = colnames(vsd.df)[3],
           partner = 'other',
           rbo.p.values = rbo.p.values))

tmp = tmp %>% 
  rbind(sapply(k.value, function(j) sapply(rbo.p.values, function(i){
    rbo(vsd.df[,3],
        vsd.df[,4],
        p = i,
        k= j)},
    USE.NAMES = TRUE
  )) %>% 
    as.data.frame() %>% 
    mutate(from = colnames(vsd.df)[3],
           to = colnames(vsd.df)[4],
           partner = 'same',
           rbo.p.values = rbo.p.values))

tmp = tmp %>% 
  rbind(sapply(k.value, function(j) sapply(rbo.p.values, function(i){
    rbo(vsd.df[,3],
        vsd.df[,2],
        p = i,
        k= j)},
    USE.NAMES = TRUE
  )) %>% 
    as.data.frame() %>% 
    mutate(from = colnames(vsd.df)[3],
           to = colnames(vsd.df)[2],
           partner = 'other',
           rbo.p.values = rbo.p.values))


tmp = tmp %>% 
  dplyr::rename('k.500' = 'V1',
                'k.1000' = 'V2')

tmp = tmp %>% 
  pivot_longer(cols = c(k.500,
                        k.1000),
               names_to = 'k.value',
               values_to = 'rbo.score')


tmp %>% 
  ggplot(aes(x = rbo.p.values,
             y = rbo.score,
             color = partner,
             shape = k.value)) +
  geom_line(aes(linetype = k.value))+
  geom_point() +
  theme_classic() +
  facet_grid(. ~ from)







# 
# 
# 
# 
# 
# 
# df = data.frame(rbo.score = numeric(),
#                 p.score = numeric(),
#                 # k.value = numeric(),
#                 type = numeric(),
#                 pair = numeric())
# for(i in rbo.values) {                                   # Head of for-loop
#   tmp = c(rbo(vsd.df[,1],
#               vsd.df[,2],
#               p = i,
#               k = 1000),
#           i,
#           # j,
#           0,
#           1)# Create new row
#   df[nrow(df) + 1, ] <- tmp                   # Append new row
#   rm(tmp)
# }
# 
# for(i in rbo.values) {                                   # Head of for-loop
#   tmp = c(rbo(MEs.df[,1],
#               MEs.df[,2],
#               p = i,
#               k = 1000),
#           i,
#           # floor(max(length(MEs.df[,'E3.2.J3.M']), length(MEs.df[,'D5.4.A2.F']))/2),
#           1,
#           1)# Create new row
#   df[nrow(df) + 1, ] <- tmp                   # Append new row
#   rm(tmp)
# }
# 
# 
# for(i in rbo.values) {                                   # Head of for-loop
#   tmp = c(rbo(vsd.df[,1],
#               vsd.df[,4],
#               p = i,
#               k = 1000),
#           i,
#           # floor(max(length(vsd.df[,'E3.2.J3.M']), length(vsd.df[,'F2.1.J2.F']))/2),
#           0,
#           0)# Create new row
#   df[nrow(df) + 1, ] <- tmp                   # Append new row
#   rm(tmp)
# }
# 
# for(i in rbo.values) {                                   # Head of for-loop
#   tmp = c(rbo(MEs.df[,1],
#               MEs.df[,4],
#               p = i,
#               k = 1000),
#           i,
#           # floor(max(length(MEs.df[,'E3.2.J3.M']), length(MEs.df[,'D5.4.A2.F']))/2),
#           1,
#           0)# Create new row
#   df[nrow(df) + 1, ] <- tmp                   # Append new row
#   rm(tmp)
# }
# 
# df %>% 
#   mutate(type = ifelse(type == 1,
#                      'MEs',
#                      'vsd'))%>% 
#   mutate(pair = ifelse(pair == 1,
#                        'partner',
#                        'other')) %>% 
#   ggplot(aes(x = p.score,
#              y = rbo.score,
#              color = type,
#              shape = pair)) +
#   geom_line(aes(linetype = pair))+
#   geom_point() +
#   theme_classic() +
#   ggtitle(colnames(vsd.df)[1])
# 
# 
# 
# 
# df = data.frame(rbo.score = numeric(),
#                 p.score = numeric(),
#                 # k.value = numeric(),
#                 type = numeric(),
#                 pair = numeric())
# for(i in rbo.values) {                                   # Head of for-loop
#   tmp = c(rbo(vsd.df[,3],
#               vsd.df[,2],
#               p = i,
#               k = 1000),
#           i,
#           # floor(max(length(vsd.df[,'E3.2.J3.M']), length(vsd.df[,'F2.1.J2.F']))/2),
#           0,
#           1)# Create new row
#   df[nrow(df) + 1, ] <- tmp                   # Append new row
#   rm(tmp)
# }
# 
# for(i in rbo.values) {                                   # Head of for-loop
#   tmp = c(rbo(MEs.df[,3],
#               MEs.df[,2],
#               p = i,
#               k = 1000),
#           i,
#           # floor(max(length(MEs.df[,'E3.2.J3.M']), length(MEs.df[,'D5.4.A2.F']))/2),
#           1,
#           1)# Create new row
#   df[nrow(df) + 1, ] <- tmp                   # Append new row
#   rm(tmp)
# }
# 
# 
# for(i in rbo.values) {                                   # Head of for-loop
#   tmp = c(rbo(vsd.df[,3],
#               vsd.df[,4],
#               p = i,
#               k = 1000),
#           i,
#           # floor(max(length(vsd.df[,'E3.2.J3.M']), length(vsd.df[,'F2.1.J2.F']))/2),
#           0,
#           0)# Create new row
#   df[nrow(df) + 1, ] <- tmp                   # Append new row
#   rm(tmp)
# }
# 
# for(i in rbo.values) {                                   # Head of for-loop
#   tmp = c(rbo(MEs.df[,3],
#               MEs.df[,4],
#               p = i,
#               k = 1000),
#           i,
#           # floor(max(length(MEs.df[,'E3.2.J3.M']), length(MEs.df[,'D5.4.A2.F']))/2),
#           1,
#           0)# Create new row
#   df[nrow(df) + 1, ] <- tmp                   # Append new row
#   rm(tmp)
# }
# 
# df %>% 
#   mutate(type = ifelse(type == 1,
#                        'MEs',
#                        'vsd'))%>% 
#   mutate(pair = ifelse(pair == 1,
#                        'other',
#                        'partner')) %>% 
#   ggplot(aes(x = p.score,
#              y = rbo.score,
#              color = type,
#              shape = pair)) +
#   geom_line(aes(linetype = pair))+
#   geom_point() +
#   theme_classic() +
#   ggtitle(colnames(vsd.df)[3])

#### WGCNA PCA ####
### create PCA for all WGCNA data
#only use numerical variables (not categorical!) 
MEs.pca = prcomp(MEs,
                      scale = TRUE)

#check PCs
summary(MEs.pca)


# PC variance
percentVar.MEs.pca <- data.frame(variance = 100*MEs.pca$sdev^2 / sum( MEs.pca$sdev^2),
                             PC = MEs.pca[["rotation"]] %>% colnames(),
                             order = rep(1:length(MEs.pca[["rotation"]] %>% colnames()))) 

# ME PCA dataframe
MEs.pca.df = MEs.pca$x 

## graphing 
##scree plot
percentVar.MEs.pca %>% 
  ggplot(aes(y=variance,
             x = fct_reorder(PC,
                             order),
             group = 1)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  ggtitle('PCA screeplot') 
ggsave('../Figures/WGCNA/PCA/Scree plot MEs.png',
       height = 10,
       width = 10)
###plot pca
## PC 1 vs PC 2
MEs.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('ID') %>%
  ggplot()   +
  geom_label(data = MEs.pca$rotation %>%
               data.frame() %>% 
               rownames_to_column('Type'),
             aes(x = PC1*max(MEs.pca$x),
                 y = PC2*max(MEs.pca$x) + 0.1,
                 label = Type ),
             size = 3,
             color = 'black') +
  geom_segment(data = MEs.pca$rotation %>% 
                 data.frame(),
               aes(x = 0, 
                   y = 0, 
                   xend = PC1*max(MEs.pca$x),
                   yend = PC2*max(MEs.pca$x)),
               arrow = arrow(length = unit(0.5, "cm")),
               color = 'red') + 
  geom_point(aes(x = PC1,
                 y = PC2),
             size=3) +
  xlab(paste0("PC1: ",round(percentVar.MEs.pca[1,1]),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar.MEs.pca[2,1]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/WGCNA/PCA/PCA MEs 1 vs 2.png',
       height = 10,
       width = 10)

# add sex
MEs.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('Sample.Name') %>% 
  mutate(sample.names.select) %>% 
  ggplot()   +
  geom_point(aes(x = PC1,
                 y = PC2,
                 color = Sex),
             size=3) +
  xlab(paste0("PC1: ",round(percentVar.MEs.pca[1,1]),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar.MEs.pca[2,1]),"% variance")) + 
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c('red',
                                'black'))
ggsave('../Figures/WGCNA/PCA/PCA MEs 1 vs 2 sex.png',
       height = 10,
       width = 10)

# add tank
MEs.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('Sample.Name') %>% 
  mutate(sample.names.select) %>% 
  ggplot()   +
  geom_line(aes(x = PC1,
                 y = PC2,
                group = Tank.Round)) +
  geom_point(aes(x = PC1,
                 y = PC2,
                 color = Sex),
             size=3) +
  xlab(paste0("PC1: ",round(percentVar.MEs.pca[1,1]),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar.MEs.pca[2,1]),"% variance")) + 
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c('red',
                                'black'))
ggsave('../Figures/WGCNA/PCA/PCA MEs 1 vs 2 sex tank.png',
       height = 10,
       width = 10)

## PC 1 vs PC 3
MEs.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('ID') %>%
  ggplot()   +
  geom_label(data = MEs.pca$rotation %>%
               data.frame() %>% 
               rownames_to_column('Type'),
             aes(x = PC1*max(MEs.pca$x),
                 y = PC3*max(MEs.pca$x) + 0.1,
                 label = Type ),
             size = 3,
             color = 'black') +
  geom_segment(data = MEs.pca$rotation %>% 
                 data.frame(),
               aes(x = 0, 
                   y = 0, 
                   xend = PC1*max(MEs.pca$x),
                   yend = PC3*max(MEs.pca$x)),
               arrow = arrow(length = unit(0.5, "cm")),
               color = 'red') + 
  geom_point(aes(x = PC1,
                 y = PC3),
             size=3) +
  xlab(paste0("PC1: ",round(percentVar.MEs.pca[1,1]),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar.MEs.pca[3,1]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/WGCNA/PCA/PCA MEs 1 vs 3.png',
       height = 10,
       width = 10)

# add sex
MEs.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('Sample.Name') %>% 
  mutate(sample.names.select) %>% 
  ggplot()   +
  geom_point(aes(x = PC1,
                 y = PC3,
                 color = Sex),
             size=3) +
  xlab(paste0("PC1: ",round(percentVar.MEs.pca[1,1]),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar.MEs.pca[3,1]),"% variance")) + 
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c('red',
                                'black'))
ggsave('../Figures/WGCNA/PCA/PCA MEs 1 vs 3 sex.png',
       height = 10,
       width = 10)

# add tank
MEs.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('Sample.Name') %>% 
  mutate(sample.names.select) %>% 
  ggplot()   +
  geom_line(aes(x = PC1,
                y = PC3,
                group = Tank.Round)) +
  geom_point(aes(x = PC1,
                 y = PC3,
                 color = Sex),
             size=3) +
  xlab(paste0("PC1: ",round(percentVar.MEs.pca[1,1]),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar.MEs.pca[3,1]),"% variance")) + 
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c('red',
                                'black'))
ggsave('../Figures/WGCNA/PCA/PCA MEs 1 vs 3 sex tank.png',
       height = 10,
       width = 10)


## PC 2 vs PC 3
MEs.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('ID') %>%
  ggplot()   +
  geom_label(data = MEs.pca$rotation %>%
               data.frame() %>% 
               rownames_to_column('Type'),
             aes(x = PC2*max(MEs.pca$x),
                 y = PC3*max(MEs.pca$x) + 0.1,
                 label = Type ),
             size = 3,
             color = 'black') +
  geom_segment(data = MEs.pca$rotation %>% 
                 data.frame(),
               aes(x = 0, 
                   y = 0, 
                   xend = PC2*max(MEs.pca$x),
                   yend = PC3*max(MEs.pca$x)),
               arrow = arrow(length = unit(0.5, "cm")),
               color = 'red') + 
  geom_point(aes(x = PC2,
                 y = PC3),
             size=3) +
  xlab(paste0("PC2: ",round(percentVar.MEs.pca[1,1]),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar.MEs.pca[3,1]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/WGCNA/PCA/PCA MEs 2 vs 3.png',
       height = 10,
       width = 10)

# add sex
MEs.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('Sample.Name') %>% 
  mutate(sample.names.select) %>% 
  ggplot()   +
  geom_point(aes(x = PC2,
                 y = PC3,
                 color = Sex),
             size=3) +
  xlab(paste0("PC2: ",round(percentVar.MEs.pca[2,1]),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar.MEs.pca[3,1]),"% variance")) + 
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c('red',
                                'black'))
ggsave('../Figures/WGCNA/PCA/PCA MEs 2 vs 3 sex.png',
       height = 10,
       width = 10)

# add tank
MEs.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('Sample.Name') %>% 
  mutate(sample.names.select) %>% 
  ggplot()   +
  geom_line(aes(x = PC2,
                y = PC3,
                group = Tank.Round)) +
  geom_point(aes(x = PC2,
                 y = PC3,
                 color = Sex),
             size=3) +
  xlab(paste0("PC2: ",round(percentVar.MEs.pca[1,1]),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar.MEs.pca[3,1]),"% variance")) + 
  coord_fixed() +
  theme_classic() +
  scale_color_manual(values = c('red',
                                'black'))
ggsave('../Figures/WGCNA/PCA/PCA MEs 2 vs 3 sex tank.png',
       height = 10,
       width = 10)



#### WGCNA PCA distance matrix ####
## create distance matrix
MEs.pca.dist = dist(MEs.pca.df, 
                    method = "euclidean",
                    upper = T,
                    diag = T) %>% 
  as.matrix()


# add names to matrix
colnames(MEs.pca.dist) = rownames(MEs.pca.df)
rownames(MEs.pca.dist) = rownames(MEs.pca.df)

# convert dataframe to long format
MEs.pca.dist = MEs.pca.dist %>% 
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
MEs.pca.dist = MEs.pca.dist %>%
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

#### graph WGCNA PCA distance matrix ####
## graph same vs. other
#poster
#females to males
MEs.pca.dist %>% 
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
  theme_classic(base_size = 30) +
  # ggtitle('Female rbo score compared to all males') +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Female') +
  ylab('Dist') +
  # labs(color = 'Comparison')
  theme(legend.position = "none")
ggsave('../Figures/WGCNA/PCA/Dist from females to all males paper.png',
       height = 5,
       width = 10)

#males to females
MEs.pca.dist %>% 
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
  theme_classic(base_size = 30) +
  # ggtitle('Male rbo score compared to all females') +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Male') +
  ylab('Dist') +
  # labs(color = 'Comparison')+ 
  theme(legend.position = "none")
ggsave('../Figures/WGCNA/PCA/Dist from males to all females paper.png',
       height = 5,
       width = 10)


#males to males
#poster
MEs.pca.dist %>% 
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
  ylab('Dist')
ggsave('../Figures/WGCNA/PCA/Dist from females to all males zscore all poster.pdf',
       height = 5.25,
       width = 5.25)

#males to females
#poster
MEs.pca.dist %>% 
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
  ylab('Dist')
ggsave('../Figures/WGCNA/PCA/Dist from males to all females zscore all poster.pdf',
       height = 5.25,
       width = 5.25)

### combine with behavior data
## male
MEs.pca.dist.male = MEs.pca.dist %>% 
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
pdf('../Figures/WGCNA/PCA/ME dist and behavior chart correlation male.pdf',
    height = 10,
    width = 10)
chart.Correlation(MEs.pca.dist.male %>% 
                    select(-c(Observation.id)),
                  histogram = TRUE,
                  pch = 19)
dev.off()


## male
MEs.pca.dist.female = MEs.pca.dist %>% 
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
pdf('../Figures/WGCNA/PCA/ME dist and behavior chart correlation female.pdf',
    height = 10,
    width = 10)
chart.Correlation(MEs.pca.dist.female %>% 
                    select(-c(Observation.id)),
                  histogram = TRUE,
                  pch = 19)
dev.off()



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
  theme_classic(base_size = 30) +
  # ggtitle('Female rbo score compared to all males') +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Female') +
  ylab('Dist') +
  # labs(color = 'Comparison')
  theme(legend.position = "none")
ggsave('../Figures/WGCNA/ME/Dist from females to all males paper.png',
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
  theme_classic(base_size = 30) +
  # ggtitle('Male rbo score compared to all females') +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Male') +
  ylab('Dist') +
  # labs(color = 'Comparison')+ 
  theme(legend.position = "none")
ggsave('../Figures/WGCNA/ME/Dist from males to all females paper.png',
       height = 5,
       width = 10)


#males to males
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
  ylab('Dist')
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
  ylab('Dist')
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






#### save zscore distance ####
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

#### gene expression PCA distance
#### compare WGCNA PCA dist, ME dist,  to rbo score ####
## load data
load('../Figures/WGCNA/Rbo.score/ME.rbo.df.sum.sync.RData')
## combine rbo score and dist

MEs.pca.dist.rbo = MEs.pca.dist %>% 
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
  full_join(ME.rbo.df.sum.sync) |> 
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
              dplyr::rename(dist.ME = dist))

# check correlation
pdf('../Figures/WGCNA/PCA/ME dist and rbo score correlation.pdf',
    height = 10,
    width = 10)
chart.Correlation(MEs.pca.dist.rbo %>% 
                    select(c(Zscore.M,
                             Zscore.F,
                             dist,
                             dist.ME)),
                  histogram = TRUE,
                  pch = 19)
dev.off()

















