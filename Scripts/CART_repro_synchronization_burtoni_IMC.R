#### Burtoni reproductive opportunity data
## CART 
# https://www.pluralsight.com/guides/explore-r-libraries:-rpart

### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_repro_opportunity/Synchronization/RNAseq/Data/")

#### load libraries ####
library(plyr)
library(readr)
library(dplyr)
library(caret)
library(rpart)
library(rpart.plot)
library(jtools)
library(ggrepel)
library(tidyverse)

#### load data ####
load("../Figures/WGCNA/Rbo.score/ME.rbo.df.sum.sync.RData")

#### CART ####
### female
set.seed(100)
# trainRowNumbers <- createDataPartition(ME.rbo.df.sum.sync$Zscore.M, p=0.5, list=FALSE)
# train <- ME.rbo.df.sum.sync[trainRowNumbers,]
# test <- ME.rbo.df.sum.sync[-trainRowNumbers,]
# dim(train); dim(test)
train.female = ME.rbo.df.sum.sync

## setup data
cols.female = c('Male.time.female.near.barrier',
                'Latency',
                'GSI.female',
                'Estradiol_pg.mL.g.female',
                'Court.total',
                'Testosterone_pg.mL.g.male')
#scale and center
pre_proc_val <- preProcess(train.female[,cols.female], method = c("center", "scale"))
# replace columns
train.female[,cols.female] = predict(pre_proc_val, train.female[,cols.female])

## run model
tree_model.female = rpart(Zscore.F ~ Male.time.female.near.barrier + Latency + GSI.female + Estradiol_pg.mL.g.female + Court.total + Testosterone_pg.mL.g.male, 
                     data = train.female,
                     minsplit = 9, 
                     minbucket= 3) 
## model summary
summary(tree_model.female)

## graph
# prp(tree_model)
## tree
png('../Figures/WGCNA/Rbo.score/CART/female similarity tree.png')
rpart.plot(tree_model.female)
title('female')
dev.off()

## variable importance
tree_model.female$variable.importance %>% 
  data.frame() %>% 
  dplyr::rename('Variable.importance' = '.') %>% 
  rownames_to_column('Variable') %>% 
  mutate(Variable.importance = 10*round(Variable.importance,
                                        digits = 2)) %>% 
  ggplot(aes(y = Variable.importance,
             x = reorder(Variable,
                         -Variable.importance))) +
  geom_col() +
  ggtitle('female similarity')+
  ylim(0,100) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab('')
ggsave('../Figures/WGCNA/Rbo.score/CART/female similarity variable importance.png',
       height = 10,
       width = 10)

### male
set.seed(100)
# trainRowNumbers <- createDataPartition(ME.rbo.df.sum.sync$Zscore.M, p=0.7, list=FALSE)
# train <- dat[trainRowNumbers,]
# test <- dat[-trainRowNumbers,]
# dim(train); dim(test) 
train.male = ME.rbo.df.sum.sync

## setup data
cols.male = c('Male.time.female.near.barrier',
              'Latency',
              'GSI.female',
              'Estradiol_pg.mL.g.female',
              'Court.total',
              'Testosterone_pg.mL.g.male')
#scale and center
pre_proc_val <- preProcess(train.male[,cols.male], method = c("center", "scale"))
# replace columns
train.male[,cols.male] = predict(pre_proc_val, train.male[,cols.male])

## run model
tree_model.male = rpart(Zscore.M ~ Male.time.female.near.barrier + Latency + GSI.female + Estradiol_pg.mL.g.female + Court.total + Testosterone_pg.mL.g.male, 
                          data = train.male,
                          minsplit = 10, 
                          minbucket=3
) 
## model summary
summary(tree_model.male)

## graph
# prp(tree_model)
## tree
png('../Figures/WGCNA/Rbo.score/CART/male similarity tree.png')
rpart.plot(tree_model.male)
title('male')
dev.off()

## variable importance
tree_model.male$variable.importance %>% 
  data.frame() %>% 
  dplyr::rename('Variable.importance' = '.') %>% 
  rownames_to_column('Variable') %>% 
  mutate(Variable.importance = 10*round(Variable.importance,
                                        digits = 2)) %>% 
  ggplot(aes(y = Variable.importance,
             x = reorder(Variable,
                         -Variable.importance))) +
  geom_col() +
  ggtitle('male similarity')+
  ylim(0,100) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab('')
ggsave('../Figures/WGCNA/Rbo.score/CART/male similarity variable importance.png',
       height = 10,
       width = 10)

## compare 
tree_model.female$variable.importance %>% 
  data.frame() %>% 
  dplyr::rename('Variable.importance.female' = '.') %>% 
  rownames_to_column('Variable') %>% 
  mutate(Variable.importance.female = 10*round(Variable.importance.female,
                                        digits = 2)) %>% 
  full_join(tree_model.male$variable.importance %>% 
              data.frame() %>% 
              dplyr::rename('Variable.importance.male' = '.') %>% 
              rownames_to_column('Variable') %>% 
              mutate(Variable.importance.male = 10*round(Variable.importance.male,
                                                    digits = 2))) %>% 
  replace(is.na(.), 0) %>% 
  ggplot(aes(x = Variable.importance.male,
             y= Variable.importance.female,
             label = Variable)) + 
  theme_classic() +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0,100)) + 
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,100)) +
  geom_abline(slope = 1,
            intercept = 0) +
  geom_point() +
  geom_text_repel() 
ggsave('../Figures/WGCNA/Rbo.score/CART/compare variable importance.png',
       height = 10,
       width = 10)  


#### CART model accuracy ####
### training data accuracy
PredictCART_train.female = predict(tree_model.female, 
                            data = train.female)
PredictCART_train.male = predict(tree_model.male, 
                                   data = train.male)


## graph training data accuracy
## graph residuals of prediction? Or 95% confidence interval of R^2?
# data.frame(female.similarity = train.female$Zscore.F, 
#       female.model.prediction = PredictCART_train.female) %>% 
#   ggplot(aes(x = female.similarity,
#              y = female.model.prediction)) +
#   geom_point() +
#   geom_smooth(method = 'lm') +
#   theme_bw() +
#   ggtitle('female model accuracy training')
# ggsave('../Figures/WGCNA/Rbo.score/CART/female model accuracy training.png',
#                                                   height = 10,
#                                                   width = 10)
# 
# 
# data.frame(male.similarity = train.male$Zscore.M, 
#            male.model.prediction = PredictCART_train.male) %>% 
#   ggplot(aes(x = male.similarity,
#              y = male.model.prediction)) +
#   geom_point() +
#   geom_smooth(method = 'lm')+
#   theme_bw()+
#   ggtitle('male model accuracy training')
# ggsave('../Figures/WGCNA/Rbo.score/CART/male model accuracy training.png',
#        height = 10,
#        width = 10)


## add equation to graph
# females
tmp = data.frame(female.similarity = train.female$Zscore.F, 
           female.model.prediction = PredictCART_train.female)

tmp2 = lm(female.model.prediction ~ female.similarity,
                          data=tmp)


summary(tmp2)
anova(tmp2)

tmp2.p = signif(summary(tmp2)[["coefficients"]]['female.similarity','Pr(>|t|)'], digits = 1)

tmp2.r = signif(sqrt(summary(tmp2)[["r.squared"]]), digits = 2)

effect_plot(tmp2,
            pred = female.similarity,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=tmp ,
             aes(x=female.similarity,
                 y=female.model.prediction),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  theme(text = element_text(size = 30))+
  annotate("text", x=1, y=-1, label= paste0('p-value = ',
                                            tmp2.p,
                                                ', R = ',
                                            tmp2.r
  )) +
  ggtitle('female model accuracy training')
ggsave('../Figures/WGCNA/Rbo.score/CART/female model accuracy training equation.png',
       height = 10,
       width = 10)


# males
tmp = data.frame(male.similarity = train.male$Zscore.M, 
                 male.model.prediction = PredictCART_train.male)

tmp2 = lm(male.model.prediction ~ male.similarity,
          data=tmp)


summary(tmp2)
anova(tmp2)

tmp2.p = signif(summary(tmp2)[["coefficients"]]['male.similarity','Pr(>|t|)'], digits = 1)

tmp2.r = signif(sqrt(summary(tmp2)[["r.squared"]]), digits = 2)

effect_plot(tmp2,
            pred = male.similarity,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=tmp ,
             aes(x=male.similarity,
                 y=male.model.prediction),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  theme(text = element_text(size = 30))+
  annotate("text", x=1, y=-1, label= paste0('p-value = ',
                                            tmp2.p,
                                            ', R = ',
                                            tmp2.r
  )) +
  ggtitle('male model accuracy training')
ggsave('../Figures/WGCNA/Rbo.score/CART/male model accuracy training equation.png',
       height = 10,
       width = 10)






### test data accuracy
PredictCART_test.female = predict(tree_model.female, 
                                   data = train.female %>% 
                                    mutate(Zscore.F = Zscore.M))

PredictCART_test.male = predict(tree_model.male, 
                                  data = train.male %>% 
                                    mutate(Zscore.M = Zscore.F))


## graph residuals of prediction? Or 95% confidence interval of R^2?
## maybe convert prediction to rank score? 
# data.frame(male.similarity = train.female$Zscore.M, 
#            female.model.prediction = PredictCART_test.female) %>% 
#   ggplot(aes(x = male.similarity,
#              y = female.model.prediction)) +
#   geom_point() +
#   geom_smooth(method = 'lm')+
#   theme_bw()+
#   ggtitle('female model accuracy test')
# ggsave('../Figures/WGCNA/Rbo.score/CART/female model accuracy test.png',
#        height = 10,
#        width = 10)
# 
# data.frame(female.similarity = train.male$Zscore.F, 
#            male.model.prediction = PredictCART_test.male) %>% 
#   ggplot(aes(x = female.similarity,
#              y = male.model.prediction)) +
#   geom_point() +
#   geom_smooth(method = 'lm')+
#   theme_bw()+
#   ggtitle('male model accuracy test')
# ggsave('../Figures/WGCNA/Rbo.score/CART/male model accuracy test.png',
#        height = 10,
#        width = 10)



## add equation to graph
# females
tmp = data.frame(male.similarity = train.female$Zscore.M, 
                 female.model.prediction = PredictCART_test.female)

tmp2 = lm(female.model.prediction ~ male.similarity,
          data=tmp)


summary(tmp2)
anova(tmp2)

tmp2.p = signif(summary(tmp2)[["coefficients"]]['male.similarity','Pr(>|t|)'], digits = 1)

tmp2.r = signif(sqrt(summary(tmp2)[["r.squared"]]), digits = 2)

effect_plot(tmp2,
            pred = male.similarity,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=tmp ,
             aes(x=male.similarity,
                 y=female.model.prediction),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  theme(text = element_text(size = 30))+
  annotate("text", x=1, y=-1, label= paste0('p-value = ',
                                            tmp2.p,
                                            ', R = ',
                                            tmp2.r
  )) +
  ggtitle('female model accuracy test')
ggsave('../Figures/WGCNA/Rbo.score/CART/female model accuracy test equation.png',
       height = 10,
       width = 10)


# males
tmp = data.frame(female.similarity = train.male$Zscore.F, 
                 male.model.prediction = PredictCART_test.male)

tmp2 = lm(male.model.prediction ~ female.similarity,
          data=tmp)


summary(tmp2)
anova(tmp2)

tmp2.p = signif(summary(tmp2)[["coefficients"]]['female.similarity','Pr(>|t|)'], digits = 1)

tmp2.r = signif(sqrt(summary(tmp2)[["r.squared"]]), digits = 2)

effect_plot(tmp2,
            pred = female.similarity,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=tmp ,
             aes(x=female.similarity,
                 y=male.model.prediction),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  theme(text = element_text(size = 30))+
  annotate("text", x=1, y=-1, label= paste0('p-value = ',
                                            tmp2.p,
                                            ', R = ',
                                            tmp2.r
  )) +
  ggtitle('male model accuracy test')
ggsave('../Figures/WGCNA/Rbo.score/CART/male model accuracy test equation.png',
       height = 10,
       width = 10)


# ## test model strength 
# summary(lm(female.model.prediction ~ male.similarity,
#            data.frame(male.similarity = train.female$Zscore.M,
#                       female.model.prediction = PredictCART_test.female)))
# 
# summary(lm(female.model.prediction ~ female.similarity,
#            data.frame(female.similarity = train.female$Zscore.F, 
#                       female.model.prediction = PredictCART_train.female)))
# 
# 
# summary(lm(male.model.prediction ~ female.similarity,
#            data.frame(female.similarity = train.female$Zscore.F, 
#                       male.model.prediction = PredictCART_test.male)))
# 
# summary(lm(male.model.prediction ~ male.similarity,
#            data.frame(male.similarity = train.male$Zscore.M, 
#                       male.model.prediction = PredictCART_train.male)))
# 
# 
# 
# 
# summary(lm(female.model.prediction ~ female.similarity,
#            data.frame(female.similarity = train.female$Zscore.F, 
#                       female.model.prediction = PredictCART_train.female)))
# 
tmp = data.frame(female.similarity = train.female$Zscore.F, 
                 female.model.prediction = PredictCART_train.female)

library(boot)
foo <- boot(tmp,function(data,indices)
  summary(lm(female.model.prediction ~ female.similarity,
             data[indices,]))$r.squared,
  R=1000)

foo <- boot(tmp,function(data,indices)
  summary(lm(female.model.prediction ~ female.similarity,
             data[indices,]))[["coefficients"]]['female.similarity','Pr(>|t|)'],
  R=1000)

foo$t0

quantile(foo$t,c(0.025,0.975))

tmp = data.frame(male.similarity = train.male$Zscore.M,
                 male.model.prediction = PredictCART_train.male)

library(boot)
foo <- boot(tmp,function(data,indices)
  summary(lm(male.model.prediction ~ male.similarity,
             data[indices,]))$r.squared,
  R=1000)

foo$t0

quantile(foo$t,c(0.025,0.975))

#### CART partition test ####
### female
set.seed(100)
trainRowNumbers <- createDataPartition(ME.rbo.df.sum.sync$Zscore.F, p=0.5, list=FALSE)
train <- ME.rbo.df.sum.sync[trainRowNumbers,]
test <- ME.rbo.df.sum.sync[-trainRowNumbers,]
dim(train); dim(test)
train.female = train

## setup data
cols.female = c('Male.time.female.near.barrier',
                'Latency',
                'GSI.female',
                'Estradiol_pg.mL.g.female',
                'Court.total',
                'Testosterone_pg.mL.g.male')
#scale and center
pre_proc_val <- preProcess(train.female[,cols.female], method = c("center", "scale"))
# replace columns
train.female[,cols.female] = predict(pre_proc_val, train.female[,cols.female])

## run model
tree_model.female = rpart(Zscore.F ~ Male.time.female.near.barrier + Latency + GSI.female + Estradiol_pg.mL.g.female + Court.total + Testosterone_pg.mL.g.male, 
                          data = train.female,
                          minsplit = 4, 
                          minbucket= 1) 
## model summary
summary(tree_model.female)

## graph
# prp(tree_model)
## tree
png('../Figures/WGCNA/Rbo.score/CART/female similarity tree.png')
rpart.plot(tree_model.female)
title('female')
dev.off()

## variable importance
tree_model.female$variable.importance %>% 
  as.data.frame() %>% 
  rename('Variable.importance' = '.') %>% 
  rownames_to_column('Variable') %>% 
  mutate(Variable.importance = 10*round(Variable.importance,
                                        digits = 2)) %>% 
  ggplot(aes(y = Variable.importance,
             x = reorder(Variable,
                         -Variable.importance))) +
  geom_col() +
  ggtitle('female similarity')+
  ylim(0,100) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab('')
ggsave('../Figures/WGCNA/Rbo.score/CART/female similarity variable importance.png',
       height = 10,
       width = 10)

### male
set.seed(100)
# # use female data
# trainRowNumbers <- createDataPartition(ME.rbo.df.sum.sync$Zscore.F, p=0.5, list=FALSE)
# train <- ME.rbo.df.sum.sync[trainRowNumbers,]
# test <- ME.rbo.df.sum.sync[-trainRowNumbers,]
# dim(train); dim(test)
train.male = train
test.male = test

## setup data
cols.male = c('Male.time.female.near.barrier',
              'Latency',
              'GSI.female',
              'Estradiol_pg.mL.g.female',
              'Court.total',
              'Testosterone_pg.mL.g.male')
#scale and center
pre_proc_val <- preProcess(train.male[,cols.male], method = c("center", "scale"))
# replace columns
train.male[,cols.male] = predict(pre_proc_val, train.male[,cols.male])

## run model
tree_model.male = rpart(Zscore.M ~ Male.time.female.near.barrier + Latency + GSI.female + Estradiol_pg.mL.g.female + Court.total + Testosterone_pg.mL.g.male, 
                        data = train.male,
                        minsplit = 4, 
                        minbucket=1
) 
## model summary
summary(tree_model.male)

## graph
# prp(tree_model)
## tree
png('../Figures/WGCNA/Rbo.score/CART/male similarity tree.png')
rpart.plot(tree_model.male)
title('male')
dev.off()

## variable importance
tree_model.male$variable.importance %>% 
  as.data.frame() %>% 
  rename('Variable.importance' = '.') %>% 
  rownames_to_column('Variable') %>% 
  mutate(Variable.importance = 10*round(Variable.importance,
                                        digits = 2)) %>% 
  ggplot(aes(y = Variable.importance,
             x = reorder(Variable,
                         -Variable.importance))) +
  geom_col() +
  ggtitle('male similarity')+
  ylim(0,100) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab('')
ggsave('../Figures/WGCNA/Rbo.score/CART/male similarity variable importance.png',
       height = 10,
       width = 10)

# replace columns
train.male[,cols.male] = predict(pre_proc_val, train.male[,cols.male])
test.male[,cols.male] = predict(pre_proc_val, test.male[,cols.male])

PredictCART_train.male = predict(tree_model.male, data = train.male)
PredictCART_test.male = predict(tree_model.male, data = test.male)

data.frame(male.similarity = train.male$Zscore.M, 
           male.model.prediction = PredictCART_train.male) %>% 
  ggplot(aes(x = male.similarity,
             y = male.model.prediction)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  theme_bw()+
  ggtitle('male model accuracy training')

data.frame(male.similarity = test.male$Zscore.M, 
           male.model.prediction = PredictCART_test.male) %>% 
  ggplot(aes(x = male.similarity,
             y = male.model.prediction)) +
  geom_point() +
  geom_smooth(method = 'lm')+
  theme_bw()+
  ggtitle('male model accuracy training')


#### party ####
library(party)
library(partykit)
library(ggparty)
#https://stackoverflow.com/questions/29131254/how-to-generate-a-prediction-interval-from-a-regression-tree-rpart-object
#https://www.statmethods.net/advstats/cart.html
#https://martinschweinberger.github.io/TreesUBonn/conditional-inference-trees.html#example-2-prepositions
### Conditional inference trees via party
#"Conditional Inference Trees (CITs) are much better at determining the true effect of a predictor, i.e. the effect of a predictor if all other effects are simultaneously considered"
# dim(train); dim(test)

## rename variables
ME.rbo.df.sum.sync.rename = ME.rbo.df.sum.sync %>% 
  dplyr::rename(Time.together = Male.time.female.near.barrier) %>% 
  dplyr::rename(Estradiol = Estradiol_pg.mL.g.female) %>% 
  dplyr::rename(Testosterone = Testosterone_pg.mL.g.male) 

# females
tree_model.female.party = ctree(Zscore.F ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                          data = ME.rbo.df.sum.sync.rename, 
                          control = ctree_control(minbucket = 3, 
                                                  minsplit = 6,
                                                  testtype = "Teststatistic",
                                                  mincriterion = 0.95))
                          
png('../Figures/WGCNA/Rbo.score/CIT/Female party model.png')
plot(tree_model.female.party, main="Female party model")
dev.off()

# plotting
png('../Figures/WGCNA/Rbo.score/CIT/Female party model present.png')
ggparty(tree_model.female.party) +
  geom_edge() +
  geom_edge_label(aes(label = paste(substr(breaks_label, start = 1, stop = 15)))) +
  geom_node_label(line_list = list(aes(label = splitvar),
                                   aes(label = "p < 0.001", 
                                       size = 10)),
                  line_gpar = list(list(size = 13), 
                                   list(size = 10)), 
                  ids = "inner") +
  geom_node_label(aes(label = paste0("N = ", nodesize)),
                  ids = "terminal", nudge_y = -0.0, nudge_x = 0.01) +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", 
                     y = Zscore.F), 
                 coef = Inf, 
                 width = 1, 
                 fill = "lightgray"),
    theme_minimal(),
    xlab(""),
    ylab("Synchronization")
  ),
  shared_axis_labels = TRUE)
dev.off()

# plotting
png('../Figures/WGCNA/Rbo.score/CIT/Female party model present space.png')
ggparty(tree_model.female.party) +
  geom_edge() +
  geom_edge_label(aes(label = paste(substr(breaks_label, start = 1, stop = 15)))) +
  geom_node_label(line_list = list(aes(label = paste0(splitvar, " (pg/mL)")),
                                   aes(label = "p < 0.001", 
                                       size = 10)),
                  line_gpar = list(list(size = 13), 
                                   list(size = 10)), 
                  ids = "inner") +
  geom_node_label(aes(label = paste0("N = ", nodesize)),
                  ids = "terminal", nudge_y = -0.0, nudge_x = 0.01) +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", 
                     y = Zscore.F), 
                 coef = Inf, 
                 width = 1, 
                 fill = "lightgray"),
    theme_minimal(),
    xlab(""),
    ylab("Synchronization")
  ),
  shared_axis_labels = TRUE)
dev.off()

# test$leaf <- predict(tree_model.female.party, type="node")
# CIleaf <- aggregate(test$Zscore.F,
#                     by=list(leaf=test$leaf),
#                     quantile,
#                     prob=c(0.025, 0.975))
# 
# 
# plot(as.factor(CIleaf$leaf), CIleaf[, 2],
#      ylab="Zscore.F", xlab="Regression tree leaf")
# legend("bottomright",
#        c(" 0.975 quantile", " 0.75 quantile", " mean", 
#          " 0.25 quantile", " 0.025 quantile"),
#        pch=c("-", "_", "_", "_", "-"),
#        pt.lwd=0.5, pt.cex=c(1, 1, 2, 1, 1), xjust=1)



# 
# 
# nodes(tree_model.female.party, 2)
# 
# nodes(tree_model.female.party, 1)[[1]]

# males
tree_model.male.party = ctree(Zscore.M ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                                data = ME.rbo.df.sum.sync.rename, 
                                control = ctree_control(minbucket = 3, 
                                                        minsplit = 6,
                                                        testtype = "Teststatistic",
                                                        mincriterion = 0.95)
) 

png('../Figures/WGCNA/Rbo.score/CIT/Male party model.png')
plot(tree_model.male.party, main="Male party model")
dev.off()

# plotting
png('../Figures/WGCNA/Rbo.score/CIT/Male party model present.png')
ggparty(tree_model.male.party) +
  geom_edge() +
  geom_edge_label(aes(label = paste(substr(breaks_label, start = 1, stop = 15)))) +
  geom_node_label(line_list = list(aes(label = splitvar),
                                   aes(label = "p < 0.001", 
                                       size = 10)),
                  line_gpar = list(list(size = 13), 
                                   list(size = 10)), 
                  ids = "inner") +
  geom_node_label(aes(label = paste0("N = ", nodesize)),
                  ids = "terminal", nudge_y = -0.0, nudge_x = 0.01) +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", 
                     y = Zscore.M), 
                 coef = Inf, 
                 width = 1, 
                 fill = "lightgray"),
    theme_minimal(),
    xlab(""),
    ylab("Synchronization")
  ),
  shared_axis_labels = TRUE)
dev.off()

# plotting
png('../Figures/WGCNA/Rbo.score/CIT/Male party model present space.png')
ggparty(tree_model.male.party) +
  geom_edge() +
  geom_edge_label(aes(label = paste(substr(breaks_label, start = 1, stop = 15)))) +
  geom_node_label(line_list = list(aes(label = paste0(splitvar, " (pg/mL)")),
                                   aes(label = "p < 0.001", 
                                       size = 10)),
                  line_gpar = list(list(size = 13), 
                                   list(size = 10)), 
                  ids = "inner") +
  geom_node_label(aes(label = paste0("N = ", nodesize)),
                  ids = "terminal", nudge_y = -0.0, nudge_x = 0.01) +
  geom_node_plot(gglist = list(
    geom_boxplot(aes(x = "", 
                     y = Zscore.M), 
                 coef = Inf, 
                 width = 1, 
                 fill = "lightgray"),
    theme_minimal(),
    xlab(""),
    ylab("Synchronization")
  ),
  shared_axis_labels = TRUE)
dev.off()

# 
# test$leaf <- predict(tree_model.male.party, type="node")
# CIleaf <- aggregate(test$Zscore.F,
#                     by=list(leaf=test$leaf),
#                     quantile,
#                     prob=c(0.025, 0.25, 0.75, 0.975))
# 
# 
# plot(as.factor(CIleaf$leaf), CIleaf[, 2],
#      ylab="Zscore.F", xlab="Regression tree leaf")
# legend("bottomright",
#        c(" 0.975 quantile", " 0.75 quantile", " mean", 
#          " 0.25 quantile", " 0.025 quantile"),
#        pch=c("-", "_", "_", "_", "-"),
#        pt.lwd=0.5, pt.cex=c(1, 1, 2, 1, 1), xjust=1)


### calculate variable importance
library(partykit)
# females
tree_model.female.partykit = partykit::ctree(Zscore.F ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                                data = ME.rbo.df.sum.sync.rename, 
                                control = ctree_control(minbucket = 3, 
                                                        minsplit = 6,
                                                        testtype = "Teststatistic",
                                                        mincriterion = 0.95))




# males
tree_model.male.partykit = partykit::ctree(Zscore.M ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                              data = ME.rbo.df.sum.sync.rename, 
                              control = ctree_control(minbucket = 3, 
                                                      minsplit = 6,
                                                      testtype = "Teststatistic",
                                                      mincriterion = 0.95)
) 

png('../Figures/WGCNA/Rbo.score/CIT/Male partykit model.png')
plot(tree_model.male.party, main="Male partykit model")
dev.off()

### calculate variable importance
library(partykit)
# females
tree_model.female.partykit = partykit::ctree(Zscore.F ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                                             data = ME.rbo.df.sum.sync.rename, 
                                             control = ctree_control(minbucket = 3, 
                                                                     minsplit = 6,
                                                                     testtype = "Teststatistic",
                                                                     mincriterion = 0.95))


png('../Figures/WGCNA/Rbo.score/CIT/Female partykit model.png')
plot(tree_model.female.partykit, main="Female partykit model")
dev.off()

# males
tree_model.male.partykit = partykit::ctree(Zscore.M ~ Time.together + Latency + GSI.male + GSI.female + Estradiol + Court.total + Testosterone, 
                                           data = ME.rbo.df.sum.sync.rename, 
                                           control = ctree_control(minbucket = 3, 
                                                                   minsplit = 6,
                                                                   testtype = "Teststatistic",
                                                                   mincriterion = 0.95)
) 

png('../Figures/WGCNA/Rbo.score/CIT/Male partykit model.png')
plot(tree_model.male.partykit, main="Male partykit model")
dev.off()

## extract variable importance
variable.importance.df = varimp(tree_model.female.partykit) %>% 
  as.data.frame() %>%
  rownames_to_column('variable') %>% 
  mutate(model = 'Female') %>% 
  full_join(varimp(tree_model.male.partykit) %>% 
              as.data.frame()%>%
              rownames_to_column('variable') %>% 
              mutate(model = 'Male')) %>% 
  dplyr::rename(variable.importance = '.') %>% 
  pivot_wider(names_from = model,
              values_from = variable.importance,
              names_prefix = "variable.importance.") %>% 
  replace(is.na(.), 0) 

# graph
variable.importance.df %>% 
  ggplot(aes(x = variable.importance.Male,
             y = variable.importance.Female,
             label = variable)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point() +
  geom_label() +
  theme_classic() +
  xlim(0,1.2) +
  ylim(0,1.2)
ggsave('../Figures/WGCNA/Rbo.score/CIT/female vs male variable importance.png',
       height = 10,
       width = 10)






#### random forest ####
# Random Forest prediction of Kyphosis data
library(randomForest)
fit <- randomForest(Zscore.F ~ Male.time.female.near.barrier + Latency + GSI.female + Estradiol_pg.mL.g.female , 
                    data = test)
print(fit) # view results
importance(fit) # importance of each predictor






