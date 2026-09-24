#### Burtoni reproductive opportunity data
# > R-4.0.3

### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_repro_opportunity/Synchronization/RNAseq/Data/")


#### load libraries ####
## installation

## load
library("DESeq2")
library("vsn")
library("pheatmap")
library("RColorBrewer")
library('WGCNA')
library(pvclust)
library(dendextend)
library(PerformanceAnalytics)
library("tidyverse")


#### load data ####
### sample IDs
sample.names = read.csv('../RNA_extraction_QC/DNA-RNA CSV Template_burtoni_repro_IMC.csv')
#replace '_' with '-'
sample.names = sample.names %>% 
  dplyr::select(Sample.Name) %>% 
  mutate(Sample.Name = str_replace_all(Sample.Name,
                                       '_',
                                       '-'),
         ID = Sample.Name) %>% 
  separate(ID,
           c('Plate.position',
             'Round',
             'Tank',
             'Sex'))


### load count files
## create empty dataframe
count.data = data.frame(gene = NA)
## run across all samples
## combine reads from files
for (i in sample.names$Sample.Name) {
  count.data = read.delim(paste('../Data/count/',
                                i,
                                '-ReadsPerGene.out.tab',
                                sep = ''),
                          header = FALSE,
                          col.names = c('gene',
                                        paste(i,
                                              '',
                                              sep = ''),
                                        'firststrand',
                                        'secondstrand'),
                          skip = 4) %>%
    select(-c(firststrand,
              secondstrand)) %>% 
    full_join(count.data)
  
}

#drop NA
count.data = count.data %>% 
  filter(!is.na(gene))

# load ensembl gene names
tilapia.gene.names.list = read.csv('Ensembl.gene.names.csv') %>% 
  dplyr::select(-c(X, # remove row numbers
                   external_gene_name)) 

  

# rename duplicates
tilapia.gene.names.list = tilapia.gene.names.list %>% 
  group_by(gene) %>%
  mutate(dup = case_when(n() == 1 ~ FALSE,
                         TRUE ~ TRUE)) %>% 
  mutate(observation = 1:n()) %>% 
  mutate(gene = ifelse(dup == TRUE, 
                       paste(gene,
                             observation,
                             sep = "_"), 
                       gene)) %>% 
  ungroup() %>% 
  select(-c(dup,
            observation))

# rename count data
count.data = count.data %>% 
  dplyr::rename(ensembl_gene_id = gene) %>% 
  full_join(tilapia.gene.names.list) %>% 
  dplyr::select(-c(ensembl_gene_id)) %>% 
  relocate(gene)

# save 
write.csv(count.data,
          '../Figures/DESEQ2/count.data.csv',
          row.names = F)

### load PCA data
## from behavior and hormones
Tank.PCA.Behavior.hormone.data = read.csv('../../IMC_Rscripts/Tank.PCA.Behavior.hormone.data.csv') %>% 
  separate(Tank,
           c('Tank',
             'Round')) %>% 
  mutate(tank = paste(Round,
                      Tank,
                      sep = '.')) %>% 
  dplyr::select(-c(X,
            Round,
            Tank))

# get behavior data
# need data.behavior.comp from 'initial_behavior_script.R'
load('../../IMC_Rscripts/data.behavior.comp.RData')
# create reduced dataframe
data.behavior.comp.reduce = data.behavior.comp %>% 
  dplyr::select(-c(Female.bower.male.not.present,
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
                Courtship = Court.total) %>% 
  mutate(Latency = 3000 - Latency)

# save as csv
write.csv(data.behavior.comp,
          '../Figures/Behavior/data.behavior.comp.csv',
          row.names = F)

## prepare sample data
sample.names.colData = sample.names %>% 
  mutate(Sample.Name = str_replace_all(Sample.Name,
                                       '-',
                                       '.')) %>% 
  column_to_rownames('Sample.Name') %>% 
  dplyr::select(c(Sex,
           Tank,
           Round)) %>% 
  mutate(Tank.ID = paste(Round,
                         Tank,
                         sep = ".")) %>% 
  dplyr::select(-c(Tank,
            Round))


#### Ensembl gene names ####
# ### get gene names for data
# ## load ensembl
# library(biomaRt)
# ###selecting biomart database
# ensembl = useMart('ensembl')
# # #list datasets
# # dataset = listDatasets(ensembl)
# ##select dataset
# #nile tilapia
# ensembl.tilapia = useDataset('oniloticus_gene_ensembl',
#                              mart=ensembl)
# 
# # # get list of all attributes
# # listAttributes(ensembl.tilapia) %>% View
# 
# #create attributes lists
# tilapia.attributes = c('external_gene_name',
#                        'ensembl_gene_id')
# 
# ##identify gene names tilapia genes
# # use gene names
# tilapia.gene.names = getBM(attributes = tilapia.attributes,
#                                     mart = ensembl.tilapia,
#                                     values = count.data$gene,
#                                     filter = 'ensembl_gene_id',
#                                     useCache = FALSE) # useCache has to do with version of R not being up to date?
# 
# ## create gene list 
# tilapia.gene.names.list = tilapia.gene.names %>% 
#   mutate(gene = ifelse(external_gene_name == '',
#                        ensembl_gene_id,
#                        external_gene_name)) %>% 
#   dplyr::select(-c(external_gene_name))
# 
# # save data
# write.csv(tilapia.gene.names.list,
#           'Ensembl.gene.names.csv')

#### Behavior Graphs ####
library(Hmisc)
library(corrplot)
#make heatmap 
#look at reduce behavior
data.behavior.comp.reduce.cor= rcorr(data.behavior.comp.reduce %>% 
                                                            dplyr::select(-c(Observation.id)) %>% 
                                                            as.matrix()) 

#graph
png('../Figures/Behavior/corrplot reduce behavior presentation.png',
    height = 5.25,
    width = 5.25)
corrplot(data.behavior.comp.reduce.cor$r, 
         type="upper", 
         order="hclust", 
         diag = FALSE,
         method = 'ellipse',
         col = colorRampPalette(c("blue", "white", "red"))(10), 
         p.mat = data.behavior.comp.reduce.cor$P, 
         sig.level = 0.05, 
         insig = "label_sig",
         tl.col = 'black',
         tl.cex = 1.5,
         cl.cex = 1)
dev.off()
# presentation
pdf('../Figures/Behavior/corrplot reduce behavior presentation.pdf',
    height = 10,
    width = 10)
corrplot(data.behavior.comp.reduce.cor$r, 
         type="upper", 
         order="hclust", 
         diag = FALSE,
         method = 'ellipse',
         col = colorRampPalette(c("blue", "white", "red"))(10), 
         p.mat = data.behavior.comp.reduce.cor$P, 
         sig.level = 0.05, 
         insig = "label_sig",
         tl.col = 'black',
         tl.cex = 1.5,
         cl.cex = 1)
dev.off()

## chart correlation
png('../Figures/Behavior/chart correlation reduce behavior presentation.png')
chart.Correlation(data.behavior.comp.reduce %>% 
                    dplyr::select(-c(Observation.id)),
                  histogram = TRUE,
                  pch = 19)
dev.off()

## court vs Time.together
time.court.lm = lm(Time.together ~ Courtship,
                           data=data.behavior.comp.reduce)


summary(time.court.lm)
anova(time.court.lm)

time.court.lm.pvalue = signif(summary(time.court.lm)[["coefficients"]]['Courtship','Pr(>|t|)'], digits = 1)

time.court.lm.rvalue = signif(sqrt(summary(time.court.lm)[["r.squared"]]), digits = 2)

# poster
library(jtools)
effect_plot(time.court.lm,
            pred = Courtship,
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=data.behavior.comp.reduce ,
             aes(x=Courtship,
                 y=Time.together),
             size = 3)+
  theme_classic()+
  theme(legend.position = 'null')+
  ylab('Time together (s)') +
  xlab('Courtship (per 10 min)')+
  theme(text = element_text(size = 15)) +
  annotate("text", x=20, y=50, label= paste0('p-value = ',
                                             time.court.lm.pvalue,
                                             ', R = ',
                                             time.court.lm.rvalue
  ))
ggsave('../Figures/Behavior/Time together vs courtship.pdf',
       height = 5.25,
       width = 5.25)

#### adjusted pvalue
library(RcmdrMisc)

#make heatmap 
#look at reduce behavior
data.behavior.comp.reduce.cor.adj= rcorr.adjust(data.behavior.comp.reduce %>% 
                                       dplyr::select(-c(Observation.id)) %>% 
                                       as.matrix(),
                                       type = 'spearman') 

# convert pvalues to numeric
class(data.behavior.comp.reduce.cor.adj$P ) <- "numeric"
  
#graph
pdf('../Figures/Behavior/corrplot reduce behavior presentation adj.pdf',
    height = 5.25,
    width = 5.25)
corrplot(data.behavior.comp.reduce.cor.adj$R$r, 
         type="upper", 
         order="hclust", 
         diag = FALSE,
         method = 'ellipse',
         col = colorRampPalette(c("blue", "white", "red"))(10), 
         p.mat = data.behavior.comp.reduce.cor.adj$P, 
         sig.level = 0.05, 
         insig = "label_sig",
         tl.col = 'black',
         tl.cex = 1.5,
         cl.cex = 1,
         pch = '*')
dev.off()

# paper
# annotate significance by hand
pdf('../Figures/Behavior/corrplot reduce behavior paper.pdf',
    height = 5.25,
    width = 5.25)
corrplot(data.behavior.comp.reduce.cor.adj$R$r, 
         type="upper", 
         order="hclust", 
         diag = FALSE,
         method = 'ellipse',
         col = colorRampPalette(c("blue", "white", "red"))(10), 
         tl.col = 'black',
         tl.cex = 1.5,
         cl.cex = 1)
dev.off()

# save results
data.behavior.comp.reduce.cor.adj.df = broom::tidy(data.behavior.comp.reduce.cor.adj$R) |> 
  mutate(FDR = p.adjust(p.value,
                        method = 'fdr'))

write.csv(data.behavior.comp.reduce.cor.adj.df,
          '../Figures/Behavior/data.behavior.comp.reduce.cor.adj.df.csv',
          row.names = F)




#### Behavior PCA ####
### create PCA for all data
## select behavior
data.behavior.comp.reduce.select = data.behavior.comp.reduce %>%
  column_to_rownames('Observation.id')

#only use numerical variables (not caegorical!) 
data.beh.pca = prcomp(data.behavior.comp.reduce.select,
                  scale = TRUE)

#check PCs
summary(data.beh.pca)


# PC variance
percentVar.beh <- data.frame(variance = 100*data.beh.pca$sdev^2 / sum( data.beh.pca$sdev^2),
                         PC = data.beh.pca[["rotation"]] %>% colnames(),
                         order = rep(1:length(data.beh.pca[["rotation"]] %>% colnames()))) 

## graphing 
##scree plot
percentVar.beh %>% 
  ggplot(aes(y=variance,
             x = fct_reorder(PC,
                             order),
             group = 1)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  ggtitle('PCA screeplot') 
ggsave('../Figures/Behavior/pca/Scree plot behavior.png',
       height = 10,
       width = 10)
#plot pca
# PC 1 vs PC 2
data.beh.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('ID') %>%
  ggplot()   +
  geom_label(data = data.beh.pca$rotation %>%
               data.frame() %>% 
               rownames_to_column('Type'),
             aes(x = PC1*max(data.beh.pca$x),
                 y = PC2*max(data.beh.pca$x) + 0.1,
                 label = Type ),
             size = 3,
             color = 'black') +
  geom_segment(data = data.beh.pca$rotation %>% 
                 data.frame(),
               aes(x = 0, 
                   y = 0, 
                   xend = PC1*max(data.beh.pca$x),
                   yend = PC2*max(data.beh.pca$x)),
               arrow = arrow(length = unit(0.5, "cm")),
               color = 'red') + 
  geom_point(aes(x = PC1,
                y = PC2),
            size=3) +
  xlab(paste0("PC1: ",round(percentVar.beh[1,1]),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar.beh[2,1]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/Behavior/pca/PCA behavior 1 vs 2.png',
       height = 10,
       width = 10)

# PC 1 vs PC 3
data.beh.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('ID') %>%
  ggplot()   +
  geom_label(data = data.beh.pca$rotation %>%
               data.frame() %>% 
               rownames_to_column('Type'),
             aes(x = PC1*max(data.beh.pca$x),
                 y = PC3*max(data.beh.pca$x) + 0.1,
                 label = Type ),
             size = 3,
             color = 'black') +
  geom_segment(data = data.beh.pca$rotation %>% 
                 data.frame(),
               aes(x = 0, 
                   y = 0, 
                   xend = PC1*max(data.beh.pca$x),
                   yend = PC3*max(data.beh.pca$x)),
               arrow = arrow(length = unit(0.5, "cm")),
               color = 'red') + 
  geom_point(aes(x = PC1,
                 y = PC3),
             size=3) +
  xlab(paste0("PC1: ",round(percentVar.beh[1,1]),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar.beh[3,1]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/Behavior/pca/PCA behavior 1 vs 3.png',
       height = 10,
       width = 10)

# PC 1 vs PC 3
data.beh.pca$x %>% 
  data.frame() %>% 
  rownames_to_column('ID') %>%
  ggplot()   +
  geom_label(data = data.beh.pca$rotation %>%
               data.frame() %>% 
               rownames_to_column('Type'),
             aes(x = PC2*max(data.beh.pca$x),
                 y = PC3*max(data.beh.pca$x) + 0.1,
                 label = Type ),
             size = 3,
             color = 'black') +
  geom_segment(data = data.beh.pca$rotation %>% 
                 data.frame(),
               aes(x = 0, 
                   y = 0, 
                   xend = PC2*max(data.beh.pca$x),
                   yend = PC3*max(data.beh.pca$x)),
               arrow = arrow(length = unit(0.5, "cm")),
               color = 'red') + 
  geom_point(aes(x = PC2,
                 y = PC3),
             size=3) +
  xlab(paste0("PC2: ",round(percentVar.beh[2,1]),"% variance")) +
  ylab(paste0("PC3: ",round(percentVar.beh[3,1]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/Behavior/pca/PCA behavior 2 vs 3.png',
       height = 10,
       width = 10)


#### DESQ2 analysis ####
#https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
### prepare data


## setup count data
count.data.deseq = count.data %>% 
  column_to_rownames('gene')

# reorder column names
count.data.deseq = count.data.deseq[, rownames(sample.names.colData)]

# check that columns and rows match
all(rownames(sample.names.colData) %in% colnames(count.data.deseq))
all(rownames(sample.names.colData) == colnames(count.data.deseq))

### create DESEQ dataset
dds <- DESeqDataSetFromMatrix(countData = count.data.deseq,
                              colData = sample.names.colData,
                              design = ~ Sex)
# check data
dds

# remove rows with low counts
#need to have a count of at least 2 and be present in at least half of samples
# 10 in more than 90%
keep = rowSums(counts(dds) >= 2) >= nrow(sample.names.colData)/2

# remove low expressed genes
# go from 33427 to 13343
dds = dds[keep,]


# #keep 50% most variable genes
# var_genes <- apply(dds@assays@data@listData[["counts"]],
#                    1, 
#                    var)
# 
# #create list of top 50% most variable genes
# num.genes.to.keep = length(var_genes)/2 
# num.genes.to.keep = ceiling(num.genes.to.keep)
# select_var = names(sort(var_genes, 
#                          decreasing=TRUE))[1:num.genes.to.keep]
# 
# dds = dds[select_var,]

# #check data
# dds

### run DESEQ
dds <- DESeq(dds)

## check results
res <- results(dds)
# res


#### DESQ2 graph ####
## MA plot
png('../Figures/DESEQ2/MA plot sex initial.png')
plotMA(res,
       ylim=c(-2,2))
dev.off()

## cooks distance
png('../Figures/DESEQ2/Cooks distance.png')
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

## volcano plot
#create dataframe
res.df = res@listData %>%
  as.data.frame() 
#add genes
rownames(res.df) = res@rownames
#add significance
#label with gene names
res.df = res.df %>%
  rownames_to_column('gene') %>%
  mutate(Sig.name = ifelse(padj <= 0.1,
                      gene,
                      NA),
         Sig = ifelse(padj <= 0.1,
                      'Significant',
                      'Not-significant'))

# graph
res.df %>%
  filter(!is.na(Sig)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             color = Sig,
             label = Sig.name)) +
  geom_hline(yintercept = -log10(0.1),
             linetype = 'dashed')+
  geom_point(size = 5)   +
  geom_label() +
  theme_classic() +
  ggtitle('Volcano plot Male vs Female, padj < 0.1') +
  xlim(-3.5,3.5)+ 
  theme(text = element_text(size = 30),
        legend.position = 'none') +
  scale_color_manual(values = c('black',
                                'red'))
ggsave('../Figures/DESEQ2/Volcano plot sex initial.png',
       height = 10,
       width = 10)

# poster
res.df %>%
  filter(!is.na(Sig)) %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             color = Sig,
             label = Sig.name)) +
  geom_hline(yintercept = -log10(0.1),
             linetype = 'dashed')+
  geom_point(size = 5)   +
  theme_classic() +
  # ggtitle('Volcano plot Male vs Female') +
  xlim(-3.5,3.5)+ 
  theme(text = element_text(size = 20),
        legend.position = 'none') +
  scale_color_manual(values = c('black',
                                'red'))
ggsave('../Figures/DESEQ2/Volcano plot sex initial poster.pdf',
       height = 5.25,
       width = 5.25)

# paper
res.df %>%
  filter(!is.na(Sig)) %>% 
  mutate(Sig.color = case_when(Sig == 'Significant' & log2FoldChange < 0~ 'Female',
                               Sig == 'Significant' & log2FoldChange > 0 ~ 'Male',
                               TRUE ~ 'Other')) |> 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             color = Sig.color,
             label = Sig.name)) +
  geom_hline(yintercept = -log10(0.1),
             linetype = 'dashed')+
  geom_point(size = 5)   +
  theme_classic() +
  # ggtitle('Volcano plot Male vs Female') +
  xlim(-3.5,3.5)+ 
  theme(text = element_text(size = 15),
        legend.position = 'none') +
  scale_color_manual(values = c('Other' = 'grey',
                                'Male' = '#E69F00',
                                'Female' = '#008080'))
ggsave('../Figures/DESEQ2/Volcano plot sex initial paper.pdf',
       height = 5.25,
       width = 5.25)

#### DESQ2 multifactor ####
# ## create new dataframe
# ## update design to include tank
# ddsMF <- DESeqDataSetFromMatrix(countData = count.data.deseq,
#                                        colData = sample.names.colData,
#                                        design = ~ Sex + Tank.ID)
# 
# # remove rows with low counts
# #need to have a count of at least 2 and be present in at least half of samples
# keep = rowSums(counts(ddsMF) >= 2) >= nrow(sample.names.colData)/2
# 
# # remove low expressed genes
# # go from 33427 to 13343
# ddsMF = ddsMF[keep,]
# 
# 
# #keep 50% most variable genes
# var_genes <- apply(ddsMF@assays@data@listData[["counts"]],
#                    1, 
#                    var)
# 
# #create list of top 50% most variable genes
# num.genes.to.keep = length(var_genes)/2 
# num.genes.to.keep = ceiling(num.genes.to.keep)
# select_var = names(sort(var_genes, 
#                         decreasing=TRUE))[1:num.genes.to.keep]
# 
# ddsMF = ddsMF[select_var,]
# 
# #check data
# ddsMF
# 
# ## run deseq
# ddsMF <- DESeq(ddsMF)
# 
# # get results of sex comparison
# resMF <- results(ddsMF,
#                      contrast=c("Sex", "M", "F"))
# 
# ### graph results
# ## MA plot
# png('../Figures/DESEQ2/MA plot sex initial multifactor.png')
# plotMA(resMF, 
#        ylim=c(-2,2))
# dev.off()
# 
# ## volcano plot
# #create dataframe
# resMF.df = resMF@listData %>% 
#   as.data.frame()
# #add genes
# rownames(resMF.df) = resMF@rownames
# #add significance
# #label with gene names
# resMF.df = resMF.df %>% 
#   rownames_to_column('gene') %>% 
#   mutate(Sig.name = ifelse(padj <= 0.1 & abs(log2FoldChange) >= 1.5,
#                            gene,
#                            NA),
#          Sig = ifelse(padj <= 0.1,
#                       'Significant',
#                       'Not-significant'))
# 
# #graph
# resMF.df %>%
#   filter(!is.na(Sig)) %>% 
#   ggplot(aes(x = log2FoldChange,
#              y = -log10(padj),
#              color = Sig)) + 
#   geom_hline(yintercept = -log10(0.1),
#              linetype = 'dashed')+ 
#   geom_vline(xintercept = 1.5,
#              linetype = 'dashed') +
#   geom_vline(xintercept = -1.5,
#              linetype = 'dashed') +
#   geom_point()  +
#   geom_label(aes(label = Sig.name),
#              vjust = -0.5,
#              hjust = 0.9,
#              label.size =  0.01) +
#   theme_classic() +
#   xlim(-max(abs(resMF.df$log2FoldChange), na.rm=T), 
#        max(abs(resMF.df$log2FoldChange), na.rm = T))+
#   ggtitle('Volcano plot Male vs Female, padj < 0.1 multifactor') 
# ggsave('../Figures/DESEQ2/Volcano plot sex initial multifactor.png',
#        height = 10,
#        width = 10)

#### data transformation ####
# ## log normal transformation
# # this gives log2(n + 1)
# ntd = normTransform(dds)
# 
# 
# #plot
# png('../Figures/DESEQ2/Log normalization.png')
# meanSdPlot(assay(ntd))
# dev.off()

#VST 
vsd = vst(dds,
          blind = F)

## save results
# save(vsd,
#      file = "vsd.RData")
# #load data
# load('vsd.RData')




#plot
# png('../Figures/DESEQ2/Vst normalization.png')
# meanSdPlot(assay(vsd))
# dev.off()
#
# ## graph normalized expression per samples
# data.frame(assay(vsd)) %>%
#   mutate(
#     Gene_id = row.names(vsd)
#   ) %>%
#   pivot_longer(-Gene_id) %>%
#   ggplot(., aes(x = name, y = value)) +
#   geom_violin() +
#   geom_point() +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text( angle = 90)
#   ) +
#   ylim(0, NA) +
#   labs(
#     title = "Normalized Expression",
#     x = "sample",
#     y = "normalized expression"
#   )
# ggsave('../Figures/DESEQ2/Normalized Expression per sample.png',
#        height = 10,
#        width = 10)
#
# # graph expression per samples
# data.frame(assay(dds)) %>%
#   mutate(
#     Gene_id = row.names(dds)
#   ) %>%
#   pivot_longer(-Gene_id) %>%
#   ggplot(., aes(x = name, y = value)) +
#   geom_violin() +
#   geom_point() +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text( angle = 90)
#   ) +
#   ylim(0, NA) +
#   labs(
#     title = "Expression",
#     x = "sample",
#     y = "expression"
#   )
# ggsave('../Figures/DESEQ2/Expression per sample.png',
#        height = 10,
#        width = 10)

# ### heatmap of samples
# ## set up data
# # get order of rows
# select = order(rowMeans(counts(dds,
#                                 normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# #create dataframe of colData
# df.coldata = as.data.frame(colData(dds)[,c("Sex","Tank.ID")])
# ## create heatmap
# png('../Figures/DESEQ2/Heatmap of genes and samples.png')
# pheatmap(assay(vsd)[select,],
#          cluster_rows=TRUE,
#          show_rownames=FALSE,
#          cluster_cols=TRUE,
#          annotation_col=df.coldata)
# dev.off()
#

### create input matrix from normalized data
vsd.mat <- assay(vsd) %>%
  t()

#### PCA ####
# calculate the variance for each gene
rv <- rowVars(assay(vsd))
#graph variance for each gene
data.frame(variance = rv) %>% 
  ggplot(aes(variance)) +
  geom_histogram() +
  theme_classic()
ggsave('../Figures/DESEQ2/PCA/Gene variance histogram.png')


# select the 1000 genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(5000, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

# create dataframe
percentVar.df = data.frame(percentVar = percentVar,
           PC = paste0("PC",
                       seq(1:length(percentVar)))) %>% 
       mutate(sum = cumsum(percentVar))

# compare PC to behavior data
# PC 1, 2, 3
data.behavior.comp.reduce.pca = data.behavior.comp.reduce %>%
  full_join(sample.names.colData %>% 
              rownames_to_column('ID') %>% 
              separate(Tank.ID,
                       into = c('round',
                                'tank'),
                       remove = F) %>% 
              mutate(Observation.id = paste(tank,
                                            round,
                                            sep = '_'))) %>% 
  full_join(pca$x %>% 
              as.data.frame() %>% 
              rownames_to_column('ID')) 

#graph PCA var
percentVar.df %>% 
  filter(sum < 0.8) |> 
  ggplot(aes(x = reorder(PC,
                         -percentVar),
             y = percentVar)) +
  geom_point() +
  geom_segment(aes(x=reorder(PC,
                             -percentVar), 
                   xend=reorder(PC,
                                -percentVar), 
                   y=0, 
                   yend=percentVar)) +
  theme_classic() +
  ggtitle(paste0(percentVar.df |> 
                   filter(sum < 0.8) |> 
                   pull(sum) |> 
                   max() |> 
                   round (digits = 2),
                 ' variance'))
ggsave('../Figures/DESEQ2/PCA/PC variance.png')


#plot pca
# PC 1 vs PC 2
pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  full_join(sample.names.colData %>% 
              rownames_to_column('ID')) %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",round(100*percentVar[1]),"% variance")) +
  ylab(paste0("PC2: ",round(100*percentVar[2]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/DESEQ2/PCA/PCA 1 vs 2 sex.png')

# PC 1 vs PC 2
#presentation
pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  full_join(sample.names.colData %>% 
              rownames_to_column('ID')) %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = Sex)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",round(100*percentVar[1]),"% variance")) +
  ylab(paste0("PC2: ",round(100*percentVar[2]),"% variance")) + 
  coord_fixed() +
  theme_classic()+ 
  theme(text = element_text(size = 30)) +
  scale_color_manual(values = c('black',
                                'red'))
ggsave('../Figures/DESEQ2/PCA/PCA 1 vs 2 sex presentation.png',
       height = 10,
       width = 10)

# PC 1 vs PC 2
#poster
pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  full_join(sample.names.colData %>% 
              rownames_to_column('ID')) %>%
  ggplot(aes(x = PC1,
             y = PC2,
             color = Sex)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",round(100*percentVar[1]),"% variance")) +
  ylab(paste0("PC2: ",round(100*percentVar[2]),"% variance")) + 
  coord_fixed() +
  theme_classic()+ 
  theme(text = element_text(size = 30)) +
  scale_color_manual(values = c('black',
                                'red'))
ggsave('../Figures/DESEQ2/PCA/PCA 1 vs 2 sex poster.pdf',
       height = 5.25,
       width = 5.25)

# PC 1 vs PC 2
# paper
pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  full_join(sample.names.colData %>% 
              rownames_to_column('ID')) %>% 
  ggplot(aes(x = PC1,
             y = PC2,
             color = Sex)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",round(100*percentVar[1]),"% variance")) + 
  ylab(paste0("PC2: ",round(100*percentVar[2]),"% variance")) + 
  coord_fixed() +
  theme_classic()+ 
  theme(text = element_text(size = 15)) +
  scale_color_manual(values = c('Other' = 'grey',
                                'M' = '#E69F00',
                                'F' = '#008080'))+
  theme( legend.position = c(0.95, 0.95),      
    legend.justification = c("right", "top"), 
    legend.background = element_rect(
      fill = "white",                   
      color = "black",                   
      linewidth = 0.5))
ggsave('../Figures/DESEQ2/PCA/PCA 1 vs 2 sex paper.pdf',
       height = 5.25,
       width = 5.25)

# PC 1 vs PC 3
pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  full_join(sample.names.colData %>% 
              rownames_to_column('ID')) %>%
  ggplot(aes(x = PC1,
             y = PC3,
             color = Sex)) +
  geom_point(size=3)  +
  xlab(paste0("PC1: ",round(100*percentVar[1]),"% variance")) +
  ylab(paste0("PC3: ",round(100*percentVar[3]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/DESEQ2/PCA/PCA 1 vs 3 sex.png')

#plot pca
# PC 2 vs PC 3
pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  full_join(sample.names.colData %>% 
              rownames_to_column('ID')) %>%
  ggplot(aes(x = PC2,
             y = PC3,
             color = Sex)) +
  geom_point(size=3) +
  xlab(paste0("PC2: ",round(100*percentVar[2]),"% variance")) +
  ylab(paste0("PC3: ",round(100*percentVar[3]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/DESEQ2/PCA/PCA 2 vs 3 sex.png')

# PC 1 vs PC 4
pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column('ID') %>% 
  full_join(sample.names.colData %>% 
              rownames_to_column('ID')) %>%
  ggplot(aes(x = PC1,
             y = PC4,
             color = Sex)) +
  geom_point(size=3)  +
  xlab(paste0("PC1: ",round(100*percentVar[1]),"% variance")) +
  ylab(paste0("PC4: ",round(100*percentVar[4]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/DESEQ2/PCA/PCA 1 vs 4 sex.png')


# PC 1 vs PC 3
data.behavior.comp.reduce.pca %>%
  ggplot(aes(x = PC1,
             y = PC3,
             color = log(Time.together))) +
  geom_point(size=3)  +
  xlab(paste0("PC1: ",round(100*percentVar[1]),"% variance")) +
  ylab(paste0("PC3: ",round(100*percentVar[3]),"% variance")) + 
  coord_fixed() +
  theme_classic()
ggsave('../Figures/DESEQ2/PCA/PCA 1 vs 3 time.together log.png',)

## chart correlation
png('../Figures/DESEQ2/PCA/chart correlation PCA and behavior.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 480)
chart.Correlation(data.behavior.comp.reduce.pca %>% 
                    select(c(Time.together,
                             Courtship,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone,
                             Estradiol,
                             PC1,
                             PC2,
                             PC3)),
                  histogram = TRUE,
                  pch = 19)
dev.off()

#males
png('../Figures/DESEQ2/PCA/chart correlation PCA and behavior males.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 480)
chart.Correlation(data.behavior.comp.reduce.pca %>% 
                    filter(Sex == 'M') %>% 
                    select(c(Time.together,
                             Courtship,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone,
                             Estradiol,
                             PC1,
                             PC2,
                             PC3)),
                  histogram = TRUE,
                  pch = 19)
mtext("Males", side=3, line=3)
dev.off()

pdf('../Figures/DESEQ2/PCA/chart correlation PCA and behavior males.pdf',
    height = 10,
    width = 10)
chart.Correlation(data.behavior.comp.reduce.pca %>% 
                    filter(Sex == 'M') %>% 
                    select(c(Time.together,
                             Courtship,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone,
                             Estradiol,
                             PC1,
                             PC2,
                             PC3)),
                  histogram = TRUE,
                  pch = 19)
mtext("Males", side=3, line=3)
dev.off()

# females
png('../Figures/DESEQ2/PCA/chart correlation PCA and behavior females.png',
    height = 10,
    width = 10,
    units = 'in',
    res = 480)
chart.Correlation(data.behavior.comp.reduce.pca %>% 
                    filter(Sex == 'F') %>% 
                    select(c(Time.together,
                             Courtship,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone,
                             Estradiol,
                             PC1,
                             PC2,
                             PC3)),
                  histogram = TRUE,
                  pch = 19)
mtext("Females", side=3, line=3)
dev.off()

pdf('../Figures/DESEQ2/PCA/chart correlation PCA and behavior females.pdf',
    height = 10,
    width = 10)
chart.Correlation(data.behavior.comp.reduce.pca %>% 
                    filter(Sex == 'F') %>% 
                    select(c(Time.together,
                             Courtship,
                             Latency,
                             GSI.male,
                             GSI.female,
                             Testosterone,
                             Estradiol,
                             PC1,
                             PC2,
                             PC3)),
                  histogram = TRUE,
                  pch = 19)
mtext("Females", side=3, line=3)
dev.off()

## compare male and female PC's
library(ggpubr)
#PC1
data.behavior.comp.reduce.pca %>% 
  pivot_wider(id_cols = 'Observation.id',
              names_from = 'Sex',
              names_prefix = 'PC1',
              values_from = 'PC1') %>%
  ggplot(aes(x = PC1M,
             y = PC1F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  xlab('Males PC1') +
  ylab('Females PC1') +
  ggtitle('Males vs Females PC1 RNAseq')+
  stat_cor(label.x.npc = 'left')
ggsave('../Figures/DESEQ2/PCA/Males vs Females PC1 RNAseq.png',
    height = 5,
    width = 5)  
#PC2
data.behavior.comp.reduce.pca %>% 
  pivot_wider(id_cols = 'Observation.id',
              names_from = 'Sex',
              names_prefix = 'PC2',
              values_from = 'PC2') %>%
  ggplot(aes(x = PC2M,
             y = PC2F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  xlab('Males PC2') +
  ylab('Females PC2') +
  ggtitle('Males vs Females PC2 RNAseq')+
  stat_cor(label.x.npc = 'left')
ggsave('../Figures/DESEQ2/PCA/Males vs Females PC2 RNAseq.png',
       height = 5,
       width = 5)  
#PC3
data.behavior.comp.reduce.pca %>% 
  pivot_wider(id_cols = 'Observation.id',
              names_from = 'Sex',
              names_prefix = 'PC3',
              values_from = 'PC3') %>%
  ggplot(aes(x = PC3M,
             y = PC3F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  xlab('Males PC3') +
  ylab('Females PC3') +
  ggtitle('Males vs Females PC3 RNAseq')+
  stat_cor(label.x.npc = 'left')
ggsave('../Figures/DESEQ2/PCA/Males vs Females PC3 RNAseq.png',
       height = 5,
       width = 5)  


## paper
#PC1
data.behavior.comp.reduce.pca %>% 
  pivot_wider(id_cols = 'Observation.id',
              names_from = 'Sex',
              names_prefix = 'PC1',
              values_from = 'PC1') %>%
  ggplot(aes(x = PC1M,
             y = PC1F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  xlab('Males PC1') +
  ylab('Females PC1') +
  ggtitle('Males vs Females PC1')+
  stat_cor(label.x.npc = 'left')
ggsave('../Figures/DESEQ2/PCA/Males vs Females PC1 RNAseq paper.pdf',
       height = 5,
       width = 5)  
#PC2
data.behavior.comp.reduce.pca %>% 
  pivot_wider(id_cols = 'Observation.id',
              names_from = 'Sex',
              names_prefix = 'PC2',
              values_from = 'PC2') %>%
  ggplot(aes(x = PC2M,
             y = PC2F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  xlab('Males PC2') +
  ylab('Females PC2') +
  ggtitle('Males vs Females PC2')+
  stat_cor(label.x.npc = 'left')
ggsave('../Figures/DESEQ2/PCA/Males vs Females PC2 RNAseq paper.pdf',
       height = 5,
       width = 5)  
#PC3
data.behavior.comp.reduce.pca %>% 
  pivot_wider(id_cols = 'Observation.id',
              names_from = 'Sex',
              names_prefix = 'PC3',
              values_from = 'PC3') %>%
  ggplot(aes(x = PC3M,
             y = PC3F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  xlab('Males PC3') +
  ylab('Females PC3') +
  ggtitle('Males vs Females PC3')+
  stat_cor(label.x.npc = 'left')
ggsave('../Figures/DESEQ2/PCA/Males vs Females PC3 RNAseq paper.pdf',
       height = 5,
       width = 5)  

#### save PCA ####
save(data.behavior.comp.reduce.pca,
     file = "data.behavior.comp.reduce.pca.RData")

save(percentVar.df,
     file = "percentVar.df.RData")

#### distance matrix ####
## check data for outliers
#Group data in a dendogram to check outliers
sampleTree = hclust(dist(vsd.mat), method = "average")

#If you want to save this plot in a pdf file, do not comment the line below:
png('../Figures/DESEQ2/Sample clustering.png')
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#Plot a line showing the cut-off
abline(h = 50, col = "red") #This value of 31000 was chosen based on my data, you need to check the best value to your data
dev.off()

## heatmap of sample distance
# create distance matrix
sampleDists = dist(t(assay(vsd)))
sampleDistMatrix = as.matrix(sampleDists)

# add row and column names
rownames(sampleDistMatrix) <- paste(vsd$Tank.ID, 
                                    vsd$Sex, 
                                    sep="-")
colnames(sampleDistMatrix) <- paste(vsd$Tank.ID, 
                                    vsd$Sex, 
                                    sep=".")



## convert distance matrix to data frame 
# remove individual
# add sample info
# add same tank
sampleDist.df = data.frame(from=colnames(sampleDistMatrix)[col(sampleDistMatrix)], 
                           to=rownames(sampleDistMatrix)[row(sampleDistMatrix)], 
                           dist=c(sampleDistMatrix)) %>% 
  filter(dist != 0) %>% 
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.')) %>% 
              dplyr::select(c(Round,
                       Tank,
                       Sex)) %>% 
              mutate(from = paste(Round,
                                  Tank,
                                  Sex,
                                  sep = '.'),
                     from.tank = paste(Round,
                                       Tank,
                                       sep = '.'),
                     from.sex = Sex) %>% 
              dplyr::select(-c(Round,
                        Tank,
                        Sex))) %>% 
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.')) %>% 
              dplyr::select(c(Round,
                       Tank,
                       Sex)) %>% 
              mutate(to.tank = paste(Round,
                                     Tank,
                                     sep = '.'),
                     to = paste(to.tank,
                                  Sex,
                                  sep = '-'),
                     to.sex = Sex) %>% 
              dplyr::select(-c(Round,
                        Tank,
                        Sex))) %>% 
    mutate(same = ifelse(from.tank == to.tank,
                         'same',
                         'other')) %>% 
  mutate(to = str_replace_all(to,
                                     '-',
                                     '.')) 

### create distance statistic
sampleDist.df.sum = sampleDist.df %>% 
  group_by(from,
           from.sex,
           to.sex) %>% 
  mutate(Mean = mean(dist),
         Sd = sd(dist)) %>% 
  ungroup() %>% 
  mutate(Zscore = (dist - Mean)/Sd)

### compare to PC data from hormones and behavior
sampleDist.df.sum.tank = sampleDist.df.sum %>% 
  filter(same == 'same') %>% 
  full_join(Tank.PCA.Behavior.hormone.data %>% 
              dplyr::select(c(tank,
                       PC1)) %>% 
              dplyr::rename(from.tank = tank))

#create dataframe of colData
df.coldata = as.data.frame(colData(dds)[,c("Sex","Tank.ID")])

#### graph distance matrix ####
### heat map
# create color list
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#graph
png('../Figures/DESEQ2/Heatmap of distance matrix samples.png')
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, 
         show_colnames = FALSE,
         annotation_col = df.coldata %>% 
           rownames_to_column('ID') %>% 
           mutate(row = substring(ID, 
                                  4)) %>% 
           column_to_rownames('row') %>% 
           select(-c(ID)))
dev.off()

## graph same vs. other
#males to females
sampleDist.df %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = dist,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Dist from males to all females')
ggsave('../Figures/DESEQ2/Dist from males to all females.png')

#females to males
sampleDist.df %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = dist,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Dist from females to all males')
ggsave('../Figures/DESEQ2/Dist from females to all males.png')

### distance stats 
### z score
## graph same vs. other
#males to females
sampleDist.df.sum %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Dist from males to all females zscore')
ggsave('../Figures/DESEQ2/Dist from males to all females zscore.png')

#females to males
sampleDist.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Dist from females to all males zscore')
ggsave('../Figures/DESEQ2/Dist from females to all males zscore.png')


#females to males
sampleDist.df.sum %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Dist from males to all females zscore') 
ggsave('../Figures/DESEQ2/Dist from females to all males zscore all.png')

#females to males
sampleDist.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Dist from females to all males zscore') 
ggsave('../Figures/DESEQ2/Dist from females to all males zscore all.png')


### compare distance zscore with PC1 of behavior and hormones
sampleDist.df.sum.tank %>% 
  ggplot(aes(x = PC1,
             y = Zscore,
             color = from.sex)) +
  geom_point() +
  geom_text(aes(label = from.tank),
            nudge_y = 0.1) +
  geom_smooth(method = 'lm') +
  theme_classic() +
  ggtitle('Dist zscore vs PC1 of behavior and hormones') +
  facet_grid(~from.sex)
ggsave('../Figures/DESEQ2/Dist zscore vs PC1 of behavior and hormones.png')




#### correlation matrix ####

### create correlation matrix 
sampleCorr = rcorr(assay(vsd))

## graph heatmap of correlations
# png('../Figures/DESEQ2/Heatmap of correlation matrix samples.png')
# pheatmap(sampleCorr$r,
#          col=colors,
#          show_colnames = FALSE,
# annotation_col = data.frame( ID = sampleCorr$r %>%
#   rownames()) %>%
#   mutate(row = substring(ID,
#                          4)) %>%
#   separate(row,
#            c('round',
#              'tank',
#              'sex')) %>%
#   mutate(tank = paste(round,
#                       tank,
#                       sep = '.')) %>%
#   column_to_rownames('ID') %>%
#   select(-c(round)))
# dev.off()

## convert distance matrix to data frame 
# remove individual
# add sample info
# add same tank
sampleCorr.df = data.frame(from=colnames(sampleCorr$r)[col(sampleCorr$r)], 
                           to=rownames(sampleCorr$r)[row(sampleCorr$r)], 
                           corr=c(sampleCorr$r)) %>% 
  mutate(from = substring(from, 
                         4),
         to = substring(to, 
                         4)) %>% 
  filter(corr != 1) %>% 
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.')) %>% 
              select(c(Round,
                       Tank,
                       Sex)) %>% 
              mutate(from = paste(Round,
                                  Tank,
                                  Sex,
                                  sep = '.'),
                     from.tank = paste(Round,
                                       Tank,
                                       sep = '.'),
                     from.sex = Sex) %>% 
              select(-c(Round,
                        Tank,
                        Sex))) %>% 
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.')) %>% 
              select(c(Round,
                       Tank,
                       Sex)) %>% 
              mutate(to.tank = paste(Round,
                                     Tank,
                                     sep = '.'),
                     to = paste(to.tank,
                                Sex,
                                sep = '.'),
                     to.sex = Sex) %>% 
              select(-c(Round,
                        Tank,
                        Sex))) %>% 
  mutate(same = ifelse(from.tank == to.tank,
                       'same',
                       'other'))

### create distance statistic
sampleCorr.df.sum = sampleCorr.df %>% 
  group_by(from,
           from.sex,
           to.sex) %>% 
  mutate(Mean = mean(corr),
         Sd = sd(corr)) %>% 
  ungroup() %>% 
  mutate(Zscore = (corr - Mean)/Sd)

### compare to PC data from hormones and behavior
sampleCorr.df.sum.tank = sampleCorr.df.sum %>% 
  filter(same == 'same') %>% 
  full_join(Tank.PCA.Behavior.hormone.data %>% 
              select(c(tank,
                       PC1)) %>% 
              dplyr::rename(from.tank = tank))



#### graph correlation matrix ####
## graph same vs. other
#males to females
sampleCorr.df %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = corr,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from males to all females')
ggsave('../Figures/DESEQ2/Corr from males to all females.png')

#females to males
sampleCorr.df %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = corr,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from females to all males')
ggsave('../Figures/DESEQ2/Corr from females to all males.png')

### correlation stats 
### z score
## graph same vs. other
#males to females
sampleCorr.df.sum %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from males to all females zscore')
ggsave('../Figures/DESEQ2/Corr from males to all females zscore.png')

#females to males
sampleCorr.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from females to all males zscore')
ggsave('../Figures/DESEQ2/Corr from females to all males zscore.png')


#males to females
sampleCorr.df.sum %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Corr from males to all females zscore') 
ggsave('../Figures/DESEQ2/Corr from females to all males zscore all.png')

#females to males
sampleCorr.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Corr from females to all males zscore') 
ggsave('../Figures/DESEQ2/Corr from females to all males zscore all.png')

#females to males
#presentation
sampleCorr.df.sum %>% 
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
  ggtitle('Correlation score females') + 
    theme(text = element_text(size = 30),
          legend.position = 'none') +
  xlab('Pair') +
  ylab('Zscore female correlation')
ggsave('../Figures/DESEQ2/Corr from females to all males zscore all presentation.png',
       height = 10,
       width = 10)


### compare correlation zscore with PC1 of behavior and hormones
sampleCorr.df.sum.tank %>% 
  ggplot(aes(x = PC1,
             y = Zscore,
             color = from.sex)) +
  geom_point() +
  geom_text(aes(label = from.tank),
            nudge_y = 0.1) +
  geom_smooth(method = 'lm') +
  theme_classic() +
  ggtitle('Corr zscore vs PC1 of behavior and hormones') +
  facet_grid(~from.sex)
ggsave('../Figures/DESEQ2/Corr zscore vs PC1 of behavior and hormones.png')

### compare correlation to distance
sample.df.sum = sampleCorr.df.sum.tank %>%
  select(c(from,
           to,
           corr,
           same,
           Zscore,
           PC1)) %>% 
  dplyr::rename(Zscore.corr = Zscore) %>% 
  full_join(sampleDist.df.sum.tank %>%
              select(c(from,
                       to,
                       dist,
                       same,
                       Zscore,
                       PC1)) %>% 
              dplyr::rename(Zscore.dist = Zscore))


## graph comparison
sample.df.sum %>% 
  ggplot(aes(x = corr,
             y = dist)) +
  geom_point() +
  theme_classic()
ggsave('../Figures/DESEQ2/Corr vs dist.png')

#zscore
sample.df.sum %>% 
  ggplot(aes(x = Zscore.corr,
             y = Zscore.dist)) +
  geom_point() +
  theme_classic()
ggsave('../Figures/DESEQ2/Corr vs dist zscores.png')


#### correlation synchronization ####
sampleCorr.df.sum.sync = sampleCorr.df.sum %>% 
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
            corr)) %>% 
  pivot_wider(names_from = from.sex,
              values_from = Zscore,
              names_prefix = 'Zscore.') %>% 
  full_join(data.behavior.comp)

## compare zscores
sampleCorr.df.sum.sync %>% 
  ggplot(aes(x = Zscore.M,
             y= Zscore.F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  ggtitle('Correlation zscore comparison')
ggsave('../Figures/DESEQ2/Correlation.score/Correlation zscore comparison sexes.png',
       height = 10,
       width = 10)

# check correlation
png('../Figures/DESEQ2/Corr zscore and behavior chart correlation.png')
chart.Correlation(sampleCorr.df.sum.sync %>% 
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

# lm
Zscore.F.court.lm = lm(Zscore.F ~ Court.total,
         data=sampleCorr.df.sum.sync)

summary(Zscore.F.court.lm)
anova(Zscore.F.court.lm)


#graph glm with outlier on graph
# presentation
library(jtools)
effect_plot(Zscore.F.court.lm, 
            pred = Court.total, 
            interval = TRUE,
            rug=TRUE) +
  geom_point(data=sampleCorr.df.sum.sync ,
             aes(x=Court.total,
                 y=Zscore.F),
             size = 5)+
  theme_bw()+
  theme(legend.position = 'null')+
  ylab('Zscore female correlation') +
  xlab('Male courtship (counts)')+ 
  theme(text = element_text(size = 30))
ggsave('../Figures/DESEQ2/Female corr zscore vs courtship with outlier presentation.png',
       height = 10,
       width = 10)

#### WGCNA ####
#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html
#https://bioinformaticsworkbook.org/tutorials/wgcna.html
#https://rstudio-pubs-static.s3.amazonaws.com/687551_ed469310d8ea4652991a2e850b0018de.html



## optional 
allowWGCNAThreads(12)
# check number of threads
WGCNAnThreads()

### pick soft power thresholding
sft <- pickSoftThreshold(vsd.mat,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

## graph soft power threshold
# setup dataframe
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

# graph scale independence
ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()
ggsave('../Figures/WGCNA/Scale independence soft threshold.png',
       height = 5,
       width = 5)

# graph mean connectivity
ggplot(sft_df, aes(x = Power, y = mean.k., label = Power)) +
  geom_text() +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Mean Connectivity") +
  ggtitle("Mean Connectivity") +
  # This adds some nicer aesthetics to our plot
  theme_classic()
ggsave('../Figures/WGCNA/Mean connectivity soft threshold.png',
       height = 5,
       width = 5)

### Run WGCNA
## use soft power threshold of 8

# bwnet <- blockwiseModules(vsd.mat,
#                           maxBlockSize = 10000, # depends on system memory
#                           TOMType = "signed", # topological overlap matrix
#                           power = 8, # soft threshold for network construction
#                           numericLabels = TRUE, # Let's use numbers instead of colors for module labels
#                           randomSeed = 1234, # there's some randomness associated with this calculation
#                           # so we should set a seed
#                           # saveTOMs = TRUE, # save TOM file
#                           # saveTOMFileBase = "BurtoniPOATOM-blockwise",
#                           verbose = 3)


## save results
# save(bwnet,
#      file = "bwnet.old.RData")
# #load data
# #load('bwnet.old.RData')


### Run WGCNA
## set min module size to 100
## use soft power threshold of 6

bwnet <- blockwiseModules(vsd.mat,
                          maxBlockSize = 15000, # depends on system memory
                          TOMType = "signed", # topological overlap matrix
                          power = 6, # soft threshold for network construction
                          numericLabels = FALSE, 
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          # so we should set a seed
                          # saveTOMs = TRUE, # save TOM file
                          # saveTOMFileBase = "BurtoniPOATOM-blockwise",
                          verbose = 3,
                          minModuleSize = 100,
                          nThreads = 12)


## save results
# save(bwnet,
#      file = "bwnet.RData")
# #load data
# load('bwnet.RData')


### graph network
## graph genes per module
library(ggrepel)
table(bwnet$colors) %>% 
  as.data.frame() %>% 
  dplyr::rename(module = Var1,
                gene.count = Freq) %>%
  ggplot(aes(x = reorder(module,
                         -gene.count),
             y = gene.count,
             label = gene.count)) +
  geom_point(size = 5) +
  geom_label_repel(point.padding = 0.25,
                   size = 10,
                   nudge_y = 10) +
  theme_classic(base_size = 30) +
  xlab('') +
  ylab('Number of genes') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggsave('../Figures/WGCNA/Gene count per module.png',
       height = 10,
       width = 10,
       units = 'in',
       dpi = 300)
# par('mar')
# par(mar = c(1, 1, 1, 1))
# Convert labels to colors for plotting
# mergedColors = labels2colors(bwnet$colors)
# Plot the dendrogram and the module colors underneath
png('../Figures/WGCNA/WGCNA dendrogram.png')
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  bwnet$colors,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 ,
  rowText = bwnet$colors,
  rowTextAlignment = 'center'
  )
dev.off()


# Plot the dendrogram and the module colors underneath
png('../Figures/WGCNA/WGCNA dendrogram paper.png',
    height = 2.6,
    width = 3.44,
    units = 'in',
    res = 300)
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  bwnet$colors,
  "",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = FALSE,
  guideHang = 0.05)
dev.off()

# #poster
# pdf('../Figures/WGCNA/WGCNA dendrogram.pdf',
#     height = 5,
#     width = 10)
# plotDendroAndColors(
#   bwnet$dendrograms[[1]],
#   bwnet$colors,
#   "Module colors",
#   dendroLabels = FALSE,
#   hang = 0.03,
#   addGuide = TRUE,
#   guideHang = 0.05 ,
#   rowText = bwnet$colors,
#   rowTextAlignment = 'center'
# )
# dev.off()



## write out module eigengenes
module_eigengenes <- bwnet$MEs

### limma results
# Create the design matrix from the `time_point` variable
des_mat <- model.matrix(~ sample.names.colData$Sex)
# Run linear model on each module. Limma wants our tests to be per row, so we also need to transpose so the eigengenes are rows

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)
# Apply multiple testing correction and obtain stats in a data frame.

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

## make module eigengene volcano plot
stats_df %>% 
  ggplot(aes(x = logFC,
             y=-log(P.Value),
             label = module)) +
  geom_label() +
  theme_bw()
ggsave('../Figures/WGCNA/DME volcano plot.png',
       height = 10,
       width = 10)

## create module specific box plot
module_eigengenes.df <- module_eigengenes %>%
  tibble::rownames_to_column("sample.id") %>%
  dplyr::inner_join(sample.names.colData %>%
                      tibble::rownames_to_column("sample.id"))

## run loop of all modules
for (i in colnames(module_eigengenes)) {
  module_eigengenes.df %>% 
  ggplot(aes(x = Sex,
      y = get(i))
  ) +
    geom_boxplot(width = 0.2, 
                 outlier.shape = NA,
                 lwd=2,
                 aes(color = Sex)) +
    ggforce::geom_sina(maxwidth = 0.3,
                       size = 5) + 
    theme_classic() +
    ylab(paste0(i)) + 
    theme_classic() + 
    theme(text = element_text(size = 30)) +
    scale_color_manual(values = c('black',
                                  'red'))
  ggsave(paste0('../Figures/WGCNA/boxplots/',
  i,
  ' boxplot.png'),
         height = 10,
         width = 10)
  
  # linked
  module_eigengenes.df %>% 
    ggplot(aes(x = Sex,
               y = get(i))
    ) +
    geom_boxplot(width = 0.2, 
                 outlier.shape = NA,
                 lwd=2,
                 aes(color = Sex)) +
    ggforce::geom_sina(maxwidth = 0.3,
                       size = 5) + 
    geom_line(aes(group = Tank.ID)) +
    theme_classic() +
    ylab(paste0(i)) + 
    theme_classic() + 
    theme(text = element_text(size = 30)) +
    scale_color_manual(values = c('black',
                                  'red'))
  ggsave(paste0('../Figures/WGCNA/boxplots/',
                i,
                ' boxplot linked.png'),
         height = 10,
         width = 10)
}

### module trait correlation
#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes = ncol(vsd.mat)
nSamples = nrow(vsd.mat)

#Recalculating MEs with label colors
MEs0 = moduleEigengenes(vsd.mat, bwnet$colors)$eigengenes
MEs = orderMEs(MEs0)

### correlation of MEs
ME.cor = cor(MEs)
ME.pvalue = corPvalueStudent(ME.cor, 
                             nrow(MEs))


#tank
# create empty boxes
ME.pvalue.empty = ME.pvalue
ME.pvalue.empty = signif(ME.pvalue.empty, 
                                      1)
ME.pvalue.empty[ME.pvalue.empty > .05] <- ''

# graph module trait correlation for presentation
pheatmap(ME.cor,
         cluster_rows = T,
         cluster_cols = T,
         scale = 'none',
         border_color = 'black',
         legend = T,
         treeheight_col = 25,
         treeheight_row = 25,
         angle_col = 315,
         fontsize = 25,
         # display_numbers = ME.pvalue.empty,
         color = blueWhiteRed(50),
         main = 'MEs',
         filename = "../Figures/WGCNA/MEs heatmap.png",
         width = 10,
         height = 10
)


## create pvclust 
# males
MEs.pvclust = pvclust(MEs, 
                                       method.hclust="average",
                                       method.dist="correlation", 
                                       nboot=1000, 
                                       parallel=T, 
                                       quiet=FALSE)

#graph
png(filename = "../Figures/WGCNA/MEs pvclust.png")
plot(MEs.pvclust)
pvrect(MEs.pvclust, alpha=0.8)
dev.off()
#dendrogram
png(filename = "../Figures/WGCNA/MEs pvclust dendro.png")
MEs.pvclust %>% 
  as.dendrogram() %>% 
  plot()
MEs.pvclust %>% 
  text
pvrect(MEs.pvclust, alpha=0.8)
dev.off()


### remove grey module
MEs = MEs %>% 
  dplyr::select(-c(MEgrey))

## save results
# save(MEs,
#      file = "MEs.RData")
# #load data
# load('MEs.RData')

### calculate KME
datKME = signedKME(vsd.mat, MEs, outputColumnName = "kME")

# get module names
module_df <- data.frame(
  gene_id = names(bwnet$colors),
  module = bwnet$colors
)

# combine lists
KMEs = datKME %>% 
  rownames_to_column('gene_id') %>% 
  full_join(module_df)

# 
# save(KMEs,
#      file = "KMEs.RData")
# #load data
# load('KMEs.RData')


#### Analysis of behavior and WGCNA sexes ####
### run on males and females seperately 

# combine all ME with behavior
data.behavior.comp.reduce.ME = MEs %>%  
  rownames_to_column('seq.id') %>% 
  separate(seq.id, 
           into = c("Position",
                    "Round", "Tank", "Sex")) %>% 
  mutate(Observation.id = paste0(Tank,
                                 "_",
                                 Round)) %>% 
  dplyr::select(-c(Position,
            Round,
            Tank,
            Sex)) %>% 
  full_join(data.behavior.comp.reduce)

# need matrix of traits that matches rows of column in order
## all
MEs.behavior.all = MEs %>%  
  rownames_to_column('seq.id') %>% 
  separate(seq.id, 
           into = c("Position",
                    "Round", "Tank", "Sex")) %>% 
  mutate(Observation.id = paste0(Tank,
                                 "_",
                                 Round)) %>% 
  dplyr::select(-c(Position,
            Round,
            Tank)) %>% 
  full_join(data.behavior.comp.reduce) %>% 
  mutate(Testosterone = ifelse(Sex == 'M',
                               Testosterone,
                               NA),
         Estradiol = ifelse(Sex == 'F',
                            Estradiol,
                               NA))
#all behavior  
behavior.all = MEs.behavior.all %>% 
  dplyr::select(c(Estradiol,
           Testosterone,
           Time.together,
           Courtship,
           Latency,
           GSI.male,
           GSI.female)) 
# module eigengenes
MEs.all = MEs.behavior.all %>% 
  dplyr::select(-c(Estradiol,
            Testosterone,
            Time.together,
            Courtship,
            Latency,
            GSI.male,
            GSI.female,
            Observation.id,
            Sex))
# all trait correlation
moduleTraitCor.all = cor(MEs.all, 
                          behavior.all, 
                          use = "p")
moduleTraitPvalue.all = corPvalueStudent(moduleTraitCor.all, 
                                          nrow(MEs.all))

# create empty boxes
moduleTraitPvalue.all.empty = moduleTraitPvalue.all
moduleTraitPvalue.all.empty = signif(moduleTraitPvalue.all.empty, 
                                               1)
# moduleTraitPvalue.all.empty[moduleTraitPvalue.all.empty <= .05] <- '*'
moduleTraitPvalue.all.empty[moduleTraitPvalue.all.empty > .1] <- ''

# #heat map of traits and modules
# pheatmap(moduleTraitCor.all,
#          cluster_rows = T,
#          cluster_cols = T,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 315,
#          fontsize = 25,
#          display_numbers = moduleTraitPvalue.all.empty,
#          color = blueWhiteRed(50),
#          breaks=seq(-0.6, 0.6, length.out=51),
#          main = 'Module trait relationship all',
#          filename = "../Figures/WGCNA/heatmaps/Module trait relationship all heatmap.png",
#          width = 10,
#          height = 10
# )


## males
MEs.males.beh = MEs %>%  
  rownames_to_column('seq.id') %>% 
  separate(seq.id, 
           into = c("Position",
                    "Round", "Tank", "Sex")) %>% 
  mutate(Observation.id = paste0(Tank,
                                 "_",
                                 Round)) %>% 
  filter(Sex == 'M') %>% 
  dplyr::select(-c(Position,
            Round,
            Tank,
            Sex)) %>% 
  full_join(data.behavior.comp.reduce)
#male behavior  
behavior.males = MEs.males.beh %>% 
  column_to_rownames('Observation.id') %>% 
  dplyr::select(c(Testosterone,
           Time.together,
           Courtship,
           Latency,
           GSI.male))
# module eigengenes
MEs.males = MEs.males.beh %>% 
  column_to_rownames('Observation.id') %>% 
  dplyr::select(-c(Estradiol,
            Testosterone,
            Time.together,
            Courtship,
            Latency,
            GSI.male,
            GSI.female))
# male trait correlation
moduleTraitCor.male = cor(MEs.males, 
                          behavior.males, 
                          use = "p",
                          method = 'spearman')
moduleTraitPvalue.male = corPvalueStudent(moduleTraitCor.male, 
                                          nrow(MEs.males))

## FDR correction
moduleTraitPvalue.male = moduleTraitPvalue.male %>% 
  as.matrix %>% 
  as.vector %>% 
  p.adjust(method='fdr') %>% 
  matrix(ncol=5)

# add names  
rownames(moduleTraitPvalue.male) = rownames(moduleTraitCor.male)
colnames(moduleTraitPvalue.male) = colnames(moduleTraitCor.male)

## females
MEs.females.beh = MEs %>%  
  rownames_to_column('seq.id') %>% 
  separate(seq.id, 
           into = c("Position",
                    "Round", "Tank", "Sex")) %>% 
  mutate(Observation.id = paste0(Tank,
                                 "_",
                                 Round)) %>% 
  filter(Sex == 'F') %>% 
  dplyr::select(-c(Position,
            Round,
            Tank,
            Sex)) %>% 
  full_join(data.behavior.comp.reduce)
#female behavior  
behavior.females = MEs.females.beh %>% 
  column_to_rownames('Observation.id') %>% 
  dplyr::select(c(Estradiol,
           Time.together,
           Latency,
           GSI.female)) 
# module eigengenes
MEs.females = MEs.females.beh %>% 
  column_to_rownames('Observation.id') %>% 
  dplyr::select(-c(Estradiol,
            Testosterone,
            Time.together,
            Courtship,
            Latency,
            GSI.male,
            GSI.female))
# female trait correlation
moduleTraitCor.female = cor(MEs.females, 
                          behavior.females, 
                          use = "p",
                          method = 'spearman')
moduleTraitPvalue.female = corPvalueStudent(moduleTraitCor.female, 
                                          nrow(MEs.females))

## FDR correction
moduleTraitPvalue.female = moduleTraitPvalue.female %>% 
  as.matrix %>% 
  as.vector %>% 
  p.adjust(method='fdr') %>% 
  matrix(ncol=4)

# add names  
rownames(moduleTraitPvalue.female) = rownames(moduleTraitCor.female)
colnames(moduleTraitPvalue.female) = colnames(moduleTraitCor.female)

## heatmaps for presentation
#male
textMatrix.male.line =  paste(signif(moduleTraitCor.male, 2), " (",
                              signif(moduleTraitPvalue.male, 1), ")", sep = "")
dim(textMatrix.male.line) = dim(moduleTraitCor.male)

# create text
#male
moduleTraitPvalue.male.empty = moduleTraitPvalue.male
moduleTraitPvalue.male.empty = signif(moduleTraitPvalue.male.empty, 
                                     2)
# moduleTraitPvalue.male.empty[moduleTraitPvalue.male.empty <= .05] <- '*'
moduleTraitPvalue.male.empty[moduleTraitPvalue.male.empty > .2] <- ''
#female
moduleTraitPvalue.female.empty = moduleTraitPvalue.female
moduleTraitPvalue.female.empty = signif(moduleTraitPvalue.female.empty, 
                                      2)
# moduleTraitPvalue.female.empty[moduleTraitPvalue.female.empty <= .05] <- '*'
moduleTraitPvalue.female.empty[moduleTraitPvalue.female.empty > .2] <- ''

# create star
#male
moduleTraitPvalue.male.empty.star = moduleTraitPvalue.male
moduleTraitPvalue.male.empty.star = signif(moduleTraitPvalue.male.empty.star, 
                                      2)
moduleTraitPvalue.male.empty.star[moduleTraitPvalue.male.empty.star <= .1] <- '+'
moduleTraitPvalue.male.empty.star[moduleTraitPvalue.male.empty.star > .1] <- ''
#female
moduleTraitPvalue.female.empty.star = moduleTraitPvalue.female
moduleTraitPvalue.female.empty.star = signif(moduleTraitPvalue.female.empty.star, 
                                        2)
moduleTraitPvalue.female.empty.star[moduleTraitPvalue.female.empty.star <= .1] <- '*'
moduleTraitPvalue.female.empty.star[moduleTraitPvalue.female.empty.star > .1] <- ''

##graph module trait correlation per sex for presentation
#male
# pheatmap(moduleTraitCor.male,
#          cluster_rows = F,
#          cluster_cols = F,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 315,
#          fontsize = 25,
#          fontsize_number = 20,
#          display_numbers = moduleTraitPvalue.male.empty,
#          color = blueWhiteRed(50),
#          main = 'Males',
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Male module trait heatmap presentation.png",
#          width = 10,
#          height = 10
# )
# #male
# pheatmap(moduleTraitCor.male,
#          cluster_rows = T,
#          cluster_cols = T,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 315,
#          fontsize = 25,
#          fontsize_number = 20,
#          display_numbers = moduleTraitPvalue.male.empty,
#          color = blueWhiteRed(50),
#          main = 'Males',
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Male module trait heatmap presentation cluster.png",
#          width = 10,
#          height = 10
# )
# 
# #female
# pheatmap(moduleTraitCor.female,
#          cluster_rows = F,
#          cluster_cols = F,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 315,
#          fontsize = 25,
#          fontsize_number = 20,
#          display_numbers = moduleTraitPvalue.female.empty,
#          color = blueWhiteRed(50),
#          main = 'Females',
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Female module trait heatmap presentation.png",
#          width = 10,
#          height = 10
# )
# #cluster
# pheatmap(moduleTraitCor.female,
#          cluster_rows = T,
#          cluster_cols = T,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 315,
#          fontsize = 25,
#          fontsize_number = 20,
#          display_numbers = moduleTraitPvalue.female.empty,
#          color = blueWhiteRed(50),
#          main = 'Females',
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Female module trait heatmap presentation cluster.png",
#          width = 10,
#          height = 10
# )
# 
## poster
#male
pheatmap(moduleTraitCor.male,
         cluster_rows = T,
         cluster_cols = T,
         scale = 'none',
         border_color = 'black',
         legend = T,
         treeheight_col = 25,
         treeheight_row = 25,
         angle_col = 315,
         fontsize = 15,
         fontsize_number = 20,
         display_numbers = moduleTraitPvalue.male.empty.star,
         color = blueWhiteRed(50),
         main = 'Males',
         breaks=seq(-0.6, 0.6, length.out=51),
         filename = "../Figures/WGCNA/heatmaps/Male module trait heatmap poster cluster fdr 0.1.pdf",
         width = 5.25,
         height = 5.25
)

#female
pheatmap(moduleTraitCor.female,
         cluster_rows = T,
         cluster_cols = T,
         scale = 'none',
         border_color = 'black',
         legend = T,
         treeheight_col = 25,
         treeheight_row = 25,
         angle_col = 315,
         fontsize = 15,
         fontsize_number = 20,
         display_numbers = moduleTraitPvalue.female.empty.star,
         color = blueWhiteRed(50),
         main = 'Females',
         breaks=seq(-0.6, 0.6, length.out=51),
         filename = "../Figures/WGCNA/heatmaps/Female module trait heatmap poster cluster fdr 0.1.pdf",
         width = 5.25,
         height = 5.25
)


## poster
#male
pheatmap(moduleTraitCor.male,
         cluster_rows = T,
         cluster_cols = T,
         scale = 'none',
         border_color = 'black',
         legend = T,
         treeheight_col = 25,
         treeheight_row = 25,
         angle_col = 315,
         fontsize = 15,
         fontsize_number = 20,
         color = blueWhiteRed(50),
         main = 'Males',
         breaks=seq(-0.6, 0.6, length.out=51),
         filename = "../Figures/WGCNA/heatmaps/Male module trait heatmap paper cluster.png",
         width = 7,
         height = 7
)

#female
pheatmap(moduleTraitCor.female,
         cluster_rows = T,
         cluster_cols = T,
         scale = 'none',
         border_color = 'black',
         legend = T,
         treeheight_col = 25,
         treeheight_row = 25,
         angle_col = 315,
         fontsize = 15,
         fontsize_number = 20,
         color = blueWhiteRed(50),
         main = 'Females',
         breaks=seq(-0.6, 0.6, length.out=51),
         filename = "../Figures/WGCNA/heatmaps/Female module trait heatmap paper cluster.png",
         width = 7,
         height = 7
)

## graph combined 
#make dataframe
moduleTraitCor.combined = moduleTraitCor.female %>% 
  as.data.frame() %>% 
  rownames_to_column('MEs') %>% 
  full_join(moduleTraitCor.male %>% 
              as.data.frame() %>% 
              rownames_to_column('MEs'),
            by = 'MEs',
            suffix = c(".F",
                       ".M")) %>% 
  column_to_rownames('MEs') %>% 
  as.matrix()
#graph
colnames(moduleTraitCor.combined)
group_df = data.frame(Sex=as.factor(rep(c("F", "M"), c(4,5))))

rownames(group_df)<-colnames(moduleTraitCor.combined)

ann_colors = list(
  Sex = c(F="black", M="red"))



pheatmap(moduleTraitCor.combined,
         annotation_col = group_df,
         annotation_colors = ann_colors,
         cluster_rows = T,
         cluster_cols = T,
         scale = 'none',
         border_color = 'black',
         legend = T,
         treeheight_col = 25,
         treeheight_row = 25,
         angle_col = 315,
         fontsize = 15,
         fontsize_number = 20,
         color = blueWhiteRed(50),
         cutree_rows = 3,
         cutree_cols = 4,
         # main = 'Females',
         breaks=seq(-0.6, 0.6, length.out=51),
         filename = "../Figures/WGCNA/heatmaps/Combined module trait heatmap paper cluster.png",
         width = 12,
         height = 5.5
)


#### Analysis of behavior and WGCNA pairs ####
## pairs
# pairs trait correlation
moduleTraitCor.tank = cor(MEs.females %>% 
                            rename_all(~ paste(., "F", sep = ".")), 
                          MEs.males %>% 
                            rename_all(~ paste(., "M", sep = ".")), 
                          use = "p",
                          method = 'spearman')
moduleTraitPvalue.tank = corPvalueStudent(moduleTraitCor.tank, 
                                          nrow(MEs.females))

# ## FDR correction
# moduleTraitPvalue.tank = moduleTraitPvalue.tank %>% 
#   as.matrix %>% 
#   as.vector %>% 
#   p.adjust(method='fdr') %>% 
#   matrix(ncol=10)

# add names  
rownames(moduleTraitPvalue.tank) = rownames(moduleTraitCor.tank)
colnames(moduleTraitPvalue.tank) = colnames(moduleTraitCor.tank)


#tank
# textMatrix.tank.line =  paste(signif(moduleTraitCor.tank, 2), " (",
#                               signif(moduleTraitPvalue.tank, 1), ")", sep = "")
# dim(textMatrix.tank.line) = dim(moduleTraitCor.tank)

#tank
# create empty boxes
moduleTraitPvalue.tank.empty = moduleTraitPvalue.tank
moduleTraitPvalue.tank.empty = signif(moduleTraitPvalue.tank.empty, 
                                      1)
# moduleTraitPvalue.tank.empty[moduleTraitPvalue.tank.empty <= .05] <- '*'
moduleTraitPvalue.tank.empty[moduleTraitPvalue.tank.empty > .1] <- ''

# create significant cutoff
moduleTraitPvalue.tank.sig = moduleTraitPvalue.tank
moduleTraitPvalue.tank.sig = round(moduleTraitPvalue.tank.sig, 
                                      3)
moduleTraitPvalue.tank.sig[moduleTraitPvalue.tank.sig > .1] <- ''


# graph module trait correlation for paper
moduleTraitCor.tank |> 
  as.data.frame() |> 
  rownames_to_column('Female.ME') |> 
  pivot_longer(cols = -Female.ME,
               names_to = 'Male.ME',
               values_to = 'Correlation') |>
  left_join(moduleTraitPvalue.tank |> 
              as.data.frame() |> 
              rownames_to_column('Female.ME') |> 
              pivot_longer(cols = -Female.ME,
                           names_to = 'Male.ME',
                           values_to = 'p.value')) |> 
  separate_wider_delim(cols = 'Female.ME',
                       delim = '.',
                       names = c('Female.ME',
                                 NA)) |> 
  separate_wider_delim(cols = 'Male.ME',
                       delim = '.',
                       names = c('Male.ME',
                                 NA)) |> 
  separate_wider_delim(cols = 'Female.ME',
                       delim = 'ME',
                       names = c(NA,
                                 'Female.ME')) |> 
  separate_wider_delim(cols = 'Male.ME',
                       delim = 'ME',
                       names = c(NA,
                                 'Male.ME')) |> 
  mutate(same = ifelse(Female.ME == Male.ME,
                       1,
                       0),
         Female.ME.fact = factor(Female.ME),
         Male.ME.fact = factor(Male.ME, 
                               levels = levels(Female.ME.fact)),
         Female.ME.fact = factor(Female.ME.fact, 
                               levels = rev(levels(Male.ME.fact)))) |>
  mutate(p.value.symbol = case_when(p.value < 0.1 & p.value > 0.05 ~ '+',
                                    p.value < 0.05 & p.value > 0.01 ~ '*',
                                    p.value < 0.01 ~ '**',
                                    TRUE ~ '')) |> 
  ggplot(aes(x = Male.ME.fact,
             y = Female.ME.fact,
             fill = Correlation,
             label = p.value.symbol)) +
  geom_tile(color = 'grey',
            linewidth = 0.5) +
  geom_tile(data = \(.x) subset(.x, same == 1), 
            color = "black", 
            linewidth = 1) + 
  theme_minimal(base_size = 15) +
  geom_text(size = 8) +
  ggtitle('Module Eigengene Correlation') +
  xlab('Male') +
  ylab('Female') +
  scale_fill_gradient2(midpoint = 0,
                       low = '#0D8CFF',
                       high = '#FF3300',
                       mid = 'white') +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  coord_equal()
ggsave(filename = "../Figures/WGCNA/heatmaps/Pairs module trait heatmap paper pvalue.pdf",
       width = 7,
       height = 7)
# 
# # graph module trait correlation for presentation
# pheatmap(moduleTraitCor.tank, 
#          cluster_rows = F,
#          cluster_cols = F,
#          scale = 'none',
#          border_color = 'black',
#          legend = F,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 315,
#          fontsize_number = 20,
#          fontsize = 25,
#          display_numbers = moduleTraitPvalue.tank.empty,
#          color = blueWhiteRed(50),
#          main = 'Pairs',
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Pairs module trait heatmap presentation pvalue fdr.png",
#          width = 7,
#          height = 7
# )
# 
# 
# #cluster
# pheatmap(moduleTraitCor.tank, 
#          cluster_rows = T,
#          cluster_cols = T,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 315,
#          fontsize_number = 20,
#          fontsize = 25,
#          display_numbers = moduleTraitPvalue.tank.empty,
#          color = blueWhiteRed(50),
#          main = 'Pairs',
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Pairs module trait heatmap presentation pvalue cluster.png",
#          width = 10,
#          height = 10
# )
# 
# 
# #order by max value
# moduleTraitCor.tank.order.df = moduleTraitCor.tank %>%
#     as.data.frame() %>%
#     summarise(across(everything(),
#                      ~ max(.x)))  %>%
#     pivot_longer(cols = everything(),
#                  names_to = 'male',
#                  values_to = 'm.value') %>%
#     cbind(moduleTraitCor.tank %>%
#     t() %>%
#     as.data.frame() %>%
#     summarise(across(everything(),
#                      ~ max(.x)))  %>%
#     pivot_longer(cols = everything(),
#                  names_to = 'female',
#                  values_to = 'f.value')) %>%
#     arrange(-m.value,
#             -f.value)
# 
# #reorder matrix
# moduleTraitCor.tank.order = moduleTraitCor.tank %>%
#   as.data.frame() %>%
#   select(moduleTraitCor.tank.order.df$male) %>%
#   t() %>%
#   as.data.frame() %>%
#   select(moduleTraitCor.tank.order.df$female) %>%
#   as.matrix() %>%
#   t()
# 
# # #get correlation
# moduleTraitCor.tank = cor(MEs.females %>%
#                             rename_all(~ paste(., "F", sep = ".")),
#                           MEs.males %>%
#                             rename_all(~ paste(., "M", sep = ".")),
#                           use = "p")
# #reorder
# moduleTraitCor.tank = moduleTraitCor.tank[,moduleTraitCor.tank.order.df$male]
# 
# moduleTraitCor.tank = moduleTraitCor.tank[moduleTraitCor.tank.order.df$female,]
# 
# 
# 
# #get pvalue
# moduleTraitPvalue.tank = corPvalueStudent(moduleTraitCor.tank,
#                                           nrow(MEs.females))
# 
# #tank
# #create empty boxes
# moduleTraitPvalue.tank.empty = moduleTraitPvalue.tank
# moduleTraitPvalue.tank.empty = signif(moduleTraitPvalue.tank.empty,
#                                       1)
# moduleTraitPvalue.tank.empty[moduleTraitPvalue.tank.empty > .1] <- ''
# 
# #graph module trait correlation for presentation
# pheatmap(moduleTraitCor.tank.order,
#          cluster_rows = F,
#          cluster_cols = F,
#          scale = 'none',
#          border_color = 'black',
#          cutree_rows = 4,
#          cutree_cols = 5,
#          legend = T,
#          angle_col = 315,
#          display_numbers = moduleTraitPvalue.tank.empty,
#          color = blueWhiteRed(50),
#          main = 'Pairs',
#          filename = "../Figures/WGCNA/heatmaps//Tanks module trait heatmap presentation pvalue order.png",
#          width = 10,
#          height = 10
# )
# dev.off()
# 
# 
# ## FDR correction
# moduleTraitCor.tank = cor(MEs.females %>% 
#                             rename_all(~ paste(., "F", sep = ".")), 
#                           MEs.males %>% 
#                             rename_all(~ paste(., "M", sep = ".")), 
#                           use = "p",
#                           method = 'spearman')
# moduleTraitPvalue.tank = corPvalueStudent(moduleTraitCor.tank, 
#                                           nrow(MEs.females))
# ## FDR correction
# moduleTraitPvalue.tank = moduleTraitPvalue.tank %>% 
#   as.matrix %>% 
#   as.vector %>% 
#   p.adjust(method='fdr') %>% 
#   matrix(ncol=10)
# 
# # add names  
# rownames(moduleTraitPvalue.tank) = rownames(moduleTraitPvalue.tank)
# colnames(moduleTraitPvalue.tank) = colnames(moduleTraitPvalue.tank)
# 
# #tank
# textMatrix.tank.line =  paste(signif(moduleTraitPvalue.tank, 2), " (",
#                                        signif(moduleTraitPvalue.tank, 1), ")", sep = "")
# dim(textMatrix.tank.line) = dim(moduleTraitCor.tank)
# 
# #tank
# #create empty boxes
# moduleTraitPvalue.tank.empty = moduleTraitPvalue.tank
# moduleTraitPvalue.tank.empty = signif(moduleTraitPvalue.tank.empty, 
#                                                2)
# # moduleTraitPvalue.tank.empty[moduleTraitPvalue.tank.empty <= .05] <- '*'
# moduleTraitPvalue.tank.empty[moduleTraitPvalue.tank.empty > .2] <- ''
# 
# #create star 
# moduleTraitPvalue.tank.empty.star = moduleTraitPvalue.tank
# moduleTraitPvalue.tank.empty.star = signif(moduleTraitPvalue.tank.empty.star, 
#                                                     2)
# moduleTraitPvalue.tank.empty.star[moduleTraitPvalue.tank.empty.star <= .1] <- '*'
# moduleTraitPvalue.tank.empty.star[moduleTraitPvalue.tank.empty.star > .1] <- ''
# 
# 
# # poster
# pheatmap(moduleTraitCor.tank, 
#          clustering_distance_rows = 'euclidean',
#          clustering_distance_cols = 'euclidean',
#          clustering_method = 'ward.D2',
#          cluster_rows = T,
#          cluster_cols = T,
#          scale = 'none',
#          border_color = 'black',
#          cutree_rows = 3,
#          cutree_cols = 4,
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 45,
#          display_numbers = moduleTraitPvalue.tank.empty.star,
#          color = blueWhiteRed(50),
#          main = 'Comparison between partners',
#          fontsize = 15,
#          fontsize_number = 20,
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Pairs module trait heatmap poster pvalue cluster fdr 0.1.png",
#          width = 10,
#          height = 10
# )
# 
# ### add behavior data to pairs 
# ## add ME to behavior
# #male behavior  
# behavior.males.me = behavior.males %>% 
#   rownames_to_column('ID') %>% 
#   full_join(MEs.males %>% 
#               rename_all(~ paste(., "M", sep = ".")) %>% 
#               rownames_to_column('ID')) %>% 
#   column_to_rownames('ID')
# # remove female behavior
# behavior.males.me = behavior.males.me %>% 
#   select(-c(Time.together,
#            Latency))
# 
# #female behavior  
# behavior.females.me = behavior.females %>% 
#   rownames_to_column('ID') %>% 
#   full_join(MEs.females %>% 
#               rename_all(~ paste(., "F", sep = ".")) %>% 
#               rownames_to_column('ID')) %>% 
#   column_to_rownames('ID')
# 
# # tank trait correlation
# moduleTraitCor.tank.behavior = cor(behavior.females.me, 
#                           behavior.males.me, 
#                           use = "p",
#                           method = 'spearman')
# moduleTraitPvalue.tank.behavior = corPvalueStudent(moduleTraitCor.tank.behavior, 
#                                           nrow(behavior.females.me))
# 
# ## FDR correction
# moduleTraitPvalue.tank.behavior = moduleTraitPvalue.tank.behavior %>% 
#   as.matrix %>% 
#   as.vector %>% 
#   p.adjust(method='fdr') %>% 
#   matrix(ncol=13)
# 
# # add names  
# rownames(moduleTraitPvalue.tank.behavior) = rownames(moduleTraitCor.tank.behavior)
# colnames(moduleTraitPvalue.tank.behavior) = colnames(moduleTraitCor.tank.behavior)
# 
# #tank
# textMatrix.tank.behavior.line =  paste(signif(moduleTraitCor.tank.behavior, 2), " (",
#                               signif(moduleTraitCor.tank.behavior, 1), ")", sep = "")
# dim(textMatrix.tank.behavior.line) = dim(moduleTraitCor.tank.behavior)
# 
# #tank
# #create empty boxes
# moduleTraitPvalue.tank.behavior.empty = moduleTraitPvalue.tank.behavior
# moduleTraitPvalue.tank.behavior.empty = signif(moduleTraitPvalue.tank.behavior.empty, 
#                                       2)
# # moduleTraitPvalue.tank.behavior.empty[moduleTraitPvalue.tank.behavior.empty <= .05] <- '*'
# moduleTraitPvalue.tank.behavior.empty[moduleTraitPvalue.tank.behavior.empty > .2] <- ''
# 
# #create star 
# moduleTraitPvalue.tank.behavior.empty.star = moduleTraitPvalue.tank.behavior
# moduleTraitPvalue.tank.behavior.empty.star = signif(moduleTraitPvalue.tank.behavior.empty.star, 
#                                                2)
# moduleTraitPvalue.tank.behavior.empty.star[moduleTraitPvalue.tank.behavior.empty.star <= .1] <- '*'
# moduleTraitPvalue.tank.behavior.empty.star[moduleTraitPvalue.tank.behavior.empty.star > .1] <- ''
# 
# 
# 
# #graph module trait correlation for presentation
# pheatmap(moduleTraitCor.tank.behavior, 
#          cluster_rows = F,
#          cluster_cols = F,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 45,
#          display_numbers = moduleTraitPvalue.tank.behavior.empty,
#          color = blueWhiteRed(50),
#          main = 'Pairs',
#          fontsize = 25,
#          fontsize_number = 20,
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Pairs module trait and behavior heatmap presentation pvalue.png",
#          width = 10,
#          height = 10
# )
# 
# 
# #cluster
# pheatmap(moduleTraitCor.tank.behavior, 
#          cluster_rows = T,
#          cluster_cols = T,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 45,
#          display_numbers = moduleTraitPvalue.tank.behavior.empty,
#          color = blueWhiteRed(50),
#          main = 'Pairs',
#          fontsize = 25,
#          fontsize_number = 20,
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Pairs module trait and behavior heatmap presentation pvalue cluster.png",
#          width = 10,
#          height = 10
# )
# 
# 
# # poster
# pheatmap(moduleTraitCor.tank.behavior, 
#          clustering_distance_rows = 'euclidean',
#          clustering_distance_cols = 'euclidean',
#          clustering_method = 'ward.D2',
#          cluster_rows = T,
#          cluster_cols = T,
#          scale = 'none',
#          border_color = 'black',
#          cutree_rows = 3,
#          cutree_cols = 4,
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 45,
#          display_numbers = moduleTraitPvalue.tank.behavior.empty.star,
#          color = blueWhiteRed(50),
#          main = 'Comparison between partners',
#          fontsize = 15,
#          fontsize_number = 20,
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Pairs module trait and behavior heatmap poster pvalue cluster fdr 0.1.png",
#          width = 10,
#          height = 10
# )
# 
# ### reorder columns and rows with pvclust
# ## create pvclust
# # males
# male.module.behavior.pvclust = pvclust(moduleTraitCor.tank.behavior,
#               method.hclust="ward.D2",
#         method.dist="euclidean",
#         use.cor="pairwise.complete.obs",
#         nboot=1000,
#         parallel=T,
#         quiet=FALSE)
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
#                                        method.hclust="ward.D2",
#                                        method.dist="euclidean",
#                                        use.cor="pairwise.complete.obs",
#                                        nboot=1000,
#                                        parallel=T,
#                                        quiet=FALSE)
# 
# #graph
# png(filename = "../Figures/WGCNA/heatmaps/Female trait and behavior pvclust.png")
# plot(female.module.behavior.pvclust)
# pvrect(female.module.behavior.pvclust, alpha=0.8)
# dev.off()
# 
# 
# ### add ALL behavior data to pairs 
# ## add ME to behavior
# #male behavior  
# MEs.males.beh = MEs.males.beh %>% 
#   rename_at(vars(starts_with('ME')), ~ paste(., "M", sep = "."))  %>% 
#   column_to_rownames('Observation.id')
# 
# #female behavior  
# MEs.females.beh = MEs.females.beh %>% 
#   rename_at(vars(starts_with('ME')), ~ paste(., "F", sep = "."))  %>% 
#   column_to_rownames('Observation.id')
# 
# # tank trait correlation
# moduleTraitCor.tank.behavior.all = cor(MEs.females.beh, 
#                                    MEs.males.beh, 
#                                    use = "p")
# moduleTraitPvalue.tank.behavior.all = corPvalueStudent(moduleTraitCor.tank.behavior.all, 
#                                                    nrow(behavior.females.me))
# 
# #tank
# #create empty boxes
# moduleTraitPvalue.tank.behavior.all.empty = moduleTraitPvalue.tank.behavior.all
# moduleTraitPvalue.tank.behavior.all.empty = signif(moduleTraitPvalue.tank.behavior.all.empty, 
#                                                1)
# # moduleTraitPvalue.tank.behavior.all.empty[moduleTraitPvalue.tank.behavior.all.empty <= .05] <- '*'
# moduleTraitPvalue.tank.behavior.all.empty[moduleTraitPvalue.tank.behavior.all.empty > .1] <- ''
# 
# #graph module trait correlation for presentation
# pheatmap(moduleTraitCor.tank.behavior.all, 
#          cluster_rows = F,
#          cluster_cols = F,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 45,
#          display_numbers = moduleTraitPvalue.tank.behavior.all.empty,
#          color = blueWhiteRed(50),
#          main = 'Pairs',
#          fontsize = 25,
#          fontsize_number = 20,
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Pairs module trait and behavior.all heatmap presentation pvalue.png",
#          width = 10,
#          height = 10
# )
# 
# 
# #cluster
# pheatmap(moduleTraitCor.tank.behavior.all, 
#          cluster_rows = T,
#          cluster_cols = T,
#          scale = 'none',
#          border_color = 'black',
#          legend = T,
#          treeheight_col = 25,
#          treeheight_row = 25,
#          angle_col = 45,
#          display_numbers = moduleTraitPvalue.tank.behavior.all.empty,
#          color = blueWhiteRed(50),
#          main = 'Pairs',
#          fontsize = 25,
#          fontsize_number = 20,
#          breaks=seq(-0.6, 0.6, length.out=51),
#          filename = "../Figures/WGCNA/heatmaps/Pairs module trait and behavior.all heatmap presentation pvalue cluster.png",
#          width = 10,
#          height = 10
# )
# 
# 
# 
# ### males vs females forest plot
# MEs.behavior.all %>% 
#   pivot_longer(cols = starts_with("ME"),
#                names_to = 'modules',
#                values_to = "ME") %>% 
#   ggplot(aes(y = modules,
#              x = ME,
#              group = Sex,
#              color = Sex,
#              fill = Sex)) +
#   geom_vline(xintercept = 0,
#              linetype="dotted") +
# stat_summary(fun.data = "mean_cl_normal",
#              geom = "errorbar",
#              color="black",
#              width=0.5,
#              position=position_dodge(width = .5))+
#   # stat_summary(fun = mean,
#   #              geom = "errorbar",
#   #              fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
#   #              fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),
#   #              color="black",
#   #              width=0.5,
#   #              position=position_dodge(width = .5)) +
#   stat_summary(fun=mean,
#                geom="point", 
#                color="black",
#                position=position_dodge(width = .5),
#                shape = 15) +
#   geom_point(position=position_dodge(width = .5)) +
#   theme_classic()
# ggsave("../Figures/WGCNA/boxplots/All modules sex comparison.png",
#        width = 10,
#        height = 10)

#presentation
MEs.behavior.all %>% 
  pivot_longer(cols = starts_with("ME"),
               names_to = 'modules',
               values_to = "ME") %>% 
  ggplot(aes(y = modules,
             x = ME,
             group = Sex,
             color = Sex,
             fill = Sex)) +
  geom_vline(xintercept = 0,
             linetype="dotted") +
  stat_summary(fun.data = "mean_cl_normal",
               geom = "errorbar",
               color="black",
               width=1,
               position=position_dodge(width = .5),
               size = 1)+
  # stat_summary(fun = mean,
  #              geom = "errorbar",
  #              fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
  #              fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),
  #              color="black",
  #              width=0.5,
  #              position=position_dodge(width = .5)) +
  stat_summary(fun=mean,
               geom="point", 
               color="black",
               position=position_dodge(width = .5),
               shape = 15,
               size = 4) +
  geom_point(position=position_dodge(width = .5),
             size = 2) +
  theme_classic(base_size = 25) +
  scale_color_manual(values = c('Other' = 'grey',
                                'M' = '#E69F00',
                                'F' = '#008080')) +
  theme(legend.position = "none") +
  ylab('')
ggsave("../Figures/WGCNA/boxplots/All modules sex comparison paper.pdf",
       width = 6,
       height = 9.5)


  

#### WGCNA correlation matrix ####
library(Hmisc)
### create correlation matrix 
ME.Corr = rcorr(MEs %>% 
                     t())


## convert correlation matrix to data frame 
# remove individual
# add sample info
# add same tank
ME.Corr.df = data.frame(from=colnames(ME.Corr$r)[col(ME.Corr$r)], 
                           to=rownames(ME.Corr$r)[row(ME.Corr$r)], 
                           corr=c(ME.Corr$r)) %>% 
  mutate(from = substring(from, 
                          4),
         to = substring(to, 
                        4)) %>% 
  filter(corr != 1) %>% 
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.')) %>% 
              select(c(Round,
                       Tank,
                       Sex)) %>% 
              mutate(from = paste(Round,
                                  Tank,
                                  Sex,
                                  sep = '.'),
                     from.tank = paste(Round,
                                       Tank,
                                       sep = '.'),
                     from.sex = Sex) %>% 
              select(-c(Round,
                        Tank,
                        Sex))) %>% 
  full_join(sample.names %>% 
              mutate(Sample.Name = str_replace_all(Sample.Name,
                                                   '-',
                                                   '.')) %>% 
              select(c(Round,
                       Tank,
                       Sex)) %>% 
              mutate(to.tank = paste(Round,
                                     Tank,
                                     sep = '.'),
                     to = paste(to.tank,
                                Sex,
                                sep = '.'),
                     to.sex = Sex) %>% 
              select(-c(Round,
                        Tank,
                        Sex))) %>% 
  mutate(same = ifelse(from.tank == to.tank,
                       'same',
                       'other'))

### create distance statistic
ME.Corr.df.sum = ME.Corr.df %>% 
  group_by(from,
           from.sex,
           to.sex) %>% 
  mutate(Mean = mean(corr),
         Sd = sd(corr)) %>% 
  ungroup() %>% 
  mutate(Zscore = (corr - Mean)/Sd)

### compare to PC data from hormones and behavior
ME.Corr.df.sum.tank = ME.Corr.df.sum %>% 
  filter(same == 'same') %>% 
  full_join(Tank.PCA.Behavior.hormone.data %>% 
              select(c(tank,
                       PC1)) %>% 
              dplyr::rename(from.tank = tank))
#### graph correlation matrix WGCNA ####
## graph same vs. other
#males to females
ME.Corr.df %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = corr,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from males to all females')
ggsave('../Figures/WGCNA/Correlation.score/Corr from males to all females.png')

#females to males
ME.Corr.df %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = corr,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from females to all males')
ggsave('../Figures/WGCNA/Correlation.score/Corr from females to all males.png')

### correlation stats 
### z score
## graph same vs. other
#males to females
ME.Corr.df.sum %>% 
  filter(from.sex == 'M') %>% 
  filter(to.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from males to all females zscore')
ggsave('../Figures/WGCNA/Correlation.score/Corr from males to all females zscore.png')

#females to males
ME.Corr.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = from,
             color = same)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle('Corr from females to all males zscore')
ggsave('../Figures/WGCNA/Correlation.score/Corr from females to all males zscore.png')


#males to females
ME.Corr.df.sum %>% 
  filter(to.sex == 'F') %>% 
  filter(from.sex == 'M') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Corr from males to all females zscore') 
ggsave('../Figures/WGCNA/Correlation.score/Corr from females to all males zscore all.png')

#females to males
ME.Corr.df.sum %>% 
  filter(to.sex == 'M') %>% 
  filter(from.sex == 'F') %>% 
  ggplot(aes(y = Zscore,
             x = same,
             color = same)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_classic() +
  ggtitle('Corr from females to all males zscore') 
ggsave('../Figures/WGCNA/Correlation.score/Corr from females to all males zscore all.png')

#females to males
#presentation
ME.Corr.df.sum %>% 
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
  ggtitle('Correlation score females') + 
  theme(text = element_text(size = 30),
        legend.position = 'none') +
  xlab('Pair') +
  ylab('Zscore female correlation')
ggsave('../Figures/WGCNA/Correlation.score/Corr from females to all males zscore all presentation.png',
       height = 10,
       width = 10)


### compare correlation zscore with PC1 of behavior and hormones
ME.Corr.df.sum.tank %>% 
  ggplot(aes(x = PC1,
             y = Zscore,
             color = from.sex)) +
  geom_point() +
  geom_text(aes(label = from.tank),
            nudge_y = 0.1) +
  geom_smooth(method = 'lm') +
  theme_classic() +
  ggtitle('Corr zscore vs PC1 of behavior and hormones') +
  facet_grid(~from.sex)
ggsave('../Figures/WGCNA/Correlation.score/Corr zscore vs PC1 of behavior and hormones.png')

#### correlation synchronization ME ####
### add behavior to ME correlation dataframe
ME.Corr.df.sum.sync = ME.Corr.df.sum %>% 
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
            corr)) %>% 
  pivot_wider(names_from = from.sex,
              values_from = Zscore,
              names_prefix = 'Zscore.') %>% 
  full_join(data.behavior.comp)

## compare zscores
ME.Corr.df.sum.sync %>% 
  ggplot(aes(x = Zscore.M,
             y= Zscore.F)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  ggtitle('ME correlation zscore comparison')
ggsave('../Figures/WGCNA/Correlation.score/ME correlation zscore comparison sexes.png',
       height = 10,
       width = 10)

# check correlation
library(PerformanceAnalytics)
pdf('../Figures/WGCNA/Correlation.score/ME corr zscore and behavior chart correlation.pdf',
    height = 10,
    width = 10)
chart.Correlation(ME.Corr.df.sum.sync %>% 
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

# # lm
# Zscore.F.court.lm = lm(Zscore.F ~ Court.total,
#                        data=sampleCorr.df.sum.sync)
# 
# summary(Zscore.F.court.lm)
# anova(Zscore.F.court.lm)
# 
# 
# #graph glm with outlier on graph
# # presentation
# library(jtools)
# effect_plot(Zscore.F.court.lm, 
#             pred = Court.total, 
#             interval = TRUE,
#             rug=TRUE) +
#   geom_point(data=sampleCorr.df.sum.sync ,
#              aes(x=Court.total,
#                  y=Zscore.F),
#              size = 5)+
#   theme_bw()+
#   theme(legend.position = 'null')+
#   ylab('Zscore female correlation') +
#   xlab('Male courtship (counts)')+ 
#   theme(text = element_text(size = 30))
# ggsave('../Figures/DESEQ2/Female corr zscore vs courtship with outlier presentation.png',
#        height = 10,
#        width = 10)



#### Create GO database ####
# making GO library
library(AnnotationForge)

## make library from NCBI
makeOrgPackageFromNCBI(version="0.1",
                       maintainer="Isaac Miller-Crews <imillerc@iu.edu>",
                       author="Isaac Miller-Crews <imillerc@iu.edu>",
                       outputDir = "./",
                       tax_id = "8128",
                       genus = "Oreochromis",
                       species = "niloticus")


# then you can call install.packages based on the return value
install.packages("/./org.Oniloticus.eg.db", 
                 repos=NULL)

#### Create GO figures per module ####
# load libraries
library(enrichplot)
library(clusterProfiler)
library(org.Oniloticus.eg.db) # need to make database in code above
library(ggnewscale)

# load all genes for GO comparison
all.genes = KMEs$gene_id

# get dataframe of genes and modules
module.genes = KMEs %>% 
  dplyr::select(gene_id,
                module)

# count number of genes
module.genes %>%
  dplyr::count(module)
# module    n
# 1      black  299
# 2       blue  921
# 3      brown  722
# 4      green  333
# 5       grey 7461
# 6    magenta  235
# 7       pink  287
# 8     purple  225
# 9        red  307
# 10 turquoise 2169
# 11    yellow  384

# convert ensembl gene name to gene name using enrichgo
tmp.enrichgo.ensembl <- enrichGO(gene = all.genes, 
                                 universe = all.genes,
                                 OrgDb = "org.Oniloticus.eg.db",
                                 keyType = "ENSEMBL",
                                 ont = "BP",
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = 'fdr',
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 readable = T)
# get gene names
tmp.enrichgo.ensembl.genes = tmp.enrichgo.ensembl@gene2Symbol %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  na.omit() %>% 
  dplyr::rename(ENSEMBL = V1) %>% 
  rownames_to_column('gene_id') %>% 
  separate_wider_delim(gene_id,
                       delim = '.',
                       names = c('gene_id',
                                 NA),
                       too_few = 'align_start') #some ID have duplicate genes

# combine with ensembl genes
module.genes.ensembl = module.genes %>% 
  left_join(tmp.enrichgo.ensembl.genes) %>% 
  mutate(gene_id = ifelse(is.na(ENSEMBL),
                            gene_id,
                            ENSEMBL)) %>% 
  dplyr::select(-c(ENSEMBL)) %>% 
  distinct()

# add ensembl genes to all genes
all.genes.ensembl = c(all.genes,
                      tmp.enrichgo.ensembl.genes$ENSEMBL) %>% 
  unique()


## get list of module colors
module.colors = module.genes %>% 
  filter(module != 'grey') %>% 
  pull(module) %>% 
  unique()

### loop through modules for enriched genes
# create empty data frame to save results
enrichGO.results.modules = data.frame()
enrichGO.results.modules.simplify = data.frame()

for (i in module.colors) {
  
  # run enrichgo  
  # gene symbol
  tmp.enrichgo <- enrichGO(gene = module.genes.ensembl %>% 
                             filter(module == i) %>% 
                             pull(gene_id) %>% 
                             unique(), 
                           universe = all.genes.ensembl,
                           OrgDb = "org.Oniloticus.eg.db",
                           keyType = "SYMBOL",
                           ont = "BP",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = 'fdr',
                           minGSSize = 10,
                           maxGSSize = 500,
                           readable = T)
  
  # save results
  enrichGO.results.modules = enrichGO.results.modules %>%
    rbind(tmp.enrichgo@result %>%
            mutate(module = i,
                   ont = 'BP'))
  
  write.csv(enrichGO.results.modules,
            '../Figures/WGCNA/GO/enrichGO.results.modules.csv',
            row.names = F)
  
  # check number of sig GO
  tmp.num = tmp.enrichgo@result %>% 
    filter(p.adjust < 0.05) %>% 
    nrow()
  

  
  if (tmp.num != 0) {
    # graph GO results
    # create plot
    tmp.enrichgo.fit <- barplot(tmp.enrichgo,
                                showCategory = 15) +
      ggtitle(paste0(i,
                     ': module enriched BP GO'))
    
    # save plot
    png(paste0('../Figures/WGCNA/GO/enrichplot/barplot/',
               i,
               ' module enriched BP GO.png'),
        res = 300,
        width = 10,
        height = 10,
        units = 'in')
    print(tmp.enrichgo.fit)
    dev.off()
    
    # create plot
    tmp.enrichgo.dot.fit <- dotplot(tmp.enrichgo,
                                    showCategory = 15) +
      ggtitle(paste0(i,
                     ': module enriched BP GO'))
    
    # save plot
    png(paste0('../Figures/WGCNA/GO/enrichplot/dotplot/',
               i,
               ' module enriched BP GO dotplot.png'),
        res = 300,
        width = 10,
        height = 10,
        units = 'in')
    print(tmp.enrichgo.dot.fit)
    dev.off()
    
  
    # png(paste0('../Figures/WGCNA/GO/enrichplot/',
    #            i,
    #            ' module enriched BP GO summary network.png'), 
    #     res = 300, 
    #     width = 10, 
    #     height = 10,
    #     units = 'in')
    # emapplot(
    #   ego,
    #   showCategory = 50,
    #   layout = "nicely"
    # ) +
    #   ggtitle(paste0(i, ": module enriched BP GO (summarized)"))
    # dev.off()
    # 
    # graph simplified semantic networks
    ego <- simplify(
      tmp.enrichgo,
      cutoff = 0.7,
      by = "p.adjust",
      select_fun = min,
      measure = "Wang"
    )
    
    ego <- pairwise_termsim(ego)
    
    # save results
    enrichGO.results.modules.simplify = enrichGO.results.modules.simplify %>%
      rbind(ego@result %>%
              mutate(module = i,
                     ont = 'BP'))
    
    write.csv(enrichGO.results.modules.simplify,
              '../Figures/WGCNA/GO/enrichGO.results.modules.simplify.csv',
              row.names = F)
  
    p.net = emapplot(
      ego,
      showCategory = 50,
      layout = "nicely"
    ) +
      ggtitle(paste0(i, ": module enriched BP GO (summarized)"))
    ggsave(paste0('../Figures/WGCNA/GO/enrichplot/network/',
               i,
               ' module enriched BP GO summary network.png'), 
        p.net,
        width = 10, 
        height = 10,
        units = 'in')
    
    
  }
  if (tmp.num == 0) {
    # save plot
    png(paste0('../Figures/WGCNA/GO/enrichplot/barplot/',
               i,
               ' module enriched BP GO.png'),
        res = 300,
        width = 10,
        height = 10,
        units = 'in')
    plot.new()
    dev.off()
    
    # save plot
    png(paste0('../Figures/WGCNA/GO/enrichplot/dotplot/',
               i,
               ' module enriched BP GO dotplot.png'),
        res = 300,
        width = 10,
        height = 10,
        units = 'in')
    plot.new()
    dev.off()
    
    # save plot
    png(paste0('../Figures/WGCNA/GO/enrichplot/network/',
               i,
               ' module enriched BP GO summary network.png'), 
        res = 300, 
        width = 10, 
        height = 10,
        units = 'in')
    plot.new()
    dev.off()
    
  
  }
}


## check enriched GO terms per module
enrichGO.results.modules %>% 
  filter(p.adjust < 0.05) %>% 
  dplyr::count(module)
# module   n
# 1     black  24
# 2      blue   8
# 3     brown  15
# 4      pink 104
# 5    purple  37
# 6       red   9
# 7 turquoise   9
# 8    yellow  38

enrichGO.results.modules.simplify %>% 
  filter(p.adjust < 0.05) %>% 
  dplyr::count(module)
# module  n
# 1     black 11
# 2      blue  3
# 3     brown  6
# 4      pink 35
# 5    purple  5
# 6       red  5
# 7 turquoise  3
# 8    yellow 24

# paper
### combined figures for paper
# make dummy enrich object that has all genes
tmp = tmp.enrichgo.ensembl

# replace results with results from modules that have less than 5 enriched terms
# tmp@result = enrichGO.results.modules  %>% 
#   filter(module %in% c(enrichGO.results.modules %>% 
#                          filter(p.adjust < 0.15) %>% 
#                          dplyr::count(module) %>% 
#                          filter(n<5) %>% 
#                          pull(module)))

tmp@result = enrichGO.results.modules.simplify  %>% 
  group_by(module) %>% 
  slice_max(n = 2,
        order_by = -p.adjust,
        with_ties = F)  %>% 
  full_join(enrichGO.results.modules.simplify  %>% 
              group_by(module) %>% 
              slice_max(n = 2,
                        order_by = Count,
                    with_ties = F)) %>% 
  full_join(
    enrichGO.results.modules.simplify %>% 
      filter(module %in% c(enrichGO.results.modules.simplify %>%
                             filter(p.adjust < 0.05) %>%
                             dplyr::count(module) %>%
                             filter(n<3) %>%
                             pull(module)))) %>% 
  distinct()

# check number of gene
tmp@result %>% 
       dplyr::count(module)

tmp@result %>% 
  dplyr::count(module) |> 
  pull(n) |> 
  sum()



# create dotplot with colors as names
p = tmp %>% 
  dotplot(showCategory = 50) + 
  aes(shape = I(21), 
      stroke = 5) + 
  aes(fill = -log10(p.adjust),
      color = module) +
  scale_color_manual(values = c("blue" = "blue",
                                "brown" = "brown",
                                "pink" = "pink",
                                "black" =  'grey25',
                                'green' = 'green',
                                'turquoise' = 'turquoise',
                                'magenta' = 'magenta',
                                'yellow' = 'yellow',
                                'red' = 'red',
                                'purple' = 'purple')) +
  scale_fill_continuous(low = "red",
                        high = "blue") +
  ggtitle(paste0(
    'Module enriched BP GO')) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.85,
                                   0.35))

# save plot
png(paste0('../Figures/WGCNA/GO/',
           'All module enriched BP GO dotplot.png'),
    res = 720,
    width = 6.5,
    height = 10,
    units = 'in')
print(p)
dev.off()

# save plot
pdf(paste0('../Figures/WGCNA/GO/',
           'All module enriched BP GO dotplot paper.pdf'),
    width = 6.5,
    height = 10)
print(p)
dev.off()


# #### Get GO terms ####
# ## load ensembl
# library(biomaRt)
# ###selecting biomart database
# ensembl = useMart('ensembl')
# # #list datasets
# # dataset = listDatasets(ensembl)
# ##select dataset
# #nile tilapia
# ensembl.tilapia = useDataset('oniloticus_gene_ensembl',
#                              mart=ensembl)
# 
# # # get list of all attributes
# # listAttributes(ensembl.tilapia) %>% View
# 
# 
# #create attributes lists
# tilapia.attributes = c('external_gene_name',
#                        'ensembl_gene_id',
#                        'go_id')
# 
# tilapia.attributes.2 = c('ensembl_gene_id',
#                          'go_id')
# 
# ##identify GO terms for WGCNA tilapia genes
# # use gene names
# tilapia.wgcna.go.terms.gene = getBM(attributes = tilapia.attributes,
#                                     mart = ensembl.tilapia,
#                                     values = KMEs$gene_id,
#                                     filter = 'external_gene_name',
#                                     useCache = FALSE) # useCache has to do with version of R not being up to date?
# # use ensembl gene IDs
# tilapia.wgcna.go.terms.ensemblID = getBM(attributes = tilapia.attributes.2,
#                                          mart = ensembl.tilapia,
#                                          values = KMEs$gene_id,
#                                          filter = 'ensembl_gene_id',
#                                          useCache = FALSE) # useCache has to do with version of R not being up to date?
# # combine gene and ensembl gene IDs
# tilapia.wgcna.go.terms = full_join(tilapia.wgcna.go.terms.gene,
#                                    tilapia.wgcna.go.terms.ensemblID)
# 
# # single gene column 
# tilapia.wgcna.go.terms = tilapia.wgcna.go.terms %>% 
#   mutate(gene_id = ifelse(is.na(external_gene_name),
#                           ensembl_gene_id,
#                           external_gene_name)) %>% 
#   dplyr::select(c(gene_id,
#                   go_id))
# 
# #check length
# #55117
# tilapia.wgcna.go.terms %>%
#   nrow()
# 
# #check number of tilapia genes?
# #12644
# tilapia.wgcna.go.terms %>%
#   pull(gene_id) %>%
#   unique() %>%
#   length()
# 
# #check number of GO IDs?
# #4909
# tilapia.wgcna.go.terms %>%
#   pull(go_id) %>%
#   unique() %>%
#   length()
# 
# ## check duplicates
# # 12644 - 54308 = 41664 duplicates
# tilapia.wgcna.go.terms %>%
#   distinct() %>% 
#   nrow()
# 
# # remove duplicates
# tilapia.wgcna.go.terms = tilapia.wgcna.go.terms %>% 
#   distinct()
# 
# ## collapse go_id into gene list
# tilapia.wgcna.go.terms = tilapia.wgcna.go.terms %>%
#   group_by(gene_id) %>%
#   summarize(go_id = str_c(go_id, collapse = ";"))
# 
# # add 'unknown' 
# tilapia.wgcna.go.terms = tilapia.wgcna.go.terms %>% 
#   mutate(go_id = ifelse(go_id == '',
#                         'unknown',
#                         go_id))
# 
# # save to csv 
# write.csv(tilapia.wgcna.go.terms,
#           file = '../Figures/WGCNA/go_terms/tilapia.wgcna.go.terms.csv')
# # save to tab delimited with no column names
# write_tsv(tilapia.wgcna.go.terms,
#           file = '../Figures/WGCNA/go_terms/tilapia.wgcna.go.terms.tsv',
#           col_names = FALSE)

# #### prepare for GO_MWU ####
# ### create data input list for each module (remove grey) 
# # want a list of genes, the kme if they are in that module, and a 0 kME if they are not in that module
# module.name.list = KMEs %>% 
#   pull(module) %>% 
#   unique() 
# # remove grey module
# module.name.list = module.name.list[! module.name.list %in% c("grey")]
# 
# ## loop through to create data frame for each module
# for (i in module.name.list) {
#   tmp.name = paste("kME",
#                    i,
#                    sep='')
#   # create temporary dataframe with module kme
#   tmp = KMEs %>% 
#     dplyr::select(c(gene_id,
#                     module,
#                     tmp.name)) %>% 
#     mutate(kme = ifelse(module != i,
#                         0,
#                         .[[tmp.name]])) %>% 
#     dplyr::select(c(gene_id,
#                     kme))
#   # save file
#   write.csv(tmp,
#             file = paste('../Figures/WGCNA/go_terms/wgcna.DMEs.neurons.kme.',
#                          i,
#                          '.csv',
#                          sep = ''),
#             row.names = FALSE,
#             quote = FALSE)
# }
# 
# 
# 
# ## tilapia.wgcna.go.terms.tsv for MWU_GO
# # need to run on desktop 
# # https://github.com/z0on/GO_MWU
# # https://github.com/schmidte10/GO-MWU-automation-

#### IEG/candidate gene analysis ####
## get list of IEGs
IEG.list = read.csv("IEG_Candidate_gene_list.csv") %>% 
  dplyr::filter(List == 'IEGS')

## subset module list
KMEs.IEG = KMEs %>% 
  dplyr::filter(gene_id %in% IEG.list$Name) %>% 
  rbind(KMEs %>% 
          dplyr::filter(gene_id %in% IEG.list$ID)) %>% 
  mutate(gene_type = 'IEG')

## get list of candidates
Candidate.list = read.csv("IEG_Candidate_gene_list.csv") %>% 
  dplyr::filter(List != 'IEGS')

## subset module list
KMEs.Candidate = KMEs %>% 
  dplyr::filter(gene_id %in% Candidate.list$Name) %>% 
  rbind(KMEs %>% 
          dplyr::filter(gene_id %in% Candidate.list$ID)) %>% 
  mutate(gene_type = 'Candidate')


## get list of candidates
Super.candidate.list = read.csv("IEG_Candidate_gene_list.csv") %>% 
  dplyr::filter(List == 'Super_Candidate_Genes')

## subset module list
KMEs.Super.candidate = KMEs %>% 
  dplyr::filter(gene_id %in% Super.candidate.list$Name) %>% 
  rbind(KMEs %>% 
          dplyr::filter(gene_id %in% Super.candidate.list$ID)) %>% 
  mutate(gene_type = 'Candidate')

## combine and save
KMEs.IEG.candidate = KMEs.IEG %>% 
  rbind(KMEs.Candidate) %>% 
  rbind(KMEs.Super.candidate) %>% 
  distinct()
# 
# write_csv(KMEs.IEG.candidate,
#           file = '../Figures/WGCNA/KMEs.IEG.candidate.csv')
# save(KMEs.IEG.candidate,
     # file = '../Figures/WGCNA/KMEs.IEG.candidate.Rdata')
# load('../Figures/WGCNA/KMEs.IEG.candidate.Rdata')
















#### WGCNA linear model ####
### load libraries 
library(WGCNA)
library(ComplexHeatmap)
library(circlize)

### load previous data
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
                                       '.'),
         ID = Sample.Name) %>% 
  separate(ID,
           c('Plate.position',
             'Round',
             'Tank',
             'Sex'))

## create Tank column 
sample.names = sample.names %>% 
  mutate(Observation.id = paste(Tank,
                                Round,
                                sep = '_'))

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

### prepare data
# already have MEs
## create list of traits
traitlist <- colnames(data.behavior.comp.reduce[,-c(1)])
## create datTraits
# need to add both males and females
# create types column for sex (1 = males, 0 = females)
# arrange by sample names
datTraits = data.behavior.comp.reduce %>% 
  full_join(sample.names) %>% 
  mutate(sex = ifelse(Sex == 'M',
                        1,
                        0)) %>% 
  arrange(Sample.Name)

## order MEs by sample names
MEs = MEs %>% 
  rownames_to_column('Sample.Name') %>% 
  arrange(Sample.Name) %>% 
  column_to_rownames('Sample.Name')

### linear model
## from Will's script
# replace "switch_time" with "sex"
# remove grey module
# generates 3 matrices corresponding to the 3 terms in the linear model: MEs[, j] ~ pop * trait
model <- as.formula("MEs[, j] ~ types * trait")
num_terms <- length(attr(terms(model), "term.labels"))



types_mat <- matrix(NA, length(traitlist), ncol(MEs))
rownames(types_mat) <- traitlist
colnames(types_mat) <- colnames(MEs)[1:ncol(MEs)]
trait_mat <- typesXtrait_mat <- types_mat

for (i in 1:length(traitlist)) {
  types <- datTraits$sex
  trait <- as.matrix(datTraits[traitlist[i]])
  print(colnames(trait))
  
  significance <- matrix(NA, ncol(MEs), num_terms)
  for (j in 1:(dim(MEs)[2])) {
    l <- lm(model)
    significance[j, ] <- car::Anova(l)[1:num_terms, 4]
  }
  rownames(significance) <- gsub("ME", "", colnames(MEs))
  colnames(significance) <- gsub("datTraits\\$", "", rownames(car::Anova(l))[1:num_terms])
  
  significance2 <- t(significance)
  
  types_mat[i, ] <- significance2[1, ]
  trait_mat[i, ] <- significance2[2, ]
  typesXtrait_mat[i, ] <- significance2[3, ]
}


# eFDR

num_perm <- 500
model_eFDR <- as.formula("newMEs[, j] ~ types * trait")
types_eFDR <- matrix(0, length(traitlist), ncol(MEs))
rownames(types_eFDR) <- traitlist
colnames(types_eFDR) <- colnames(MEs)[1:ncol(MEs)]
trait_eFDR <- typesXtrait_eFDR <- types_eFDR

for (z in 1:num_perm) {
  
  newMEs <- MEs[sample(nrow(MEs)), ]
  
  for (i in 1:length(traitlist)) {
    
    types <- datTraits$sex
    trait <- as.matrix(datTraits[traitlist[i]])
    print(z)
    print(colnames(trait))
    
    significance <- matrix(NA, ncol(newMEs), num_terms)
    for (j in 1:(dim(newMEs)[2])) {
      l <- lm(model_eFDR)
      significance[j, ] <- car::Anova(l)[1:num_terms, 4]
    }
    rownames(significance) <- gsub("ME", "", colnames(newMEs))
    colnames(significance) <- gsub("datTraits\\$", "", rownames(car::Anova(l))[1:num_terms])
    
    significance2 <- t(significance)
    
    types_eFDR[i, ] <- types_eFDR[i, ] + (significance2[1, ] <= types_mat[i, ])
    trait_eFDR[i, ] <- trait_eFDR[i, ] + (significance2[2, ] <= trait_mat[i, ])
    typesXtrait_eFDR[i, ] <- typesXtrait_eFDR[i, ] + (significance2[3, ] <= typesXtrait_mat[i, ])
  }
}

types_eFDR <- types_eFDR / num_perm
trait_eFDR <- trait_eFDR / num_perm
typesXtrait_eFDR <- typesXtrait_eFDR / num_perm

new_mat <- matrix(NA, nrow = dim(trait_eFDR)[2], ncol = dim(trait_eFDR)[1] * 2)
new_mat <- as.data.frame(new_mat)

for (i in 1:(dim(trait_eFDR)[1])) {
  new_mat[, i * 2 - 1] <- trait_eFDR[i, ]
  new_mat[, i * 2] <- typesXtrait_eFDR[i, ]
  colnames(new_mat)[i * 2 - 1] <- rownames(trait_eFDR)[i]
  colnames(new_mat)[i * 2] <- rownames(typesXtrait_eFDR)[i]
}
rownames(new_mat) <- colnames(trait_eFDR)


# dendrogram of modules
#MEs = MEs.dl.dom
l <- t(MEs)
# l <- l[rownames(l) != "MEgrey", ]
MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")


# dendrogram of scaled traits
traits_dend <- hclust(dist(t(scale(datTraits[traitlist]))), method = "average")
# save for later use in Illustrator
pdf("../Figures/WGCNA/linear_model/trait_dend.pdf")
plot(traits_dend)
dev.off()


# color functions for heatmap
# for trait and pop x trait p-values; I made the limit slighly above the desired value (0.1)
col_fun = colorRamp2(c(1e-10, 0.105), c("red", "white"))
# for population p-value summary
col_fun2 = colorRamp2(c(25, 0), c("blue", "white"))#the range is up to the trait terms number


# for my figure, population p-values were largely similar across linear models, so we summarized
# them in a single column indicating how many models had a signfificant population effect
types_sumStar <- apply(types_eFDR[-24,] < 0.05, 2, sum)


side_ha = rowAnnotation(Module = gsub("ME", "", rownames(new_mat)), types.Term = types_sumStar,
                        # color coding of modules for annotating module egiengene dendrogram; need to "automate" this in future
                        col = list(Module = c("darkred" = "darkred", "darkgreen" = "darkgreen", "magenta" = "magenta",
                                              "lightyellow" = "lightyellow", "cyan" = "cyan", "greenyellow" = "greenyellow",      
                                              "grey60" = "grey60", "turquoise" = "turquoise", 
                                              "blue" = "blue", "midnightblue" = "midnightblue",      
                                              "red" = "red", "black" = "black", "green" = "green",
                                              "lightgreen" = "lightgreen", "brown" = "brown", "pink" = "pink",
                                              "yellow" = "yellow", "royalblue" = "royalblue", "salmon" = "salmon", 
                                              "lightcyan" = "lightcyan", "purple" = "purple", "tan" = "tan", "grey" = "grey"), 
                                   # I did the Population.Term this way to make a discrete scale legend
                                   types.Term = c("0" = col_fun2(0), "1" = col_fun2(1), "2" = col_fun2(2), "3" = col_fun2(3),
                                                  "4" = col_fun2(4), "5" = col_fun2(5), "6" = col_fun2(6),"7" = col_fun2(7), "8" = col_fun2(8), "9" = col_fun2(9), "10" = col_fun2(10),
                                                  "11" = col_fun2(11), "12" = col_fun2(12), "13" = col_fun2(13),"14" = col_fun2(14), "15" = col_fun2(15), 
                                                  "16" = col_fun2(16), "17" = col_fun2(17),
                                                  "18" = col_fun2(18), "19" = col_fun2(19),"20" = col_fun2(20),"21" = col_fun2(21),"22" = col_fun2(22), "23" = col_fun2(23), "24" = col_fun2(24))),
                        gp = gpar(col = "black"),
                        show_annotation_name = TRUE
)


# heatmap
new_mat.test <- new_mat[,-c(47,48)]
ht <- Heatmap(as.matrix(new_mat), col = col_fun,
              # this outputs significance stars
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (is.na(new_mat[i, j]))
                  grid.text("NA", x, y, gp = gpar(fontsize = 4))
                else if (new_mat[i, j] < 0.001)
                  grid.text("***", x, y, gp = gpar(fontsize = 8))
                else if (new_mat[i, j] < 0.01 & new_mat[i, j] > 0.001)
                  grid.text("**", x, y, gp = gpar(fontsize = 8))
                else if (new_mat[i, j] < 0.05 & new_mat[i, j] > 0.01)
                  grid.text("*", x, y, gp = gpar(fontsize = 8))
              },
              # rows = module relationships
              cluster_rows = MEs_dend,
              show_row_dend = TRUE,
              row_names_side = "left",
              show_row_names = TRUE,
              # columns = traits, dendrogram added later in Illustrator
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              show_column_names = FALSE,
              column_names_side = "top",
              column_title_rot = 90,
              # titles
              name = "eFDR",
              row_title = "Module Eigengenes",
              # add annotations
              left_annotation = side_ha,
              # formatting of cells
              rect_gp = gpar(col = "black", lwd = 0.5),
              height = unit(6, "cm"),
              # split model terms by traits and label the column "pairs"
              # this is also where the columns are reordered to correspond to trait dendrogram
              column_split = factor(gsub("\\.typesXtrait", "", gsub("\\.trait", "", names(new_mat))), 
                                    levels = traits_dend$labels[traits_dend$order]),
              column_gap = unit(2, "mm")
)


# save as pdf
pdf("../Figures/WGCNA/linear_model/heatmap_example.pdf", 
    width = 6, height = 8)
print(ht)
dev.off()




#### WGCNA linear model males ####
### load libraries 
library(WGCNA)
library(ComplexHeatmap)
library(circlize)

### load previous data
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
                                       '.'),
         ID = Sample.Name) %>% 
  separate(ID,
           c('Plate.position',
             'Round',
             'Tank',
             'Sex'))

## create Tank column 
sample.names = sample.names %>% 
  mutate(Observation.id = paste(Tank,
                                Round,
                                sep = '_'))

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

### prepare data
# already have MEs
## create list of traits
traitlist <- colnames(data.behavior.comp.reduce[,-c(1)])
## create datTraits
# need to add both males and females
# create types column for sex (1 = males, 0 = females)
# arrange by sample names
datTraits = data.behavior.comp.reduce %>% 
  full_join(sample.names) %>% 
  mutate(sex = ifelse(Sex == 'M',
                      1,
                      0)) %>% 
  arrange(Sample.Name)

## order MEs by sample names
MEs = MEs %>% 
  rownames_to_column('Sample.Name') %>% 
  arrange(Sample.Name) %>% 
  column_to_rownames('Sample.Name')

## subset to males
datTraits = datTraits %>% 
  filter(Sex == 'M')

MEs = MEs %>% 
  rownames_to_column('Sample.Name') %>% 
  filter(Sample.Name %in% datTraits$Sample.Name) %>% 
  column_to_rownames('Sample.Name')


### linear model
## from Will's script
# replace "switch_time" with "sex"
# remove grey module
# generates 3 matrices corresponding to the 3 terms in the linear model: MEs[, j] ~ pop * trait
model <- as.formula("MEs[, j] ~  trait")
num_terms <- length(attr(terms(model), "term.labels"))



types_mat <- matrix(NA, length(traitlist), ncol(MEs))
rownames(types_mat) <- traitlist
colnames(types_mat) <- colnames(MEs)[1:ncol(MEs)]
trait_mat <- typesXtrait_mat <- types_mat

for (i in 1:length(traitlist)) {
  types <- datTraits$sex
  trait <- as.matrix(datTraits[traitlist[i]])
  print(colnames(trait))
  
  significance <- matrix(NA, ncol(MEs), num_terms)
  for (j in 1:(dim(MEs)[2])) {
    l <- lm(model)
    significance[j, ] <- car::Anova(l)[1:num_terms, 4]
  }
  rownames(significance) <- gsub("ME", "", colnames(MEs))
  colnames(significance) <- gsub("datTraits\\$", "", rownames(car::Anova(l))[1:num_terms])
  
  significance2 <- t(significance)
  
  # types_mat[i, ] <- significance2[1, ]
  trait_mat[i, ] <- significance2[1, ]
  # typesXtrait_mat[i, ] <- significance2[3, ]
}


# eFDR

num_perm <- 500
model_eFDR <- as.formula("newMEs[, j] ~  trait")
types_eFDR <- matrix(0, length(traitlist), ncol(MEs))
rownames(types_eFDR) <- traitlist
colnames(types_eFDR) <- colnames(MEs)[1:ncol(MEs)]
trait_eFDR <- typesXtrait_eFDR <- types_eFDR

for (z in 1:num_perm) {
  
  newMEs <- MEs[sample(nrow(MEs)), ]
  
  for (i in 1:length(traitlist)) {
    
    types <- datTraits$sex
    trait <- as.matrix(datTraits[traitlist[i]])
    print(z)
    print(colnames(trait))
    
    significance <- matrix(NA, ncol(newMEs), num_terms)
    for (j in 1:(dim(newMEs)[2])) {
      l <- lm(model_eFDR)
      significance[j, ] <- car::Anova(l)[1:num_terms, 4]
    }
    rownames(significance) <- gsub("ME", "", colnames(newMEs))
    colnames(significance) <- gsub("datTraits\\$", "", rownames(car::Anova(l))[1:num_terms])
    
    significance2 <- t(significance)
    
    # types_eFDR[i, ] <- types_eFDR[i, ] + (significance2[1, ] <= types_mat[i, ])
    trait_eFDR[i, ] <- trait_eFDR[i, ] + (significance2[1, ] <= trait_mat[i, ])
    # typesXtrait_eFDR[i, ] <- typesXtrait_eFDR[i, ] + (significance2[3, ] <= typesXtrait_mat[i, ])
  }
}

# types_eFDR <- types_eFDR / num_perm
trait_eFDR <- trait_eFDR / num_perm
# typesXtrait_eFDR <- typesXtrait_eFDR / num_perm

new_mat <- matrix(NA, nrow = dim(trait_eFDR)[2], ncol = dim(trait_eFDR)[1])
new_mat <- as.data.frame(new_mat)

for (i in 1:(dim(trait_eFDR)[1])) {
  new_mat[, i  ] <- trait_eFDR[i, ]
  # new_mat[, i * 2] <- typesXtrait_eFDR[i, ]
  colnames(new_mat)[i  ] <- rownames(trait_eFDR)[i]
  # colnames(new_mat)[i * 2] <- rownames(typesXtrait_eFDR)[i]
}
rownames(new_mat) <- colnames(trait_eFDR)


# dendrogram of modules
#MEs = MEs.dl.dom
l <- t(MEs)
# l <- l[rownames(l) != "MEgrey", ]
MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")


# dendrogram of scaled traits
traits_dend <- hclust(dist(t(scale(datTraits[traitlist]))), method = "average")
# save for later use in Illustrator
pdf("../Figures/WGCNA/linear_model/trait_dend males.pdf")
plot(traits_dend)
dev.off()


# color functions for heatmap
# for trait and pop x trait p-values; I made the limit slighly above the desired value (0.1)
col_fun = colorRamp2(c(1e-10, 0.105), c("red", "white"))
# for population p-value summary
col_fun2 = colorRamp2(c(25, 0), c("blue", "white"))#the range is up to the trait terms number


# for my figure, population p-values were largely similar across linear models, so we summarized
# them in a single column indicating how many models had a signfificant population effect
types_sumStar <- apply(types_eFDR[-24,] < 0.05, 2, sum)


side_ha = rowAnnotation(Module = gsub("ME", "", rownames(new_mat)), types.Term = types_sumStar,
                        # color coding of modules for annotating module egiengene dendrogram; need to "automate" this in future
                        col = list(Module = c("darkred" = "darkred", "darkgreen" = "darkgreen", "magenta" = "magenta",
                                              "lightyellow" = "lightyellow", "cyan" = "cyan", "greenyellow" = "greenyellow",      
                                              "grey60" = "grey60", "turquoise" = "turquoise", 
                                              "blue" = "blue", "midnightblue" = "midnightblue",      
                                              "red" = "red", "black" = "black", "green" = "green",
                                              "lightgreen" = "lightgreen", "brown" = "brown", "pink" = "pink",
                                              "yellow" = "yellow", "royalblue" = "royalblue", "salmon" = "salmon", 
                                              "lightcyan" = "lightcyan", "purple" = "purple", "tan" = "tan", "grey" = "grey"), 
                                   # I did the Population.Term this way to make a discrete scale legend
                                   types.Term = c("0" = col_fun2(0), "1" = col_fun2(1), "2" = col_fun2(2), "3" = col_fun2(3),
                                                  "4" = col_fun2(4), "5" = col_fun2(5), "6" = col_fun2(6),"7" = col_fun2(7), "8" = col_fun2(8), "9" = col_fun2(9), "10" = col_fun2(10),
                                                  "11" = col_fun2(11), "12" = col_fun2(12), "13" = col_fun2(13),"14" = col_fun2(14), "15" = col_fun2(15), 
                                                  "16" = col_fun2(16), "17" = col_fun2(17),
                                                  "18" = col_fun2(18), "19" = col_fun2(19),"20" = col_fun2(20),"21" = col_fun2(21),"22" = col_fun2(22), "23" = col_fun2(23), "24" = col_fun2(24))),
                        gp = gpar(col = "black"),
                        show_annotation_name = TRUE
)


# heatmap
new_mat.test <- new_mat[,-c(47,48)]
ht <- Heatmap(as.matrix(new_mat), col = col_fun,
              # this outputs significance stars
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (is.na(new_mat[i, j]))
                  grid.text("NA", x, y, gp = gpar(fontsize = 4))
                else if (new_mat[i, j] < 0.001)
                  grid.text("***", x, y, gp = gpar(fontsize = 8))
                else if (new_mat[i, j] < 0.01 & new_mat[i, j] > 0.001)
                  grid.text("**", x, y, gp = gpar(fontsize = 8))
                else if (new_mat[i, j] < 0.05 & new_mat[i, j] > 0.01)
                  grid.text("*", x, y, gp = gpar(fontsize = 8))
              },
              # rows = module relationships
              cluster_rows = MEs_dend,
              show_row_dend = TRUE,
              row_names_side = "left",
              show_row_names = TRUE,
              # columns = traits, dendrogram added later in Illustrator
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              show_column_names = FALSE,
              column_names_side = "top",
              column_title_rot = 90,
              # titles
              name = "eFDR",
              row_title = "Module Eigengenes",
              # add annotations
              left_annotation = side_ha,
              # formatting of cells
              rect_gp = gpar(col = "black", lwd = 0.5),
              height = unit(6, "cm"),
              # split model terms by traits and label the column "pairs"
              # this is also where the columns are reordered to correspond to trait dendrogram
              column_split = factor(gsub("\\.typesXtrait", "", gsub("\\.trait", "", names(new_mat))), 
                                    levels = traits_dend$labels[traits_dend$order]),
              column_gap = unit(2, "mm")
)


# save as pdf
pdf("../Figures/WGCNA/linear_model/heatmap_example males.pdf", 
    width = 6, height = 8)
print(ht)
dev.off()




#### WGCNA linear model females ####
### load libraries 
library(WGCNA)
library(ComplexHeatmap)
library(circlize)

### load previous data
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
                                       '.'),
         ID = Sample.Name) %>% 
  separate(ID,
           c('Plate.position',
             'Round',
             'Tank',
             'Sex'))

## create Tank column 
sample.names = sample.names %>% 
  mutate(Observation.id = paste(Tank,
                                Round,
                                sep = '_'))

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

### prepare data
# already have MEs
## create list of traits
traitlist <- colnames(data.behavior.comp.reduce[,-c(1)])
## create datTraits
# need to add both males and females
# create types column for sex (1 = males, 0 = females)
# arrange by sample names
datTraits = data.behavior.comp.reduce %>% 
  full_join(sample.names) %>% 
  mutate(sex = ifelse(Sex == 'M',
                      1,
                      0)) %>% 
  arrange(Sample.Name)

## order MEs by sample names
MEs = MEs %>% 
  rownames_to_column('Sample.Name') %>% 
  arrange(Sample.Name) %>% 
  column_to_rownames('Sample.Name')

## subset to males
datTraits = datTraits %>% 
  filter(Sex == 'F')

MEs = MEs %>% 
  rownames_to_column('Sample.Name') %>% 
  filter(Sample.Name %in% datTraits$Sample.Name) %>% 
  column_to_rownames('Sample.Name')


### linear model
## from Will's script
# replace "switch_time" with "sex"
# remove grey module
# generates 3 matrices corresponding to the 3 terms in the linear model: MEs[, j] ~ pop * trait
model <- as.formula("MEs[, j] ~  trait")
num_terms <- length(attr(terms(model), "term.labels"))



types_mat <- matrix(NA, length(traitlist), ncol(MEs))
rownames(types_mat) <- traitlist
colnames(types_mat) <- colnames(MEs)[1:ncol(MEs)]
trait_mat <- typesXtrait_mat <- types_mat

for (i in 1:length(traitlist)) {
  types <- datTraits$sex
  trait <- as.matrix(datTraits[traitlist[i]])
  print(colnames(trait))
  
  significance <- matrix(NA, ncol(MEs), num_terms)
  for (j in 1:(dim(MEs)[2])) {
    l <- lm(model)
    significance[j, ] <- car::Anova(l)[1:num_terms, 4]
  }
  rownames(significance) <- gsub("ME", "", colnames(MEs))
  colnames(significance) <- gsub("datTraits\\$", "", rownames(car::Anova(l))[1:num_terms])
  
  significance2 <- t(significance)
  
  # types_mat[i, ] <- significance2[1, ]
  trait_mat[i, ] <- significance2[1, ]
  # typesXtrait_mat[i, ] <- significance2[3, ]
}


# eFDR

num_perm <- 500
model_eFDR <- as.formula("newMEs[, j] ~  trait")
types_eFDR <- matrix(0, length(traitlist), ncol(MEs))
rownames(types_eFDR) <- traitlist
colnames(types_eFDR) <- colnames(MEs)[1:ncol(MEs)]
trait_eFDR <- typesXtrait_eFDR <- types_eFDR

for (z in 1:num_perm) {
  
  newMEs <- MEs[sample(nrow(MEs)), ]
  
  for (i in 1:length(traitlist)) {
    
    types <- datTraits$sex
    trait <- as.matrix(datTraits[traitlist[i]])
    print(z)
    print(colnames(trait))
    
    significance <- matrix(NA, ncol(newMEs), num_terms)
    for (j in 1:(dim(newMEs)[2])) {
      l <- lm(model_eFDR)
      significance[j, ] <- car::Anova(l)[1:num_terms, 4]
    }
    rownames(significance) <- gsub("ME", "", colnames(newMEs))
    colnames(significance) <- gsub("datTraits\\$", "", rownames(car::Anova(l))[1:num_terms])
    
    significance2 <- t(significance)
    
    # types_eFDR[i, ] <- types_eFDR[i, ] + (significance2[1, ] <= types_mat[i, ])
    trait_eFDR[i, ] <- trait_eFDR[i, ] + (significance2[1, ] <= trait_mat[i, ])
    # typesXtrait_eFDR[i, ] <- typesXtrait_eFDR[i, ] + (significance2[3, ] <= typesXtrait_mat[i, ])
  }
}

# types_eFDR <- types_eFDR / num_perm
trait_eFDR <- trait_eFDR / num_perm
# typesXtrait_eFDR <- typesXtrait_eFDR / num_perm

new_mat <- matrix(NA, nrow = dim(trait_eFDR)[2], ncol = dim(trait_eFDR)[1])
new_mat <- as.data.frame(new_mat)

for (i in 1:(dim(trait_eFDR)[1])) {
  new_mat[, i  ] <- trait_eFDR[i, ]
  # new_mat[, i * 2] <- typesXtrait_eFDR[i, ]
  colnames(new_mat)[i  ] <- rownames(trait_eFDR)[i]
  # colnames(new_mat)[i * 2] <- rownames(typesXtrait_eFDR)[i]
}
rownames(new_mat) <- colnames(trait_eFDR)


# dendrogram of modules
#MEs = MEs.dl.dom
l <- t(MEs)
# l <- l[rownames(l) != "MEgrey", ]
MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")


# dendrogram of scaled traits
traits_dend <- hclust(dist(t(scale(datTraits[traitlist]))), method = "average")
# save for later use in Illustrator
pdf("../Figures/WGCNA/linear_model/trait_dend females.pdf")
plot(traits_dend)
dev.off()


# color functions for heatmap
# for trait and pop x trait p-values; I made the limit slighly above the desired value (0.1)
col_fun = colorRamp2(c(1e-10, 0.105), c("red", "white"))
# for population p-value summary
col_fun2 = colorRamp2(c(25, 0), c("blue", "white"))#the range is up to the trait terms number


# for my figure, population p-values were largely similar across linear models, so we summarized
# them in a single column indicating how many models had a signfificant population effect
types_sumStar <- apply(types_eFDR[-24,] < 0.05, 2, sum)


side_ha = rowAnnotation(Module = gsub("ME", "", rownames(new_mat)), types.Term = types_sumStar,
                        # color coding of modules for annotating module egiengene dendrogram; need to "automate" this in future
                        col = list(Module = c("darkred" = "darkred", "darkgreen" = "darkgreen", "magenta" = "magenta",
                                              "lightyellow" = "lightyellow", "cyan" = "cyan", "greenyellow" = "greenyellow",      
                                              "grey60" = "grey60", "turquoise" = "turquoise", 
                                              "blue" = "blue", "midnightblue" = "midnightblue",      
                                              "red" = "red", "black" = "black", "green" = "green",
                                              "lightgreen" = "lightgreen", "brown" = "brown", "pink" = "pink",
                                              "yellow" = "yellow", "royalblue" = "royalblue", "salmon" = "salmon", 
                                              "lightcyan" = "lightcyan", "purple" = "purple", "tan" = "tan", "grey" = "grey"), 
                                   # I did the Population.Term this way to make a discrete scale legend
                                   types.Term = c("0" = col_fun2(0), "1" = col_fun2(1), "2" = col_fun2(2), "3" = col_fun2(3),
                                                  "4" = col_fun2(4), "5" = col_fun2(5), "6" = col_fun2(6),"7" = col_fun2(7), "8" = col_fun2(8), "9" = col_fun2(9), "10" = col_fun2(10),
                                                  "11" = col_fun2(11), "12" = col_fun2(12), "13" = col_fun2(13),"14" = col_fun2(14), "15" = col_fun2(15), 
                                                  "16" = col_fun2(16), "17" = col_fun2(17),
                                                  "18" = col_fun2(18), "19" = col_fun2(19),"20" = col_fun2(20),"21" = col_fun2(21),"22" = col_fun2(22), "23" = col_fun2(23), "24" = col_fun2(24))),
                        gp = gpar(col = "black"),
                        show_annotation_name = TRUE
)


# heatmap
new_mat.test <- new_mat[,-c(47,48)]
ht <- Heatmap(as.matrix(new_mat), col = col_fun,
              # this outputs significance stars
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (is.na(new_mat[i, j]))
                  grid.text("NA", x, y, gp = gpar(fontsize = 4))
                else if (new_mat[i, j] < 0.001)
                  grid.text("***", x, y, gp = gpar(fontsize = 8))
                else if (new_mat[i, j] < 0.01 & new_mat[i, j] > 0.001)
                  grid.text("**", x, y, gp = gpar(fontsize = 8))
                else if (new_mat[i, j] < 0.05 & new_mat[i, j] > 0.01)
                  grid.text("*", x, y, gp = gpar(fontsize = 8))
              },
              # rows = module relationships
              cluster_rows = MEs_dend,
              show_row_dend = TRUE,
              row_names_side = "left",
              show_row_names = TRUE,
              # columns = traits, dendrogram added later in Illustrator
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              show_column_names = FALSE,
              column_names_side = "top",
              column_title_rot = 90,
              # titles
              name = "eFDR",
              row_title = "Module Eigengenes",
              # add annotations
              left_annotation = side_ha,
              # formatting of cells
              rect_gp = gpar(col = "black", lwd = 0.5),
              height = unit(6, "cm"),
              # split model terms by traits and label the column "pairs"
              # this is also where the columns are reordered to correspond to trait dendrogram
              column_split = factor(gsub("\\.typesXtrait", "", gsub("\\.trait", "", names(new_mat))), 
                                    levels = traits_dend$labels[traits_dend$order]),
              column_gap = unit(2, "mm")
)


# save as pdf
pdf("../Figures/WGCNA/linear_model/heatmap_example females.pdf", 
    width = 6, height = 8)
print(ht)
dev.off()





#### WGCNA linear model control tank doesn't work ####
### load libraries 
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(lme4)

### load previous data
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
                                       '.'),
         ID = Sample.Name) %>% 
  separate(ID,
           c('Plate.position',
             'Round',
             'Tank',
             'Sex'))

## create Tank column 
sample.names = sample.names %>% 
  mutate(Observation.id = paste(Tank,
                                Round,
                                sep = '_'))

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

### prepare data
# already have MEs
## create list of traits
traitlist <- colnames(data.behavior.comp.reduce[,-c(1)])
## create datTraits
# need to add both males and females
# create types column for sex (1 = males, 0 = females)
# arrange by sample names
datTraits = data.behavior.comp.reduce %>% 
  full_join(sample.names) %>% 
  mutate(sex = ifelse(Sex == 'M',
                      1,
                      0)) %>% 
  arrange(Sample.Name)

## order MEs by sample names
MEs = MEs %>% 
  rownames_to_column('Sample.Name') %>% 
  arrange(Sample.Name) %>% 
  column_to_rownames('Sample.Name')

### linear model
## from Will's script
# replace "switch_time" with "sex"
# remove grey module
# generates 3 matrices corresponding to the 3 terms in the linear model: MEs[, j] ~ pop * trait
model <- as.formula("MEs[, j] ~ types * trait + (1 | trait:tank)")
num_terms <- length(attr(terms(model), "term.labels"))



types_mat <- matrix(NA, length(traitlist), ncol(MEs))
rownames(types_mat) <- traitlist
colnames(types_mat) <- colnames(MEs)[1:ncol(MEs)]
trait_mat <- typesXtrait_mat <- types_mat

for (i in 1:length(traitlist)) {
  types <- datTraits$sex
  trait <- as.matrix(datTraits[traitlist[i]])
  tank <- datTraits$Tank 
  print(colnames(trait))
  
  significance <- matrix(NA, ncol(MEs), num_terms)
  for (j in 1:(dim(MEs)[2])) {
    l <- lmer(model)
    significance[j, ] <- car::Anova(l)[1:num_terms, 4]
  }
  rownames(significance) <- gsub("ME", "", colnames(MEs))
  colnames(significance) <- gsub("datTraits\\$", "", rownames(car::Anova(l))[1:num_terms])
  
  significance2 <- t(significance)
  
  types_mat[i, ] <- significance2[1, ]
  trait_mat[i, ] <- significance2[2, ]
  typesXtrait_mat[i, ] <- significance2[3, ]
}


# eFDR

num_perm <- 500
model_eFDR <- as.formula("newMEs[, j] ~ types * trait")
types_eFDR <- matrix(0, length(traitlist), ncol(MEs))
rownames(types_eFDR) <- traitlist
colnames(types_eFDR) <- colnames(MEs)[1:ncol(MEs)]
trait_eFDR <- typesXtrait_eFDR <- types_eFDR

for (z in 1:num_perm) {
  
  newMEs <- MEs[sample(nrow(MEs)), ]
  
  for (i in 1:length(traitlist)) {
    
    types <- datTraits$sex
    trait <- as.matrix(datTraits[traitlist[i]])
    print(z)
    print(colnames(trait))
    
    significance <- matrix(NA, ncol(newMEs), num_terms)
    for (j in 1:(dim(newMEs)[2])) {
      l <- lm(model_eFDR)
      significance[j, ] <- car::Anova(l)[1:num_terms, 4]
    }
    rownames(significance) <- gsub("ME", "", colnames(newMEs))
    colnames(significance) <- gsub("datTraits\\$", "", rownames(car::Anova(l))[1:num_terms])
    
    significance2 <- t(significance)
    
    types_eFDR[i, ] <- types_eFDR[i, ] + (significance2[1, ] <= types_mat[i, ])
    trait_eFDR[i, ] <- trait_eFDR[i, ] + (significance2[2, ] <= trait_mat[i, ])
    typesXtrait_eFDR[i, ] <- typesXtrait_eFDR[i, ] + (significance2[3, ] <= typesXtrait_mat[i, ])
  }
}

types_eFDR <- types_eFDR / num_perm
trait_eFDR <- trait_eFDR / num_perm
typesXtrait_eFDR <- typesXtrait_eFDR / num_perm

new_mat <- matrix(NA, nrow = dim(trait_eFDR)[2], ncol = dim(trait_eFDR)[1] * 2)
new_mat <- as.data.frame(new_mat)

for (i in 1:(dim(trait_eFDR)[1])) {
  new_mat[, i * 2 - 1] <- trait_eFDR[i, ]
  new_mat[, i * 2] <- typesXtrait_eFDR[i, ]
  colnames(new_mat)[i * 2 - 1] <- rownames(trait_eFDR)[i]
  colnames(new_mat)[i * 2] <- rownames(typesXtrait_eFDR)[i]
}
rownames(new_mat) <- colnames(trait_eFDR)


# dendrogram of modules
#MEs = MEs.dl.dom
l <- t(MEs)
# l <- l[rownames(l) != "MEgrey", ]
MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")


# dendrogram of scaled traits
traits_dend <- hclust(dist(t(scale(datTraits[traitlist]))), method = "average")
# save for later use in Illustrator
pdf("../Figures/WGCNA/linear_model/trait_dend.pdf")
plot(traits_dend)
dev.off()


# color functions for heatmap
# for trait and pop x trait p-values; I made the limit slighly above the desired value (0.1)
col_fun = colorRamp2(c(1e-10, 0.105), c("red", "white"))
# for population p-value summary
col_fun2 = colorRamp2(c(25, 0), c("blue", "white"))#the range is up to the trait terms number


# for my figure, population p-values were largely similar across linear models, so we summarized
# them in a single column indicating how many models had a signfificant population effect
types_sumStar <- apply(types_eFDR[-24,] < 0.05, 2, sum)


side_ha = rowAnnotation(Module = gsub("ME", "", rownames(new_mat)), types.Term = types_sumStar,
                        # color coding of modules for annotating module egiengene dendrogram; need to "automate" this in future
                        col = list(Module = c("darkred" = "darkred", "darkgreen" = "darkgreen", "magenta" = "magenta",
                                              "lightyellow" = "lightyellow", "cyan" = "cyan", "greenyellow" = "greenyellow",      
                                              "grey60" = "grey60", "turquoise" = "turquoise", 
                                              "blue" = "blue", "midnightblue" = "midnightblue",      
                                              "red" = "red", "black" = "black", "green" = "green",
                                              "lightgreen" = "lightgreen", "brown" = "brown", "pink" = "pink",
                                              "yellow" = "yellow", "royalblue" = "royalblue", "salmon" = "salmon", 
                                              "lightcyan" = "lightcyan", "purple" = "purple", "tan" = "tan", "grey" = "grey"), 
                                   # I did the Population.Term this way to make a discrete scale legend
                                   types.Term = c("0" = col_fun2(0), "1" = col_fun2(1), "2" = col_fun2(2), "3" = col_fun2(3),
                                                  "4" = col_fun2(4), "5" = col_fun2(5), "6" = col_fun2(6),"7" = col_fun2(7), "8" = col_fun2(8), "9" = col_fun2(9), "10" = col_fun2(10),
                                                  "11" = col_fun2(11), "12" = col_fun2(12), "13" = col_fun2(13),"14" = col_fun2(14), "15" = col_fun2(15), 
                                                  "16" = col_fun2(16), "17" = col_fun2(17),
                                                  "18" = col_fun2(18), "19" = col_fun2(19),"20" = col_fun2(20),"21" = col_fun2(21),"22" = col_fun2(22), "23" = col_fun2(23), "24" = col_fun2(24))),
                        gp = gpar(col = "black"),
                        show_annotation_name = TRUE
)


# heatmap
new_mat.test <- new_mat[,-c(47,48)]
ht <- Heatmap(as.matrix(new_mat), col = col_fun,
              # this outputs significance stars
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (is.na(new_mat[i, j]))
                  grid.text("NA", x, y, gp = gpar(fontsize = 4))
                else if (new_mat[i, j] < 0.001)
                  grid.text("***", x, y, gp = gpar(fontsize = 8))
                else if (new_mat[i, j] < 0.01 & new_mat[i, j] > 0.001)
                  grid.text("**", x, y, gp = gpar(fontsize = 8))
                else if (new_mat[i, j] < 0.05 & new_mat[i, j] > 0.01)
                  grid.text("*", x, y, gp = gpar(fontsize = 8))
              },
              # rows = module relationships
              cluster_rows = MEs_dend,
              show_row_dend = TRUE,
              row_names_side = "left",
              show_row_names = TRUE,
              # columns = traits, dendrogram added later in Illustrator
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              show_column_names = FALSE,
              column_names_side = "top",
              column_title_rot = 90,
              # titles
              name = "eFDR",
              row_title = "Module Eigengenes",
              # add annotations
              left_annotation = side_ha,
              # formatting of cells
              rect_gp = gpar(col = "black", lwd = 0.5),
              height = unit(6, "cm"),
              # split model terms by traits and label the column "pairs"
              # this is also where the columns are reordered to correspond to trait dendrogram
              column_split = factor(gsub("\\.typesXtrait", "", gsub("\\.trait", "", names(new_mat))), 
                                    levels = traits_dend$labels[traits_dend$order]),
              column_gap = unit(2, "mm")
)


# save as pdf
pdf("../Figures/WGCNA/linear_model/heatmap_example.pdf", 
    width = 6, height = 8)
print(ht)
dev.off()



