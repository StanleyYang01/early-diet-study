## ---- load cleaned RNA-seq count file and design file ready
load("~/Documents/code/short-term-diet/01_EdgeR_glmQLFit/data/01_data_cleaning.RData")
source("~/Documents/code/01_function/my_edgeR.R") # tidyverse library is loaded in my_edgeR.R 

## ---- data for cortical tissue

out_prefix<-"Cortex.RMoutlier"
out_analysis <- "./analysis/"
out_figure <- "./figure/"

data.raw=data.raw.C
data.design=data.design.C
group=group.C

## ---- remove outlier sample CD.2weeks.Sed.C_17989C and get rid of all running groups
library(tidyverse)
data.raw <- data.raw %>% select(-contains("17989"),-contains("Run"))
data.design <- data.design %>% filter(!Sample_Name %in% c("17989C"), TERM_3 == "Sed")
group <- factor(data.design$Group)

## ---- edgeR object, filtering, normalization

# build edgeR object
library(edgeR)
y <- DGEList(data.raw, group=group, genes=row.names(data.raw)) # must specify
options(digits=3)
y$samples

# before filtering and normalization
log.cpm <- cpm(y, log=TRUE, prior.count=2) 
boxplot(log.cpm)
hist(apply(log.cpm, 2, median), xlab="median of log2(cpm)", main="")

# attach gene symbols to edgeR object 
y<- symbols(y)
detach("package:org.Mm.eg.db", unload=TRUE) # conflict with select function, so to detach
detach("package:AnnotationDbi", unload=TRUE)

head(y$genes)
## dropNAsymbols
y <- y[!is.na(y$genes$Symbol),] 

# keep gene more than 1cpm in at least 2 samples 
dim(y)
keep <- rowSums(cpm(y) > 1) >= 2 ## filter
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

# after filtering and before normalization
log.cpm <- cpm(y, log=TRUE, prior.count=2) 
boxplot(log.cpm)
hist(apply(log.cpm, 2, median), xlab="median of log2(cpm)", main="")

# normalization
y <- calcNormFactors(y)
y$samples
dim(y)

# after normalization
log.cpm <- cpm(y, log=TRUE, prior.count=2, normalized.lib.size=T) 
boxplot(log.cpm)
hist(apply(log.cpm, 2, median), xlab="median of log2(cpm)", main="")

## Explore the RNA seq data

# person corrleation
cpm <- cpm(y, normalized.lib.size=T) 
cpm.cor <- cor(cpm, method = "pearson")
hist(cpm.cor)

library(corrplot)
corrplot(cpm.cor, method="square", order="hclust", cl.lim = c(0.85, 1),  tl.col="purple", tl.cex = 0.75, type = "full", is.corr = FALSE)
colorlegend(colbar = grey(1:100 / 100), 1:10, col = "red", align = "l",
            xlim = c(0, 6), ylim = c(-0.5,-0.1), vertical = FALSE)

# MDS plot
  # no label version
par(xpd = T, mar = par()$mar + c(0,0,0,7))  # to make legends outside of the plot
colors <- rep(c("darkgreen", "red", "blue", "purple"), 2)
pch <- c(0,1,2,5,15,16,17,18)
plotMDS(y, top = 500, cex = 1, pch=pch[group], dim.plot = c(1,2), ndim = 2, gene.selection = "pairwise", col=colors[group]) # col = as.numeric(group) differentiate colors between groups
legend(1.2, 0.7,levels(group), pch=pch, col=colors)
par(mar=c(5, 4, 4, 2.3) + 0.1)

  # with label version, to identify outliers
par(xpd = T, mar = par()$mar + c(0,0,0,7))  # to make legends outside of the plot
colors <- rep(c("darkgreen", "red", "blue", "purple"), 2)
  # pch <- c(0,1,2,5,15,16,17,18)
plotMDS(y, top = 500, cex = 1, dim.plot = c(1,2), ndim = 2, gene.selection = "pairwise", col=colors[group]) # col = as.numeric(group) differentiate colors between groups
legend(1.2, 0.7,levels(group), pch=pch, col=colors)
par(mar=c(5, 4, 4, 2.3) + 0.1)

## ---- Pincicpal component analysis ---
logCPM.PCA<-log.cpm 
rownames(logCPM.PCA) <- y$genes$Symbol
#colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
colnames(logCPM.PCA) <- data.design$sample_short # get it into the y project

pca_original = prcomp(t(logCPM.PCA),scale=T, center=T)
pca_x <- pca_original$x
pca_table <- data.frame(pca_x, data.design)
x <- pca_original$sdev^2/sum(pca_original$sdev^2) # Proportion of Variance Explained for all components

## Scree plot
plot(x, xlab="Principal Component", ylab="Proportion of Variance Explained", type="b")
plot(cumsum(x), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", type="b")

## PCA plot 
library(ggplot2)
PCA_plot <- function(pca_table, PC_x, PC_y, color, shape){
  #PC_x,PC_y are type of interger
  #color, shape, are type of string
  g <- ggplot(pca_table, aes_string(x=names(pca_table[PC_x]), y=names(pca_table[PC_y]), color=color, shape=shape)) 
  g <- g + geom_point(alpha=0.7, size=3) 
  # g <- g + labs(color = "Group", shape="Tissue")
  g + labs(x = paste(names(pca_table[PC_x]), scales::percent(x[PC_x]),"variance explained", sep=" "), y=paste(names(pca_table[PC_y]), scales::percent(x[PC_y]),"variance explained", sep=" "))
  #filename <- paste()
  #ggsave(filename, width=7, height=7, units="in")
}

PCA_plot(pca_table, 1, 2, "TERM_1", "TERM_2")
PCA_plot(pca_table, 2, 3, "TERM_1", "TERM_2")

## ---- estimate dispersion

# design--------------------------------------------------------------
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design
# estimateDisp--------------------------------------------------------
y <- estimateDisp(y, design, robust=TRUE)
# plotBCV, width="3.8in", fig.cap="Scatterplot of the biological coefficient of variation (BCV) against the average abundance of each gene. The plot shows the square-root estimates of the common, trended and tagwise NB dispersions."----
plotBCV(y)
# glmQLFit------------------------------------------------------------
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
# QLDisp, out.width="3.8in", fig.cap="A plot of the quarter-root QL dispersion against the average abundance of each gene. Estimates are shown for the raw (before EB moderation), trended and squeezed (after EB moderation) dispersions. Note that the QL dispersions and trend shown here are relative to the NB dispersion trend shown in Figure~\ref{fig:plotBCV}."----
plotQLDisp(fit)
# df.prior------------------------------------------------------------
summary(fit$df.prior)


## ---- data manipulation for GLM

# transpose cpm, observation (row) is sample, variable (column) is each gene. 
cpm <- cpm(y, normalized.lib.size=T) 
cpm.t<-t(cpm)
colnames(cpm)
all(data.design$sample_short == rownames(cpm.t))
# add design file to transposed cpm, data manipulation using filter based on design file TERMs (TERM_1 - TERM_4)
cpm.t.meta <- cbind(cpm.t, data.design)
# select all the WD and Sed group
cpm.WD.t <-cpm.t.meta %>% filter(TERM_1 == "WD", TERM_3 == "Sed") 
# remove metadata and transpose table back, rows = genes, columns = samples
gene_number <- dim(cpm)[1]
gene_number
cpm.WD <- cpm.WD.t[, c(1: gene_number)] %>% t()  
# check whether metadata (which contains characters has been removed, values should all be doulbe type now)
all(map_lgl(cpm.WD, is.double))
# assign sample_name_short it to the colnames after transposition.
colnames(cpm.WD) <- as.character(cpm.WD.t$sample_short)

# select all the CD and Sed group, and calculate the mean of each group based on weeks
cpm.ave.CD.meta <- cpm.t.meta %>% filter(TERM_1 == "CD") %>% group_by(TERM_2) %>% summarise_all(mean)
# check the column boundray of cpm
cpm.ave.CD.meta[,1:2]
cpm.ave.CD.meta[,(gene_number+1):(gene_number+2)]
# remove meta column and transform
cpm.ave.CD <- cpm.ave.CD.meta[, 2:(gene_number+1)] %>% t() 
colnames(cpm.ave.CD) <- cpm.ave.CD.meta[,1] %>% unlist() %>% paste("CD.ave.", ., sep="")
all(map_lgl(cpm.ave.CD, is.double))

# The following is to subtract the mean of CD from WD of each week, and glm model the difference between WD and CD group
dim(cpm.WD)
dim(cpm.ave.CD)
cpm.clean <- data.frame(cpm.WD, cpm.ave.CD) ## column 1:14395 are WD count read, 14396:14399 are CD count read.
dim(cpm.clean)
head(cpm.clean)

# WD-CD.ave, result will be combined into D_all value
D_2wk <- cpm.clean %>% select(contains("2weeks")) %>% mutate_all(funs(Dif=.- CD.ave.2weeks)) %>% select(contains("C_Dif"))
D_4wk <- cpm.clean %>% select(contains("4weeks")) %>% mutate_all(funs(Dif=.- CD.ave.4weeks)) %>% select(contains("C_Dif"))
D_6wk <- cpm.clean %>% select(contains("6weeks")) %>% mutate_all(funs(Dif=.- CD.ave.6weeks)) %>% select(contains("C_Dif"))
D_8wk <- cpm.clean %>% select(contains("8weeks")) %>% mutate_all(funs(Dif=.- CD.ave.8weeks)) %>% select(contains("C_Dif"))

D_all <- data.frame(D_2wk, D_4wk, D_6wk, D_8wk) 
rownames(D_all) <- rownames(cpm.ave.CD)
head(D_all)

## ---- GLM, z ~ u + duration, z= y(WD)-y(CD)

# use durataion as predictor, duration matches the column of D_all
duration <- c(rep("2weeks", 6), rep("4weeks", 6), rep("6weeks", 6), rep("8weeks", 6))

# for each gene, generate GLM with duration as prediction, collect the summary of the model for each gene
glm.fit <- D_all %>% t() %>% as.data.frame() %>%
  map(~ lm(. ~  duration)) %>% 
  map(summary)

colnames(D_all)
glm.fit %>% names() %>% head()

# look at two examples "Apoe", "Pecam1" in glm model
y$genes[y$genes$Symbol=="Apoe",]
glm.fit[["ENSMUSG00000002985"]]
y$genes[y$genes$Symbol=="Pecam1",]
glm.fit[["ENSMUSG00000020717"]]

## ---- obtain the statistics from the GLM model

# Use home-made pipeline stat_table_all to: 
source("~/Documents/code/01_function/my_GLM.R") 

# 1. retrieve the R-square(r2) and the corrected p value (q) from the fit summary of each gene
# 2. filter genes based on qval<0.05, r2>0.5, arrange based on q (ascending).
# 3. from filtered genes, retrieve coefficient of estimate and p val for each gene for down stream analysis.
stat_table_all <- stat_table(glm.fit)
dim(stat_table_all)
head(stat_table_all)

# add the gene names to stat_table_all 
stat_table_all<-left_join(stat_table_all, y$genes, by=c("ENSMUSG_ID"="genes"))
dim(stat_table_all)
head(stat_table_all)

# export the gene list table for Kegg and GO search
write.table(stat_table_all, paste(out_analysis, out_prefix, "_gene_list.txt", sep=""), sep="\t", row=F, quote=F)

# export the background gene list for Kegg and GO search, subtract genes of stat_table_all from y$genes
bg_gene_list <- y$genes 
dim(bg_gene_list)
write.table(bg_gene_list, paste(out_analysis, out_prefix, "_bg_gene.txt", sep=""), sep="\t", row=F, quote=F)

## ---- Explore the GLM model

# principal component analysis of estimate of the coefficient
stat_est <- stat_table_all %>% select(contains("est"))

pca_original = prcomp(t(stat_est),scale=T, center=T)
pca_x <- data.frame(pca_original$x, estimate = rownames(pca_original$x))
pca_table <- pca_x
x <- pca_original$sdev^2/sum(pca_original$sdev^2) # Proportion of Variance Explained for all components

# Scree plot
plot(x, xlab="Principal Component", ylab="Proportion of Variance Explained", type="b")
plot(cumsum(x), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", type="b")
# PCA plot
g <- ggplot(pca_table, aes(x=pca_table[1], y=pca_table[2], color= estimate))
  g <- g + geom_point(alpha=0.7, size=3) 
  # g <- g + labs(color = "Group", shape="Tissue")
  g + labs(x = paste(names(pca_table[1]), scales::percent(x[1]),"variance explained", sep=" "), y=paste(names(pca_table[2]), scales::percent(x[2]),"variance explained", sep=" "))

# hirachical analysis of estimate of the coefficient
stat_est_h<- t(scale(t(stat_est))) %>% as.data.frame() %>% as.matrix()
  
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(stat_est_h, col=col.pan, Rowv=TRUE, scale="none", 
          trace="none", dendrogram="both", cexRow=0.5, cexCol=0.7, density.info="none",
          margin=c(9,4), lhei=c(2,10), lwid=c(2,6))  


# plot the estimate of coefficient of each gene
stat_table_est.gather <- stat_table_all %>% 
  select("ENSMUSG_ID", contains("_est")) %>% 
  gather(key=coeff_group, value=estimate, -ENSMUSG_ID)
str(stat_table_est.gather)
stat_table_est.gather$coeff_group <- factor(stat_table_est.gather$coeff_group,
                                            levels =c("X.Intercept._est", "duration4weeks_est", "duration6weeks_est","duration8weeks_est"),
                                            ordered = TRUE)
ggplot(stat_table_est.gather, aes(x=coeff_group,y=estimate, group=ENSMUSG_ID)) + 
  geom_line(alpha=0.1, color="red") 
ggplot(stat_table_est.gather, aes(x=coeff_group,y=estimate, group=ENSMUSG_ID))+ 
  geom_line(alpha=0.1, color="red") + scale_y_continuous(limits= c(-50, 75))

# To identify each gene’s critical predictor. binarized the estimate and significance matrices by hard cutoffs.
stat_pval <- stat_table_all %>% select(contains("_pval"))
head(stat_pval)
head(stat_est)

# 
stat_binary <- (apply(stat_pval, 2, function(x) x<0.05) & apply(stat_est, 2, function(x) abs(x)>0.2)) %>% as.data.frame()
rowSums(stat_binary) %>% table()

# summary of the numbers of genes passes the filter
apply(stat_binary, 2, table)

# sepertate out genes by their coefficients
col <- as.list(names(stat_binary))
genes_by_coeff <- map(col, ~ stat_table_all[["Symbol"]][stat_binary[[.]]])
names(genes_by_coeff) <- c("Intercept", "duration4wks", "duration6wks", "duration8wks")
glimpse(genes_by_coeff)


# to check the overlap of genes in different coefficient catagories
library(VennDiagram)
grid.newpage()
within(genes_by_coeff,
       draw.quad.venn(area1 = length(Intercept), 
                      area2 = length(duration4wks), 
                      area3 = length(duration6wks), 
                      area4 = length(duration8wks),
                      n12 = length(intersect(Intercept, duration4wks)), 
                      n13 = length(intersect(Intercept, duration6wks)),
                      n14 = length(intersect(Intercept, duration8wks)),
                      n23 = length(intersect(duration4wks, duration6wks)),
                      n24 = length(intersect(duration4wks, duration8wks)),
                      n34 = length(intersect(duration6wks, duration8wks)),
                      n123 = Intercept %>% intersect(duration4wks) %>% intersect(duration6wks) %>% length(),
                      n124 = Intercept %>% intersect(duration4wks) %>% intersect(duration8wks) %>% length(),
                      n134 = Intercept %>% intersect(duration6wks) %>% intersect(duration8wks) %>% length(),
                      n234 = duration4wks %>% intersect(duration6wks) %>% intersect(duration8wks) %>% length(),
                      n1234 = Intercept %>% intersect(duration4wks) %>% intersect(duration6wks) %>% intersect(duration8wks) %>% length(),
                      category = names(genes_by_coeff),
                      fill = c("orange", "red", "green", "blue"))
        )


# teasing about age specific difference
attach(genes_by_coeff)
area1 = Intercept 
area2 = duration4wks  
area3 = duration6wks  
area4 = duration8wks 
n12 = intersect(Intercept, duration4wks)  
n13 = intersect(Intercept, duration6wks) 
n14 = intersect(Intercept, duration8wks) 
n23 = intersect(duration4wks, duration6wks) 
n24 = intersect(duration4wks, duration8wks) 
n34 = intersect(duration6wks, duration8wks) 
n123 = Intercept %>% intersect(duration4wks) %>% intersect(duration6wks) 
n124 = Intercept %>% intersect(duration4wks) %>% intersect(duration8wks) 
n134 = Intercept %>% intersect(duration6wks) %>% intersect(duration8wks) 
n234 = duration4wks %>% intersect(duration6wks) %>% intersect(duration8wks) 
n1234 = Intercept %>% intersect(duration4wks) %>% intersect(duration6wks) %>% intersect(duration8wks)

D4 <- area2 %>% setdiff(n12) %>% setdiff(n13) %>% setdiff(n14) %>% setdiff(n23) %>% setdiff(n14) %>% setdiff(n24) %>% setdiff(n34) %>%
  setdiff(n123) %>% setdiff(n124) %>% setdiff(n234)%>% setdiff(n1234)

D6 <- area3 %>% setdiff(n12) %>% setdiff(n13) %>% setdiff(n14) %>% setdiff(n23) %>% setdiff(n14) %>% setdiff(n24) %>% setdiff(n34) %>%
  setdiff(n123) %>% setdiff(n124) %>% setdiff(n234)%>% setdiff(n1234)

D8 <- area4 %>% setdiff(n12) %>% setdiff(n13) %>% setdiff(n14) %>% setdiff(n23) %>% setdiff(n14) %>% setdiff(n24) %>% setdiff(n34) %>%
  setdiff(n123) %>% setdiff(n124) %>% setdiff(n234)%>% setdiff(n1234)

all_intersect <- n12 %>% union(n13) %>% union(n13) %>% union(n14) %>%
  union(n23) %>% union(n24) %>% union(n34) %>% union(n123) %>% union(n124) %>%
  union(n134) %>% union(n234) %>% union(n1234)

gene_list <-list(D4, D6, D8, all_intersect)
names(gene_list) <- c("D4", "D6", "D8", "all_inter")
path <- paste(out_analysis, out_prefix, names(gene_list), ".txt", sep="")
walk2(gene_list, path, write.table, sep="\t", row=F, quote=F, col.names=F)


# To classify genes in more detail, we collapsed the binarized matrix gene by gene and grouped genes from the unique collapsed patterns.
head(stat_binary)
colnames(stat_binary) <- c("Intercept", "duration4wks", "duration6wks", "duration8wks")
stat_binary <- stat_binary %>% add_column(Symbol=stat_table_all$Symbol, .before=1)
stat_binary <- stat_binary %>% mutate(pattern= paste(Intercept,duration4wks,duration6wks, duration8wks ,sep = "_"))
pattern_TF <- table(stat_binary$pattern) # %>% as.data.frame()
stat_binary.gather <- stat_binary %>% gather(key=coeff_group, value=value, -pattern, -Symbol)
stat_binary.gather$coeff_group <- factor(stat_binary.gather$coeff_group, 
                                         levels = c("Intercept","duration4wks", "duration6wks", "duration8wks"),
                                         ordered = TRUE)

ggplot(stat_binary.gather, aes(x = coeff_group, y = pattern, fill = value)) + geom_tile(colour = "white") +
  theme_bw() + xlab("") + ylab("") + coord_flip() +
  # scale_x_discrete(labels = c("Intercept", "duration4wks", "duration6wks", "duration8wks")) +
  scale_y_discrete(labels = pattern_TF) +
  scale_fill_manual(values = c("grey80", "firebrick1")) 

## find the genes that coeffient is remained on, following pattern:
# TRUE_TRUE_TRUE_TRUE   -> All_On :  genes whose coefficients of all category are significant   
# FALSE_TRUE_TRUE_TRUE  -> 4wk_to_8wk : genes whose coefficients of 4wk, 6wk and 8wk are significant
# FALSE_FALSE_TRUE_TRUE -> 6wk_to_8wk : genes whose coefficient of 6wk and 8wk are significant
# FALSE_FALSE_TRUE_TRUE -> 8wk:  genes whose coefficient of 8wk is significant
gene_pattern <- c("All_On", "4wk_to_8wk", "6wk_to_8wk", "8wk")
gene_pattern_binary <- c("TRUE_TRUE_TRUE_TRUE", "FALSE_TRUE_TRUE_TRUE", "FALSE_FALSE_TRUE_TRUE", "FALSE_FALSE_FALSE_TRUE") %>% as.list()
pattern_filter <- function(pattern_binary, gene_pattern_table){
          gene_pattern_table %>% 
                  filter(pattern==pattern_binary) %>% 
                  select(Symbol)
}

gene_list <- gene_pattern_binary %>% map(pattern_filter, gene_pattern_table=stat_binary)
names(gene_list) <- gene_pattern
path <- paste(out_analysis, out_prefix, names(gene_list), ".txt", sep="")
walk2(gene_list, path, write.table, sep="\t", row=F, quote=F, col.names=F)
str(gene_list)

save.image("~/Documents/code/short-term-diet/03_modified_GLM/data/Cortex_outlier_removed.RData")

# heatmap 
logCPM <- log.cpm
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- data.design$sample_short # get it into the y project
logCPM<- logCPM[stat_table_all$Symbol[1:100],]
      # get rid of the running group
logCPM<- t(scale(t(logCPM))) %>% as.data.frame() %>% as.matrix()
logCPM.order <- logCPM[, colnames(logCPM) %>% order()]

library(gplots)
#heatmap, columns ordered by the original sample order 
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM.order, col=col.pan, Rowv=TRUE, Colv= FALSE, scale="none", 
          trace="none", dendrogram="row", cexRow=0.5, cexCol=0.7, density.info="none",
          margin=c(9,4), lhei=c(2,10), lwid=c(2,6))

#heatmap, columns ordered by the clustering
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM.order, col=col.pan, Rowv=TRUE, scale="none", 
          trace="none", dendrogram="both", cexRow=0.5, cexCol=0.7, density.info="none",
          margin=c(9,4), lhei=c(2,10), lwid=c(2,6))

# ---- GO and Kegg pathway analysis
# Xulong's function
source("~/Documents/code/01_function/Kegg_function.R")
gk.diet <- myGK(stat_table_all$Symbol)
head(gk.diet$BP[, c("Term", "Pvalue")])
head(gk.diet$MF[, c("Term", "Pvalue")])
head(gk.diet$CC[, c("Term", "Pvalue")])
gk.diet$KEGG[, c("Term", "Pvalue")]

file_name<- list("./analysis/Cortex_BP.txt", "./analysis/Cortex_MF.txt", "./analysis/Cortext_CC.txt", "./analysis/Cortext_KEGG.txt")
walk2(gk.diet, file_name, write.table,  sep="\t", row=F, quote=F)
write.table(stat_table_all, "./analysis/diet_gene.txt", sep="\t", row=F, quote=F)

sessionInfo()

