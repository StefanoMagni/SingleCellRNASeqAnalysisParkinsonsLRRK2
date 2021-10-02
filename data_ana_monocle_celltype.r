# loadsetwd('/Volumes/Macintosh HD 2/projects/science_projects/dropseq/jonas_iPSC/')
# Put yours here
setwd('/Users/stefano.magni/ownCloud/Documents/AlexGroup/SingleCell/JonasProject/Code_Cleaned_For_Publication/Pipeline_R_Alex/') 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("Biobase")

# packages
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)

# BiocManager::install("HSMMSingleCell")
# BiocManager::install("matrix")
# BiocManager::install("monocle")
# install.packages("fastICA", type="binary")
library(HSMMSingleCell)
library(monocle) # Restart R Studio if it gives you errors.

library(RColorBrewer) # for color definition
library(corrgram)
library(Hmisc)

source('./heatmap.3.R') # modified heatmap function


#read data
mu_0 <-read.csv('./Day0_10and42/A_0_filt.csv',header=T, row.names=1)
mu_10<-read.csv('./Day0_10and42/A_10_filt.csv',header=T, row.names=1)
mu_14<-read.csv('./Day0_10and42/A_14_filt.csv',header=T, row.names=1)
mu_42<-read.csv('./Day0_10and42/A_42_filt.csv',header=T, row.names=1)

wt_0 <-read.csv('./Day0_10and42/B_0_filt.csv',header=T, row.names=1)
wt_10<-read.csv('./Day0_10and42/B_10_filt.csv',header=T, row.names=1)
wt_14<-read.csv('./Day0_10and42/B_14_filt.csv',header=T, row.names=1)
wt_42<-read.csv('./Day0_10and42/B_42_filt.csv',header=T, row.names=1)


#for dimessions ... also for colors
nwt0  <- dim(wt_0)
nwt10 <- dim(wt_10)
nwt14 <- dim(wt_14)
nwt42 <- dim(wt_42)

nmu0  <- dim(mu_0)
nmu10 <- dim(mu_10)
nmu14 <- dim(mu_14)
nmu42 <- dim(mu_42)

#checksum of cell numbers
nwt0[2]  +
nwt10[2]  +
nwt14[2]  +
nwt42[2]  +
nmu0[2]   +
nmu10[2]  +
nmu14[2]  +
nmu42[2]


#read gene lists
#phenotype lists
neuron <-read.delim('./gene_lists/list_ph_Neuron.txt',header=F)
astro  <-read.delim('./gene_lists/list_ph_Astrocyte.txt',header=F)
micro  <-read.delim('./gene_lists/list_ph_Microglia.txt',header=F)
ips    <-read.delim('./gene_lists/list_stem.txt',header=F)
oligo  <-read.delim('./gene_lists/list_ph_Oligodendrocyte.txt',header=F)
dopneu <-read.delim('./gene_lists/list_ph_mDANeurons.txt',header=F)
endo   <-read.delim('./gene_lists/list_ph_Endothelial.txt',header=F)

#for expression ploitting all
neuron_all <-read.delim('./gene_lists/list_ph_Neuron_all.txt',header=F)
astro_all  <-read.delim('./gene_lists/list_ph_Astrocyte_all.txt',header=F)
micro_all  <-read.delim('./gene_lists/list_ph_Microglia_all.txt',header=F)
ips_all    <-read.delim('./gene_lists/list_stem_all.txt',header=F)
oligo_all  <-read.delim('./gene_lists/list_ph_Oligodendrocyte_all.txt',header=F)
dopneu_all <-read.delim('./gene_lists/list_ph_mDANeurons_all.txt',header=F)
endo_all   <-read.delim('./gene_lists/list_ph_Endothelial_all.txt',header=F)

????????????
# color definitions

wtcodef <- colorRampPalette(c("royalblue1","blue","darkblue"))(4) #"blue for WT"
mucodef <- colorRampPalette(c("salmon","red","darkred"))(4) # for mutant

hmcols<-colorRampPalette(c("grey","orange","red"))(256) #for expression
cola <- c('grey75', rev(brewer.pal(11, 'RdYlBu')), 'black') # best for dropseq visualization


# to merge sets
dat1           <- merge(wt_0, wt_10, by= "row.names", all.x= T, all.y= T)
rownames(dat1) <- dat1$Row.names #reset rownames
dat1$Row.names <- NULL  #remove added rownames col

dat2 <- merge(dat1, wt_14, by= "row.names", all.x= T, all.y= T)
rownames(dat2) <- dat2$Row.names #reset rownames
dat2$Row.names <- NULL  #remove added rownames col

allwt <- merge(dat2, wt_42, by= "row.names", all.x= T, all.y= T)
rownames(allwt) <- allwt$Row.names #reset rownames
allwt$Row.names <- NULL  #remove added rownames col

#mutant
dat1           <- merge(mu_0, mu_10, by= "row.names", all.x= T, all.y= T)
rownames(dat1) <- dat1$Row.names #reset rownames
dat1$Row.names <- NULL  #remove added rownames col

dat2 <- merge(dat1, mu_14, by= "row.names", all.x= T, all.y= T)
rownames(dat2) <- dat2$Row.names #reset rownames
dat2$Row.names <- NULL  #remove added rownames col

allmu <- merge(dat2, mu_42, by= "row.names", all.x= T, all.y= T)
rownames(allmu) <- allmu$Row.names #reset rownames
allmu$Row.names <- NULL  #remove added rownames col

allcell <- merge(allwt, allmu, by= "row.names", all.x= T, all.y= T)
rownames(allcell) <- allcell$Row.names #reset rownames
allcell$Row.names <- NULL  #remove added rownames col

dat1 <- 0
dat2 <- 1

# to merge sets in a day-wise
dat1           <- merge(wt_0, mu_0, by= "row.names", all.x= T, all.y= T)
rownames(dat1) <- dat1$Row.names #reset rownames
dat1$Row.names <- NULL  #remove added rownames col

dat2 <- merge(dat1, wt_10, by= "row.names", all.x= T, all.y= T)
rownames(dat2) <- dat2$Row.names #reset rownames
dat2$Row.names <- NULL  #remove added rownames col

dat3 <- merge(dat2, mu_10, by= "row.names", all.x= T, all.y= T)
rownames(dat3) <- dat3$Row.names #reset rownames
dat3$Row.names <- NULL  #remove added rownames col

dat4 <- merge(dat3, wt_14, by= "row.names", all.x= T, all.y= T)
rownames(dat4) <- dat4$Row.names #reset rownames
dat4$Row.names <- NULL  #remove added rownames col

dat5 <- merge(dat4, mu_14, by= "row.names", all.x= T, all.y= T)
rownames(dat5) <- dat5$Row.names #reset rownames
dat5$Row.names <- NULL  #remove added rownames col

dat6 <- merge(dat5, wt_42, by= "row.names", all.x= T, all.y= T)
rownames(dat6) <- dat6$Row.names #reset rownames
dat6$Row.names <- NULL  #remove added rownames col

dat7 <- merge(dat6, mu_42, by= "row.names", all.x= T, all.y= T)
rownames(dat7) <- dat7$Row.names #reset rownames
dat7$Row.names <- NULL  #remove added rownames col

dat1 <- 0
dat2 <- 1
dat3 <- 0
dat4 <- 1
dat5 <- 0
dat6 <- 1

allcell_day <- dat7
dat7 <-0
#replace NA by 0
allcell[is.na(allcell)] <- 0
allcell_day[is.na(allcell_day)] <- 0

#backuped

#FROM here for expression matrix # STE: NEXT 6 COMMANDS RETURN ERRORS, SKIP
cellcol <-c(rep(wtcodef[1],length(colnames(wt_0))),rep(wtcodef[2],length(colnames(wt_10))),rep(wtcodef[3],length(colnames(wt_14))),rep(wtcodef[4],length(colnames(wt_42))),rep(mucodef[1],length(colnames(mu_0))),rep(mucodef[2],length(colnames(mu_10))),rep(mucodef[3],length(colnames(mu_14))),rep(mucodef[4],length(colnames(mu_42))))

cellcol_day <-c(rep(wtcodef[1],length(colnames(wt_0))), rep(mucodef[1],length(colnames(mu_0))), rep(wtcodef[2],length(colnames(wt_10))), rep(mucodef[2],length(colnames(mu_10))), rep(wtcodef[3],length(colnames(wt_14))), rep(mucodef[3],length(colnames(mu_14))), rep(wtcodef[4],length(colnames(wt_42))),rep(mucodef[4],length(colnames(mu_42))))

### for plotting expression
column_annotation <- cellcol
column_annotation <- as.matrix(column_annotation)

column_annotation_day <- cellcol_day
column_annotation_day <- as.matrix(column_annotation_day)


#pdf('expression_viz_no_cluster.pdf')
#DO NOT RUN IF NOT ENOUGHT TIME!   STE: SKIP THIS AND NEXT BLOCK
dim_matrix <- dim(allcell)
jpeg("expression_viz_no_cluster.jpg", width=1200, height=900)
heatmap.3((as.matrix(log(log(allcell[,]+1)+1))),col=cola,trace="none",dendrogram="none",Rowv=F,Colv=F,key=T,KeyValueName="Expression",ColSideColors=column_annotation, keysize=0.6, distfun=function (y) dist(y,method = "euclidean"))
dev.off()

jpeg("expression_viz_no_cluster_day.jpg", width=1200, height=900)
heatmap.3((as.matrix(log(log(allcell_day[,]+1)+1))),col=cola,trace="none",dendrogram="none",Rowv=F,Colv=F,key=T,KeyValueName="Expression",ColSideColors=column_annotation_day, keysize=0.7, distfun=function (y) dist(y,method = "euclidean"))
dev.off()



#for reducing number of genes
siggenes <- rowSums(allcell_day)>50 #only genes that are expressed in 50 0r more cells
allcell_day_red <- allcell_day[siggenes,]

#for estimate dendogramm e.g. always the first 100 most expressed and 50 fewest expressed cells of day 0 and 42 of WT and Mut
#for Kamil ... maybe only use a subset of 10000 cellls
todendocells <- allcell_day_red[c(1:100,450:500, 501:600,950:1000, 2355:2455,2506:2655, 2656:2755,3005:3055)]

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x,method="euclidean")

#then here ordering genes for a subset of cells
hh<-heatmap.3(as.matrix(log(log(todendocells+1)+1)),dendrogram="row", trace="none", margin=c(8,9),  hclust=hclustfunc,distfun=distfunc);
#row indices are then in hh$rowInd
#here backuped

#DONT RUN IF IN A REAL HURRY (;   STE: SKIP THIS

jpeg("expression_viz_estimated_clustered.jpg", width=1200, height=900)
heatmap.3((as.matrix(log(log(allcell_day_red[hh$rowInd,]+1)+1))),col=cola,trace="none",dendrogram="none",Rowv=F,Colv=F,key=T,KeyValueName="Expression",ColSideColors=column_annotation_day, keysize=0.7, distfun=function (y) dist(y,method = "euclidean"))
dev.off()


#to put this into CellDataSet

genenames_all <- rownames(allcell)
phenoano_all <- data.frame(c(colnames(wt_0),colnames(wt_10),colnames(wt_14),colnames(wt_42),colnames(mu_0),colnames(mu_10),colnames(mu_14),colnames(mu_42)), c(rep("WT",length(colnames(wt_0))),rep("WT",length(colnames(wt_10))),rep("WT",length(colnames(wt_14))),rep("WT",length(colnames(wt_42))),rep("MUT",length(colnames(mu_0))),rep("MUT",length(colnames(mu_10))),rep("MUT",length(colnames(mu_14))),rep("MUT",length(colnames(mu_42)))),c(rep(0,length(colnames(wt_0))),rep(10,length(colnames(wt_10))),rep(14,length(colnames(wt_14))),rep(42,length(colnames(wt_42))),rep(0,length(colnames(mu_0))),rep(10,length(colnames(mu_10))),rep(14,length(colnames(mu_14))),rep(42,length(colnames(mu_42)))),c(rep("0",length(colnames(wt_0))),rep("10",length(colnames(wt_10))),rep("14",length(colnames(wt_14))),rep("42",length(colnames(wt_42))),rep("0",length(colnames(mu_0))),rep("10",length(colnames(mu_10))),rep("14",length(colnames(mu_14))),rep("42",length(colnames(mu_42)))))

geneano_all <- data.frame(genenames_all, rep("coding",length(genenames_all)), rep(0,length(genenames_all)))
colnames(geneano_all) <- c("gene_short_name", "cell_line", "specs") # !!!ESSENTIAL TO HAVE "gene_short_name" !!!

colnames(phenoano_all) <- c("cell_sample","cell_line", "days", "Days" ) #

aHSMM_expr_matrix <- as.matrix(allcell)
rownames(aHSMM_expr_matrix) <- c(1:dim(aHSMM_expr_matrix)[1]) # this has to match with afd

samplenames <- paste(phenoano_all[,1],phenoano_all[,2],phenoano_all[,4],sep="_")

colnames(aHSMM_expr_matrix) <- samplenames # this has to match with mfd

aHSMM_sample_sheet <- phenoano_all
aHSMM_gene_annotation <- geneano_all

afd <- new("AnnotatedDataFrame", data = aHSMM_gene_annotation)
apd <- new("AnnotatedDataFrame", data = aHSMM_sample_sheet)
rownames(apd)<- colnames(aHSMM_expr_matrix) # this has to match !
colnames(apd)<- colnames(phenoano_all)#c("cell_sample","cell_line", "cell_type", "days", "nsc_score", "endo_score", "oligo_score", "micro_score", "astro_score", "neuro_score","DAneuro_score" ) # this has to match !

# without sparese matrix
aHSMM <- newCellDataSet(as.matrix(aHSMM_expr_matrix), phenoData = apd, featureData = afd, expressionFamily=negbinomial.size())

#aHSMM <- newCellDataSet(as(as.matrix(aHSMM_expr_matrix), "sparseMatrix"),phenoData = apd, featureData = afd,lowerDetectionLimit=0.5, expressionFamily=negbinomial.size())

#naHSMM <- newCellDataSet(as(as.matrix(aHSMM_expr_matrix), "sparseMatrix"),phenoData = apd, featureData = afd,lowerDetectionLimit=0.5, expressionFamily=negbinomial.size())

# preparing for monocle- for this aHSMM has to be initialized to get gene names etc
#
# NOTICE NAME MATCH: The expression value matrix must have the same number of columns as the phenoData has rows, and it must have the same number of rows as the featureData data frame has rows. Row names of the phenoData object should match the column names of the expression matrix. Row names of the featureData object should match row names of the expression matrix. Also, one of the columns of the featureData must be named ???gene short name???.
#
#
#for (1) MEAN OF RELATED EXPRESSED GENES
score_neuron_avexp <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(neuron), ])
score_astro_avexp  <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(astro), ])
score_micro_avexp  <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(micro), ])
score_ips_avexp    <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(ips), ])
score_oligo_avexp  <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(oligo), ])
score_dopneu_avexp <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(dopneu), ])
score_endo_avexp   <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(endo), ])


# for (2) NUMBER OF RELEVANT GENES with 1 as treshhold
score_neuron <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(neuron), ]>=1)
score_astro  <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(astro), ]>=1)
score_micro  <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(micro), ]>=1)
score_ips    <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(ips), ]>=1)
score_oligo  <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(oligo), ]>=1)
score_dopneu <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(dopneu), ]>=1)
score_endo   <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(endo), ]>=1)


# for (2) NUMBER OF RELEVANT GENES with tre as treshhold
tre = 2
score_ips    <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(ips), ]>=tre)
score_endo   <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(endo), ]>=tre)
score_oligo  <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(oligo), ]>=tre)
score_micro  <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(micro), ]>=tre)
score_astro  <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(astro), ]>=tre)
score_neuron <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(neuron), ]>=tre)
score_dopneu <- colMeans(allcell[fData(aHSMM)$gene_short_name %in% unlist(dopneu), ]>=tre)

score_celltype <- rbind(score_ips,score_endo, score_oligo, score_micro ,score_astro,score_neuron, score_dopneu)


#pdf('score_celltype_continous_score.pdf')
#plot(as.data.frame(t(score_celltype)))
#dev.off()


#pdf('score_celltype_correlaiton_continous_score.pdf')
#corrgram(t(score_celltype),lower.panel=panel.shade, upper.panel=panel.pie,diag.panel=panel.density,gap=.2,col.regions=colorRampPalette(c("navy","royalblue","white","salmon","red")))
#dev.off()


#normalizing cell type
score_celltype_norm <- score_celltype/rowSums(score_celltype )

#library(matrixStats)

tomycelltype_norm <- max.col(t(score_celltype_norm),ties.method ="last")
tomycelltype <- max.col(t(score_celltype),ties.method="last")

posscelltype <- c("IPS", "Endothelia" ,"Oligodendrocyte", "Microglia", "Astrocyte", "Neuron" ,"DA_Neuron")

mycelltype <- (posscelltype[tomycelltype])
mycelltype_norm <- (posscelltype[tomycelltype_norm])

#to check for ambigous cells ... with same max for different cell scores! check to define undefined cells after this!
masa <-(lapply(1:dim(score_celltype)[2], function(i) which(score_celltype[,i] == max(score_celltype[,i]))))
ambicell <- unlist(lapply(1:length(masa), function(i) length(unlist(masa[c(i)][]))>1))
#same for the normalized version
masa_norm <-(lapply(1:dim(score_celltype_norm)[2], function(i) which(score_celltype_norm[,i] == max(score_celltype_norm[,i]))))
ambicell_norm <- unlist(lapply(1:length(masa_norm), function(i) length(unlist(masa_norm[c(i)][]))>1))

mycelltype[ambicell] <- "Ambiguous"
mycelltype_norm[ambicell_norm] <- "Ambiguous"

#to check for empty all scores -> undefined - should be the same for both norm and raw
undefcell <-(colSums(score_celltype==0)==7) #cells that have no signature at all
undefcell_norm <-(colSums(score_celltype_norm==0)==7) #cells that have no signature at all

mycelltype[undefcell] <- "Undefined"
mycelltype_norm[undefcell_norm] <- "Undefined"

score_ambi <- c(rep(0,length(score_ips)))
score_undef <- c(rep(0,length(score_ips)))

score_ambi_norm <- c(rep(0,length(score_ips)))
score_undef_norm <- c(rep(0,length(score_ips)))

score_ambi[ambicell] <- 2* max(score_celltype)
score_undef[undefcell] <- 4* max(score_celltype)

score_ambi_norm[ambicell_norm] <- 2* max(score_celltype_norm)
score_undef_norm[undefcell_norm] <- 4* max(score_celltype_norm)

score_celltype_all <- rbind(score_celltype,score_ambi,score_undef)
score_celltype_all_norm <- rbind(score_celltype_norm,score_ambi_norm,score_undef_norm)

summary(as.data.frame(mycelltype),maxsum=10)
summary(as.data.frame(mycelltype_norm),maxsum=10)
#backuped


pdf('score_celltype_correlaiton_binary_score.pdf')
corrgram(t(score_celltype),lower.panel=panel.shade, upper.panel=panel.pie,diag.panel=panel.density,gap=.2,col.regions=colorRampPalette(c("navy","royalblue","white","salmon","red")))
dev.off()


#for cell type plotting
pie <- ggplot(as.data.frame(mycelltype), aes(x = factor(1), fill = factor(mycelltype))) +
geom_bar(width = 1)

pdf('celltypes.pdf')
pie + coord_polar(theta = "y") + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
dev.off()

summary(as.data.frame(mycelltype_norm))

pie_norm <- ggplot(as.data.frame(mycelltype_norm), aes(x = factor(1), fill = factor(mycelltype_norm))) +
geom_bar(width = 1)

pdf('celltypes_norm.pdf')
pie_norm + coord_polar(theta = "y") +  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
dev.off()
#rownames(allcell), gene_short_name %in% unlist(neuron)))

#pathway lists -not needed here!
mito      <-read.delim('./gene_lists/list_Mito.txt',header=F)
cellcyc   <-read.delim('./gene_lists/list_CellCycle.txt',header=F)

#### END DUMMY CELL THING
########### RESTART FROM HERE ###########
### HERE OBJECT WITH CELL TYPE
phenoano_all_ct <- data.frame(c(colnames(wt_0),colnames(wt_10),colnames(wt_14),colnames(wt_42),colnames(mu_0),colnames(mu_10),colnames(mu_14),colnames(mu_42)), c(rep("WT",length(colnames(wt_0))),rep("WT",length(colnames(wt_10))),rep("WT",length(colnames(wt_14))),rep("WT",length(colnames(wt_42))),rep("MUT",length(colnames(mu_0))),rep("MUT",length(colnames(mu_10))),rep("MUT",length(colnames(mu_14))),rep("MUT",length(colnames(mu_42)))),mycelltype_norm,c(rep(0,length(colnames(wt_0))),rep(10,length(colnames(wt_10))),rep(14,length(colnames(wt_14))),rep(42,length(colnames(wt_42))),rep(0,length(colnames(mu_0))),rep(10,length(colnames(mu_10))),rep(14,length(colnames(mu_14))),rep(42,length(colnames(mu_42)))),c(rep("0",length(colnames(wt_0))),rep("10",length(colnames(wt_10))),rep("14",length(colnames(wt_14))),rep("42",length(colnames(wt_42))),rep("0",length(colnames(mu_0))),rep("10",length(colnames(mu_10))),rep("14",length(colnames(mu_14))),rep("42",length(colnames(mu_42)))),score_ips_avexp, score_endo_avexp, score_oligo_avexp, score_micro_avexp, score_astro_avexp, score_neuron_avexp,score_dopneu_avexp)

geneano_all_ct <- data.frame(genenames_all, rep("coding",length(genenames_all)), rep(0,length(genenames_all)))
colnames(geneano_all_ct) <- c("gene_short_name", "cell_line", "specs") # !!!ESSENTIAL TO HAVE "gene_short_name" !!!

#colnames(phenoano_all_ct) <- c("cell_sample","cell_line", "days", "Days" ) #

colnames(phenoano_all_ct) <- c("cell_sample","cell_line", "cell_type", "days", "DAYS", "nsc_score", "endo_score", "oligo_score", "micro_score", "astro_score", "neuro_score","DAneuro_score" ) #

ctHSMM_expr_matrix <- as.matrix(allcell)
rownames(ctHSMM_expr_matrix) <- c(1:dim(aHSMM_expr_matrix)[1]) # this has to match with afd

samplenames_ct <- paste(phenoano_all_ct[,1],phenoano_all_ct[,2],phenoano_all_ct[,4],sep="_")

colnames(ctHSMM_expr_matrix) <- samplenames_ct # this has to match with mfd

ctHSMM_sample_sheet <- phenoano_all_ct
ctHSMM_gene_annotation <- geneano_all_ct

ctfd <- new("AnnotatedDataFrame", data = ctHSMM_gene_annotation)
ctpd <- new("AnnotatedDataFrame", data = ctHSMM_sample_sheet)
rownames(ctpd)<- colnames(ctHSMM_expr_matrix) # this has to match !
colnames(ctpd)<- colnames(phenoano_all_ct)#c("cell_sample","cell_line", "cell_type", "days", "nsc_score", "endo_score", "oligo_score", "micro_score", "astro_score", "neuro_score","DAneuro_score" ) # this has to match !

# without sparese matrixaHSMM <- newCellDataSet(as.matrix(mHSMM_expr_matrix), phenoData = mpd, featureData = mfd, expressionFamily=negbinomial.size())

#ctHSMM <- newCellDataSet(as(as.matrix(ctHSMM_expr_matrix), "sparseMatrix"),phenoData = ctpd, featureData = ctfd,lowerDetectionLimit=0.5, expressionFamily=negbinomial.size())

ctHSMM <- newCellDataSet(as(as.matrix(ctHSMM_expr_matrix), "sparseMatrix"),phenoData = ctpd, featureData = ctfd, expressionFamily=negbinomial.size())

### END OF CELLT PYE OBJECT DEFINITIONS

#for this aHSMM has to be initiallized!!!
neuron_id <- row.names(subset(fData(aHSMM), gene_short_name %in% unlist(neuron)))
astro_id <- row.names(subset(fData(aHSMM), gene_short_name %in% unlist(astro)))
micro_id <- row.names(subset(fData(aHSMM), gene_short_name %in% unlist(micro)))
ips_id <- row.names(subset(fData(aHSMM), gene_short_name %in% unlist(ips)))
oligo_id <-  row.names(subset(fData(aHSMM), gene_short_name %in% unlist(oligo)))
dopneu_id <-  row.names(subset(fData(aHSMM), gene_short_name %in% unlist(dopneu)))
endo_id <-  row.names(subset(fData(aHSMM), gene_short_name %in% unlist(endo)))



aHSMM <- detectGenes(aHSMM, min_expr = 0.1)
print(head(fData(aHSMM)))
expressed_genes_all <- row.names(subset(fData(aHSMM), num_cells_expressed >= 10))

ctHSMM <- detectGenes(ctHSMM, min_expr = 0.1)
print(head(fData(ctHSMM)))
expressed_genes_all_ct <- row.names(subset(fData(ctHSMM), num_cells_expressed >= 10))

#for toatal mRNAs
pData(aHSMM)$Total_mRNAs <- Matrix::colSums(exprs(aHSMM))
pData(ctHSMM)$Total_mRNAs <- Matrix::colSums(exprs(ctHSMM))


aHSMM <- aHSMM[,pData(aHSMM)$Total_mRNAs < 1e6]
ctHSMM <- ctHSMM[,pData(ctHSMM)$Total_mRNAs < 1e6]

#this then for dipersion!
aHSMM <- estimateSizeFactors(aHSMM)
aHSMM <- estimateDispersions(aHSMM)

ctHSMM <- estimateSizeFactors(ctHSMM)
ctHSMM <- estimateDispersions(ctHSMM)

upper_bound <- 10^(mean(log10(pData(ctHSMM)$Total_mRNAs)) + 2*sd(log10(pData(ctHSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(ctHSMM)$Total_mRNAs)) - 2*sd(log10(pData(ctHSMM)$Total_mRNAs)))


pdf('total_mRNAs_celltype.pdf')
qplot(Total_mRNAs, data=pData(ctHSMM), color=cell_type,geom="density") +
geom_vline(xintercept=lower_bound) +
geom_vline(xintercept=upper_bound)
dev.off()

pdf('total_mRNAs_DAYS.pdf')
qplot(Total_mRNAs, data=pData(ctHSMM), color=DAYS,geom="density") +
geom_vline(xintercept=lower_bound) +
geom_vline(xintercept=upper_bound)
dev.off()

# filtering the data
ctHSMM <- ctHSMM[,pData(ctHSMM)$Total_mRNAs > lower_bound &
pData(ctHSMM)$Total_mRNAs < upper_bound]
ctHSMM <- detectGenes(ctHSMM, min_expr = 0.1)

# test of roughly lognormal
L <- log(exprs(ctHSMM[expressed_genes_all_ct,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

#this takes some time (; the large number of removed entries corresponds to the 0 imputation
pdf('FPKM_all_cells_merged.pdf')
qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') + xlab("Standardized log(FPKM)") + ylab("Density")
dev.off()


#dispersion analysis
disp_table <- dispersionTable(ctHSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ctHSMM <- setOrderingFilter(ctHSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(ctHSMM)
pdf('dispersion_all_samples.pdf')
plot_ordering_genes(ctHSMM)
dev.off()

pdf('PC_variance_explained_all_samples.pdf')
plot_pc_variance_explained(ctHSMM, return_all = F)
dev.off()

#dimension reduction - tsne only with monocle 2.4.0 otherwise DDRTree
#aHSMM <- reduceDimension(aHSMM, max_components=2, num_dim = 6, reduction_method = 'tSNE', verbose = T)

ctHSMM <- reduceDimension(ctHSMM, max_components=2, num_dim=6, reduction_method = 'tSNE', verbose = T,check_duplicates=FALSE)

ctHSMM <- clusterCells(ctHSMM,num_clusters=12)
#mHSMM <- clusterCells(mHSMM,num_clusters=10)

#mHSMM <- reduceDimension(mHSMM, max_components=2, num_dim = 6, reduction_method = 'DDRTree', verbose = T)

pdf('Cluster_tSNE_11clusters_all_samples.pdf')
plot_cell_clusters(ctHSMM, 1, 2)
dev.off()

pdf('Cluster_tSNE_11clusters_all_samples_by_cell_Line.pdf')
plot_cell_clusters(ctHSMM, 1, 2, color="cell_line")
dev.off()

pdf('Cluster_tSNE_8clusters_all_samples_by_day_new.pdf')
plot_cell_clusters(ctHSMM, 1, 2, color="days")
dev.off()

pdf('Cluster_tSNE_clusters_cell_line_by_day.pdf')
plot_cell_clusters(ctHSMM, 1, 2, color="cell_line") + facet_wrap(~days)
dev.off()

pdf('Cluster_tSNE_clusters_cell_line_by_day.pdf')
plot_cell_clusters(ctHSMM, 1, 2, color="cell_type") + facet_wrap(~days)
dev.off()

#to exclude batch effect by num_expressed_genes
ctHSMM <- reduceDimension(ctHSMM, max_components=2, num_dim = 6, reduction_method = 'tSNE', residualModelFormulaStr="~num_genes_expressed", verbose = T, check_duplicates=FALSE)

ctHSMM <- clusterCells(ctHSMM,num_clusters=12)

pdf('Cluster_tSNE_clusters_by_day_batcheffect_by_num_genes_removed_cell_type.pdf')
 plot_cell_clusters(ctHSMM, 1, 2, color="cell_type") + facet_wrap(~days)
dev.off()

pdf('Cluster_tSNE_clusters_by_day_batcheffect_by_num_genes_removed_cell_line.pdf')
plot_cell_clusters(ctHSMM, 1, 2, color="cell_line") + facet_wrap(~days)
dev.off()

pdf('Cluster_tSNE_clusters_by_cell_line_batcheffect_by_num_genes_removed_cell_type.pdf')
plot_cell_clusters(ctHSMM, 1, 2, color="cell_type") + facet_wrap(~cell_line)
dev.off()

pdf('Cluster_tSNE_clusters_batcheffect_by_num_genes_removed_cell_line.pdf')
plot_cell_clusters(ctHSMM, 1, 2, color="cell_line")
dev.off()

pdf('Cluster_tSNE_clusters_batcheffect_by_num_genes_removed_Days.pdf')
plot_cell_clusters(ctHSMM, 1, 2, color="DAYS")
dev.off()

#Pseudotimes
disp_table_all <- dispersionTable(ctHSMM)
ordering_genes_all <- subset(disp_table_all, mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit)$gene_id

ctHSMM_all <- setOrderingFilter(ctHSMM, ordering_genes_all)
pdf('dispersion_ALL_PT_new.pdf')
plot_ordering_genes(ctHSMM_all)
dev.off()

ctHSMM_all <- reduceDimension(ctHSMM_all, max_components=11) # before 3 comp
ctHSMM_all <- orderCells(ctHSMM_all)

pdf('cell_trajectory_col_day_ALL_days_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="days")
dev.off()

pdf('cell_trajectory_col_day_ALL_cell_type_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="cell_type")
dev.off()

pdf('cell_trajectory_col_day_ALL_cell_line_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="cell_line")
dev.off()

pdf('cell_trajectory_col_state_ALL_state_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="State")
dev.off()


pdf('cell_trajectory_col_pseudotime_ALL_new_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="Pseudotime")
dev.off()

plot_cell_trajectory(ctHSMM_all, color_by="State") + facet_wrap(~State, nrow=1)

pdf('cell_trajectory_col_day_separated_cell_line_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="days") + facet_wrap(~cell_line, nrow=1)
dev.off()

pdf('cell_trajectory_col_pseudotime_separated_cell_line_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="Pseudotime") + facet_wrap(~cell_line, nrow=1)
dev.off()


pdf('cell_trajectory_col_day_separated_cell_type_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="days") + facet_wrap(~cell_type, nrow=1)
dev.off()

pdf('cell_trajectory_col_pseudotime_separated_days_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="Pseudotime") + facet_wrap(~days, nrow=1)
dev.off()

pdf('cell_trajectory_col_cellline_separated_days_11comp.pdf')
plot_cell_trajectory(ctHSMM_all, color_by="cell_line") + facet_wrap(~days, nrow=1)
dev.off()


#another filtering - and to get things out of monocle structures!
nsiggenes <- as.numeric(row.names(subset(fData(ctHSMM_all), num_cells_expressed >= 15))) # gives geenes that are at least expressed in 15 cells
filt_express_ma <-(allcell[nsiggenes,colnames(allcell) %in% pData(ctHSMM_all)$cell_sample]) # gives expression matrix of siggenes and filtered

nsiggenes50 <- as.numeric(row.names(subset(fData(ctHSMM_all), num_cells_expressed >= 50))) # gives geenes that are at least expressed in 15 cells
filt_express_50 <-(allcell[nsiggenes50,colnames(allcell) %in% pData(ctHSMM_all)$cell_sample]) # gives expression matrix of siggenes and filtered


# to make a more hand-made comparison
ctHSMM_wt_forcomp <- ctHSMM_all[,pData(ctHSMM_all)$cell_type == "WT"]
ctHSMM_mu_forcomp <- ctHSMM_all[,pData(ctHSMM_all)$cell_type == "MUT"]


pdf("densities_of_pseudotimes.pdf")
qplot(Pseudotime, data=pData(ctHSMM_all), color=cell_type,geom="density")
dev.off()

#here for cell lines WT and MUT
ctHSMM_wt <- ctHSMM_all[,pData(ctHSMM_all)$cell_line == "WT"]
ctHSMM_mu <- ctHSMM_all[,pData(ctHSMM_all)$cell_line == "MUT"]

ctHSMM_d0 <- ctHSMM_all[,pData(ctHSMM_all)$DAYS == "0"]
ctHSMM_d10 <- ctHSMM_all[,pData(ctHSMM_all)$DAYS == "10"]
ctHSMM_d14 <- ctHSMM_all[,pData(ctHSMM_all)$DAYS == "14"]
ctHSMM_d42 <- ctHSMM_all[,pData(ctHSMM_all)$DAYS == "42"]

ctHSMM_d0_wt  <- ctHSMM_wt[,pData(ctHSMM_wt)$DAYS == "0"]
ctHSMM_d10_wt <- ctHSMM_wt[,pData(ctHSMM_wt)$DAYS == "10"]
ctHSMM_d14_wt <- ctHSMM_wt[,pData(ctHSMM_wt)$DAYS == "14"]
ctHSMM_d42_wt <- ctHSMM_wt[,pData(ctHSMM_wt)$DAYS == "42"]

ctHSMM_d0_mu  <- ctHSMM_mu[,pData(ctHSMM_mu)$DAYS == "0"]
ctHSMM_d10_mu <- ctHSMM_mu[,pData(ctHSMM_mu)$DAYS == "10"]
ctHSMM_d14_mu <- ctHSMM_mu[,pData(ctHSMM_mu)$DAYS == "14"]
ctHSMM_d42_mu <- ctHSMM_mu[,pData(ctHSMM_mu)$DAYS == "42"]

q1 <- qplot(Pseudotime, data=pData(ctHSMM_d0), color=cell_line,geom="density")
 q2 <- qplot(Pseudotime, data=pData(ctHSMM_d10), color=cell_line,geom="density")
q3 <- qplot(Pseudotime, data=pData(ctHSMM_d14), color=cell_line,geom="density")
q4 <- qplot(Pseudotime, data=pData(ctHSMM_d42), color=cell_line,geom="density")




plotas_wt_0 <- as.data.frame(pData(ctHSMM_d0_wt))
plotas_mu_0 <- as.data.frame(pData(ctHSMM_d0_mu))
plotas_wt_10 <- as.data.frame(pData(ctHSMM_d10_wt))
plotas_mu_10 <- as.data.frame(pData(ctHSMM_d10_mu))
plotas_wt_14 <- as.data.frame(pData(ctHSMM_d14_wt))
plotas_mu_14 <- as.data.frame(pData(ctHSMM_d14_mu))
plotas_wt_42 <- as.data.frame(pData(ctHSMM_d42_wt))
plotas_mu_42 <- as.data.frame(pData(ctHSMM_d42_mu))


pdf('./pseudotime_distributions_new2.pdf',height=4, width=8)
par(mfrow = c(4, 1))
ggplot(plotas_mu_0, aes(x=Pseudotime, group=cell_line)) + geom_density(alpha = 0.5, fill="red") + geom_density(data=plotas_wt_0,alpha = 0.5, fill="blue")+ xlim(0,25) + theme(aspect.ratio = .3)

ggplot(plotas_mu_10, aes(x=Pseudotime, group=cell_line)) + geom_density(alpha = 0.5, fill="red") + geom_density(data=plotas_wt_10,alpha = 0.5, fill="blue")+ xlim(0,25)+ ylim(0,1.5) + theme(aspect.ratio = .3)

ggplot(plotas_mu_14, aes(x=Pseudotime, group=cell_line)) + geom_density(alpha = 0.5, fill="red") + geom_density(data=plotas_wt_14,alpha = 0.5, fill="blue") + xlim(0,25)+ ylim(0,1.5)+ theme(aspect.ratio = .3)

ggplot(plotas_mu_42, aes(x=Pseudotime, group=cell_line)) + geom_density(alpha = 0.5, fill="red") + geom_density(data=plotas_wt_42,alpha = 0.5, fill="blue") + xlim(0,25)+ theme(aspect.ratio = .3)

dev.off()

ggplot(plotas_mu_0, aes(x=Pseudotime, group=cell_line)) + geom_density(alpha = 0.5, fill="red") + geom_density(data=plotas_wt_0,alpha = 0.5, fill="blue")+
geom_density(data=plotas_mu_10,alpha = 0.5, fill="red") + geom_density(data=plotas_wt_10,alpha = 0.5, fill="blue") +xlim(0,25) + theme(aspect.ratio = 1)


pdf('./pseudotime_distributions.pdf',height=12, width=6)
#colordefs: wtcodef[1],length(colnames(wt_0))), rep(mucodef[1]
par(mfrow = c(4, 1))
 plot(density(pData(ctHSMM_d0_wt)$Pseudotime),col=rgb(0,0,1,0.2),main="Day 0", xlab="",xlim=c(0,25)) +#xlim=c(0.05,.2),ylim=c(0.,120)
#line(density(tocoav_p3na))
#line(density(tocoav_p8na))
polygon(density(pData(ctHSMM_d0_wt)$Pseudotime),density=-1,col=rgb(0,0,1,0.2)) +
polygon(density(pData(ctHSMM_d0_mu)$Pseudotime),density=-1,col=rgb(1,0,0,0.2)) +
lines(c(0,25),c(0,0)) +
lines(density(pData(ctHSMM_d0_mu)$Pseudotime),col=rgb(1,0,0,0.2))

plot(density(pData(ctHSMM_d10_mu)$Pseudotime),col=rgb(1,0,0,0.2),main="Day 10", xlab="",xlim=c(0,25)) #xlim=c(0.05,.2),ylim=c(0.,120)
#line(density(tocoav_p3na))
#line(density(tocoav_p8na))
polygon(density(pData(ctHSMM_d10_wt)$Pseudotime),density=-1,col=rgb(0,0,1,0.2))
polygon(density(pData(ctHSMM_d10_mu)$Pseudotime),density=-1,col=rgb(1,0,0,0.2))
lines(c(0,25),c(0,0))
lines(density(pData(ctHSMM_d10_wt)$Pseudotime),col=rgb(0,0,1,0.2))

plot(density(pData(ctHSMM_d14_mu)$Pseudotime),col=rgb(1,0,0,0.2),main="Day 14", xlab="",xlim=c(0,25)) #xlim=c(0.05,.2),ylim=c(0.,120)
#line(density(tocoav_p3na))
#line(density(tocoav_p8na))
polygon(density(pData(ctHSMM_d14_wt)$Pseudotime),density=-1,col=rgb(0,0,1,0.2))
polygon(density(pData(ctHSMM_d14_mu)$Pseudotime),density=-1,col=rgb(1,0,0,0.2))
lines(c(0,25),c(0,0))
lines(density(pData(ctHSMM_d14_wt)$Pseudotime),col=rgb(0,0,1,0.2))

plot(density(pData(ctHSMM_d42_mu)$Pseudotime),col=rgb(1,0,0,0.2),main="Day 42", xlab="Pseudotime",xlim=c(0,25)) #xlim=c(0.05,.2),ylim=c(0.,120)
#line(density(tocoav_p3na))
#line(density(tocoav_p8na))
polygon(density(pData(ctHSMM_d42_wt)$Pseudotime),density=-1,col=rgb(0,0,1,0.2))
polygon(density(pData(ctHSMM_d42_mu)$Pseudotime),density=-1,col=rgb(1,0,0,0.2))
lines(c(0,25),c(0,0))
lines(density(pData(ctHSMM_d42_wt)$Pseudotime),col=rgb(0,0,1,0.2))
dev.off()

polygon(density(pData(ctHSMM_d10_wt)$Pseudotime),density=-1,col=rgb(1,0,0,0.2))
polygon(density(pData(ctHSMM_d10_mu)$Pseudotime),density=-1,col=rgb(0,0,1,0.2))




### for Index approach

#SHOULD THE DATA BE SPLITED WRT TO CELL TYPE?

#number of cells that need to express a sig gene
ncells=15
fisiggenes <- as.numeric(row.names(subset(fData(ctHSMM_all), num_cells_expressed >= ncells)))

mu_0_na <- mu_0[fisiggenes,]
is.na(mu_0_na) <- !mu_0_na

mu_10_na <- mu_10[fisiggenes,]
is.na(mu_10_na) <- !mu_10_na

mu_14_na <- mu_14[fisiggenes,]
is.na(mu_14_na) <- !mu_14_na

mu_42_na <- mu_42[fisiggenes,]
is.na(mu_42_na) <- !mu_42_na


wt_0_na <-wt_0[fisiggenes,]
is.na(wt_0_na) <- !wt_0_na

wt_10_na <-wt_10[fisiggenes,]
is.na(wt_10_na) <- !wt_10_na

wt_14_na <-wt_14[fisiggenes,]
is.na(wt_14_na) <- !wt_14_na

wt_42_na <-wt_42[fisiggenes,]
is.na(wt_42_na) <- !wt_42_na

#pearson cell-cell correaltion: HERE DONE FOR mu_0_na ... for Index approach to be done for all days and cell lines - mainly copy paste (;

alpha <- 0.01 # correlation significance check
### MOD BY STE FROM HERE ###
############################################################
#NOW WITH PEARSON CORRELATION
############################################################
##### DAY 0 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="pearson")
cyx_mu_0_na <- rcorr(as.matrix(mu_0_na),type="pearson") #calculates pearson correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_mu_0_na$P < alpha # this checks for significant level
csigcor_mu_0_na <- cyx_mu_0_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_mu_0_na  <- mean(csigcor_mu_0_na,na.rm=T)

#SD of cell cell correlation for NA sets
###sd_corr_mu_0_na  <- sd(corr_mu_0_na,na.rm=T)
sd_corr_mu_0_na  <- sd(csigcor_mu_0_na,na.rm=T)

##### DAY 10 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="pearson")
cyx_mu_10_na <- rcorr(as.matrix(mu_10_na),type="pearson") #calculates pearson correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_mu_10_na$P < alpha # this checks for significant level
csigcor_mu_10_na <- cyx_mu_10_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_mu_10_na  <- mean(csigcor_mu_10_na,na.rm=T)

#SD of cell cell correlation for NA sets
sd_corr_mu_10_na  <- sd(csigcor_mu_10_na,na.rm=T)

##### DAY 14 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="pearson")
cyx_mu_14_na <- rcorr(as.matrix(mu_14_na),type="pearson") #calculates pearson correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_mu_14_na$P < alpha # this checks for significant level
csigcor_mu_14_na <- cyx_mu_14_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_mu_14_na  <- mean(csigcor_mu_14_na,na.rm=T)

#SD of cell cell correlation for NA sets
sd_corr_mu_14_na  <- sd(csigcor_mu_14_na,na.rm=T)

##### DAY 42 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="pearson")
cyx_mu_42_na <- rcorr(as.matrix(mu_42_na),type="pearson") #calculates pearson correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_mu_42_na$P < alpha # this checks for significant level
csigcor_mu_42_na <- cyx_mu_42_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_mu_42_na  <- mean(csigcor_mu_42_na,na.rm=T)

#SD of cell cell correlation for NA sets
sd_corr_mu_42_na  <- sd(csigcor_mu_42_na,na.rm=T)
####################
Ave_CC_Cor_mu_ALL <- c(ave_corr_mu_0_na,ave_corr_mu_10_na,ave_corr_mu_14_na,ave_corr_mu_42_na)
SD_CC_mu_ALL <- c(sd_corr_mu_0_na, sd_corr_mu_10_na, sd_corr_mu_14_na, sd_corr_mu_42_na)
plot(Ave_CC_Cor_mu_ALL)
####################

##### DAY 0 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="pearson")
cyx_wt_0_na <- rcorr(as.matrix(wt_0_na),type="pearson") #calculates pearson correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_wt_0_na$P < alpha # this checks for significant level
csigcor_wt_0_na <- cyx_wt_0_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_wt_0_na  <- mean(csigcor_wt_0_na,na.rm=T)

#SD of cell cell correlation for NA sets
sd_corr_mu_0_na  <- sd(corr_mu_0_na,na.rm=T)

##### DAY 10 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="pearson")
cyx_wt_10_na <- rcorr(as.matrix(wt_10_na),type="pearson") #calculates pearson correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_wt_10_na$P < alpha # this checks for significant level
csigcor_wt_10_na <- cyx_wt_10_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_wt_10_na  <- mean(csigcor_wt_10_na,na.rm=T)

#SD of cell cell correlation for NA sets
###sd_corr_mu_10_na  <- sd(corr_mu_10_na,na.rm=T)

##### DAY 14 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="pearson")
cyx_wt_14_na <- rcorr(as.matrix(wt_14_na),type="pearson") #calculates pearson correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_wt_14_na$P < alpha # this checks for significant level
csigcor_wt_14_na <- cyx_wt_14_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_wt_14_na  <- mean(csigcor_wt_14_na,na.rm=T)

#SD of cell cell correlation for NA sets
###sd_corr_mu_10_na  <- sd(corr_mu_10_na,na.rm=T)

##### DAY 42 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="pearson")
cyx_wt_42_na <- rcorr(as.matrix(wt_42_na),type="pearson") #calculates pearson correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_wt_42_na$P < alpha # this checks for significant level
csigcor_wt_42_na <- cyx_wt_42_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_wt_42_na  <- mean(csigcor_wt_42_na,na.rm=T)

#SD of cell cell correlation for NA sets
#sd_corr_mu_42_na  <- sd(corr_mu_42_na,na.rm=T)
####################
Ave_CC_Cor_wt_ALL <- c(ave_corr_wt_0_na,ave_corr_wt_10_na,ave_corr_wt_14_na,ave_corr_wt_42_na)
plot(Ave_CC_Cor_wt_ALL)
###################
plot(c(0,10,14,42), Ave_CC_Cor_wt_ALL, xlab = "", ylab = "", type="o",  pch=19)
par(new = TRUE)
plot(c(0,10,14,42), Ave_CC_Cor_mu_ALL, axes = FALSE, col="red", type="o",  pch=19, xlab="time (days)", ylab="Average Cell-Cell Correlation")
####################
############################################################
#NOW WITH SPEARMAN
############################################################
############################################################
##### DAY 0 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="spearman")
cyx_mu_0_na <- rcorr(as.matrix(mu_0_na),type="spearman") #calculates spearman correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_mu_0_na$P < alpha # this checks for significant level
csigcor_mu_0_na <- cyx_mu_0_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_mu_0_na  <- mean(csigcor_mu_0_na,na.rm=T)

#SD of cell cell correlation for NA sets
###sd_corr_mu_0_na  <- sd(corr_mu_0_na,na.rm=T)

##### DAY 10 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="spearman")
cyx_mu_10_na <- rcorr(as.matrix(mu_10_na),type="spearman") #calculates spearman correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_mu_10_na$P < alpha # this checks for significant level
csigcor_mu_10_na <- cyx_mu_10_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_mu_10_na  <- mean(csigcor_mu_10_na,na.rm=T)

#SD of cell cell correlation for NA sets
###sd_corr_mu_10_na  <- sd(corr_mu_10_na,na.rm=T)

##### DAY 14 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="spearman")
cyx_mu_14_na <- rcorr(as.matrix(mu_14_na),type="spearman") #calculates spearman correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_mu_14_na$P < alpha # this checks for significant level
csigcor_mu_14_na <- cyx_mu_14_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_mu_14_na  <- mean(csigcor_mu_14_na,na.rm=T)

#SD of cell cell correlation for NA sets
###sd_corr_mu_10_na  <- sd(corr_mu_10_na,na.rm=T)

##### DAY 42 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="spearman")
cyx_mu_42_na <- rcorr(as.matrix(mu_42_na),type="spearman") #calculates spearman correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_mu_42_na$P < alpha # this checks for significant level
csigcor_mu_42_na <- cyx_mu_42_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_mu_42_na  <- mean(csigcor_mu_42_na,na.rm=T)

#SD of cell cell correlation for NA sets
#sd_corr_mu_42_na  <- sd(corr_mu_42_na,na.rm=T)
####################
Ave_CC_Cor_mu_ALL_Spearman <- c(ave_corr_mu_0_na,ave_corr_mu_10_na,ave_corr_mu_14_na,ave_corr_mu_42_na)
plot(Ave_CC_Cor_mu_ALL_Spearman)
####################

##### DAY 0 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="spearman")
cyx_wt_0_na <- rcorr(as.matrix(wt_0_na),type="spearman") #calculates spearman correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_wt_0_na$P < alpha # this checks for significant level
csigcor_wt_0_na <- cyx_wt_0_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_wt_0_na  <- mean(csigcor_wt_0_na,na.rm=T)

#SD of cell cell correlation for NA sets
###sd_corr_mu_0_na  <- sd(corr_mu_0_na,na.rm=T)

##### DAY 10 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="spearman")
cyx_wt_10_na <- rcorr(as.matrix(wt_10_na),type="spearman") #calculates spearman correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_wt_10_na$P < alpha # this checks for significant level
csigcor_wt_10_na <- cyx_wt_10_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_wt_10_na  <- mean(csigcor_wt_10_na,na.rm=T)

#SD of cell cell correlation for NA sets
###sd_corr_mu_10_na  <- sd(corr_mu_10_na,na.rm=T)

##### DAY 14 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="spearman")
cyx_wt_14_na <- rcorr(as.matrix(wt_14_na),type="spearman") #calculates spearman correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_wt_14_na$P < alpha # this checks for significant level
csigcor_wt_14_na <- cyx_wt_14_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_wt_14_na  <- mean(csigcor_wt_14_na,na.rm=T)

#SD of cell cell correlation for NA sets
###sd_corr_mu_10_na  <- sd(corr_mu_10_na,na.rm=T)

##### DAY 42 #####
#cyx_mu_0 <-rcorr(as.matrix(mu_0),type="spearman")
cyx_wt_42_na <- rcorr(as.matrix(wt_42_na),type="spearman") #calculates spearman correlation between cells and omits NA

#to filter for correlations with significance smaller than alpha only
#filta <-cyx_mu_0$P < alpha
#csigcor_mu_0 <- cyx_mu_0$r*filta

filta <-cyx_wt_42_na$P < alpha # this checks for significant level
csigcor_wt_42_na <- cyx_wt_42_na$r*filta # and this filters the correlation entries for significance alpha

#average over cell cell correlation for NA sets
ave_corr_wt_42_na  <- mean(csigcor_wt_42_na,na.rm=T)

#SD of cell cell correlation for NA sets
#sd_corr_mu_42_na  <- sd(corr_mu_42_na,na.rm=T)
####################
Ave_CC_Cor_wt_ALL_Spearman <- c(ave_corr_wt_0_na,ave_corr_wt_10_na,ave_corr_wt_14_na,ave_corr_wt_42_na)
plot(Ave_CC_Cor_wt_ALL_Spearman)
###################
plot(c(0,10,14,42), Ave_CC_Cor_wt_ALL_Spearman, xlab = "", ylab = "", type="o",  pch=19)
par(new = TRUE)
plot(c(0,10,14,42), Ave_CC_Cor_mu_ALL_Spearman, axes = FALSE, col="red", type="o",  pch=19, xlab="time (days)", ylab="Average Cell-Cell Correlation")
title(main="Spearman") 
####################
#pearson gene-gene correaltion: HERE DONE FOR mu_0_na ... for Index approach to be done for all days and cell lines - mainly copy paste (;

# THIS CAN TAKE A LONG TIME DUE TO LARGE GENE NUMBERS ... so bootstrapping one approach
ns <- 2000 #sampling size of genes
nboot <- 200 #number of bootstrapping steps
alpha <- 0.01 # correlation significance check

#######################################################################
##### DAY0 MU ######
nma <- dim(mu_0_na)[1] #largest row number of matrix
tocoav_mu_0_na <- 0
SDtocoav_mu_0_na <- 0
#this over loop
for (i in c(1:nboot)){
    rangen <- sample(1:nma, ns, replace=TRUE) #random genes
    
    gyx_mu_0_na<-rcorr(as.matrix(t(mu_0_na[rangen,])),type="pearson")
    filta<-gyx_mu_0_na$P<alpha
    sigcor_mu_0_na<-gyx_mu_0_na$r*filta
    
    tocoav_mu_0_na[i] <- mean(abs(sigcor_mu_0_na),na.rm=T)
    #SDtocoav_mu_0_na[i] <- sd(abs(sigcor_mu_0_na),na.rm=T)
}
hist(tocoav_mu_0_na,breaks=20)


##### DAY10 MU ######

nma <- dim(mu_10_na)[1] #largest row number of matrix
tocoav_mu_10_na <- 0
SDtocoav_mu_10_na <- 0
#this over loop
for (i in c(1:nboot)){
  rangen <- sample(1:nma, ns, replace=TRUE) #random genes
  
  gyx_mu_10_na<-rcorr(as.matrix(t(mu_10_na[rangen,])),type="pearson")
  filta<-gyx_mu_10_na$P<alpha
  sigcor_mu_10_na<-gyx_mu_10_na$r*filta
  
  tocoav_mu_10_na[i] <- mean(abs(sigcor_mu_10_na),na.rm=T)
  #SDtocoav_mu_10_na[i] <- sd(abs(sigcor_mu_10_na),na.rm=T)
}
hist(tocoav_mu_10_na,breaks=20)

##### DAY14 MU ######

nma <- dim(mu_14_na)[1] #largest row number of matrix
tocoav_mu_14_na <- 0
SDtocoav_mu_14_na <- 0
#this over loop
for (i in c(1:nboot)){
  rangen <- sample(1:nma, ns, replace=TRUE) #random genes
  
  gyx_mu_14_na<-rcorr(as.matrix(t(mu_14_na[rangen,])),type="pearson")
  filta<-gyx_mu_14_na$P<alpha
  sigcor_mu_14_na<-gyx_mu_14_na$r*filta
  
  tocoav_mu_14_na[i] <- mean(abs(sigcor_mu_14_na),na.rm=T)
  #SDtocoav_mu_14_na[i] <- sd(abs(sigcor_mu_14_na),na.rm=T)
}
hist(tocoav_mu_14_na,breaks=20)

##### DAY42 MU ######

nma <- dim(mu_42_na)[1] #largest row number of matrix
tocoav_mu_42_na <- 0
SDtocoav_mu_42_na <- 0
#this over loop
for (i in c(1:nboot)){
  rangen <- sample(1:nma, ns, replace=TRUE) #random genes
  
  gyx_mu_42_na<-rcorr(as.matrix(t(mu_42_na[rangen,])),type="pearson")
  filta<-gyx_mu_42_na$P<alpha
  sigcor_mu_42_na<-gyx_mu_42_na$r*filta
  
  tocoav_mu_42_na[i] <- mean(abs(sigcor_mu_42_na),na.rm=T)
  #SDtocoav_mu_42_na[i] <- sd(abs(sigcor_mu_42_na),na.rm=T)
}
hist(tocoav_mu_42_na,breaks=20)

####################################################################

##### DAY0 MU ######
nma <- dim(wt_0_na)[1] #largest row number of matrix
tocoav_wt_0_na <- 0
SDtocoav_wt_0_na <- 0
#this over loop
for (i in c(1:nboot)){
  rangen <- sample(1:nma, ns, replace=TRUE) #random genes
  
  gyx_wt_0_na<-rcorr(as.matrix(t(wt_0_na[rangen,])),type="pearson")
  filta<-gyx_wt_0_na$P<alpha
  sigcor_wt_0_na<-gyx_wt_0_na$r*filta
  
  tocoav_wt_0_na[i] <- mean(abs(sigcor_wt_0_na),na.rm=T)
  #SDtocoav_mu_0_na[i] <- sd(abs(sigcor_mu_0_na),na.rm=T)
}
hist(tocoav_wt_0_na,breaks=20)

##### DAY10 MU ######

nma <- dim(wt_10_na)[1] #largest row number of matrix
tocoav_wt_10_na <- 0
SDtocoav_wt_10_na <- 0
#this over loop
for (i in c(1:nboot)){
  rangen <- sample(1:nma, ns, replace=TRUE) #random genes
  
  gyx_wt_10_na<-rcorr(as.matrix(t(wt_10_na[rangen,])),type="pearson")
  filta<-gyx_wt_10_na$P<alpha
  sigcor_wt_10_na<-gyx_wt_10_na$r*filta
  
  tocoav_wt_10_na[i] <- mean(abs(sigcor_wt_10_na),na.rm=T)
  #SDtocoav_mu_10_na[i] <- sd(abs(sigcor_mu_10_na),na.rm=T)
}
hist(tocoav_wt_10_na,breaks=20)

##### DAY14 MU ######

nma <- dim(wt_14_na)[1] #largest row number of matrix
tocoav_wt_14_na <- 0
SDtocoav_wt_14_na <- 0
#this over loop
for (i in c(1:nboot)){
  rangen <- sample(1:nma, ns, replace=TRUE) #random genes
  
  gyx_wt_14_na<-rcorr(as.matrix(t(wt_14_na[rangen,])),type="pearson")
  filta<-gyx_wt_14_na$P<alpha
  sigcor_wt_14_na<-gyx_wt_14_na$r*filta
  
  tocoav_wt_14_na[i] <- mean(abs(sigcor_wt_14_na),na.rm=T)
  #SDtocoav_mu_14_na[i] <- sd(abs(sigcor_mu_14_na),na.rm=T)
}
hist(tocoav_wt_14_na,breaks=20)

##### DAY42 MU ######

nma <- dim(wt_42_na)[1] #largest row number of matrix
tocoav_wt_42_na <- 0
SDtocoav_wt_42_na <- 0
#this over loop
for (i in c(1:nboot)){
  rangen <- sample(1:nma, ns, replace=TRUE) #random genes
  
  gyx_wt_42_na<-rcorr(as.matrix(t(wt_42_na[rangen,])),type="pearson")
  filta<-gyx_wt_42_na$P<alpha
  sigcor_wt_42_na<-gyx_wt_42_na$r*filta
  
  tocoav_wt_42_na[i] <- mean(abs(sigcor_wt_42_na),na.rm=T)
  #SDtocoav_mu_42_na[i] <- sd(abs(sigcor_mu_42_na),na.rm=T)
}
hist(tocoav_wt_42_na,breaks=20)
####################################################################
MeansOfGGMeanCorrWTall <- c(mean(tocoav_wt_0_na), mean(tocoav_wt_10_na), mean(tocoav_wt_14_na), mean(tocoav_wt_42_na))
MeansOfGGMeanCorrMUall <- c(mean(tocoav_mu_0_na), mean(tocoav_mu_10_na), mean(tocoav_mu_14_na), mean(tocoav_mu_42_na))
plot(c(0,10,14,42), MeansOfGGMeanCorrWTall, xlab = "", ylab = "", type="o",  pch=19)
par(new = TRUE)
plot(c(0,10,14,42), MeansOfGGMeanCorrMUall, axes = FALSE, col="red", type="o",  pch=19, xlab="time (days)", ylab="Average Gene-Gene Correlation")
title(main="") 
####################################################################
CritTransINDEX_WTall <- MeansOfGGMeanCorrWTall / Ave_CC_Cor_wt_ALL
CritTransINDEX_MUall <- MeansOfGGMeanCorrMUall / Ave_CC_Cor_mu_ALL
plot(c(0,10,14,42), CritTransINDEX_WTall, xlab = "", ylab = "", type="o",  pch=19)
par(new = TRUE)
plot(c(0,10,14,42), CritTransINDEX_MUall, axes = FALSE, col="red", type="o",  pch=19, xlab="time (days)", ylab="Critical Transition Index")
title(main="") 
####################################################################
# FOR SPEARMAN correlation simply exchange type="pearson" with type="spearman" in rcorr function - again copy paste


#THEN AVERAGING OVER CORRELATIONS AND CALCULATING INDEX


# FOR GETTING THINGS FROM THE BRANCHING (MONOCLE ANALYSIS) CHECK ON

head(pData(ctHSMM_all))

# FROM THE COLUMN NAMES U CAN SEE WHAT IS DEFINED AND GET THIS OUT BY e.g.

pseuti <- pData(ctHSMM_all)$Pseudotime #OR

plot(pData(ctHSMM_all)$Pseudotime)

# sort the pseudotimes
pseuti_ordered <- sort(pseuti)
plot(pseuti_ordered)

#OR ONLY for WT
pseuti_wt <- sort(ctHSMM_wt$Pseudotime)

#OR ONLY for MUT
pseuti_mu <- sort(pData(ctHSMM_mu)$Pseudotime)

plot(pseuti_wt)
points(pseuti_mu,col="red")

# DOES THE PSEUDOTIME AS BIAS FOR DEVELOPMENTAL STATE CORRELATES WITH ENTROPY OF THE CELL?

# HOW CAN THIS BE "NORMALIZED" WRT TO DIFFERENT CELL NUMBERS FOR MUTANT AND WT TO COMPARE THE DYNAMICS?

# CAN U INFER THE ROUGHNESS OF THE LANDSCAPE?




#For Violin plots ... summarizing distributions see also these here:

#Violin plots
library('vioplot')
library('ggplot2')
library('randomcoloR')
library('infotheo')
library('vioplot')
library('caroline')# for other violineplots
library('gplots')

#and a !!!NOT WORKING!!! template to plot these - taken from another data analysis script!!!
#for making violins work
vp3 <- !is.na(p3na_exp[vargene_p3[126],])
vp3li <-(as.list(p3na_exp[vargene_p3[126],vp3]))
violins(as.vector(unlist(vp3li)))
violins(as.vector(unlist(p3na_exp[vargene_p3[126],vp3])))

#this is compact version for violies
violins(as.vector(unlist(p3na[vargene_p3_exp[126],!is.na(p3na_exp[vargene_p3[126],])])))



pdf("overview_Geneexpression_p3_100.pdf",height=10	,width=10)
par(mfrow=c(5,2))
par(las=0,mar=c(5,4,2,2),par(yaxp = c(0, 100, 4)))
violins(genelist[c(1:10)],connect=F,drawRect=F,horizontal=F,col=cop3[c(1:10)])
violins(genelist[c(11:20)],connect=F,drawRect=F,horizontal=F,col=cop3[c(11:20)])
violins(genelist[c(21:30)],connect=F,drawRect=F,horizontal=F,col=cop3[c(21:30)])
violins(genelist[c(31:40)],connect=F,drawRect=F,horizontal=F,col=cop3[c(31:40)])
violins(genelist[c(41:50)],connect=F,drawRect=F,horizontal=F,col=cop3[c(41:50)])
violins(genelist[c(51:60)],connect=F,drawRect=F,horizontal=F,col=cop3[c(51:60)])
violins(genelist[c(61:70)],connect=F,drawRect=F,horizontal=F,col=cop3[c(61:70)])
violins(genelist[c(71:80)],connect=F,drawRect=F,horizontal=F,col=cop3[c(71:80)])
violins(genelist[c(81:90)],connect=F,drawRect=F,horizontal=F,col=cop3[c(81:90)])
violins(genelist[c(91:100)],connect=F,drawRect=F,horizontal=F,col=cop3[c(91:100)])
dev.off()

