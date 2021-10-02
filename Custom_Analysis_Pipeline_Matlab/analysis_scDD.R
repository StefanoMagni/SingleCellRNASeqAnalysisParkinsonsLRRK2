
library(SummarizedExperiment)
library(scde)
thres_zero <- 0.5

# rm(list=ls())
mydataA <- read.table('./jonas/A_14.txt',header = TRUE, row.names = 1)
mydataB <- read.table('./jonas/B_14.txt',header = TRUE, row.names = 1)

##########

num_sampA <- dim(mydataA)[2]
num_sampB <- dim(mydataB)[2]
condition <- c(rep(1, num_sampA), rep(2, num_sampB))

genes_inter <- intersect(rownames(mydataA),rownames(mydataB))
mydataA_inter <- mydataA[genes_inter,]
mydataB_inter <- mydataB[genes_inter,]
rm(mydataA)
rm(mydataB)

mydata <- cbind(mydataA_inter,mydataB_inter)
val_dup <- duplicated(rownames(mydata))
ind_true <- which(val_dup, TRUE)

# It seems that there were repeated colnames that were caussing problems
colnames(mydata) <- names(condition) <- paste0("Sample",1:ncol(mydata), sep="")

# We need a matrix to create the correct object
mydata_mat <- data.matrix(mydata)

mydata_obj <- SummarizedExperiment(assays=list("NormCounts"=mydata_mat),colData=data.frame(condition))

condition.names <- names(mydata_obj)

mydata_preproc <- preprocess(mydata_obj, ConditionNames=condition.names, zero.thresh = thres_zero)

###########

num_sampA <- dim(mydataA)[2]
num_sampB <- dim(mydataB)[2]

conditionA <- c(rep(1, num_sampA))
conditionB <- c(rep(2, num_sampB))
condition <- c(rep(1, num_sampA), rep(2, num_sampB))

colnames(mydataA) <- names(conditionA) <- paste0("C1.Sample",1:ncol(mydataA), sep="")
colnames(mydataB) <- names(conditionB) <- paste0("C2.Sample",1:ncol(mydataB), sep="")

mydata_matA <- data.matrix(mydataA_inter)
mydata_matB <- data.matrix(mydataB_inter)
mydata_mat_list <- c(mydata_matA, mydata_matB)
mydata_obj_list <- SummarizedExperiment(assays=list("NormCounts"=mydata_mat_list),colData=data.frame(conditionB))
mydata_obj_list <- SummarizedExperiment(assays=list("C1"=mydata_matA, "C2"=mydata_matB), colData=data.frame(conditionA,conditionB))
mydata_obj_list <- SummarizedExperiment(assays=list("C1"=mydataA, "C2"=mydataB), colData=data.frame(condition))

mydata_obj_list <- c(mydata_objA, mydata_objB)

condition.names <- c('C1', 'C2')
mydata_preproc <- preprocess(mydata_obj_list, ConditionNames=condition.names, zero.thresh = thres_zero)


#############

min_read_det = 4
mydataA_cd <- clean.counts(mydataA, min.lib.size=1000, min.reads = min_read_det, min.detected = min_read_det)
mydataB_cd <- clean.counts(mydataB, min.lib.size=1000, min.reads = min_read_det, min.detected = min_read_det)
genes_inter <- intersect(rownames(mydataA_cd),rownames(mydataB_cd))
mydataA_inter <- mydataA_cd[genes_inter,]
mydataB_inter <- mydataB_cd[genes_inter,]
mydata <- cbind(mydataA_inter,mydataB_inter)

num_sampA <- dim(mydataA_cd)[2]
num_sampB <- dim(mydataB_cd)[2]
condition <- factor(c(rep(1, num_sampA), rep(2, num_sampB)))

mydata_sub <- mydata[1:10000,]
o.ifm <- scde.error.models(counts = mydata_sub, groups = condition, threshold.segmentation = TRUE, n.cores = 1)

names(sg) <- colnames(cd)
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 10, threshold.segmentation = TRUE)
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
# make sure groups corresponds to the models (o.ifm)
groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels = c("ESC", "MEF"))
names(groups) <- row.names(o.ifm)
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100)
