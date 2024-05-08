# The script can be run with this command.
# Windows: to run the R-script write in the R Console window: source ("EdgeR-GLM_Bi3013.r")                                 #


library(edgeR);                               # require that edgeR is installed

setwd("C:/TempBi3013");                       # working directory


# Generalized linear models (GLM)
# "GLM likelihood ratio test is recommended for experiments with multiple factors"
# Make sure you have created folders to store and analyse the data

ROT <- ("C:/TempBi3013");
setwd (ROT);

IMG <- ("C:/TempBi3013/IMG");                 # folder for images

# folders for counttab and results
CountTAB <- ("C:/TempBi3013/CountTAB");       # Count table is placed here 
Results <- ("C:/TempBi3013/Results");         # output of DE genes result tables

# Describe the targets
targets <- read.delim("Samples.txt", row.names="Samples"); # design file, define what samples you are analysing (in ROT directory)

# read Count table                                   # contain the information of sequencing reads that maps to the individual genes

setwd(CountTAB);

x <- read.delim("CountTable_Bi3013.txt", row.names="Phatr3_ID");      # the Count tab is loaded into x

# Filtering

countsPerMillion <- cpm(x);                                # Calculate counts pr. million reads
countCheck <- countsPerMillion > 1                         # TRUE FALSE matrix, remove genes with to low expression 
keep <- which(rowSums(countCheck) >= 2);                   # Must have at least 2 samples above threshold (~ 20 reads)
List_OkExp <- x[keep,];                                    # List where those genes with lowest expression are removed
summary(cpm(List_OkExp));                                  # Summary of the samples


Group <- factor(paste(targets$group));                     # make group used for design matrix

y <- DGEList(counts=List_OkExp, group=Group);              # y contain DGEList data, here using only genes above threshold

y <- calcNormFactors(y);                                   # Trimmed Mean of M-values

design <- model.matrix(~0+Group);                          # here it is a simple design, P-deplete versus P-replete

colnames(design) <- levels(Group);


y <- estimateGLMCommonDisp(y,design);                      # mean dispersion across all genes

y <- estimateGLMTrendedDisp(y,design);                     # mean dispersion across all genes with similar abundance

y <- estimateGLMTagwiseDisp(y,design);                     # calculates gene-specific dispersions

# fit <- glmFit(y, design);                                # a likelihood ratio test to which produces an object of class DGEGLM
                                                           # uses chisquare approximation to the likelihood ratio statistic

fit <- glmQLFit(y, design);                                # fits the negative binomial GLM for each tag and produces an object of class DGEGLM
                                                           # uses "quasi-likelihood" F-test in the likelihood ratio statistic

my.contrasts <- makeContrasts(WTPd_vs_WTPr=WT_Pd-WT_Pr, levels=design);


setwd(Results);

qlfWTPd_vs_WTPr <- glmQLFTest(fit, contrast=my.contrasts[,"WTPd_vs_WTPr"]);
TopScore <- topTags(qlfWTPd_vs_WTPr, n=nrow(x));                                # include all genes
write.table(TopScore, file="qlfWTPd_vs_WTPr.txt", sep="\t", quote=FALSE);       # write out all genes

keep <- TopScore$table$FDR <=0.0001 & abs(TopScore$table$logFC) >=2;            # only write out genes which has FDR <= 0.0001 and abs(log2)= +1
write.table(TopScore[keep,], file="qlfWTPd_vs_WTPr_FDR.0.0001_Log2=2.txt", sep="\t", quote=FALSE);


# some plotting

setwd (IMG);
deGenes <- decideTestsDGE(qlfWTPd_vs_WTPr, p=0.0001);
deGenes <- rownames(qlfWTPd_vs_WTPr)[as.logical(deGenes)];
plotSmear(qlfWTPd_vs_WTPr, de.tags=deGenes);
abline(h=c(-2, 2), col=2)
savePlot(filename = "logFC avg-logCPM", type = c("png"), device = dev.cur(), restoreConsole = TRUE);  # does not work for Mac


# Plot genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
# The BCV is the square root of the negative binomial dispersion. Displays the common, trended and tagwise BCV estimates.

# setwd (IMG);
plotBCV(y, xlab="Average log CPM", ylab="Biological coefficient of variation");
savePlot(filename = "BCVplot", type = c("png"), device = dev.cur(), restoreConsole = TRUE);


setwd (ROT);
sessionInfo();           # give info about installed packages
