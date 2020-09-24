rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
setwd("../")

suppressPackageStartupMessages(library(sgPLS))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(pheatmap))

source("Scripts/pls_functions.R")

X_pooled=readRDS("Data/Proteins_denoised.rds")
covars=readRDS("Data/Covariates.rds")
Y_pooled=covars$type

MyPLSDA_pooled <- plsda(X_pooled, Y_pooled, ncomp = 1,mode = "regression")

MyPLSDA_pooled$loadings$X
MyPLSDA_pooled$loadings$Y
MyPLSDA_pooled$explained_variance

MysPLSDA_pooled <- splsda(X_pooled, Y_pooled, ncomp=1, mode="regression", keepX=5)

MysPLSDA_pooled$loadings$X
MysPLSDA_pooled$loadings$X[MysPLSDA_pooled$loadings$X!=0,]

Groups=c("EGF", "FGF2", "GCSF", "VEGF", "GMSCF", "TGFa",
         "Eotaxin", "Fractalkine", "GRO", "MCP1", "MCP3", "MDC", "MIP1a", "MIP1b", "IP10", "IL8", 
         "IL1b", "IL2", "IL4", "IL5", "IL6", "IL7", "IL10", "IL13", "INFa", "INFg", "TNFa", "sCD40L")
X_pooled=X_pooled[,Groups]
Xgroups=c(6, 16)

MygPLSDA_pooled <- gPLSda(X_pooled, Y_pooled, ncomp = 1,ind.block.x = Xgroups, keepX = 1)
MygPLSDA_pooled$loadings$X

MysgPLSDA_pooled <- sgPLSda(X_pooled, Y_pooled, ncomp = 1,ind.block.x = Xgroups, keepX = 1, alpha.x = 0.1)
MysgPLSDA_pooled$loadings$X

MysgPLSDA_pooled <- sgPLSda(X_pooled, Y_pooled, ncomp = 1,ind.block.x = Xgroups, keepX = 1, alpha.x = 0.9)
MysgPLSDA_pooled$loadings$X


set.seed(1)
res_splsda=CalibratesPLSDA(dataX=X_pooled, dataY=Y_pooled, ncomp=1, Nrepeat=5)
PlotCalib(res=res_splsda)

set.seed(1)
res_sgplsda=CalibratesgPLSDA(dataX = X_pooled, dataY = Y_pooled, ncomp = 1, Nrepeat = 5, Xgroups = Xgroups)
PlotCalib(res=res_sgplsda, type="sgPLSDA")


for (subtype in c("BCLL", "DLBL", "FL", "MM")){
  sets=covars$egm_set[covars$LY_subtype==subtype]
  X=X_pooled[covars$egm_set%in%sets,]
  Y=Y_pooled[covars$egm_set%in%sets]
  assign(paste0("X_", subtype), X)
  assign(paste0("Y_", subtype), Y)
}

Y_diagnostic=covars$Time_to_diag[as.character(covars$type)=="1"] ## time to diagnostic in days
X_diagnostic=X_pooled[as.character(covars$type)=="1", ]
names(Y_diagnostic)=rownames(X_diagnostic)

for (subtype in c("BCLL", "DLBL", "FL", "MM")){
  ids=covars$egm_id[covars$LY_subtype==subtype]
  X=X_diagnostic[ids,]
  Y=Y_diagnostic[ids]
  assign(paste0("X_diag_", subtype), X)
  assign(paste0("Y_diag_", subtype), Y)
}

MyPal=brewer.pal("Paired", n = 12)

MyParameters=readRDS("Data/MyParameters.rds")
str(MyParameters)

MyPLSDA_pooled <- plsda(X_pooled, Y_pooled, ncomp=1, mode='regression')
MysPLSDA_pooled <- splsda(X_pooled, Y_pooled, keepX=MyParameters$Pooled$sPLSDA$NVar, ncomp=1, mode='regression')
MygPLSDA_pooled <- gPLSda(X_pooled, Y_pooled, ncomp = 1, ind.block.x = Xgroups, keepX = MyParameters$Pooled$gPLSDA$NGroups)
MysgPLSDA_pooled <- sgPLSda(X_pooled, Y_pooled, ncomp = 1, ind.block.x = Xgroups, 
                            keepX = MyParameters$Pooled$sgPLSDA$NGroups, alpha.x = MyParameters$Pooled$sgPLSDA$alpha)

Loadings=cbind(MysPLSDA_pooled$loadings$X, MysgPLSDA_pooled$loadings$X, rep(NA, 28), rep(NA, 28))
Loadings=as.vector(t(Loadings))
Loadings=Loadings[-c(length(Loadings)-1, length(Loadings))]

par(mar=c(10,5,3,3))
plot(Loadings, col=c(MyPal[6], MyPal[10], NA, NA), xaxt="n", ylab="Loadings Coefficients",  type="h", lwd=3, xlab="")
axis(1, at=seq(1.5, 28*4, by=4), labels = colnames(MyPLSDA_pooled$X), las=2)
axis(1, at=c(0, Xgroups, 28)*4, line = 6, labels = NA)
axis(1, at=c(3, 10.5, 21.5)*4, labels = c("Growth Factors", "Chemokines", "Cytokines"), line=6, tick = FALSE)
abline(v=c(0, Xgroups, 28)*4, lty=3, col="black")
abline(h=0, lty=2)
legend("bottomright", legend = c("sPLS-DA", "sgPLS-DA"), lty=1, lwd=3,
       col=c(MyPal[6], MyPal[10]), cex=0.75)

MyPredict=predict(MysPLSDA_pooled, newdata = X_pooled)
fitted=MyPredict$class$max.dist
table(fitted)

set.seed(1)
Stability_results=StabilityPlot(X = X_pooled, Y = Y_pooled, NIter = 100)
pheatmap(Stability_results, cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = TRUE, filename = "Figures/Fig2_B.pdf", height=5, width=10)

MyParameters$MM$sPLSDA$NVar
MyParameters$MM$sgPLSDA$NGroups
MyParameters$MM$sgPLSDA$alpha

dim(X_diagnostic)
length(Y_diagnostic)
summary(Y_diagnostic)

MyParameters$Time$sPLS$NVar
