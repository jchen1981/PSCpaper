###########################################################
setwd('~/Dropbox/Workspace/MayoClinic/Temp/Zhuyao/2017_03_02_Analysis/PenileCancerGitHub')

source('Code/Functions.R')
temp <- load('Data/data.RData')

temp

################################################################
# Distribution of surrvial time
sort(meta.dat$Survival)

# Methylation profiles after removing obvious outliers
pdf('Result/IndividualMethylationProfile.pdf')
plot(x=c(0, 1), y=c(0, 3.5), xlim=c(0, 1), ylim=c(0, 3.5), xlab='Methylation', ylab='Density', type='n')
for (i in 1:ncol(data)) {
	lines(density(data[, i]), col='gray')
}
dev.off()

################################################################
# PERMANOVA test
distM <- fastDist(B2M(data))
meta.dat$SurvivalB <- factor(meta.dat$SurvivalB)

adonis.output <- adonis(as.dist(distM) ~ meta.dat$SurvivalB)
adonis.output

# PCoA plot / equivalent to PCA using Euclidean distance
obj <- cmdscale(as.dist(distM), k=2, eig=T)
pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
y <- cbind(obj$points[, 1], obj$points[, 2])

xlab <- paste0('PC1(', pve[1], '%)')
ylab <- paste0('PC2(', pve[2], '%)')

colnames(y) <- c("PC1", "PC2")		
xlim <- c(-max(range(abs(y[, 1]))) * 1.25, max(range(abs(y[, 1]))) * 1.25)
ylim <- c(-max(range(abs(y[, 2]))) * 1.25, max(range(abs(y[, 2]))) * 1.25)

pdf('Result/PCA.plot.pdf', width = 5, height = 5)
plot(y[, 1], y[, 2], type = 'n', xlim = xlim, ylim = ylim, 
		xlab = xlab, ylab = ylab)

s.class(y, 
		fac = meta.dat$SurvivalB,
		axesell = F,
		col = c('darkred', 'darkblue'),
		grid = TRUE,
		add.plot=T
)

title(main = paste("Euclidean distance"), sub =  paste0('(P=', adonis.output$aov.tab[1, 6], ')'))
dev.off()

################################################################
# Overall methylation 
meth.mean <- colMeans(data)
pdf('Result/MeanMethylationLevelBoxplot.pdf', width=3, height=5)
for (covar in c('SurvivalB')) {
	ind <- !is.na(meta.dat[, covar]) 
	# Parametric
	if (covar != 'Age') {
		pv <- summary(aov(meth.mean[ind] ~ meta.dat[ind, covar]))[[1]][1, 5]
		boxplot(meth.mean[ind] ~ meta.dat[ind, covar], ylab='Mean methylation over all CpGs', 
				xlab=covar, col='steelblue',
				sub=paste0('P=', round(pv, 3)))
	} else {
		pv <- cor.test(meth.mean[ind], meta.dat[ind, covar])$p.value
		plot(meta.dat[ind, covar], meth.mean[ind], ylab='Mean methylation over all CpGs', 
				xlab=covar,
				sub=paste0('P=', round(pv, 3)))
		abline(lm(meth.mean[ind] ~ meta.dat[ind, covar]), col='red')
	}
	
}
dev.off()

#######################################################################
# Genome-wide association analysis using the R package CpgAssoc
assoc.obj.outcome <- cpg.assoc(beta.val=data, meta.dat$SurvivalB, fdr.cutoff=0.20)
assoc.obj.age <- cpg.assoc(beta.val=data, meta.dat$Age, fdr.cutoff=0.20)

ind <- !is.na(meta.dat[, 'Age']) 
assoc.obj.outcome.adj <- cpg.assoc(beta.val=data[, ind], meta.dat$SurvivalB[ind], data.frame(factor(meta.dat$Age[ind])),  fdr.cutoff=0.20)

# use storey's q value procedure instead of BH procedure for better power
pv <- assoc.obj.outcome$results$P.value
qv.obj <- qvalue(pv)
qv.obj$pi0   # 0.539
qv <- qv.obj$qvalues
assoc.obj.outcome$results[, 'FDR'] <- qv

pv <- assoc.obj.age$results$P.value
qv.obj <- qvalue(pv)
qv.obj$pi0  # 0.848
qv <- qv.obj$qvalues
assoc.obj.age$results[, 'FDR'] <- qv

qv <- assoc.obj.outcome$results[, 'FDR']
pv <- assoc.obj.outcome$results[, 'P.value']

sig.res <- assoc.obj.outcome$results[qv <= 0.125, ]
sig.age <- assoc.obj.outcome.adj$results[qv <= 0.125, ]
sig.res <- cbind(sig.res, age.adj.P=sig.age$P.value)

sig.res$CPG.Labels <- as.character(sig.res$CPG.Labels)

data.sig <- data[(sig.res$CPG.Labels), ]
data.sig <- data.frame(t(data.sig), Outcome=meta.dat$SurvivalB)
obj <- aggregate(.~ Outcome, data.sig, mean)
obj <- t(obj[, -1])
colnames(obj) <- c('Good', 'Poor')
obj <- cbind(obj, difference=abs(obj[, 1] - obj[, 2]))
sig.res 
sig.res <- cbind(sig.res, obj)
write.csv(sig.res, 'Result/TopHits.csv')

# Plotting
CpG.sel <- sig.res$CPG.Labels[sig.res$difference >= 0.05 & sig.res$age.adj.P <= 0.05]
data.sig.df <- melt(data.sig)
colnames(data.sig.df) <- c('Prognosis', 'CpG', 'Methylation')
data.sig.df$CpG <- factor(data.sig.df$CpG, levels = sig.res$CPG.Labels[order(sig.res$Good, decreasing = TRUE)])
data.sig.df2 <- subset(data.sig.df, CpG %in% CpG.sel)


pdf('Result/TopHits(FDR<0.124;AjdustP<0.05;Difference>0.05.pdf', width = 8, height = 4)
dodge <- position_dodge(width=0.88)
obj <- ggplot(data.sig.df2, aes(x = CpG, y = Methylation, fill = Prognosis)) +
		geom_boxplot(position=dodge, alpha = 0.75, outlier.alpha = 0,  lwd = 0.25, fatten = 1) +
		geom_point(position=position_jitterdodge(dodge.width=0.88), size = 1, alpha = 0.3) +
		xlab('') +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 60, hjust = 1, size=8)) +
		theme(legend.position="right", legend.title = element_blank()) 
print(obj)	
dev.off()

#######################################################################
# Generate heatmap
data.obj <- list()
data.obj$data <- t(scale(t(data[rownames(data) %in% sig.res$CPG.Labels, ])))
data.obj$meta.dat <- meta.dat
generate_heatmap(data.obj,  meta.info=c('SurvivalB',  'HPV_risk', 'N.Stage', 'T.Stage', 'Grade'), 
		sam.ord=NULL, data.type='P', 
		sepwidth=0.01, colsep=NULL, rowsep=NULL,
		colFnsC=NULL, colFnsF=NULL, Rowv=T, Colv=T, dendrogram='both', margins=c(5, 5),
		in.grid=F, sepcolor='black', is.labCol=T, cexCol=1, cexRow=NULL,
		omas=c(1, 1, 1, 8), width=6, height=9, file.name='Result/Heatmap.PrognosisCpGFDR<0.15(n=101)') 

###################################################################




