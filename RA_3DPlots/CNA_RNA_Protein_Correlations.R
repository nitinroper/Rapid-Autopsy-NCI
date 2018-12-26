#correlating protein/rna/copy number heterogeneity
#RA000
t<-read.csv("RA00*_3DPlots.csv")

cor.test(t$protein_correlation, t$cnv_purity, method="pearson")
cor.test(t$protein_correlation, t$rna_correlation, method="pearson")
cor.test(t$rna_correlation, t$facets_cnv_purity_adjusted_correlation, method="pearson")



