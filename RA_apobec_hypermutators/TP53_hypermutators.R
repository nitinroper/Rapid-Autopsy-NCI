#identify apobec outliers in LUAD TCGA and correlate with TP53 status
t<-read.csv("apobec_hypermutators_TP53_tcga.csv")
apobec<-t$signature2_13
apobec1<-na.omit(apobec)
lowerq = quantile(apobec1)[2]
upperq = quantile(apobec1)[4]
iqr = upperq - lowerq 
#IQR
mild.threshold.upper = (iqr * 1.5) + upperq
mild.threshold.upper
extreme.threshold.upper = (iqr * 3) + upperq
extreme.threshold.upper

tbl = table(t$TP53status, t$hypermutator) 
tbl
chisq.test(tbl)

#identify apobec outliers in Broad luad dataset  and correlate with TP53 status
t<-read.csv("apobec_hypermutators_TP53_broad.csv")
apobec<-t$signature2_13
apobec1<-na.omit(apobec)
lowerq = quantile(apobec1)[2]
upperq = quantile(apobec1)[4]
iqr = upperq - lowerq 
#IQR
mild.threshold.upper = (iqr * 1.5) + upperq
mild.threshold.upper
extreme.threshold.upper = (iqr * 3) + upperq
extreme.threshold.upper

tbl = table(t$TP53status, t$hypermutator) 
tbl
chisq.test(tbl)

#identify apobec outliers in tracerx dataset and correlate with TP53 status
t<-read.csv("apobec_hypermutators_TP53_tracerx.csv")
apobec<-t$apobec_mutations
apobec1<-na.omit(apobec)
lowerq = quantile(apobec1)[2]
upperq = quantile(apobec1)[4]
iqr = upperq - lowerq 
#IQR
mild.threshold.upper = (iqr * 1.5) + upperq
mild.threshold.upper
extreme.threshold.upper = (iqr * 3) + upperq
extreme.threshold.upper

tbl = table(t$TP53, t$hypermutator) 
tbl
chisq.test(tbl)

