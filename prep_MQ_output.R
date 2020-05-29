p <-read.delim("C:/covid_proteomics/peptides.txt", stringsAsFactors = F)

library(reshape)
??melt
grep(colnames(p), pattern="LFQ")
colnames(p)[1:50]

idcol = c("Sequence", "Leading.razor.protein")
pmelt <-melt(p, 
             id.vars=idcol,
             measure.vars=grep(colnames(p), pattern="LFQ"))
head(pmelt)

levels(as.factor(pmelt$variable))

pmelt$Condition = rep(0, nrow(pmelt))
pmelt$Condition[grep(pmelt$variable, pattern="_COVID19")]<-"COVID"
head(pmelt)
pmelt$Condition[grep(pmelt$variable, pattern="NON.COVID19")]<-"NONCO"
pmelt$Condition[grep(pmelt$variable, pattern="control")]<-"CNTRL"
levels(as.factor(pmelt$Condition))

splt_last = function(x){tail(unlist(strsplit(x, "_")), n=1)}
charnames<-as.character(pmelt$variable)
charnames[1]
splt_last(charnames[1])
pmelt$BioReplicate<-lapply(FUN=splt_last, charnames)
head(pmelt)
?write.table
df <- apply(pmelt,2,as.character)
head(df)
write.table(df, file="pecora.txt", sep="\t")
