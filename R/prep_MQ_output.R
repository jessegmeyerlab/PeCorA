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
df[,1]
colnames(df)
colnames(df)[1]<-"Peptide.Modified.Sequence"
colnames(df)[2]<-"Protein"

write.table(df, file="pecora.txt", sep=",", row.names = F)


#### drop the plasma pools

pmelt_nocontrol<-pmelt[-grep(pmelt$variable, pattern="control"),]
nrow(pmelt)
nrow(pmelt_nocontrol)
levels(as.factor(pmelt_nocontrol$Condition))
head(pmelt_nocontrol)
df_noc <- apply(pmelt_nocontrol,2,as.character)
colnames(df_noc)
colnames(df_noc)[1]<-"Peptide.Modified.Sequence"
colnames(df_noc)[2]<-"Protein"
write.table(df_noc, file="pecora_nocontrol.txt", sep=",", row.names = F)

