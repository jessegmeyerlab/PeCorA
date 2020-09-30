if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("UniProt.ws")

library(UniProt.ws)


getwd()
setwd("C:/covid_proteomics/pecora/no_control/")

o<-read.delim("out.txt", stringsAsFactors = F)

head(o)


split_sp <- function(x)unlist(strsplit(x, " "))[1]

up <- UniProt.ws(taxId=9606)  ### 9606 is human, 10090 is mouse
taxId(up) <- 9606
uniprot <- o$protein
genemap<-select(up,unique(uniprot), "GENES")


genes <- unlist(lapply(genemap[,2], split_sp))
geneid.reg <- genes[match(uniprot, genemap[,"UNIPROTKB"])]



protein.description.map<-select(up,keys=unique(uniprot),columns="PROTEIN-NAMES")

prot.descrs<-unlist(protein.description.map[,2])
descr.id.reg <- prot.descrs[match(uniprot, protein.description.map[,"UNIPROTKB"])]

geneid.reg
reg.out<-cbind(geneid.reg,descr.id.reg, o)
head(reg.out)
getwd()
write.table(file="out.annotated.txt",reg.out,
            quote=F, 
            sep="\t",
            row.names = F)


