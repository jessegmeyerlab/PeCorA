t_init=Sys.time()

if (!require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
pacman::p_load(phia, ggplot2, standardize, optparse)


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=getwd(), 
              help="output directory", metavar="character"),
  make_option(c("-c", "--cont"), type="character", default="cntrl", 
              help="control group name", metavar="character"),
  make_option(c("-a", "--area"), type="character", default="Normalized.Area", 
              help="area column name", metavar="character"),
  make_option(c("-p", "--pval"), type="numeric", default=0.01, 
              help="significant p-value, [default = %default]", metavar="numeric")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# define variables
input_file = opt$file  # input file is first
output_file = opt$out # output file name
control_name = opt$cont # control name is second
area_column_name = opt$area # area column name
working_directory = opt$dir  # output directory is last
sig_pval_lvl = opt$pval  # minimum pvalue

print(sig_pval_lvl)

#print(paste("file name=", opt$file))
#print(paste("value column=", area_column_name))
#print(paste("value column=", opt$area))
#print(paste("output dir=", opt$dir))


if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}


# make sure the working directory ends in /
if(substr(working_directory,nchar(working_directory), nchar(working_directory))!="/"){
  working_directory = paste(working_directory,"/", sep="")
}
suppressWarnings(try(dir.create(working_directory)))



# read the file
t<-read.csv(input_file, stringsAsFactors = F)

#print(colnames(t))
### CHECKS FOR CORRECT COLUMNS
if(!"Peptide.Modified.Sequence" %in% colnames(t)){
  stop("missing column named Peptide.Modified.Sequence", call.=FALSE)
}

if(!"Condition" %in% colnames(t)){
  stop("missing column named Condition", call.=FALSE)
}
if(!"BioReplicate" %in% colnames(t)){
  stop("missing column named BioReplicate", call.=FALSE)
}
if(!"Protein" %in% colnames(t)){
  stop("missing column named Protein", call.=FALSE)
}

### CHECK for at least 2 conditions
if(length(unique(t$Condition))<2){
  stop("must have at least 2 conditions", call.=FALSE)
}



# make new column that combines petide and charge
t["modpep_z"] = paste(t$Peptide.Modified.Sequence, "all", sep = "_")
# make column that combines condition and replicate
t["condition_rep"] = paste(t$Condition, t$BioReplicate, sep = "_")
n_reps_total<-length(unique(t$condition_rep))


## filter the data to remove NA, <100 area, and then reps without all measures
if(length(which(is.na(as.numeric(t[,area_column_name]))==TRUE))>0){
  t <- t[-suppressWarnings(which(is.na(as.numeric(t[,area_column_name]))==TRUE)),]
}
t <- t[-which(as.numeric(t[,area_column_name])<=100),]
t <- t[t$modpep_z %in% names(which(table(t$modpep_z)==n_reps_total)),]


## scale the data and then plot
t$ms1log2<-log2(as.numeric(t[,area_column_name]))
t$ms1scaled <- scale_by(ms1log2 ~ Condition*BioReplicate, t)

png(paste(working_directory,"gobal_data_scaling.png",sep=""),
    width = 8, height = 6, units="in", res=300)
par(las=2, mfcol=c(2,1))
boxplot(ms1log2 ~ Condition*BioReplicate, t, main="raw data")
boxplot(ms1scaled ~ Condition*BioReplicate, t, main="globally scaled")
quiet<-dev.off()

########### data now scaled to 0 center ###################
###  next subtract the mean value of control from each peptide
### takes 3 minutes for 31k peptides, 10k/min
### change to vectors
ms1adj <- rep(0, nrow(t))
i=1
t_cntrl<-t[t$Condition==control_name,]
ms1scaled_cntrl<-t_cntrl[,"ms1scaled"]
ms1scaled_full<-t[,"ms1scaled"]
t0<-Sys.time()
print("scaling peptides to control == 0")
pb <- txtProgressBar(min = 0, max = length(unique(t$modpep_z)), style = 3)
for(x in unique(t$modpep_z) ){
  ms1adj[which(t$modpep_z==x)] <- ms1scaled_full[which(t$modpep_z==x)] - mean(ms1scaled_cntrl[which(t_cntrl$modpep_z==x)])  # subtract the ave of control
  setTxtProgressBar(pb, i)
  i=i+1
}
t$ms1adj<-ms1adj
print(" ")
print(paste("peptide scaling completed in", round(Sys.time()-t0,2), "minutes", sep=" "))

# get all unique protein groups
print("checking which proteins still have at least 2 peptides")
pgs<-levels(as.factor(t$Protein)) 
pgs_morethan2<-c()
for(x in pgs){
  if(length(unique(t[t$Protein %in% x, "modpep_z"]))>1){
    pgs_morethan2<-c(pgs_morethan2, x)
  }
}

### loop through proteins with >2 pep, get subset df
###### look through peptide precursors in that protein
######### check linear model, record p value
############ adjust pvals within each protein ### could do that differently?
allp<-list() # empty list to store pvalues for each protein
j=1 # progress bar iterator
t0<-Sys.time()  # initial system time
print("computing the interaction p-values")
pb <- txtProgressBar(min = 0, max = length(pgs_morethan2), style = 3) 
for( x in pgs_morethan2){ #loop through protein groups
  tmpdf <- t[t$Protein ==x,]  ## get the subset dataframe
  tmpdf["allothers"]<-rep("allothers", times=nrow(tmpdf))
  pvalues<-c(rep(0, length(unique(tmpdf$modpep_z))))
  i=1
  for( y in unique(tmpdf$modpep_z)){
    subtmpdf <- tmpdf
    subtmpdf[which(tmpdf$modpep_z ==y),"allothers"] <- y
    tmplm<-lm(subtmpdf$ms1adj ~ subtmpdf$Condition*subtmpdf$allothers)
    tmpanova<-Anova(tmplm)
    pvalues[i]<-tmpanova$`Pr(>F)`[3]
    i=i+1
  }
  allp[[x]]<-pvalues # record p-values
  setTxtProgressBar(pb, j)
  j=j+1
}

print(" ")
print(paste("PeCorA finished in ", round(Sys.time()-t0,2), " minutes", sep=""))


### post testing analysis
print(  paste("number of proteins tested =", length(allp), sep=" ")  )
print(  paste("number of peptides tested =", length(unlist(allp)), sep=" ")  )



######## make table of the peptides that disagree   ##########
### dataframe with all values
print("started making data table")
alldf= data.frame()
x<-names(allp)[1]
for(x in names(allp)){
  #print(x)
  tmpdf <- t[t$Protein ==x,]
  #tmpdf["allothers"]<-rep("allothers", times=nrow(tmpdf))
  tmp_peps <- unique(tmpdf$modpep_z)
  if(length(tmp_peps)>0){
    tmp_pval <- allp[[x]]
    tmpout= cbind.data.frame(protein=rep(x, length(allp[[x]])), tmp_peps, tmp_pval=as.numeric(tmp_pval))
    alldf=rbind( alldf,  tmpout)
  }
}
print("correcting p-values")
alldf$adj_pval<-p.adjust(alldf$tmp_pval, method="BH")     ## adjust p-values
alldf_ordered<-alldf[order(alldf$adj_pval),]              ## sort 

## print summary results
print(paste("number of uncorrelated peptides =", nrow(alldf[alldf$adj_pval<=sig_pval_lvl,]),  sep=" "))
print(paste("number of proteins with uncorrelated peptides =", 
            length(unique(alldf[alldf$adj_pval<=sig_pval_lvl,]$protein)),
            sep=" "
))

colnames(alldf_ordered)[2]<-"peptide"
colnames(alldf_ordered)[3]<-"pvalue"

write.table(alldf_ordered, paste(working_directory, opt$out, sep=""), 
            sep="\t", quote=F, row.names = F)

print("output table written")


########### print images from table #############3

############ visualize those that are supposedly different   #############
### loop through all the peptides that are different
## make new directory with protein name for them
# put plots of the peptide in there to look at


sign_prots<-as.character(unique(alldf[alldf$adj_pval<=sig_pval_lvl,]$protein))
print("started printing peptide figures")
suppressWarnings(dir.create(paste(working_directory, "/allplots", sep="")))
setwd(paste(working_directory, "/allplots", sep=""))

for(x in sign_prots){
  tmpdf <- t[t$Protein ==x,]
  tmpdf["allothers"]<-rep("allothers", times=nrow(tmpdf))
  tmpalldf<-alldf[alldf$protein %in% x,]
  tmp_peps<-unique(tmpalldf$tmp_peps)[which(tmpalldf$adj_pval<sig_pval_lvl)]
  if(length(tmp_peps)>0){
    for(y in tmp_peps){
      subtmpdf <- tmpdf
      subtmpdf[which(tmpdf$modpep_z ==y),"allothers"] <- y
      #boxplot(ms1adj ~ allothers*Condition, subtmpdf)
      p=ggplot(subtmpdf, aes(x=Condition, y=ms1adj, fill=allothers)) + 
        geom_boxplot(size=1.5) +
        geom_point(pch = 21, position = position_jitterdodge(), size=3)+ 
        theme(text = element_text(size=20))  +
        ylab("Log2(intensity)-Log2(control)") +
        xlab("Group")+
        ggtitle(x)+
        scale_fill_manual(values=c("grey", "#00BA38")) +theme(legend.position = "bottom")
      ### part for protein string formatted with prefix sp|
      if(substr(x, 3,3)=="|"){
        png(paste(substr(x, start=4,9), y, ".png", sep="_"))
        print(p)
        dev.off()
      }
      if(substr(x, 3,3)!="|"){
        png(paste(x, y, ".png", sep="_"))
        print(p)
        dev.off()
      }
    }
  }
}

print("printing images finshed")
print(paste("total time = ", round(Sys.time()-t_init,2), " minutes", sep="" ))

