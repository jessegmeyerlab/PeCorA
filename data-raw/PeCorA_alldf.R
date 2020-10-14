## code to prepare `PeCorA_alldf` dataset goes here
PeCorA_alldf <- function(t){
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
  print("output alldf")
  alldf
}

#alldf <-  PeCorA (scaled_peptides)
usethis::use_data(PeCorA_alldf, overwrite = TRUE)
