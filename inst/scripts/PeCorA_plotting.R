## code to prepare `PeCorA_plotting` dataset goes here

PeCorA_plotting <- function (w,u,v){
  sign_prots<-as.character(unique(u[u$adj_pval<=0.01,]$protein))
  #suppressWarnings(dir.create(paste(getwd(), "/test_allplots", sep="")))
  #setwd(paste(paste(getwd(), "/test_allplots", sep="")))

  for(x in sign_prots){
    tmpdf <- v[v$Protein ==x,]
    tmpdf["allothers"]<-rep("allothers", times=nrow(tmpdf))
    tmpalldf<-w[w$protein %in% x,]
    tmp_peps<-unique(tmpalldf$peptide)[which(tmpalldf$adj_pval<0.01)]
    if(length(tmp_peps)>0){
      for(y in tmp_peps){
        subtmpdf <- tmpdf
        subtmpdf[which(tmpdf$modpep_z ==y),"allothers"] <- y
        #boxplot(ms1adj ~ allothers*Condition, subtmpdf)
        p=ggplot(subtmpdf, aes(x=Condition, y=ms1adj, fill=allothers)) +
          geom_boxplot(size=1.5) +
          geom_point(pch = 21, position = position_jitterdodge(), size=3)+
          theme(text = element_text(size=20)) +
          ylab("Log2(intensity)-Log2(control)") +
          xlab("Group")+
          ggtitle(x)+
          scale_fill_manual(values=c("grey", "#00BA38")) +theme(legend.position = "bottom")
        ### part for protein string formatted with prefix sp|
        #if(substr(x, 3,3)=="|"){
         # png(paste(substr(x, start=4,9), y, ".png", sep="_"))
          #print(p)
          #dev.off()
        #}
       # if(substr(x, 3,3)!="|"){
          #png(paste(x, y, ".png", sep="_"))
          #print(p)
          #dev.off()
       # }
      }
    }
  }
  p
}

#PeCorA_plotting_plot<-PeCorA_plotting(disagree_peptides,disagree_peptides[1,],scaled_peptides)
usethis::use_data(PeCorA_plotting, overwrite = TRUE)
