################################################################
 #### Script by Marion Talbi and Milan Malinsky
 #### This is used to analyse the correlations in 5Mb windows
 #### for  simulated and empirical data and to obtain
 #### and analyse The Population Recombination Index (PRDI)
################################################################
 #### R library
library(dplyr)
library(tidyr)
library(ggplot2)
################################################################
 #### Functions
### Return the plot with the correlation between and within ecotypes
ViolinPlot<-function(recMaps2kb,method_name){
  loop_benabenb=NULL;loop_litalitb=NULL;loop_benalitb=NULL;loop_benalita=NULL;loop_benblitb=NULL;loop_benblita=NULL
  i=1
  while (i<=length(recMaps2kb$Chr)) {
    data_subset=recMaps2kb[i:(i+2499),]
    cor_litalitb=cor.test(data_subset$mean_r_lit_a,data_subset$mean_r_lit_b,method = method_name)$estimate
    cor_benabenb=cor.test(data_subset$mean_r_ben_a,data_subset$mean_r_ben_b,method = method_name)$estimate
    cor_benalitb=cor.test(data_subset$mean_r_ben_a,data_subset$mean_r_lit_b,method = method_name)$estimate
    cor_benalita=cor.test(data_subset$mean_r_ben_a,data_subset$mean_r_lit_a,method = method_name)$estimate
    cor_benblitb=cor.test(data_subset$mean_r_ben_b,data_subset$mean_r_lit_b,method = method_name)$estimate
    cor_benblita=cor.test(data_subset$mean_r_ben_b,data_subset$mean_r_lit_a,method = method_name)$estimate
    
    loop_litalitb=c(loop_litalitb,cor_litalitb)
    loop_benabenb=c(loop_benabenb,cor_benabenb)
    loop_benalitb=c(loop_benalitb,cor_benalitb)
    loop_benalita=c(loop_benalita,cor_benalita)
    loop_benblitb=c(loop_benblitb,cor_benblitb)
    loop_benblita=c(loop_benblita,cor_benblita)
    i=i+2500
  }
  
  all_cor<-as.data.frame(cbind(loop_benabenb,loop_litalitb,loop_benalitb,loop_benalita,loop_benblitb,loop_benblita))
  all_cor_violin<-pivot_longer(all_cor,cols=1:6,names_to="comparison",values_to="spearman_cor")
  print(ggplot(all_cor_violin)+
          geom_violin(aes(x=reorder(comparison,-spearman_cor),y=spearman_cor))+
          geom_boxplot(aes(x=comparison,y=spearman_cor),alpha=0.3,size=0.2)+
          theme_bw()+ylim(-0,1)+ggtitle("Correlation r - 2Kb windows "))
}
### Return the Median Distance (Dms or Dme)
CorrelationDistanceMedian<-function(recMaps2kb,method_name){
  loop_benabenb=NULL;loop_litalitb=NULL;loop_benalitb=NULL;loop_benalita=NULL;loop_benblitb=NULL;loop_benblita=NULL
  
  i=1
  while (i<=length(recMaps2kb$Chr)) {
    data_subset=recMaps2kb[i:(i+2499),]
    cor_litalitb=cor.test(data_subset$mean_r_lit_a,data_subset$mean_r_lit_b,method = method_name)$estimate
    cor_benabenb=cor.test(data_subset$mean_r_ben_a,data_subset$mean_r_ben_b,method = method_name)$estimate
    cor_benalitb=cor.test(data_subset$mean_r_ben_a,data_subset$mean_r_lit_b,method = method_name)$estimate
    cor_benalita=cor.test(data_subset$mean_r_ben_a,data_subset$mean_r_lit_a,method = method_name)$estimate
    cor_benblitb=cor.test(data_subset$mean_r_ben_b,data_subset$mean_r_lit_b,method = method_name)$estimate
    cor_benblita=cor.test(data_subset$mean_r_ben_b,data_subset$mean_r_lit_a,method = method_name)$estimate
    
    loop_litalitb=c(loop_litalitb,cor_litalitb)
    loop_benabenb=c(loop_benabenb,cor_benabenb)
    loop_benalitb=c(loop_benalitb,cor_benalitb)
    loop_benalita=c(loop_benalita,cor_benalita)
    loop_benblitb=c(loop_benblitb,cor_benblitb)
    loop_benblita=c(loop_benblita,cor_benblita)
    i=i+2500
  }
  return(median(min(loop_litalitb,loop_benabenb))-median(c(loop_benalitb,loop_benalita,loop_benblita,loop_benblitb)))
}
### Return the plot of the PRDI values function of the local PCA outlier
ReturnPRDIlocalPCAoutlier<-function(DataDmeDms){
  par(mfrow=c(1,2))
  plot((DataDmeDms$PRDI[c(1,4,5,7)])~DataDmeDms$SplitTime[c(1,4,5,7)], cex=3, pch= 8,lwd=2, type = "b",lty = 2,ylab=c("PRDI"),
       xlab=c("Split time simulated"),ylim=c(0,0.2),main="High Migration Rate")
  lines((DataDmeDms$PRDIinMDS[c(1,4,5,7)])~DataDmeDms$SplitTime[c(1,4,5,7)],col="#FF7256", cex=3, pch= 8,lwd=2, type = "b",lty = 2,ylab=c("PRDI"),
        xlab=c("Split time simulated"),ylim=c(0,0.2))
  lines((DataDmeDms$PRDIoutMDS[c(1,4,5,7)])~DataDmeDms$SplitTime[c(1,4,5,7)],col="slateblue1", cex=3, pch= 8,lwd=2, type = "b",lty = 2,ylab=c("PRDI"),
        xlab=c("Split time simulated"),ylim=c(0,0.2))
  lines((DataDmeDms$PRDIinMDSsign[c(1,4,5,7)])~DataDmeDms$SplitTime[c(1,4,5,7)],col="#FFF68F", cex=3, pch= 8,lwd=2, type = "b",lty = 2,ylab=c("PRDI"),
        xlab=c("Split time simulated"),ylim=c(0,0.2))
  
  plot(DataDmeDms$PRDI[c(2,3,6,8)]~ DataDmeDms$SplitTime[c(2,3,6,8)], cex=3, pch= 9,lwd=2, type = "b",lty = 2,ylab=c("PRDI"),
       xlab=c("Split time simulated"),ylim=c(0,0.2),main="Low Migration Rate")
  lines(DataDmeDms$PRDIinMDS[c(2,3,6,8)]~ DataDmeDms$SplitTime[c(2,3,6,8)],col="#FF7256", cex=3,pch = 9,lwd=2, type = "b", lty = 2)
  lines(DataDmeDms$PRDIoutMDS[c(2,3,6,8)]~ DataDmeDms$SplitTime[c(2,3,6,8)],col="slateblue1", cex=3,pch = 9,lwd=2, type = "b", lty = 2)
  lines(DataDmeDms$PRDIinMDSsign[c(2,3,6,8)]~ DataDmeDms$SplitTime[c(2,3,6,8)],col="#FFF68F", cex=3,pch = 9,lwd=2, type = "b", lty = 2)
  legend(2000, 0.05, legend=c("All genome", "Inside local PCA outlier","Outside local PCA outlier","Inside significant local PCA outlier"),
         col=c("black","#FF7256","slateblue1","#FFF68F"), lty=1, cex=1,lwd=2)
}

################################################################
 #### Preparation of the dataset
 #### Load the dataset from Dryad
data_2kb=read.table("All_genome_replicate_pyrho_2kb.optimize")[,1:7]
colnames(data_2kb)<-c("Chr","Start","End","mean_r_ben_a","mean_r_ben_b","mean_r_lit_a","mean_r_lit_b")

 #### Output the Dm values used for the calcul of PRDI
CorrelationDistanceMedian(data_2kb,"spearman") ## We run this function for each simulated dataset

##We saved the file and used it here to calculate PRDI in different parts of the genome
##Calcul of PRDI - inside, outside and inside significant local PCA outlier
DataDmeDms=read.table("MigrRateSplitTimeDmsDmeData.txt",h=T)
DataDmeDms<-cbind(DataDmeDms,(as.numeric(DataDmeDms[9,3])-as.numeric(DataDmeDms$Dm)), ##AllGenome
                  (as.numeric(DataDmeDms[10,3])-as.numeric(DataDmeDms$Dm)), ##Inside local PCA outliers
                  (as.numeric(DataDmeDms[11,3])-as.numeric(DataDmeDms$Dm)), ##Outside local PCA outliers
                  (as.numeric(DataDmeDms[12,3])-as.numeric(DataDmeDms$Dm))) ##Inside significant local PCA outliers
colnames(DataDmeDms)=c("MigrRate","SplitTime","Dm","PRDI","PRDIinMDS","PRDIoutMDS","PRDIinMDSsign")
 

################################################################
#### Plot the Figure
 ####Figure 2 - Panel A and B
 ### Output the Violin plot of the correlation between 5Mb bin of 2kb windows
ViolinPlot(data_2kb,"spearman")
 ####Figure 2 - Panel C
 ##Return the plot with the values of PRDI function of the split time and the migration rates 
plot((DataDmeDms$PRDI[c(1,4,5,7)])~DataDmeDms$SplitTime[c(1,4,5,7)], cex=2, pch= 8, type = "b",lty = 2,ylab=c("PRDI"),
     xlab=c("Split time simulated"),ylim=c(0,0.2))
lines(DataDmeDms$PRDI[c(2,3,6,8)]~ DataDmeDms$SplitTime[c(2,3,6,8)], cex=2,pch = 9, type = "b", lty = 2)
abline(h=DataDmeDms$Dm[9],lty=2)

 ####Supplementary Figure 12
 ### Supp figure in out MDS
ReturnPRDIlocalPCAoutlier(DataDmeDms)






