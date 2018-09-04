#Simulate Data
#For this example we first simulate WGBS data with known regions of conserved methylation and differential methylation between treatment groups.  For each of 6 individuals we construct a methylome consisting of 800 cytosine nucleotides.  The seuqencing depth for every nucleotide is 10 in this example, and we use a random binomial distriobution with different probabilities to segment the genome into 6 different chunks.  3 of these segments have the same probability distributrions for all 6 samples (segments 1,3,and 5), 2 have consistent differences in the probability of methylation between group C and D (segments 2 and 4), and 1 segment has highly variable differences in methylation between individuals that leads to a difference in mean methylation between groups, but also a lot of noise, making it an inconsistent response.

set.seed(400)

#Using binomial draws to generate the number of methylated nucleotides out of 10 for each of the 800 basepairs for each of the 6 individuals.
c(rbinom(100,10,0.5), rbinom(200,10,0.1), rbinom(50,10,0.9), rbinom(30,10,0.7), rbinom(300,10,0.4), rbinom(120,10,0.2))->C1
c(rbinom(100,10,0.5), rbinom(200,10,0.1), rbinom(50,10,0.9), rbinom(30,10,0.7), rbinom(300,10,0.4), rbinom(120,10,0.35))->C2
c(rbinom(100,10,0.5), rbinom(200,10,0.1), rbinom(50,10,0.9), rbinom(30,10,0.7), rbinom(300,10,0.4), rbinom(120,10,0.3))->C3
c(rbinom(100,10,0.5), rbinom(200,10,0.35), rbinom(50,10,0.9), rbinom(30,10,0.1), rbinom(300,10,0.4), rbinom(120,10,0.3))->D1
c(rbinom(100,10,0.5), rbinom(200,10,0.35), rbinom(50,10,0.9), rbinom(30,10,0.1), rbinom(300,10,0.4), rbinom(120,10,0.1))->D2
c(rbinom(100,10,0.5), rbinom(200,10,0.35), rbinom(50,10,0.9), rbinom(30,10,0.1), rbinom(300,10,0.4), rbinom(120,10,0.95))->D3

#Making the data.frame that contains the # of methylated and total reads at each site
data.frame(C1, 10)->C1
data.frame(C2, 10)->C2
data.frame(C3, 10)->C3
data.frame(D1, 10)->D1
data.frame(D2, 10)->D2
data.frame(D3, 10)->D3
colnames(C1)<-c("Methylated", "Total")
colnames(C2)<-c("Methylated", "Total")
colnames(C3)<-c("Methylated", "Total")
colnames(D1)<-c("Methylated", "Total")
colnames(D2)<-c("Methylated", "Total")
colnames(D3)<-c("Methylated", "Total")

#Calculating percent methylation at each site for each nucleotide
C1$Percent<-C1$Methylated/C1$Total
C2$Percent<-C2$Methylated/C2$Total
C3$Percent<-C3$Methylated/C3$Total
D1$Percent<-D1$Methylated/C1$Total
D2$Percent<-D2$Methylated/D2$Total
D3$Percent<-D3$Methylated/D3$Total

#Everything up to this point would be done with actual data through mapping reads to your refrence genome, followed by using methimpute to impute methylation on a nucleotide by nucleotide level. Methimpute is an important step in that it allows for more accurate imputation of the precent methylation levels at each nucleotide. This would come from the rc.meth.lvl value for each individual at each nucleotide.

#Calculating the average percent methylation for each of our two groups. and combing them into a single data.frame
C_Ave<-((C1$Percent+ C2$Percent+ C3$Percent)/3)
D_Ave<-((D1$Percent + D2$Percent+ D3$Percent)/3)
cbind(C_Ave, D_Ave)->Mean_Groups
data.frame(Mean_Groups)->MG

#Calculate differences in methylation between groups for each nucleotide
MG$Dif<-MG$D_Ave-MG$C_Ave

#load ChangepointLibrary
library(changepoint)

#Identifying changepoints in the data, scanning across the genome to find regions where there is a shift in the difference in methylation between our two groups. With real data this would be done for a single type of methylation (CG, CHG, CHH) for a single chromosome at a time.
#We use a PELT changepoint method with a manual penalty of 1.4.  Decreasing penalty values will be less conservative, and increasing will be more.
cpt.mean(MG$Dif, method="PELT", penalty = "Manual", pen.value = 1.4)->BreakPoints

#Plots all breakpoints
plot(BreakPoints)

#construct a file (BPs) that has the breakdown of all of your genome segments by the difference in methylation between groups
c(1,cpts(BreakPoints))->start
c(cpts(BreakPoints)-1, nrow(MG))->stop
data.frame(cbind(start,stop))->BPs
BPs$DifMean<-param.est(BreakPoints)$mean

#construct a file that filters segments only to get those with at least a 4% difference in methylation.  This can be modified to only include more extreme changes in methylation. 
BPs[abs(BPs$DifMean)>0.04,]->DSegs
#Setting some parameters to make this code applicable to real data as well as simulated data

#Number of Individuals=a
a<-6
#Number of individuals in Group 1
b<-3
#Number of individuals in Group 2
c<-3

#Constructing the RegionDMRs file which goes back to the original binomial style data rtather than the percent data (or smoothed percent when used in conmjunction with methimpute) to explicity test DMRs in these regions
RegionDMRs<- data.frame(Region=rep(0,nrow(DSegs)*a),Individual=rep(0,nrow(DSegs)*a), Group=rep(0,nrow(DSegs)*a), Methylated=rep(0,nrow(DSegs)*a), Total=rep(0,nrow(DSegs)*a))
RegionDMRs[,1]<-rep(1:nrow(DSegs), each=a)
RegionDMRs[,2]<-rep(c("C1", "C2", "C3", "D1", "D2", "D3"), nrow(DSegs))
RegionDMRs[,3]<-rep(rep(c("C","D"), c(b,c)), nrow(DSegs))

for(i in 1:nrow(DSegs)){
  DSegs$start[i]->begin
  DSegs$stop[i]->end
  sum(C1[begin:end,1])->RegionDMRs[(((i-1)*6)+1),4]
  sum(C1[begin:end,2])->RegionDMRs[(((i-1)*6)+1),5]
  sum(C2[begin:end,1])->RegionDMRs[(((i-1)*6)+2),4]
  sum(C2[begin:end,2])->RegionDMRs[(((i-1)*6)+2),5]  
  sum(C3[begin:end,1])->RegionDMRs[(((i-1)*6)+3),4]
  sum(C3[begin:end,2])->RegionDMRs[(((i-1)*6)+3),5]  
  sum(D1[begin:end,1])->RegionDMRs[(((i-1)*6)+4),4]
  sum(D1[begin:end,2])->RegionDMRs[(((i-1)*6)+4),5]  
  sum(D2[begin:end,1])->RegionDMRs[(((i-1)*6)+5),4]
  sum(D2[begin:end,2])->RegionDMRs[(((i-1)*6)+5),5]  
  sum(D3[begin:end,1])->RegionDMRs[(((i-1)*6)+6),4]
  sum(D3[begin:end,2])->RegionDMRs[(((i-1)*6)+6),5] 
}


#Loading stats package
library(lme4)

#Adding PValue and FDRPvalue columns to the data frame
DSegs$Pvals<-rep(1)
DSegs$FDR_PVals<-rep(1)

#Testing for consistency of differential methylation in putative DMR regions
for(i in 1:nrow(DSegs)){
  RegionDMRs[RegionDMRs$Region==i,]->x1
  jj3<-glmer(cbind(Methylated,Total-Methylated)~Group+(1|Individual),data=x1,family=binomial,glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
  summary(jj3)->k
  k$coefficients->coefs
  coefs[2,4]->DSegs$Pvals[i]
}

#Adding size and FDR data
DSegs$Size_Region<-DSegs$stop-DSegs$start
p.adjust(DSegs$Pvals)->DSegs$FDR_PVals

#In this simulated data set we identify nearly the full DMRs in the two consistent DMRs (199 and 29 bp when real was 200 and 30) as well as the variable DMR of 120 basepairs (found here as 122 bps).  The first two DMRs are found to be extremely significant in terms of consistency of the DMR, while the 3rd putative DMR is found not to be a consistent DMR (p=0.49) as would be expected.