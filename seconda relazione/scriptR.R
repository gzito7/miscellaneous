t=read.csv2("~/Università/StatSup/Seconda relazione/tabella.csv", row.names = 1)
tab=scale(t[,-9]) #standardizzo tabella
tab.pca=princomp(tab) #breve analisi delle componenti principali
summary(tab.pca)
loadings(tab.pca)
biplot(tab.pca,col=c("darkgray","red"),cex=0.65) #piano principale
    #metodi a prototipo
library(cluster)
  #kmeans
layout(t(1:2))
wss=rep(0,10) #valutazione somma delle mutue distanze tra elementi di uno stesso cluster
for(k in 2:10){
  wss[k]=kmeans(tab,k,nstart=25)$tot.withinss
}
plot(2:10,wss[2:10],type="b",pch=16)
as=rep(0,10) #valutazione andamento silhouette kmeans
for(k in 2:10){ 
  cl=kmeans(tab,k,nstart=25)$cluster
  as[k]=mean(silhouette(cl,dist(tab))[,3])
}
plot(2:10,as[2:10],type="b",pch=16)
layout(matrix(c(1:4),2,2)) #confronto silhouette interessanti
plot(silhouette(kmeans(tab,2,nstart=25)$cluster,dist(tab)))
plot(silhouette(kmeans(tab,3,nstart=25)$cluster,dist(tab)))
plot(silhouette(kmeans(tab,4,nstart=25)$cluster,dist(tab)))
plot(silhouette(kmeans(tab,5,nstart=25)$cluster,dist(tab)))
#le silhouette più interessanti si hanno per k=4,5
  #pam
layout(t(1:2))
c=rep(0,10) #valutazione andamento silhouette pam euclidea
for(i in 2:10){
  c[i]=pam(tab,i)$silinfo$avg.width
}
plot(2:10,c[2:10],type="b",pch=19)
c=rep(0,10)#valutazione andamento silhouette pam manhattan
for(i in 2:10){
  c[i]=pam(tab,i,metric="manhattan")$silinfo$avg.width
}
plot(2:10,c[2:10],type="b",pch=19)
layout(matrix(1:4,2,2)) #confronto silhouette interessanti
plot(silhouette(pam(tab,2)))
plot(silhouette(pam(tab,5)))
plot(silhouette(pam(tab,2,metric="manhattan")))
plot(silhouette(pam(tab,4,metric="manhattan")))
#le silhouette più interessante sono k=5 euclidea e k=4 manhattan
    #metodi gerarchici
  #distanza euclidea
layout(t(1:3))
tab.hcc.eu=hclust(dist(tab),"complete") 
plot(tab.hcc.eu,hang=-1,cex=0.3) #sembra interessante k=5
tab.hca.eu=hclust(dist(tab),"average")
plot(tab.hca.eu,hang=-1,cex=0.3) #non interessante
tab.hcs.eu=hclust(dist(tab),"single")
plot(tab.hcs.eu,hang=-1,cex=0.3) #non interessante
  #distanza del massimo
layout(t(1:3))
tab.hcc.max=hclust(dist(tab,method="maximum"),"complete") 
plot(tab.hcc.max,hang=-1,cex=0.3) #k=2 sembra interessante
tab.hca.max=hclust(dist(tab,method="maximum"),"average")
plot(tab.hca.max,hang=-1,cex=0.3) #non interessante
tab.hcs.max=hclust(dist(tab,method="maximum"),"single")
plot(tab.hcs.max,hang=-1,cex=0.3) #non interessante
layout(t(1:2)) #confronto silhouette interessanti
plot(silhouette(cutree(tab.hcc.eu,5),dist(tab))) #interessante
plot(silhouette(cutree(tab.hcc.max,2),dist(tab,method="maximum"))) #non interessante
    #confronto tra i cluster con silhouette più interessante
tab.km4=kmeans(tab,4,nstart=25)
tab.km5=kmeans(tab,5,nstart=25)
tab.pam5=pam(tab,5)
tab.pamman4=pam(tab,4,metric="manhattan")
tab.hcc.eu5=cutree(hclust(dist(tab),"complete"),5)
layout(t(1:2))
plot(tab.pca$scores,col=2+tab.km4$cluster,pch=20) #kmeans k=4
text(tab.pca$scores,labels=as.character(row.names(tab)),col=2+tab.km4$cluster,pos=3,cex=0.8)
abline(0,0,col="red")
plot(tab.pca$scores,col=2+tab.km5$cluster,pch=20) #kmeans k=5
text(tab.pca$scores,labels=as.character(row.names(tab)),col=2+tab.km5$cluster,pos=3,cex=0.8)
abline(0,0,col="red")
layout(t(1:2))
plot(tab.pca$scores,col=tab.pam5$clustering,pch=20) #pam
text(tab.pca$scores,labels=as.character(row.names(tab)),col=tab.pam5$clustering,pos=3,cex=0.8)
abline(0,0,col="red")
plot(tab.pca$scores,col=tab.pamman4$clustering,pch=20) #pam manhattan
text(tab.pca$scores,labels=as.character(row.names(tab)),col=tab.pamman4$clustering,pos=3,cex=0.8)
abline(0,0,col="red")
layout(1)
plot(tab.pca$scores,col=tab.hcc.eu5,pch=20) #gerarchico complete linkage k=5
text(tab.pca$scores,labels=as.character(row.names(tab)),col=tab.hcc.eu5,pos=3,cex=0.8)
abline(0,0,col="red")
    #analisi
tab.pam=pam(tab,5)
layout(t(1:2))
clusplot(tab,tab.pam$clustering,stand=F,shade=T,labels=2,col.p=tab.pam$clustering,cex.txt=0.7,col.txt=tab.pam$clustering,col.clus="darkgrey")
plot(silhouette(tab.pam))
layout(1)
silhouette(tab.pam5) #chi sono i tennisti con silhouette negativa?
library(MASS)
parcoord(tab,col=as.numeric(tab.pam$cluster),lty=5)