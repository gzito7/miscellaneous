tab <- read.csv2("~/Università/StatSup/Prima relazione/tabella.csv",row.names = 1)
#analisi preliminare della tabella
head(tab,10)
str(tab)
summary(tab)
library(corrplot)
corrplot(cor(tab),"number",number.cex = 0.8)
tab.pca=princomp(tab,cor=T) #calcolo le componenti principali
summary(tab.pca)
#visualizzo l'andamento della varianza cumulata
pv=(tab.pca$sdev^2)/(sum(tab.pca$sdev^2))
plot(cumsum(pv), type = "b", col= "red", ylim=c(0,1), lwd=2, pch=20)
abline(0.8,0,col="blue")
#visualizzo i piani principali relativi alle prime tre componenti
biplot(tab.pca,col=c("darkgray","red"),cex=0.6)
biplot(tab.pca,col=c("darkgray","red"),choice=c(1,3),cex=0.6)
biplot(tab.pca,col=c("darkgray","red"),choice=c(2,3),cex=0.6)
#studio la matrice dei loadings per l'interpretazione quantitativa delle componenti
loadings(tab.pca)
varimax(loadings(tab.pca)[,1:3])$loadings
#valuto la stabilità del risultato
res=rep(0,50)
for(i in 1:50){
  tab.r=data.frame(scale(tab))[-i,]
  tab.r.pca=princomp(tab.r)
  tabp=predict(tab.r.pca,newdata=data.frame(scale(tab))[i,])[1:3]
  res[i]=mean((tabp-predict(tab.pca)[i,1:3])^2)
}
sqrt(mean(res))
hist(res)
round(res,2)
tab[c(6,9,43,50),] #tennisti che più incidono nell'analisi
res[c(6,9,43,50)]
#cosa succede alle componenti se rimuovo uno di questi tennisti?
loadings(princomp(tab[-50,],cor=T))[,1:3] #cambia il verso della seconda componente
#verifichiamo se ciò accade anche per gli altri
res=rep(0,50)
for(i in 1:50){
   tab.r=data.frame(scale(tab))[-i,]
   tab.r.pca=princomp(tab.r)
   tabp1=predict(tab.r.pca,newdata=data.frame(scale(tab))[i,])[1:3]
   tabp2=c(tabp1[1],-tabp1[2],tabp1[3])
   res[i]=min(mean((tabp1-predict(tab.pca)[i,1:3])^2),mean((tabp2-predict(tab.pca)[i,1:3])^2))
}
sqrt(mean(res))
hist(res) #adesso è stabile
#assegno i pesi dei nuovi fattori ai tennisti
tab.p=predict(tab.pca)
head(tab.p[,1:3],10)
summary(tab.p[,1:3])
tab[c(8,22),] #confronto i fattori iniziali relativi a Rune e Mannarino
rev(sort(tab.p[,1]))[1:10] #i dieci migliori tennisti per la prima componente
sort(tab.p[,2])[1:17] #tennisti molto forti al servizio
sort(tab.p[,2])[18:33] #tennisti equilibrati
sort(tab.p[,2])[34:50] #tennisti molto forti in risposta
sort(tab.p[,3]) #tennisti tendenti all'errore in ordine crescente