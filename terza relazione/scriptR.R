#funzione utile
#restituisce la data dopo l'unità successive a s nella forma (maggiore,minore)
ts_data <- function(s,t,l){
  #s: data di inizio come (maggiore,minore)
  #t: periodo
  #l: unità  successive
  e=s+c(floor(l/t),l%%t)
  e=c(e[1]+floor(e[2]/t),e[2]%%t)
  return(e)
}
tab=read.csv("~/Università/StatSup/Terza relazione/tabella.csv") #carico la tabella
sgs=ts(tab[,2],frequency=12,start=c(1992,1))
sgs
plot(sgs)
n=length(sgs)
idt=start(sgs)
fdt=end(sgs)
pdt=frequency(sgs)
    #decomposizone serie storica
layout(t(1:2))
acf(sgs,30,main="")
acf(diff(sgs),main="")
layout(t(1:2)) #confronto anni
par(bg="black")
m_sgs=matrix(ts(c(tab[,2],rep(NA,3)),frequency=12,start=c(1992,1)),12,33)
ts.plot(m_sgs,col=heat.colors(33))
ts.plot(scale(m_sgs,scale=F),col=heat.colors(33))
par(bg="white")
layout(1)
  #decomposizone additiva
sgs.da=decompose(sgs)
plot(sgs.da)
layout(1:2) #analisi dei residui escluso 2020
sgs.da.r=c(na.omit(window(sgs.da$random,end=c(2019,12))),na.omit(window(sgs.da$random,c(2021,1))))
plot(sgs.da.r,pch=20)
acf(sgs.da.r,main="")
sd(acf(sgs.da.r,plot=F)$acf)
layout(t(1:2)) #confronto gaussiana
hist(sgs.da.r,20,freq=F,ylim=c(0,0.0033),main="")
lines(density(sgs.da.r),col="blue")
lines(sort(sgs.da.r),dnorm(sort(sgs.da.r),mean(sgs.da.r),sd(sgs.da.r)),col="red")
mean((sgs.da.r-mean(sgs.da.r))^3)/(mean((sgs.da.r-mean(sgs.da.r))^2)^(3/2)) #skewness
mean((sgs.da.r-mean(sgs.da.r))^4)/(mean((sgs.da.r-mean(sgs.da.r))^2)^2)-3 #kurtosi
qqnorm(sgs.da.r,main="")
qqline(sgs.da.r)
shapiro.test(sgs.da.r)
layout(1)
  #decomposizione moltiplicativa
sgs.dm=decompose(sgs,type="multiplicative")
plot(sgs.dm)
layout(1:2) #analisi dei residui escluso 2020
sgs.dm.r=mean(sgs.dm$trend,na.rm=T)*(c(na.omit(window(sgs.dm$random,end=c(2019,12))),na.omit(window(sgs.dm$random,c(2021,1))))-1) #rumore portato sulla stessa scala del caso additivo
plot(sgs.dm.r,pch=20)
acf(sgs.dm.r,main="")
sd(acf(sgs.dm.r,plot=F)$acf)
layout(t(1:2)) #confronto gaussiana
hist(sgs.dm.r,20,freq=F,ylim=c(0,0.0045),main="")
lines(density(sgs.dm.r),col="blue")
lines(sort(sgs.dm.r),dnorm(sort(sgs.dm.r),mean(sgs.dm.r),sd(sgs.dm.r)),col="red")
mean((sgs.dm.r-mean(sgs.dm.r))^3)/(mean((sgs.dm.r-mean(sgs.dm.r))^2)^(3/2)) #skewness
mean((sgs.dm.r-mean(sgs.dm.r))^4)/(mean((sgs.dm.r-mean(sgs.dm.r))^2)^2)-3 #kurtosi
qqnorm(sgs.dm.r,main="")
qqline(sgs.dm.r)
shapiro.test(sgs.dm.r)
layout(1)
    #analisi
  #metodi di smorzamento esponenziale (Holt-Winters)
#coefficienti ottimizzati automaticamente
sgs.hwm=HoltWinters(sgs,seasonal="multiplicative")
plot(sgs.hwm,main="")
print(paste("alpha =",round(sgs.hwm$alpha,3),"beta =",sgs.hwm$beta,"gamma =",round(sgs.hwm$gamma,3)))
layout(t(1:2)) #analisi residui
sgs.hwm.r=resid(sgs.hwm)
plot(sgs.hwm.r,type="p",pch=20)
acf(sgs.hwm.r,main="")
layout(1)
var(sgs.hwm.r)/var(window(sgs,start(sgs.hwm.r),end(sgs.hwm.r))) #varianza non spiegata
#scelta coefficienti con grid-search
nt=15 #numero di test set
ft=12 #unità di tempo nel futuro su cui valutare la previsione
min=c(0,0,0)
err_gsm=1e+10
err_hwm=0
a=c(sgs.hwm$alpha,1:9/10) #includiamo anche i coefficienti ottimizzati automaticamente
b=c(sgs.hwm$beta,1:9/10)
c=c(sgs.hwm$gamma,1:9/10)
for(i in a){
  for(j in b){
    for(k in c){
      err=rep(0,nt)
      for(l in (n-nt-ft):(n-1-ft)){
        train=window(sgs,idt,ts_data(idt,pdt,l))
        test=window(sgs,ts_data(idt,pdt,l+1),ts_data(idt,pdt,l+ft))
        train.hw=HoltWinters(train,alpha=i,beta=j,gamma=k,seasonal='multiplicative')
        err[l-(n-nt-ft)+1]=mean((as.numeric(test)-as.numeric(predict(train.hw,ft)))^2)
      }
      if(err_gsm>sum(err)){
        err_gsm=sum(err)
        min[1]=i
        min[2]=j
        min[3]=k
      }
      if(i==sgs.hwm$alpha && j==sgs.hwm$beta && k==sgs.hwm$gamma){
        err_hwm=sum(err)
      }
    }
  }
}
print(paste("alpha =",min[1],"beta =",min[2],"gamma =",min[3]))
sgs.gsm=HoltWinters(sgs,alpha=min[1],beta=min[2],gamma=min[3],seasonal='multiplicative')
plot(sgs.gsm,main="")
print(paste("Modello Holt-Winters automatico - errore:",round(err_hwm)))
print(paste("Modello Holt-Winters grid-search - errore:",round(err_gsm)))
layout(t(1:2)) #analisi residui
sgs.gsm.r=resid(sgs.gsm)
plot(sgs.gsm.r,type="p",pch=20)
acf(sgs.gsm.r,main="")
layout(1)
var(sgs.gsm.r)/var(window(sgs,start(sgs.gsm.r),end(sgs.gsm.r))) #varianza non spiegata
  #metodi autoregressivi
#metodo diretto ridotto
pacf(sgs,main="") #l'ultimo valore di lag rilevante è 13
lg=13	#creazione modello
msgs=matrix(nrow=n-lg,ncol=lg+1)
for(i in 1:(lg+1)) {
  msgs[,i]=sgs[i:(n-lg-1+i)]
}
msgs<-data.frame(msgs)
msgs.lm<-lm(X14~.,data=msgs) #riduzione modello
summary(msgs.lm) #R^2=0.9691, X12 p-value=0.44569
msgs.lm<-lm(X14~.-X12,data=msgs)
summary(msgs.lm) #R^2=0.969, X6 p-value=0.57174
msgs.lm<-lm(X14~.-X12-X6,data=msgs)  
summary(msgs.lm) #R^2=0.969, X5 p-value=0.19685
msgs.lm<-lm(X14~.-X12-X6-X5,data=msgs)  
summary(msgs.lm) #R^2=0.9689, X7 p-value=0.09986
msgs.lm<-lm(X14~.-X12-X6-X5-X7,data=msgs)  
summary(msgs.lm) #R^2=0.9686, X4 p-value=0.0581
msgs.lm<-lm(X14~.-X12-X6-X5-X7-X4,data=msgs)  
summary(msgs.lm) #R^2=0.9683, X8 p-value=0.1016
msgs.lm<-lm(X14~.-X12-X6-X5-X7-X4-X8,data=msgs)  
summary(msgs.lm) #R^2=0.9681, X10 p-value=0.1766
msgs.lm<-lm(X14~.-X12-X6-X5-X7-X4-X8-X10,data=msgs)  
summary(msgs.lm) #R^2=0.9679, X3 p-value=0.05690
msgs.lm<-lm(X14~.-X12-X6-X5-X7-X4-X8-X10-X3,data=msgs)  
summary(msgs.lm) #R^2=0.9676, X9 p-value=0.176
msgs.lm<-lm(X14~.-X12-X6-X5-X7-X4-X8-X10-X3-X9,data=msgs)  
summary(msgs.lm) #R^2=0.9675, X11 p-value=0.161
msgs.lm<-lm(X14~.-X12-X6-X5-X7-X4-X8-X10-X3-X9-X11,data=msgs)  
summary(msgs.lm) #R^2=0.967
sgs.lm=window(sgs,ts_data(idt,pdt,lg))-resid(msgs.lm)
ts.plot(sgs.lm,sgs,col=c("red","black"))
layout(t(1:3)) #analisi residui
sgs.lm.r=ts(resid(msgs.lm),frequency=12,start=ts_data(idt,pdt,lg))
plot(sgs.lm.r,type="p",pch=20)
acf(sgs.lm.r,main="")
pacf(sgs.lm.r,main="")
layout(1)
var(sgs.lm.r)/var(window(sgs,ts_data(idt,pdt,lg))) #varianza non spiegata
#OLS
ols=ar(sgs,method="ols")
ols$order #25 lag
sgs.ols=window(sgs,ts_data(idt,pdt,ols$order))-na.omit(ols$resid)
ts.plot(sgs.ols,sgs,col=c("red","black"))
layout(t(1:3)) #analisi residui
sgs.ols.r=na.omit(ols$resid)
plot(sgs.ols.r,type="p",pch=20)
acf(sgs.ols.r,main="")
pacf(sgs.ols.r,30,main="")
layout(1)
var(sgs.ols.r)/var(window(sgs,ts_data(idt,pdt,ols$order))) #varianza non spiegata
#confronto modelli autoregressivi per autovalidazione
err_lm=rep(0,nt)
err_ols=rep(0,nt)
for(l in (n-nt-ft):(n-1-ft)){
  train=window(sgs,idt,ts_data(idt,pdt,l))
  test=window(sgs,ts_data(idt,pdt,l+1),ts_data(idt,pdt,l+ft))
  #metodo diretto ridotto
  L=length(train)
  mtrain=matrix(nrow=L-lg,ncol=lg+1)
  for(i in 1:(lg+1)){
    mtrain[,i]=train[i:(L-lg-1+i)]
  }
  mtrain<-data.frame(mtrain)
  train.lm<-lm(X14~X1+X2+X13,data=mtrain)
  train.lm.p=rep(0,L+ft)
  train.lm.p[1:L]=train
  for(i in 1:ft){
    train.lm.p[L+i]=coef(train.lm)%*%c(1,train.lm.p[L+i-13],train.lm.p[L+i-12],train.lm.p[L+i-1])
  }
  err_lm[l-(n-nt-ft)+1]=mean((as.numeric(test)-train.lm.p[(L+1):(L+ft)])^2)
  #metodo OLS
  train.ols=ar(train,method="ols")
  err_ols[l-(n-nt-ft)+1]=mean((as.numeric(test)-as.numeric(predict(train.ols,n.ahead=ft,se.fit=F)))^2)
}
print(paste("Modello autoregressivo ridotto - errore:",round(sum(err_lm))))
print(paste("Modello autoregressivo OLS - errore:",round(sum(err_ols))))
    #previsione finale con HW con coefficienti grid-search
sgs.gsm.p=predict(sgs.gsm,12)
ts.plot(sgs,sgs.gsm.p,col=c("black","red"))
lines(sgs.gsm.p+quantile(sgs.gsm.r,0.025),col="green4") #intervalli non parametrici
lines(sgs.gsm.p+quantile(sgs.gsm.r,0.975),col="green4")
  #zoom del periodo finale
wsgs=window(sgs,start=c(2023,10))
ts.plot(wsgs,sgs.gsm.p,ylim=c(3200,7300),col=c("black","red"))
lines(sgs.gsm.p+quantile(sgs.gsm.r,0.025),col="green4")
lines(sgs.gsm.p+quantile(sgs.gsm.r,0.975),col="green4")