###################################################################
make.anboot<-function(csdata,cldata,species,catchCat,strataDesc,lower,upper,method,nboot=null){
  
est.list <- dbeObject(desc = "Simulated Landings", species = species, catchCat =catchCat, strataDesc =strataDesc ,
                                methodDesc =method)
if(method=="analytical")est.list<- RaiseLgth(est.list,csDataCons(csDataVal(csdata)),clDataCons((clDataVal(cldata))))
if(method=="bootstrap")est.list<-RaiseLgthBoot(est.list,csDataCons(csDataVal(csdata)),clDataCons((clDataVal(cldata))),B=nboot)
est.list<- dbeCalc(est.list, type = "CI", probs = c(lower,upper), update = TRUE)
est.list<- dbeCalc(est.list, type = "CV", update =TRUE)
if(method=="analytical")est.list<- RaiseAge(est.list,csDataCons(csDataVal(csdata)) , type="fixed")
if(method=="bootstrap")est.list<- RaiseAgeBoot(est.list,csDataCons(csDataVal(csdata)) , type="fixed")
est.list<- dbeCalc(est.list, type = "CI", vrbl = "a", probs = c(lower,upper), update = TRUE)
est.list<- dbeCalc(est.list, type = "CV", vrbl = "a", update = TRUE)
list<-est.list@ageNum$ci
yr.list<- aggregate(list(value =list$value, inf = list$inf, sup =list$sup),list(age = list$age), sum)

yr.list
}
###################################################################

anboot.summary<-function(est,minAge,maxAge){
  nages<-1+maxAge-minAge
a<-as.numeric(as.character(est$age))+1-minAge
x<-l<-u<-rep(NA,nages)
x[a]<-est$value
l[a]<-est$inf
u[a]<-est$sup
list(mean=x,lower=l,upper=u)
}
###################################################################
bay.summary<-function(est,lower,upper,mult=1){
caa<-mult*apply(est$totcatch.land,c(2,3),sum)
q<-NULL
for(iage in 1:nrow(caa))q<-rbind(q,quantile(caa[iage,],c(lower,upper)))
list(mean=rowMeans(caa),lower=q[,1],upper=q[,2])
}
###################################################################

make.sum.list<-function(meanmat,upper,lower,include,true){
mean<-rowMeans(meanmat,na.rm=T)
bias<-rowMeans(meanmat-true,na.rm=T)
mse<-sqrt(rowMeans((meanmat-true)^2,na.rm=T))
meanwidth<-rowMeans(upper-lower,na.rm=T)
include<-rowMeans(include,na.rm=T)
list(overall.mean=mean,mean=meanmat,lower=lower,upper=upper,bias=bias,mse=mse,meanwidth=meanwidth,include=include)
}
###################################################################

plot.sim.intervals<-function(true,an.out,bo.out,bay.out,ageMin,ageMax,iplot){
ymax<-max(c(true,an.out$upper[,iplot],bo.out$upper[,iplot],bay.out$upper[,iplot]),na.rm=T)
plot(ageMin:ageMax,true,pch = 16, col = 2,ylim=c(0,ymax),xlab="age",ylab="numbers")
error.bar(ageMin:ageMax-0.1,an.out$mean[,iplot],lower=an.out$lower[,iplot],upper=an.out$upper[,iplot],lwd=1,incr=F,
   log='',gap=F,add=T)
points(ageMin:ageMax-0.1,an.out$mean[,iplot],col=3,pch=3)
error.bar(ageMin:ageMax+0.1,bo.out$mean[,iplot],lower=bo.out$lower[,iplot],upper=bo.out$upper[,iplot],lwd=1,incr=F,
   log='',gap=F,add=T)
points(ageMin:ageMax+0.1,bo.out$mean[,iplot],pch=4,col=4)
error.bar(ageMin:ageMax+0.3,bay.out$mean[,iplot],lower=bay.out$lower[,iplot],upper=bay.out$upper[,iplot],lwd=1,incr=F,
   log='',gap=F,add=T)
points(ageMin:ageMax+0.3,bay.out$mean[,iplot],pch=5,col=5)
legend("topright",c("true","analytical","bootstrap","Bayesian"),pch=c(16,3,4,5),col=c(2,3,4,5))
title("Interval estimation by different methods")
}
###################################################################

method.plot4<-function(true,v1,v2,v3,ageMin,ageMax,xlab,ylab,title){
  ymax<-max(c(true,v1,v2,v3))
plot(ageMin:ageMax,true,pch = 16, col = 2,ylim=c(0,ymax),xlab=xlab,ylab=ylab)
points(ageMin:ageMax,rowMeans(an.meancaa),col=3,pch=3)
points(ageMin:ageMax,rowMeans(bo.meancaa),col=4,pch=4)
points(ageMin:ageMax,rowMeans(bay.meancaa),col=5,pch=5)
legend("topright",c("true","analytical","bootstrap","Bayesian"),pch=c(16,3,4,5),col=c(2,3,4,5))
title(title)
}
###################################################################

method.plot3<-function(v1,v2,v3,ageMin,ageMax,xlab,ylab,title){
  ymax<-max(c(v1,v2,v3),na.rm=T)
  ymin<-min(c(0,v1,v2,v3),na.rm=T)
plot(ageMin:ageMax,v1,ylim=c(ymin,ymax),xlab=xlab,ylab=ylab)
points(ageMin:ageMax,v2,col=4,pch=2)
points(ageMin:ageMax,v3,col=2,pch=3)
legend("topright",c("bootstrap","analytical","Bayesian"),pch=1:3,col=c(1,4,2))
title(title)
}
###################################################################

make.sim.data<-function(params,setup.data,nlsamp.land,nasamp.land,nlsamp.disc,nasamp.disc,length.list){
ml.data<-obs.data<-NULL
sim.data<-cost.simulate(params,setup.data,nlsamp.land,nasamp.land,nlsamp.disc,nasamp.disc,length.list)
sim.data<<-sim.data
setup.data$structure$n.size<-sim.data$nsize
setup.data$structure$n.haul<-sim.data$nhaul
ml<-setup.data$structure$dis==0
ml.structure<-obs.structure<-setup.data$structure
ml.use<-(1:length(sim.data$sample.land))[ml]
obs.use<-(1:length(sim.data$sample.land))[!ml]
for(i in 1:length(setup.data$structure)){
	ml.structure[[i]]<- setup.data$structure[[i]][ml]
	obs.structure[[i]]<- setup.data$structure[[i]][!ml]}
if(sum(setup.data$structure$dis==0)>0){
ml.data<-sim.to.cost(sim.data$sample.land[ml.use],sim.data$sample.disc[ml.use],
	ml.structure, setup.data$species,"ML")}
if(sum(setup.data$structure$dis==1)>0){
obs.data<-sim.to.cost(sim.data$sample.land[obs.use],sim.data$sample.disc[obs.use],
	obs.structure, setup.data$species,"OBS")}
COSTml<-COSTobs<-new("csData")
COSTml@tr<-as.data.frame(ml.data$tr)
COSTml@hl<-as.data.frame(ml.data$hl)
COSTml@sl<-as.data.frame(ml.data$sl)
COSTml@ca<-as.data.frame(ml.data$ca)
COSTml@hh<-as.data.frame(ml.data$hh)
COSTobs@tr<-as.data.frame(obs.data$tr)
COSTobs@hl<-as.data.frame(obs.data$hl)
COSTobs@sl<-as.data.frame(obs.data$sl)
COSTobs@ca<-as.data.frame(obs.data$ca)
COSTobs@hh<-as.data.frame(obs.data$hh)


list(obs = COSTobs, ml = COSTml,paa.land=sim.data$paa.land)
}

###################################################################

make.predict.cov<-function(seas.list,gear.list,area.list){
  x<-expand.grid(seas.list,gear.list,area.list)
  seas<-as.integer(x[,1])
  gear<-as.integer(x[,2])
  area<-as.integer(x[,3])
  list(seas=seas,area=area,gear=gear)
}
  
#######################################################################
new.make.structure<-function(sampmat,seaslist,gearlist,arealist){
if(is.null(seaslist))seaslist<-1
if(is.null(gearlist))gearlist<-1
if(is.null(arealist))arealist<-1

  landratio<-1.93
trip<-sampmat[,4]+sampmat[,5]
ntrip<-sum(trip)
ncell<-nrow(sampmat)
ncelluse<-sum(trip>0)
rowno<-c(rep(1:ncell,sampmat[,4]),rep(1:ncell,sampmat[,5]))
newmat<-sampmat[rowno,]
cell<-sampmat[rowno,6]
obs.cell<-sampmat[rowno,7]
discard<-rep(0:1,colSums(sampmat[,4:5]))
seas<-as.factor(seaslist[newmat[,1]])
gear<-as.factor(gearlist[newmat[,2]])
area<-as.factor(arealist[newmat[,3]])

n.land<-rnorm(length(gear),90,70)
n.land<-n.land*n.land

is<-rep(1,length(n.land))
div<-c(4,9)
for(i in 1:2)is<-is+(log(n.land)>div[i])
n.size<-n.haul<-NULL
for(i in 1:length(n.land)){
  if(is[i]==1){n.haul<-c(n.haul,sample(c(2,5),1,prob=c(0.5,0.5)))
            n.size<-c(n.size,1) }
  if(is[i]==2){n.haul<-c(n.haul,sample(c(2,3,4,5,6,7,10,11,12,13,14,15,17,18,19,21),1,
             prob=c(0.08,0.12,0.04,0.04,0.04,0.04,0.08,0.08,0.08,0.04,0.04,0.12,0.08,0.04,0.04,0.04)))
	n.size<-c(n.size,sample(1:6,1,prob=c(0.2684825,0.307393,0.2140078,0.1322957,0.07003891,0.007782101)))
       }
  if(is[i]==3){n.haul<-c(n.haul,sample(c(4,5,8,9,11,12,14,16,20,21,23,26,36,37),1,
         prob=c(0.05,0.09,0.09,0.05,0.09,0.05,0.05,0.14,0.09,0.05,0.14,0.05,0.05,0.05)))
         n.size<-c(n.size,sample(1:7,1,prob=c(0.009345794,0.04361371,0.2834891,0.2523364,0.2834891,0.1152648,0.01246106)))      
	}}
n.land[n.land<1000]<-1000
list(seas=seas,gear=gear,area=area,cell=cell,obs.cell=obs.cell,
	dis=discard,n.size=n.size,
	n.land=n.land,n.haul=n.haul,landratio=landratio)
}
###################################################################
new.setup<-function(params,seaslist=NULL,arealist=NULL,gearlist=NULL,species,age.covariates,weight.covariates,
                    nmland=20,nobs=20,ageMin,ageMax){
if(is.null(seaslist))seaslist<-1
if(is.null(arealist))arealist<-1
if(is.null(gearlist))gearlist<-1
  models<-make.reg.models(age.covariates,weight.covariates)

sampmat<-get.sampled.cells(params,seaslist,gearlist,arealist,nmland=nmland,nobs=nobs,
	force=NULL)
structure<-new.make.structure(sampmat,seaslist,gearlist,arealist)

list(structure=structure,agemodel=models$agemodel,lgamodel=models$lgamodel,wglmodel=models$wglmodel,
	seaslist=seaslist,arealist=arealist,gearlist=gearlist,
       ageMin=ageMin,ageMax=ageMax,species=species)
}
###################################################################
getparams<-function(newparams,CAA,simno){

Const<-newparams$age$Int$eff$Const[,simno]
HZ<-newparams$age$Hsz$eff$Const[,simno]
seas<-newparams$age$Int$eff$seas[,,simno]
gear<-newparams$age$Int$eff$gear[,,simno]
area<-newparams$age$Int$eff$area[,,simno]
if(!is.null(dim(newparams$age$Int$tau)))tau<-newparams$age$Int$tau[,simno]
else tau<-newparams$age$Int$tau
nseas<-dim(newparams$age$Int$eff$seas)[2]
p<-0

for(i in 1:nseas){
x<-Const+newparams$age$Int$eff$seas[,i,simno]
px<-t(exp(x))
px<-t(px/rowSums(px))
p<-p+px
}
p<-p/nseas
if(dim(p)[2]==1)p<-as.vector(p)
ageparams<-list(p=p,Const=Const,HZ=HZ,seas=seas,area=area,gear=gear,tau=tau)

Const<-newparams$lga$Int$eff$Const[simno]
HZ<-newparams$lga$Hsz$eff$Const[simno]
if(!is.null(dim(newparams$lga$Int$eff$seas)))seas<-newparams$lga$Int$eff$seas[,simno]
else seas<-newparams$lga$Int$eff$seas
if(!is.null(dim(newparams$lga$Int$eff$gear)))gear<-newparams$lga$Int$eff$gear[,simno]
else gear<-newparams$lga$Int$eff$gear
if(!is.null(dim(newparams$lga$Int$eff$area)))area<-newparams$lga$Int$eff$area[,simno]
else area<-newparams$lga$Int$eff$area
if(!is.null(dim(newparams$lga$Int$tau)))tau<-newparams$lga$Int$tau[,simno]
else tau<-newparams$lga$Int$tau
tau.obs<-newparams$lga$tau.obs[simno]
slope<-newparams$lga$Slp$eff$Const[,simno]
lgaparams<-list(Const=Const,HZ=HZ,seas=seas,area=area,gear=gear,tau=tau,tau.obs=tau.obs,slope=slope)

Const<-newparams$wgl$Int$eff$Const[simno]
HZ<-newparams$wgl$Hsz$eff$Const[simno]
if(!is.null(dim(newparams$wgl$Int$eff$seas)))seas<-newparams$wgl$Int$eff$seas[,simno]
else seas<-newparams$wgl$Int$eff$seas
if(!is.null(dim(newparams$wgl$Int$eff$gear)))gear<-newparams$wgl$Int$eff$gear[,simno]
else gear<-newparams$wgl$Int$eff$gear
if(!is.null(dim(newparams$wgl$Int$eff$area)))area<-newparams$wgl$Int$eff$area[,simno]
else area<-newparams$wgl$Int$eff$area
if(!is.null(dim(newparams$wgl$Int$tau)))tau<-newparams$wgl$Int$tau[,simno]
else tau<-newparams$wgl$Int$tau
tau.obs<-newparams$wgl$tau.obs[simno]
slope<-newparams$wgl$Slp$eff$Const[,simno]
wglparams<-list(Const=Const,HZ=HZ,seas=seas,area=area,gear=gear,tau=tau,tau.obs=tau.obs,slope=slope)

truecaa.land<-wga.land<-truecaa.disc<-wga.disc<-planded<-lga.land<-lga.disc<-NULL

if(length(simno)==1){wga.land<-CAA$mean.wga.land[,simno]
             truecaa.land<-colSums(CAA$totcatch.land[,,simno])
             if(!is.null(CAA$mean.wga.disc))wga.disc<-CAA$mean.wga.disc[,simno]
             if(!is.null(CAA$totcatch.disc))truecaa.disc<-colSums(CAA$totcatch.disc[,,simno])
             planded<-CAA$planded[,simno]                     
                   }
else {wga.land<-CAA$mean.wga.land[,simno]
      truecaa.land<-apply(CAA$totcatch.land[,,simno],c(2,3),sum)
      if(!is.null(CAA$mean.wga.disc))wga.disc<-CAA$mean.wga.disc[,simno]
      if(!is.null(CAA$totcatch.disc))truecaa.disc<-apply(CAA$totcatch.disc[,,simno],c(2,3),sum)
      planded<-CAA$planded[,simno]
    }
lga.land<-CAA$mean.lga.land
lga.disc<-CAA$mean.lga.disc
caaparams<-list(wga.land=wga.land,CAA.land=truecaa.land,wga.disc=wga.disc,
              CAA.disc=truecaa.disc,planded=planded,lga.land=lga.land,lga.disc=lga.disc)

list(age=ageparams,lga=lgaparams,wgl=wglparams,caa=caaparams)
}
###################################################################
jointrue<-function(alltrue,true){
alltrue$age$p<-cbind(alltrue$age$p,true$age$p)
alltrue$age$Const<-cbind(alltrue$age$Const,true$age$Const)
alltrue$age$HZ<-cbind(alltrue$age$HZ,true$age$HZ)
alltrue$age$seas<-makearray(alltrue$age$seas,true$age$seas)
alltrue$age$gear<-makearray(alltrue$age$gear,true$age$gear)
alltrue$age$area<-makearray(alltrue$age$area,true$age$area)
alltrue$age$tau<-makearray(alltrue$age$tau,true$age$tau)

alltrue$lga$Const<-c(alltrue$lga$Const,true$lga$Const)
alltrue$lga$HZ<-c(alltrue$lga$HZ,true$lga$HZ)
alltrue$lga$seas<-cbind(alltrue$lga$seas,true$lga$seas)
alltrue$lga$gear<-cbind(alltrue$lga$gear,true$lga$gear)
alltrue$lga$area<-cbind(alltrue$lga$area,true$lga$area)
alltrue$lga$tau<-cbind(alltrue$lga$tau,true$lga$tau)
alltrue$lga$tau.obs<-c(alltrue$lga$tau.obs,true$lga$tau.obs)
alltrue$lga$slope<-c(alltrue$lga$slope,true$lga$slope)

alltrue$wgl$Const<-c(alltrue$wgl$Const,true$wgl$Const)
alltrue$wgl$HZ<-c(alltrue$wgl$HZ,true$wgl$HZ)
alltrue$wgl$seas<-cbind(alltrue$wgl$seas,true$wgl$seas)
alltrue$wgl$gear<-cbind(alltrue$wgl$gear,true$wgl$gear)
alltrue$wgl$area<-cbind(alltrue$wgl$area,true$wgl$area)
alltrue$wgl$tau<-cbind(alltrue$wgl$tau,true$wgl$tau)
alltrue$wgl$tau.obs<-c(alltrue$wgl$tau.obs,true$wgl$tau.obs)
alltrue$wgl$slope<-c(alltrue$wgl$slope,true$wgl$slope)

alltrue$caa$wga.land<-cbind(alltrue$caa$wga.land,true$caa$wga.land)
alltrue$caa$CAA.land<-cbind(alltrue$caa$CAA.land,true$caa$CAA.land)
alltrue$caa$wga.disc<-cbind(alltrue$caa$wga.disc,true$caa$wga.disc)
alltrue$caa$CAA.disc<-cbind(alltrue$caa$CAA.disc,true$caa$CAA.disc)
alltrue$caa$planded<-cbind(alltrue$caa$planded,true$caa$planded)
alltrue$caa$lga.land<-cbind(alltrue$caa$lga.land,true$caa$lga.land)
alltrue$caa$lga.disc<-cbind(alltrue$caa$lga.disc,true$caa$lga.disc)

alltrue
}
###################################################################
makearray<-function(a,m){
  n<-length(c(a,m))
  d1<-nrow(m)
  d2<-ncol(m)
  d3<-n/(d1*d2)
  out<-array(c(a,m),dim=c(d1,d2,d3))
  out
}
###################################################################
get.obs.cell<-function(sampmat,params){
model<-params$age$model$Int
  mat<-NULL
  if(model$seas)mat<-cbind(mat,sampmat[,1])
  if(model$gear)mat<-cbind(mat,sampmat[,2])
  if(model$area)mat<-cbind(mat,sampmat[,3])
  ncov<-model$seas+model$gear+model$area
  refmat<-unique(params$age$data$mcov$Int$cov[,2:(ncov+2)])
  obs.cell<-rep(0,nrow(mat))
  for(i in 1:length(obs.cell)){for(j in 1:nrow(refmat)){
    if(sum(mat[i,]==refmat[j,1:ncov])==ncov){obs.cell[i]<-1
    break}}}
  obs.cell

}
###################################################################
###################################################################
###################################################################
###################################################################
