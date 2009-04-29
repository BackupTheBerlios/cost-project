
#######################################################################
cost.simloop <- function(params,setup.data, burnin,nmcmc,l.int,Int,Slp,
                     landings,nHaul,nseas, fit = TRUE){
                     
  sim.data <- cost.simulate(params,setup.data,simno= 1)

    setup.data$structure$n.size<-sim.data$nsize
    setup.data$structure$n.haul<-sim.data$nhaul
    ml<-setup.data$structure$dis==0
    ml.structure<-obs.structure<-setup.data$structure
    ml.use<-(1:length(sim.data$sample.land))[ml]
    obs.use<-(1:length(sim.data$sample.land))[!ml]
    
    for(i in 1:length(setup.data$structure)){
	   ml.structure[[i]]<- setup.data$structure[[i]][ml]
       obs.structure[[i]]<- setup.data$structure[[i]][!ml]}
    
    ml.data <- sim.to.cost(sim.data$sample.land[ml.use],sim.data$sample.disc[ml.use],
                    ml.structure, setup.data$species,"ML")
    obs.data <- sim.to.cost(sim.data$sample.land[obs.use],sim.data$sample.disc[obs.use],
                    obs.structure, setup.data$species,"OBS")
    
    tempdata<<-sim.data$sample.land[ml.use]
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
     ml.data<<-ml.data
    obs.data<<-obs.data
    input.data<-read.cost.data(COSTobs,COSTml, setup.data$species)
    if(nseas==4){input.data$season_obs<-1+floor((input.data$season_obs-1)/3)    
    input.data$season_mland<-1+floor((input.data$season_mland-1)/3)
    input.data$seas<-1+floor((input.data$seas-1)/3) }
  
    input.data$sampsize_disc[is.na(input.data$sampsize_disc)]<-0
    input.data$haulsize_disc[is.na(input.data$haulsize_disc)]<-0
    input.data$cens.mu=c(1,30.0,3.4)
    input.data$cens.tau=c(1.0,1.0,1.0)
    input.data$cens.pri=c(3.4,0.001,0.0001,0.001)
  
    simulation.data<<-input.data

    aobs = NULL
    acov = list(year=input.data$year,seas=input.data$seas,
	   gear=input.data$gear,area=input.data$area)
        neigh = list(num=c(1,1),adj=c(2,1)) 

    simobj=list(aobs=aobs,acov=acov,neigh=neigh)
    class(simobj) <- "caa.data"
    dump(c("simobj","input.data","setup.data","nHaul","Int","Slp"))
    mbe.fit <- NULL
    if(fit){
    fit = cost.fit(simobj,input.data,burnin=burnin,numit.inner=1,numit.outer=nmcmc,constr=1,seed=213421,
               ageMin= setup.data$ageMin,ageMax= setup.data$ageMax,nSeason=nseas,
    agemodel= setup.data$agemodel,
               lgamodel= setup.data$lgamodel,lgarel="log-linear",model1=T,model2=F)
    fit$aobs=list(data=list(nBoats=nHaul))
    res2<-insert.wgl.param(fit,Int,Slp) 

    nage<-1+setup.data$ageMax-setup.data$ageMin
    s<-as.integer(setup.data$structure$seas)
    t.year <-t.gear <-t.area <-rep(1,nseas)
    t.seas <-1:nseas
    print(dim(fit$lga$Int$eff$seas))

    res.pred <-  predict.fit.COST(res2,fit$COST.list,t.year,t.seas,t.gear,t.area,landings,1:nseas,
          t2.year=NULL,t2.seas=NULL,t2.gear=NULL,t2.area=NULL,
          burnin=0,nMC=100,l.int=l.int,
          par.haulsize=NULL)
    
        mbe.fit <- list(fit=fit,caa=res.pred)
    }
    #  fit<-res.pred<-NULL

    return(list(obs = COSTobs, ml = COSTml, input.data = input.data, mbe.fit = mbe.fit))
}
#######################################################################
cost.simulate<-function(params, setup.data,simno){
seas<-code.cov(setup.data$structure$seas, setup.data$seaslist)
gear<-code.cov(setup.data$structure$gear, setup.data$gearlist)
area<-code.cov(setup.data$structure$area, setup.data$arealist)
realseas<-setup.data$structure$seas
ageMin<-setup.data$ageMin
ageMax<- setup.data$ageMax
simno<-1
ntrip<-length(seas)
cell.effects<-make.cost.cell(setup.data$structure$cell,params,ageMin:ageMax,simno)
cell.effects<<-cell.effects
p<-cost.make.cell.p(params,cell.effects,seas,gear,area,
	setup.data$structure$cell,simno, setup.data$species)
p<<-p
cost.cell.lga<-cost.make.cell.lga(p, setup.data$structure,realseas)
cost.cell.lga<<-cost.cell.lga
nlsamp.land<-1000
nasamp.land<-2
nlsamp.disc<-30
nasamp.disc<-2
sample.land<-sample.disc<-vector("list",length=ntrip)
nsize<-nhaul<-NULL
for(itrip in 1:ntrip){
#for(itrip in 1){
simfish<-cost.sim.fish(itrip,
        setup.data$structure$landratio*setup.data$structure$n.land[itrip],p,
	cost.cell.lga,params,itrip,simno)
if(sum(simfish$laland)>0){
sample.land[[itrip]]<-make.size.sample(simfish$laland,simfish$l,
	simfish$wglSlp,simfish$wglInt,simfish$wglsd,ageMin:ageMax,
	setup.data$structure$n.size[itrip],nlsamp.land,nasamp.land)
nsize<-c(nsize,sample.land[[itrip]]$nsize)}
else nsize<-c(nsize,0)
if(sum(simfish$ladisc)>0){ 
sample.disc[[itrip]]<-make.disc.sample(simfish$ladisc,simfish$l,
	simfish$wglSlp,simfish$wglInt,simfish$wglsd,ageMin:ageMax,
	setup.data$sim.models$phaul,nlsamp.disc,nasamp.disc,
        setup.data$div)
nhaul<-c(nhaul,sample.disc[[itrip]]$nhaul)}
else nhaul<-c(nhaul,0)
}
sample.land<<-sample.land
sample.disc<<-sample.disc
list(sample.land=sample.land,sample.disc=sample.disc,nsize=nsize,nhaul=nhaul)
}
#######################################################################
get.sampled.cells<-function(seasons,gears,areas,ntrip=0,ndisc=0,force=NULL){
nseason<-length(seasons)
ngear<-length(gears)
narea<-length(areas)
ncell<-narea*ngear*nseason
sampmat<-matrix(0,nrow=ncell,ncol=5)
sampmat[,1]<-rep(1:nseason,times=ngear*narea)
sampmat[,2]<-rep(1:ngear,each=nseason,times=narea)
sampmat[,3]<-rep(1:narea,each=ngear*nseason)
useas<-ugear<-uarea<-NULL
in.force<-F
if(!is.null(force)){
for(i in 1:nrow(force)){
use<-sampmat[,1]==force[i,1]&sampmat[,2]==force[i,2]&sampmat[,3]==force[i,3]
sampmat[use,4:5]<-force[i,4:5]
in.force<-in.force|use
}

useas<-unique(force[,1])
ugear<-unique(force[,2])
uarea<-unique(force[,3])
}
seaslist<-setdiff(1:nseason,useas)
gearlist<-setdiff(1:ngear,ugear)
arealist<-setdiff(1:narea,uarea)
n<-max(c(length(seaslist),length(gearlist),length(arealist)))

v1<-c(sample(seaslist),sample(1:nseason,n-length(seaslist),replace=T))
v2<-c(sample(gearlist),sample(1:ngear,n-length(gearlist),replace=T))
v3<-c(sample(arealist,2),sample(1:narea,n-2,replace=T))
v1<-sample(v1)
v2<-sample(v2)
v3<-sample(v3)
for(i in 1:n){suse<-sampmat[,1]==v1[i]&sampmat[,2]==v2[i]&
	sampmat[,3]==v3[i]
sampmat[suse,4]<-max(1,sampmat[suse,4])}
nalloc<-sum(sampmat[,4])
if(ntrip>nalloc){
extra<-sample((1:ncell)[!in.force],ntrip-nalloc,replace=T)
for(e in extra)sampmat[e,4]<-sampmat[e,4]+1}
nsamp<-sum(sampmat[,4])
if(nsamp>ntrip)cat("ntrip increased to",nsamp,"\n")

ndisc.rem<-ndisc-sum(sampmat[in.force,5])
nsamp.rem<-nsamp-sum(sampmat[in.force,4])
if(ndisc.rem>0){
for(i in sample((1:ncell)[!in.force])){
 pdisc<-ndisc.rem/nsamp.rem
 ns<-sampmat[i,4]
 if(ns>0){nd<-min(rbinom(ns,1,pdisc),ndisc.rem)
   ndisc.rem<-ndisc.rem-nd
   nsamp.rem<-nsamp.rem-ns
   if(ndisc.rem>nsamp.rem){nd<-nd+ndisc.rem-nsamp.rem
                 ndisc.rem<-nsamp.rem}
   sampmat[i,5]<-nd}
 }
}
sampmat
}
#######################################################################
model.n.land<-function(season,gear,area,n.land){
lognland<-log(n.land)
seas<-as.factor(season)
gear<-as.factor(gear)
area<-as.factor(area)
seas<-seas[n.land>0]
gear<-gear[n.land>0]
area<-area[n.land>0]
lognland<-lognland[n.land>0]
#lm<-lm(lognland~seas+gear+area)
lm<-lm(lognland~gear)
lm
}
#######################################################################
model.n.haul<-function(nhaul,n.land,season,gear,area){
lognhaul<-log(nhaul)
seas<-as.factor(season)
gear<-as.factor(gear)
area<-as.factor(area)
seas<-seas[nhaul>0&n.land>0]
gear<-gear[nhaul>0&n.land>0]
area<-area[nhaul>0&n.land>0]
lognhaul<-lognhaul[nhaul>0&n.land>0]
n.land<-n.land[nhaul>0&n.land>0]
#lm<-lm(lognhaul~log(n.land)+seas+gear+area,na.action=na.omit)
lm<-lm(lognhaul~log(n.land),na.action=na.omit)

lm
}
#######################################################################
model.n.size<-function(nsize,size,n.land,div){
logn.land<-log(n.land)
logn.land<-logn.land[n.land>0]
size<-size[n.land>0]
psize<-matrix(0,nrow=length(div)+1,ncol=nsize)
for(i in 1:nsize){
 psize[1,i]<-mean(logn.land<div[1]&size==i,na.rm=T)
 for(j in 2:length(div)){
   psize[j,i]<-mean(logn.land>=div[j-1]&logn.land<div[j]&size==i,na.rm=T)}
 psize[length(div)+1,i]<-mean(logn.land>=div[length(div)]&size==i,na.rm=T)
   }
psize<-psize/rowSums(psize)
psize[is.na(psize)]<-0
psize
}
#######################################################################
get.n.size<-function(psize,totland,div){
n<-rep(NA,length(totland))
ndiv<-length(div)
cat<-1
for(i in 1:ndiv)cat<-cat+(log(totland)>div[i])
for(i in 1:(ndiv+1)){nsize<-sample(1:ncol(psize),sum(cat==i),T,psize[i,])
 n[cat==i]<-nsize}
n
}
#######################################################################
allocate.size<-function(l,lfreq,ns){
  if(length(l)==0)scat<-matrix(c(1,rep(0,ns-1)),ncol=1)
if(length(l)==1&ns>1)scat<-matrix(c(1,rep(0,ns-1)),ncol=1)
if(length(l)==1&ns==1)scat<-matrix(1,ncol=1)
if(ns<1)print("warning, ns<1")
if(length(l)>1){
minl<-min(l)
maxl<-max(l)
sd<-0.5*(maxl-minl)/ns
mid<-minl-sd+2*sd*(1:ns)
ps<-NULL
for(i in 1:ns)ps<-cbind(ps,dnorm(l,mid[i],sd))
if(is.matrix(ps))ps<-ps/rowSums(ps)
else ps<-ps/sum(ps)
scat<-NULL
for(i in 1:length(l)){
  if(is.matrix(ps))prob<-ps[i,]
  else  prob<-ps[i]
  scat<-cbind(scat,rmultinom(1,lfreq[i],prob))}}
t(scat)
}
#######################################################################
allocate.disc.haul<-function(nhauls,l,lfreq,nph.dist){
newdisc<-sample(nph.dist,nhauls,replace=T)
haul<-NULL
for(i in 1:length(l))haul<-cbind(haul,rmultinom(1,lfreq[i],newdisc))
haul
}
#######################################################################
make.structure<-function(sampmat,sim.models,seaslist,gearlist,arealist,model.data,nseas){
landratio<-(mean(model.data$n.land)+mean(model.data$n.disc))/mean(model.data$n.land)
ntrip<-sum(sampmat[,4])
ncell<-nrow(sampmat)
foNum<-paste("trip",1:ntrip,sep='')
ncell<-nrow(sampmat)
cell<-rep(0,ncell)
ncelluse<-sum(sampmat[,4]>0)
cell[sampmat[,4]>0]<-1:ncelluse
rowno<-rep(1:ncell,sampmat[,4])
newmat<-sampmat[rowno,]
cell<-cell[rowno]
d<-rep(1:ncell,sampmat[,5])
nod<-rep(1:ncell,sampmat[,4]-sampmat[,5])
dd<-rep(1:0,c(length(d),length(nod)))
order<-order(c(d,nod))
dd<-dd[order]
seas<-as.factor(seaslist[newmat[,1]])
gear<-as.factor(gearlist[newmat[,2]])
area<-as.factor(arealist[newmat[,3]])
dis<-dd
if(0){
temparea<-subvar(area,sim.models$lm.tot.land$xlevels$area)
tempgear<-subvar(gear,sim.models$lm.tot.land$xlevels$gear)
tempseas<-subvar(seas,sim.models$lm.tot.land$xlevels$seas)
tot.land<-predict(sim.models$lm.tot.land,
	data.frame(seas=tempseas,gear=tempgear,area=temparea))
sd<-sd(sim.models$lm.tot.land$res)
tot.land<-tot.land+rnorm(length(tot.land),0,sd)
}
temparea<-subvar(area,sim.models$lm.n.land$xlevels$area)
tempgear<-subvar(gear,sim.models$lm.n.land$xlevels$gear)
tempseas<-subvar(seas,sim.models$lm.n.land$xlevels$seas)
n.land<-rep(NA,length(tempgear))
for(g in as.character(tempgear)){
ng<-sum(tempgear==g)
n.land[tempgear==g]<-sample(model.data$n.land[model.data$gear==g&
	model.data$n.land>0],ng,replace=T)}
if(1){
#n.land<-predict(sim.models$lm.n.land,
#	data.frame(seas=tempseas,gear=tempgear,area=temparea))
#res<-sample(sim.models$lm.n.land$res,length(n.land),replace=T)
#n.land<-n.land+res

temparea<-subvar(area,sim.models$lm.n.haul$xlevels$area)
tempgear<-subvar(gear,sim.models$lm.n.haul$xlevels$gear)
tempseas<-subvar(seas,sim.models$lm.n.haul$xlevels$seas)
n.haul<-predict(sim.models$lm.n.haul,
	data.frame(n.land,seas=tempseas,gear=tempgear,area=temparea))
sd<-sd(sim.models$lm.n.haul$res)
n.haul<-n.haul+rnorm(length(n.haul),0,sd)
}

is<-rep(1,length(n.land))
for(i in 1:length(sim.models$div))is<-is+(log(n.land)>sim.models$div[i])
n.size<-n.haul<-NULL
nh<-ncol(sim.models$phaul)
ns<-ncol(sim.models$psize)
for(i in 1:length(n.land)){
	n.size<-c(n.size,sample(1:ns,1,prob=sim.models$psize[is[i],]))
	n.haul<-c(n.haul,sample(1:nh,1,prob=sim.models$phaul[is[i],]))}
n.land[n.land<1000]<-1000
list(seas=seas,gear=gear,area=area,cell=cell,
	dis=dis,n.size=n.size,
	n.land=n.land,n.haul=n.haul,landratio=landratio)
}
#######################################################################
get.sim.model.data<-function(datalist,species,usequarter=T){
tot.trip.land<-n.size.class<-n.hauls<-seas<-gear<-area<-n.land<-
	tot.disc<-n.disc.dist<-NULL
for(data in datalist){
ca<-data@ca[data@ca$spp==species,]
hl<-data@hl[data@hl$spp==species,]
hl.land<-hl[hl$catchCat=="LAN"&hl$spp==species,]
hl.disc<-hl[hl$catchCat=="DIS"&hl$spp==species,]
sl<-data@sl[data@sl$spp==species,]
sl.land<-sl[sl$catchCat=="LAN"&sl$spp==species,]
sl.disc<-sl[sl$catchCat=="DIS"&sl$spp==species,]
hh<-data@hh
tr<-data@tr

trips<-tr$trpCode
ntrip<-length(trips)

nsize<-as.integer(tapply(sl$commCat,
	factor(sl$trpCode,levels=trips),max,na.rm=T))
nsize[is.na(nsize)]<-0

nsamp<-tapply(hl.land$lenNum,
	list(factor(hl.land$trpCode,levels=trips),
	hl.land$commCat),sum,na.rm=T)
sampwt<-tapply(sl.land$subSampWt,
	list(factor(sl.land$trpCode,levels=trips),
	sl.land$commCat),sum,na.rm=T)
totwt<-tapply(sl.land$wt,
	list(factor(sl.land$trpCode,levels=trips),
	sl.land$commCat),sum,na.rm=T)
nland<-nsamp*totwt/sampwt
nland<-round(rowSums(nland,na.rm=T))

totland<-as.integer(tapply(sl.land$wt,
	factor(sl.land$trpCode,levels=trips),sum,na.rm=T))
totland[is.na(totland)]<-0

m.gear<-tapply(hh$foCatEu5,
	factor(hh$trpCode,levels=trips),getmode)
m.area<-tapply(hh$area,
	factor(hh$trpCode,levels=trips),getmode)
month<-as.integer(matrix(unlist(strsplit(hh$date,'-')),nrow=3)[2,])
if(usequarter)month<-1+floor((month-1)/3)
m.seas<-tapply(month,
	factor(hh$trpCode,levels=trips),getmode)

tot.trip.land<-c(tot.trip.land,totland)
n.size.class<-c(n.size.class,nsize)
n.hauls<-c(n.hauls,tr$foNum)
n.land<-c(n.land,nland)
seas<-c(seas,m.seas)
gear<-c(gear,m.gear)
area<-c(area,m.area)
if(sum(sl.disc$staNum)>0){
nsamp.disc<-tapply(hl.disc$lenNum,
	list(factor(hl.disc$trpCode,levels=trips),
	factor(hl.disc$staNum,levels=1:max(sl.disc$staNum))),sum,na.rm=T)
sampwt.disc<-tapply(sl.disc$subSampWt,
	list(factor(sl.disc$trpCode,levels=trips),
	factor(sl.disc$staNum,levels=1:max(sl.disc$staNum))),sum,na.rm=T)
totwt.disc<-tapply(sl.disc$wt,
	list(factor(sl.disc$trpCode,levels=trips),
	factor(sl.disc$staNum,levels=1:max(sl.disc$staNum))),sum,na.rm=T)
ndisc<-nsamp.disc*totwt.disc/sampwt.disc
is.haul<-matrix(0,nrow=length(trips),ncol=max(sl.disc$staNum))
for(i in 1:length(trips)){
 hn<-sl.disc$staNum[sl.disc$trpCode==trips[i]]
 if(sum(hn)>0){
 use<-1:max(hn)
 is.haul[i,use]<-1}}
ndisc[is.haul==0]<-NA
ndisc[is.haul==1&is.na(ndisc)]<-0
ndiscdist<-round(as.vector(ndisc[!is.na(ndisc)]))
totdisc<-round(rowSums(ndisc,na.rm=T))
tot.disc<-c(tot.disc,totdisc)
n.disc.dist<-c(n.disc.dist,ndiscdist)}
}
return(list(tot.trip.land=tot.trip.land,n.size.class=n.size.class,
	n.hauls=n.hauls,n.land=n.land,seas=seas,gear=gear,area=area,
	n.disc.dist=n.disc.dist,n.disc=tot.disc))
}
#######################################################################
make.sim.models<-function(data,div){
nsize<-max(data$n.size.class,na.rm=T)
psize<-model.n.size(nsize,data$n.size.class,
	data$n.land,div=div)
nhaul<-max(data$n.haul,na.rm=T)
phaul<-model.n.size(nhaul,data$n.haul[1:53],
	data$n.disc,div=div)
lm.n.haul<-model.n.haul(data$n.hauls,data$n.land,data$seas,
	data$gear,data$area)
#lm.tot.land<-model.tot.land(data$seas,data$gear,data$area,
#	data$tot.trip.land)
lm.n.land<-model.n.land(data$seas,data$gear,data$area,
	data$n.land)

return(list(psize=psize,div=div,lm.n.haul=lm.n.haul,phaul=phaul,
	lm.n.land=lm.n.land))
}
#######################################################################
subvar<-function(var,xlevels){
tempvar<-var
if(!is.null(xlevels)){
tempvar[!is.element(var,xlevels)]<-xlevels[1]
tempvar<-as.factor(tempvar)}
tempvar
}
#######################################################################
code.cov<-function(cov,covlist){
cov<-as.character(cov)
cov.lev<-unique(cov)
new.cov<-rep(0,length(cov))
for (c in cov.lev)new.cov[cov==c]<-(1:length(covlist))[covlist==c]
new.cov
}
###################################################################
zerovar1<-function(n,param){
x<-rep(0,n)
if(!is.null(param))x<-rnorm(n,0,1/(sqrt(param)))
x
}
###################################################################
zerovar2<-function(param,dim){
x<-0
if(!is.null(param)){
d<-dim(param)
if(length(d)>length(dim))d<-d[d!=1]
if(length(dim)>length(d))dim<-dim[dim!=1]
if(is.null(d))d<-length(param)
checkdim<-0
for(i in 1:length(d))if(dim[i]>d[i])checkdim<-1
if(length(d)!=length(dim))checkdim<-1
if(checkdim==0){
	if(length(d)==1)x<-param[dim]
	if(length(d)==2)x<-param[dim[1],dim[2]]
	if(length(d)==3)x<-param[dim[1],dim[2],dim[3]]}
if(length(x)==2){print("d")
}}
x
}
###################################################################
make.cost.cell<-function(cell,params,agerange,simno){
ncell<-max(cell)
agecell<-matrix(zerovar1(length(agerange)*ncell,
	params$age$Int$cell.tau[simno]),nrow=ncell)
lgaSlpcell<-zerovar1(ncell,params$lga$Slp$cell.tau[simno])
lgaIntcell<-zerovar1(ncell,params$lga$Int$cell.tau[simno])
wglSlpcell<-zerovar1(ncell,params$wgl$Slp$cell.tau[simno])
wglIntcell<-zerovar1(ncell,params$wgl$Int$cell.tau[simno])
list(agecell=agecell,lgaSlpcell=lgaSlpcell,lgaIntcell=lgaIntcell,
	wglSlpcell=wglSlpcell,wglIntcell=wglIntcell)
}
###################################################################
cost.make.cell.p<-function(params,cell.effects,seas,gear,area,cell,simno,species){
n<-nrow(params$age$Int$eff$Const)

p<-make.null.sum(list(params$age$Int$eff$Const[,simno],
	params$age$Int$eff$seas[,seas,simno],
	params$age$Int$eff$gear[,gear,simno],
	params$age$Int$eff$area[,area,simno],
	t(cell.effects$agecell[cell,])),n)
p<-t(p)
age.p<-p
#age.p<-exp(p)/rowSums(exp(p))
n<-1
if(is.null(dim(params$lga$Int$eff$seas)))vs<-params$lga$Int$eff$seas[seas]
else vs<-params$lga$Int$eff$seas[seas,simno]
if(is.null(dim(params$lga$Int$eff$gear)))vg<-params$lga$Int$eff$gear[gear]
else vg<-params$lga$Int$eff$gear[gear,simno]
if(is.null(dim(params$lga$Int$eff$area)))va<-params$lga$Int$eff$area[area]
else va<-params$lga$Int$eff$area[area,simno]
lga.int<-make.null.sum(list(params$lga$Int$eff$Const[simno],
	vs,vg,va,t(cell.effects$lgaIntcell[cell])),n)

lga.slp<-params$lga$Slp$eff$Const[simno]

if(!is.null(params$wgl$Int$eff$Const)){
n<-1
if(is.null(dim(params$wgl$Int$eff$seas)))vs<-params$wgl$Int$eff$seas[seas]
else vs<-params$wgl$Int$eff$seas[seas,simno]
if(is.null(dim(params$wgl$Int$eff$gear)))vg<-params$wgl$Int$eff$gear[gear]
else vg<-params$wgl$Int$eff$gear[gear,simno]
if(is.null(dim(params$wgl$Int$eff$area)))va<-params$wgl$Int$eff$area[area]
else va<-params$wgl$Int$eff$area[area,simno]
wgl.int<-make.null.sum(list(params$wgl$Int$eff$Const[simno],
	vs,vg,va,t(cell.effects$wglIntcell[cell])),n)

wgl.slp<-params$wgl$Slp$eff$Const[simno]}
else {
  if(species=="Melanogrammus aeglefinus")wgl<-getwgl("Haddock")
  if(species!="Melanogrammus aeglefinus")wgl<-getwgl("Cod")
  wgl.int<-wgl$const+wgl$seas[seas]
  wgl.slp<-wgl$slope}

list(age.p=age.p,lga.int=lga.int,lga.slp=lga.slp,
	wgl.int=wgl.int,wgl.slp=wgl.slp)
}
###################################################################
cost.make.cell.lga<-function(p,structure,realseas){
nseas<-4
if(max(as.integer(realseas))>4)nseas<-12
age<-0:(dim(p$age.p)[2]-1)
logage<-matrix(age,ncol=length(p$lga.int),nrow=length(age))
logage<-logage+matrix((as.integer(realseas))/nseas,
	ncol=length(p$lga.int),nrow=length(age),byrow=T)
logage<-log(logage)
maxage<-log(max(age)+1)
minage<-log(min(age)+1/nseas)
logage<-(logage-minage)/(maxage-minage)
lga<-logage*p$lga.slp+matrix(p$lga.int,ncol=ncol(logage),nrow=nrow(logage),byrow=T)
lga
}
###################################################################
cost.sim.fish<-function(trip,n.land,p,cost.cell.lga,params,itrip,simno){
  page<-p$age.p[trip,]

haul.age<-zerovar1(length(page),params$age$Int$haul.tau[simno])
page<-page+haul.age
page<-exp(page)/sum(exp(page))
  page2<<-page
plga<-cost.cell.lga[,trip]
haul.lga<-zerovar1(1,params$lga$Int$haul.tau[simno])

plga<-plga+haul.lga

wglSlp<-p$wgl.slp
wglInt<-p$wgl.int[trip]
  if(!is.null(params$wgl$Int$tau)){
haul.wgl<-zerovar1(1,params$wgl$Int$haul.tau[simno])
wglInt<-wglInt+haul.wgl
wglsd<-1/sqrt(params$wgl$tau.obs[simno])
}
  else wglsd<-0
    

l<-log(1:100)
lint1<-log((1:100)-0.5)
lint2<-log((1:100)+0.5)
lmid<-0.5*(lint1+lint2)
lsd<-1/sqrt(params$lga$tau.obs[simno])
pla<-NULL
for(i in 1:length(plga))pla<-cbind(pla,page[i]*(pnorm(lint2,plga[i],lsd)-pnorm(lint1,plga[i],lsd)))
pla<-pla/sum(pla)
#la<-matrix(rmultinom(1,n.land,pla),nrow=nrow(pla),ncol=ncol(pla))
#r<-3.4
#k<-1
#m<-10
  k<-params$COST.list$cens.mcmc$k[simno]
  m<-params$COST.list$cens.mcmc$m[simno]
  if(is.matrix(params$COST.list$cens.mcmc$mu)){
  rmean<-params$COST.list$cens.mcmc$mu[3,simno]
  rsd<-1/sqrt(params$COST.list$cens.mcmc$tau[3,simno])}
  else {rmean<-params$COST.list$cens.mcmc$mu[3]
  rsd<-1/sqrt(params$COST.list$cens.mcmc$tau[3])}
if(rsd==1)r<-rmean
else  r<-rnorm(n=1,mean=rmean,sd=rsd)
pland<-k*pnorm(m*(lmid-r))

plaland<-pla*matrix(pland,nrow=nrow(pla),ncol=ncol(pla))
pladisc<-pla*matrix(1-pland,nrow=nrow(pla),ncol=ncol(pla))
sim<-matrix(rmultinom(1,n.land,cbind(plaland,pladisc)),nrow=nrow(pla),ncol=2*ncol(pla))

laland<-sim[,1:ncol(pla)]
ladisc<-sim[,-(1:ncol(pla))]
list(laland=laland,ladisc=ladisc,l=lmid,wglSlp=wglSlp,wglInt=wglInt,wglsd=wglsd)
}
###################################################################
make.size.sample<-function(lamat,l.in,wglSlp,wglInt,wglsd,
	a.in,nsize,nlsamp,nasamp){
  l.in<-round(exp(l.in))
nage<-length(a.in)
l.all<-rep(l.in,nage)
age.all<-rep(a.in,each=length(l.in))
l.all.freq<-as.vector(lamat)
use<-l.all.freq>0
l.all<-l.all[use]
l.all.freq<-l.all.freq[use]
age.all<-age.all[use]

size<-allocate.size(l.all,l.all.freq,nsize)
nofish<-colSums(size)==0
nsize<-sum(!nofish)
size<-size[,!nofish]
if(!is.matrix(size))size<-matrix(size,nrow=length(size))
l.age.out<-age.out<-l.out<-l.freq.out<-wt.samp<-wt.land<-vector("list",length=nsize)
for(i in 1:nsize){
  if(sum(size[,i])>0){
  lsamp<-rep(l.all,size[,i])
agesamp<-rep(age.all,size[,i])
wt.land[[i]]<-NA
if(!is.null(wglsd))wt.land[[i]]<-sum(wglInt+wglSlp*lsamp+rnorm(length(lsamp),0,wglsd))
n<-length(lsamp)
if(nlsamp>=n){
long.l.samp<-lsamp
long.age.samp<-agesamp
}
else{

samp<-sample(1:length(lsamp),nlsamp)
long.l.samp<-lsamp[samp]
long.age.samp<-agesamp[samp]
}
wt.samp[[i]]<-NA
if(!is.null(wglsd))wt.samp[[i]]<-sum(wglInt+wglSlp*long.l.samp+rnorm(length(long.l.samp),0,wglsd))

la<-table(long.l.samp,long.age.samp)
 
lamat2<-matrix(la,ncol=ncol(la))
l.freq.out[[i]]<-rowSums(la)
l.out[[i]]<-as.numeric(rownames(la))
anames<-as.integer(colnames(la))
l.age.temp<-age.temp<-NULL
for(j in 1:length(l.out[[i]])){
	if(l.freq.out[[i]][j]<nasamp){
		l.age.temp<-c(l.age.temp,rep(l.out[[i]][j],l.freq.out[[i]][j]))
				age.temp<-c(age.temp,rep(anames,lamat2[j,]))}
	else
		{l.age.temp<-c(l.age.temp,rep(l.out[[i]][j],nasamp))
			age.temp<-c(age.temp,sample(rep(anames,lamat2[j,]),nasamp))}
				}
l.age.out[[i]]<-l.age.temp
age.out[[i]]<-age.temp
}}
list(l.age=l.age.out,age=age.out,l=l.out,lfreq=l.freq.out,
	wt.land=wt.land,wt.samp=wt.samp,nsize=nsize)
}
###################################################################
make.disc.sample<-function(lamat,l.in,wglSlp,wglInt,wglsd,
	a.in,phaul,nlsamp,nasamp,div){
nfish<-sum(lamat)
l.age.out<-age.out<-l.out<-l.freq.out<-nhaul<-wt.land<-wt.samp<-NULL
if(nfish>2){
nh<-ncol(phaul)
is<-1
for(i in 1:length(div))is<-is+(log(nfish)>div[i])
nhaul<-sample(1:nh,1,prob=phaul[is,])
if(nhaul>=nfish)nhaul<-nfish-1

p<-c(0,sort(sample(1:(nfish-1),nhaul-1)),nfish)
n<-diff(p)

l.in<-round(exp(l.in))
nage<-length(a.in)
l.all<-rep(l.in,nage)
age.all<-rep(a.in,each=length(l.in))
l.all.freq<-as.vector(lamat)
use<-l.all.freq>0
l.all<-l.all[use]
l.all.freq<-l.all.freq[use]
age.all<-age.all[use]
l.age.out<-age.out<-l.out<-l.freq.out<-
	vector("list",length=nhaul)
ltot<-rep(l.all,l.all.freq)
atot<-rep(age.all,l.all.freq)
for(i in 1:nhaul){
  sampuse<-sample(1:length(ltot),n[i])
	lsamp<-ltot[sampuse]
	agesamp<-atot[sampuse]
	if(nlsamp>=n[i]){
	long.l.samp<-lsamp
	long.age.samp<-agesamp
	}
	else{
	samp<-sample(1:length(lsamp),nlsamp)
	long.l.samp<-lsamp[samp]
	long.age.samp<-agesamp[samp]
	}
wt.land[[i]]<-wt.samp[[i]]<-NA
if(!is.null(wglsd)){
wt.land[[i]]<-sum(wglInt+wglSlp*lsamp+rnorm(length(lsamp),0,wglsd))
wt.samp[[i]]<-sum(wglInt+wglSlp*long.l.samp+rnorm(length(long.l.samp),0,wglsd))}

la<-table(long.l.samp,long.age.samp)

lamat2<-matrix(la,ncol=ncol(la))
l.freq.out[[i]]<-rowSums(la)
l.out[[i]]<-as.numeric(rownames(la))
anames<-as.integer(colnames(la))
l.age.temp<-age.temp<-NULL
for(j in 1:length(l.out[[i]])){
	if(l.freq.out[[i]][j]<nasamp){
		l.age.temp<-c(l.age.temp,rep(l.out[[i]][j],l.freq.out[[i]][j]))
				age.temp<-c(age.temp,rep(anames,lamat2[j,]))}
	else
		{l.age.temp<-c(l.age.temp,rep(l.out[[i]][j],nasamp))
			age.temp<-c(age.temp,sample(rep(anames,lamat2[j,]),nasamp))}
				}
l.age.out[[i]]<-l.age.temp
age.out[[i]]<-age.temp

ltot<-ltot[-sampuse]
atot<-atot[-sampuse]
}}
list(l.age=l.age.out,age=age.out,l=l.out,lfreq=l.freq.out,nhaul=nhaul,
	wt.land=wt.land,wt.samp=wt.samp)
}
###################################################################
sim.to.cost<-function(sample.land,sample.disc,structure,species,type){
#-------------------------------------------------------------------
ntrip<-length(structure$n.haul)
trpCode<-paste(type,"trip",1:ntrip,sep="")
#-------------------------------------------------------------------
#tr
foNum<-structure$n.haul # no hauls per trip
foNum[structure$dis==0]<-0
obs.first<-order(-structure$dis)
tr<-list(foNum=foNum,trpCode=trpCode)
#-------------------------------------------------------------------
#hl
land.lenCls<-land.lenNum<-disc.lenCls<-disc.lenNum<-land.commCat<-
	land.trpCode<-disc.staNum<-disc.trpCode<-NULL
for(i in obs.first){
	land.lenCls<-c(land.lenCls,unlist(sample.land[[i]]$l))
	land.lenNum<-c(land.lenNum,unlist(sample.land[[i]]$lfreq))
	disc.lenCls<-c(disc.lenCls,unlist(sample.disc[[i]]$l))
	disc.lenNum<-c(disc.lenNum,unlist(sample.disc[[i]]$lfreq))
	land.trpCode<-c(land.trpCode,rep(trpCode[i],
		length(unlist(sample.land[[i]]$l))))
	nh<-sample.disc[[i]]$nhaul
	if(!is.null(nh)){
	for(j in 1:nh)disc.staNum<-c(disc.staNum,
		rep(j,length(sample.disc[[i]]$l[[j]])))}

	disc.trpCode<-c(disc.trpCode,rep(trpCode[i],
		length(unlist(sample.disc[[i]]$l))))
	for(j in 1:length(sample.land[[i]]$l)){
		land.commCat<-c(land.commCat,rep(j,length(sample.land[[i]]$l[[j]])))}
}
disc.commCat<-rep(NA,length(disc.lenNum))
disc.catchCat<-rep("DIS",length(disc.trpCode))
land.catchCat<-rep("LAN",length(land.trpCode))
land.staNum<-rep(NA,length(land.trpCode))
hl.lenCls<-c(land.lenCls,disc.lenCls)
hl.lenNum<-c(land.lenNum,disc.lenNum)
hl.commCat<-c(land.commCat,disc.commCat)
hl.trpCode<-c(land.trpCode,disc.trpCode)
hl.catchCat<-c(land.catchCat,disc.catchCat)
hl.spp<-rep(species,length(hl.catchCat))
hl.staNum<-c(land.staNum,disc.staNum)

hl<-list(lenCls=hl.lenCls,lenNum=hl.lenNum,commCat=hl.commCat,trpCode=hl.trpCode,
	catchCat=hl.catchCat,spp=hl.spp,staNum=hl.staNum)
#-------------------------------------------------------------------
#sl
disc.wt<-disc.subSampWt<-land.wt<-land.subSampWt<-disc.trpCode<-
	land.commCat<-disc.staNum<-NULL
for(i in obs.first){
	disc.wt<-c(disc.wt,unlist(sample.disc[[i]]$wt.land))
	disc.subSampWt<-c(disc.subSampWt,unlist(sample.disc[[i]]$wt.samp))
	land.wt<-c(land.wt,unlist(sample.land[[i]]$wt.land))
	land.subSampWt<-c(land.subSampWt,unlist(sample.land[[i]]$wt.samp))
	nh<-sample.disc[[i]]$nhaul
	if(is.null(nh))nh<-0
	disc.trpCode<-c(disc.trpCode,rep(trpCode[i],nh))
	if(nh>0){for(j in 1:nh)disc.staNum<-c(disc.staNum,j)}
        if(length(sample.land[[i]]$wt.land)>0){
	for(j in 1:length(sample.land[[i]]$wt.land))land.commCat<-c(land.commCat,j)}
}
disc.commCat<-rep(NA,length(disc.wt))
land.staNum<-rep(NA,length(land.wt))
land.trpCode<-rep(trpCode,structure$n.size[obs.first])
disc.catchCat<-rep("DIS",length(disc.wt))
land.catchCat<-rep("LAN",length(land.wt))
sl.wt<-c(land.wt,disc.wt)
sl.subSampWt<-c(land.subSampWt,disc.subSampWt)
sl.catchCat<-c(land.catchCat,disc.catchCat)
sl.trpCode<-c(land.trpCode,disc.trpCode)
sl.spp<-rep(species,length(sl.catchCat))
sl.commCat<-c(land.commCat,disc.commCat)
sl.staNum<-c(land.staNum,disc.staNum)
sl<-list(wt=sl.wt,subSampWt=sl.subSampWt,catchCat=sl.catchCat,
	spp=sl.spp,trpCode=sl.trpCode,staNum=sl.staNum,commCat=sl.commCat)
#-------------------------------------------------------------------
#ca
disc.lenCls<-disc.age<-land.lenCls<-land.age<-land.trpCode<-disc.trpCode<-NULL

for(i in obs.first){
	disc.lenCls<-c(disc.lenCls,unlist(sample.disc[[i]]$l.age))
	disc.age<-c(disc.age,unlist(sample.disc[[i]]$age))
	land.lenCls<-c(land.lenCls,unlist(sample.land[[i]]$l.age))
	land.age<-c(land.age,unlist(sample.land[[i]]$age))
	land.trpCode<-c(land.trpCode,rep(trpCode[i],length(unlist(sample.land[[i]]$age))))
	disc.trpCode<-c(disc.trpCode,rep(trpCode[i],length(unlist(sample.disc[[i]]$age))))
}
catchCat<-c(rep("LAN",length(land.age)),rep("DIS",length(disc.age)))

ca.lenCls<-c(land.lenCls,disc.lenCls)
ca.age<-c(land.age,disc.age)
ca.trpCode<-c(land.trpCode,disc.trpCode)
ca.spp<-rep(species,length(catchCat))
ca<-list(lenCls=ca.lenCls,age=ca.age,trpCode=ca.trpCode,
	catchCat=catchCat,spp=ca.spp)

#-------------------------------------------------------------------
#hh
#n.seas<-4
#if(max(as.integer(structure$seas))>4)n.seas<-12
#obs.seas<-as.integer(structure$seas[obs.first])
#if(n.seas==4)obs.seas<-obs.seas*3-1
#hh.date<-paste(2000,"-",obs.seas,"-",15,sep='')
#hh.foCatEu5<-structure$gear[obs.first]
#hh.area<-structure$area[obs.first]
#trpCode<-trpCode[obs.first]
#hh<-list(date=hh.date,foCatEu5=hh.foCatEu5,area=hh.area,trpCode=trpCode)

n.seas<-4
if(max(as.integer(structure$seas))>4)n.seas<-12
obs.seas<-as.integer(structure$seas[obs.first])
if(n.seas==4)obs.seas<-obs.seas*3-1
hh.date<-paste(2000,"-",obs.seas,"-",15,sep='')
hh.foCatEu5<-structure$gear[obs.first]
hh.area<-structure$area[obs.first]
trpCode<-trpCode[obs.first]

haulrep<-haulno<-NULL
if(sum(structure$dis==1)>0){
haulrep<-structure$n.haul[structure$dis==1]
haulno<-NULL
for(i in 1:length(haulrep))haulno<-c(haulno,1:haulrep[i])}

hh.date<-c(hh.date[structure$dis==0],rep(hh.date[structure$dis==1],haulrep))
hh.foCatEu5<-c(hh.foCatEu5[structure$dis==0],rep(hh.foCatEu5[structure$dis==1],haulrep))
hh.area<-c(hh.area[structure$dis==0],rep(hh.area[structure$dis==1],haulrep))
trpCode<-c(trpCode[structure$dis==0],rep(trpCode[structure$dis==1],haulrep))

staNum<-c(rep(NA,sum(structure$dis==0)),haulno)
hh<-list(date=hh.date,foCatEu5=hh.foCatEu5,area=hh.area,trpCode=trpCode,staNum=staNum)

list(tr=tr,hl=hl,sl=sl,ca=ca,hh=hh)
}

###################################################################
getmode<-function(x){
if(is.factor(x))x<-as.character(x)
mode<-NA
if(sum(!is.na(x))>0){
tab<-as.integer(table(x))
m<-max(tab)
mode<-unique(x)[tab==m]
mode<-mode[1]}
mode
}
###################################################################
make.null.sum<-function(list,n){
x<-rep(0,n)
  for(i in 1:length(list)){
    if(!is.null(list[[i]]))x<-x+list[[i]]}
x
}
###################################################################

quadplot2<-function(params,loc="topleft",isim=NULL,title=NULL){
nage<-dim(params)[1]
q<-NULL
for(i in 1:nage)q<-rbind(q,quantile(params[i,],c(0.05,0.95)))
m<-rowMeans(params)
plot(1:nage,m,ylim=range(cbind(m,q)),col=2,xlab='age group',ylab='parameter')
points(1:nage,q[,1],pch=2,col=3)
points(1:nage,q[,2],pch=3,col=3)
if(loc!='false'){
legend(loc,pch=1:4,col=c(2,3,3,4),legend=c("mean","5%","95%","true"))}
if(!is.null(isim))title(paste('sim',isim,title))
{}
}
###################################################################
# dga: setup arguments change to a list
setup<-function(data.args){
	
    datalist <- data.args$datalist
    species  <- data.args$species
    ageseas  <- data.args$ageseas; agegear <- data.args$agegear; agearea <- data.args$agearea
    lgaseas  <- data.args$lgaseas; lgagear <- data.args$lgagear; lgaarea <- data.args$lgaarea
    wglseas  <- data.args$wglseas; wglgear <- data.args$wglgear; wglarea <- data.args$wglarea
    ntrip    <- data.args$ntrip
    ndisc    <- data.args$ndisc
    div      <- data.args$div
    ageMin   <- data.args$ageMin
    ageMax   <- data.args$ageMax
    nseas    <- data.args$nseas
    	
agecell<-ageseas+agegear+agearea>1
lgacell<-lgaseas+lgagear+lgaarea>1
wglcell<-wglseas+wglgear+wglarea>1
model.data<-get.sim.model.data(datalist,species,nseas==4)
seaslist<-sort(unique(model.data$seas))
arealist<-sort(unique(model.data$area))
gearlist<-sort(unique(model.data$gear))
agemodel=list(Int=list(year=FALSE,seas=ageseas,gear=agegear,area=agearea,
                         cell=agecell,haul=TRUE),
                Hsz=NULL)  
  lgamodel=list(Int=list(year=FALSE,seas=lgaseas,gear=lgagear,area=lgaarea,
                         cell=lgacell,haul=FALSE),
                Slp=list(year=FALSE,seas=FALSE,gear=FALSE,area=FALSE,
                         cell=FALSE,haul=FALSE),
                Hsz=NULL)        
  wglmodel=list(Int=list(year=FALSE,seas=wglseas,gear=wglgear,area=wglarea,
                         cell=wglcell,haul=TRUE),
                Slp=list(year=FALSE,seas=FALSE,gear=FALSE,area=FALSE,
                         cell=FALSE,haul=FALSE),
                Hsz=NULL)           
  prior.par=list(age=list(Int=NULL,Hsz=NULL),
                 lga=list(Int=NULL,Slp=NULL,Hsz=NULL),
                 wgl=list(Int=NULL,Slp=NULL,Hsz=NULL))


sampmat<-get.sampled.cells(seaslist,gearlist,arealist,ntrip=ntrip,ndisc=ndisc,
	force=NULL)
sim.models<-make.sim.models(model.data,div)
structure<-make.structure(sampmat,sim.models,seaslist,gearlist,arealist,model.data,nseas)

list(structure=structure,agemodel=agemodel,lgamodel=lgamodel,wglmodel=wglmodel,
	prior.par=prior.par,seaslist=seaslist,arealist=arealist,gearlist=gearlist,
       ageMin=ageMin,ageMax=ageMax,div=div,species=species,sim.models=sim.models)
}

###################################################################
make.tau<-function(params,model){
area.tau<-cell.tau<-haul.tau<-NULL
taulist<-c(model$area,model$cell,model$haul)

if(sum(taulist)==1){
  if(model$area)area.tau<-params$tau
  if(model$cell)cell.tau<-params$tau
  if(model$haul)haul.tau<-params$tau}
if(sum(taulist)==3){
  area.tau<-params$tau[1,]
  cell.tau<-params$tau[2,]
  haul.tau<-params$tau[3,]}
if(sum(taulist)==2){
  id<-cumsum(taulist)
 if(model$area)area.tau<-params$tau[1,]
 if(model$cell) cell.tau<-params$tau[id[2],]
 if(model$haul)haul.tau<-params$tau[id[3],]}
params$area.tau<-area.tau
params$cell.tau<-cell.tau
params$haul.tau<-haul.tau
params
}
###################################################################
getwgl<-function(species)
{
       if(species == "Cod") {
               slope <- 2.8571
               a <- log(c(0.0176, 0.0169, 0.0168, 0.0163, 0.0174, 0.0172,
                       0.017, 0.0185, 0.018, 0.0181, 0.0182, 0.0177))
               const <- mean(a)
               seas <- a - const
       }
       if(species == "Haddock") {
               slope <- 2.8268
               a <- log(c(0.0155, 0.0153, 0.0145, 0.0148, 0.015, 0.0151,
                       0.016, 0.0163, 0.0164, 0.017, 0.0164, 0.0159))
               const <- mean(a)
               seas <- a - const
       }
list(const=const,slope=slope,seas=seas)
}
