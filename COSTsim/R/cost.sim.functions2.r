#######################################################################
cost.simulate<-function(params, setup.data,nlsamp.land,nasamp.land,nlsamp.disc,nasamp.disc,length.list){
seas<-code.cov(setup.data$structure$seas, setup.data$seaslist)
gear<-code.cov(setup.data$structure$gear, setup.data$gearlist)
area<-code.cov(setup.data$structure$area, setup.data$arealist)
realseas<-setup.data$structure$seas
ageMin<-setup.data$ageMin
ageMax<- setup.data$ageMax
simno<-1
ntrip<-length(seas)
cell.effects<-make.cost.cell(setup.data$structure,params,ageMin:ageMax,simno)
cell.effects<<-cell.effects
p<-cost.make.cell.p(params,cell.effects,seas,gear,area,
	simno, setup.data$species)
p<<-p
cost.cell.lga<-cost.make.cell.lga(p, setup.data$structure,realseas)

sample.land<-sample.disc<-vector("list",length=ntrip)
nsize<-nhaul<-NULL
paa.land<-NULL
for(itrip in 1:ntrip){
#for(itrip in 1){
simfish<-cost.sim.fish(itrip,
        setup.data$structure$landratio*setup.data$structure$n.land[itrip],p,
	cost.cell.lga,params,itrip,simno,length.list)
paa.land<-rbind(paa.land,colSums(simfish$laland)/sum(simfish$laland))
if(sum(simfish$laland)>0){
  sample.land[[itrip]]<-make.size.sample(simfish$laland,simfish$l,
	simfish$wglSlp,simfish$wglInt,simfish$wglsd,ageMin:ageMax,
	setup.data$structure$n.size[itrip],nlsamp.land,nasamp.land)
nsize<-c(nsize,sample.land[[itrip]]$nsize)}
else nsize<-c(nsize,0)
if(sum(simfish$ladisc)>0){ 
sample.disc[[itrip]]<-make.disc.sample(simfish$ladisc,simfish$l,
	simfish$wglSlp,simfish$wglInt,simfish$wglsd,ageMin:ageMax,
	nlsamp.disc,nasamp.disc)
nhaul<-c(nhaul,sample.disc[[itrip]]$nhaul)}
else nhaul<-c(nhaul,0)
}
list(sample.land=sample.land,sample.disc=sample.disc,nsize=nsize,nhaul=nhaul,paa.land=paa.land)
}
#######################################################################
get.sampled.cells<-function(params,seasons=NULL,gears=NULL,areas=NULL,nmland=0,nobs=0,force=NULL){
ntrip<-nmland+nobs
if(is.null(seasons))nseason<-1
else nseason<-length(seasons)
if(is.null(gears))ngear<-1
else ngear<-length(gears)
if(is.null(areas))narea<-1
else narea<-length(areas)
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
trip.vec<-sampmat[,4]+sampmat[,5]
seaslist<-setdiff(1:nseason,useas)
gearlist<-setdiff(1:ngear,ugear)
arealist<-setdiff(1:narea,uarea)
n<-max(c(length(seaslist),length(gearlist),length(arealist)))

v1<-c(sample(seaslist),sample(1:nseason,n-length(seaslist),replace=T))
v2<-c(sample(gearlist),sample(1:ngear,n-length(gearlist),replace=T))
if(length(arealist)>1)v3<-c(sample(arealist,2),sample(1:narea,n-2,replace=T))
else v3<-c(arealist,sample(1:narea,n-1,replace=T))
if(length(v1)>1)v1<-sample(v1)
if(length(v2)>1)v2<-sample(v2)
if(length(v2)>1)v3<-sample(v3)
for(i in 1:n){suse<-sampmat[,1]==v1[i]&sampmat[,2]==v2[i]&
	sampmat[,3]==v3[i]
trip.vec[suse]<-max(1,trip.vec[suse])}
nalloc<-sum(trip.vec)
if(ntrip>nalloc){
extra<-sample((1:ncell)[!in.force],ntrip-nalloc,replace=T)
for(e in extra)trip.vec[e]<-trip.vec[e]+1}
nsamp<-sum(trip.vec)
if(nsamp>ntrip)cat("ntrip increased to",nsamp,"\n")

nobs.rem<-nobs-sum(sampmat[,5])
nml.rem<-nmland-sum(sampmat[,4])

cells.left<-sample(rep(1:ncell,trip.vec))
if(nobs.rem>0){obs<-cells.left[1:nobs.rem]
               for(i in obs)sampmat[i,5]<-sampmat[i,5]+1}
if(nml.rem>0){ml<-cells.left[(nobs.rem+1):(nobs.rem+nml.rem)]
               for(i in ml)sampmat[i,4]<-sampmat[i,4]+1}
full.cell<-sampmat[,3]+(sampmat[,2]-1)*narea+(sampmat[,1]-1)*ngear*narea
obs.cell<-get.obs.cell(sampmat,params)
cell<-rank(full.cell-100000*obs.cell)
nobs<-sum(obs.cell)
cell[obs.cell==0]<-cell[obs.cell==0]-nobs
sampmat<-cbind(sampmat,cell,obs.cell)
sampmat
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
if(length(nph.dist)>1)newdisc<-sample(nph.dist,nhauls,replace=T)
else newdisc<-rep(nph.dist,nhauls)
haul<-NULL
for(i in 1:length(l))haul<-cbind(haul,rmultinom(1,lfreq[i],newdisc))
haul
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
make.cost.cell<-function(structure,params,agerange,simno){
  ncell<-length(structure$dis)
obscell<-structure$obs.cell==1
  agecell<-matrix(0,nrow=ncell,ncol=length(agerange))
if(!is.null(params$age$Int$eff$cell))agecell[obscell,]<-params$age$Int$eff$cell[,structure$cell[obscell],simno]
lgaIntcell<-rep(0,ncell)
  if(!is.null(params$lga$Int$eff$cell)){
if(!is.null(dim(params$lga$Int$eff$cell)))
  lgaIntcell[obscell]<-params$lga$Int$eff$cell[structure$cell[obscell,simno]]
else lgaIntcell[obscell]<-params$lga$Int$eff$cell[structure$cell[obscell]]}
wglIntcell<-rep(0,ncell)
  if(!is.null(params$wgl$Int$eff$cell)){  
if(!is.null(dim(params$wgl$Int$eff$cell)))
  wglIntcell[obscell]<-params$wgl$Int$eff$cell[structure$cell[obscell,simno]]
else if(!is.null(params$wgl$Int$eff$cell))wglIntcell[obscell]<-params$wgl$Int$eff$cell[structure$cell[obscell]]}

if(sum(obscell==0)>0){
  narea<-length(params$age$data$neigh$num)
if(!is.null(params$lga$Int$eff$cell)) lga.cell<-get.unobs.cell.params(params$lga,1,narea,simno) 
if(!is.null(params$age$Int$eff$cell)) age.cell<-get.unobs.cell.params(params$age,length(agerange),narea,simno) 
if(!is.null(params$wgl$Int$eff$cell))if(!is.null(params$wgl)) wgl.cell<-get.unobs.cell.params(params$wgl,1,narea,simno) 
if(!is.null(params$lga$Int$eff$cell))lgaIntcell[!obscell]<-lga.cell[structure$cell[!obscell]]
if(!is.null(params$wgl$Int$eff$cell))if(!is.null(params$wgl)&!is.null(wgl.cell))wglIntcell[!obscell]<-wgl.cell[structure$cell[!obscell]]
if(!is.null(params$age$Int$eff$cell))agecell[!obscell,]<-age.cell[structure$cell[!obscell],]
}

lgaSlpcell<-wglSlpcell<-rep(0,ncell)
list(agecell=agecell,lgaSlpcell=lgaSlpcell,lgaIntcell=lgaIntcell,
	wglSlpcell=wglSlpcell,wglIntcell=wglIntcell)
}
###################################################################
cost.make.cell.p<-function(params,cell.effects,seas,gear,area,simno,species){
n<-nrow(params$age$Int$eff$Const)
nseas<-dim(params$age$Int$eff$seas)[2]
p<-make.null.sum(list(params$age$Int$eff$Const[,simno],
	params$age$Int$eff$seas[,seas,simno],
	params$age$Int$eff$gear[,gear,simno],
	params$age$Int$eff$area[,area,simno],
	t(cell.effects$agecell)),n)
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
	vs,vg,va,t(cell.effects$lgaIntcell)),n)

lga.slp<-params$lga$Slp$eff$Const[simno]

if(!is.null(params$wgl$Int$eff$Const)){
n<-1
if(is.null(dim(params$wgl$Int$eff$seas)))vs<-params$wgl$Int$eff$seas[seas]
else vs<-params$wgl$Int$eff$seas[seas,simno]
if(is.null(dim(params$wgl$Int$eff$gear)))vg<-params$wgl$Int$eff$gear[gear]
else vg<-params$wgl$Int$eff$gear[gear,simno]
if(is.null(dim(params$wgl$Int$eff$area)))va<-params$wgl$Int$eff$area[area]
else va<-params$wgl$Int$eff$area[area,simno]
cw<-t(cell.effects$wglIntcell)
if(sum(is.na(cw))>0)cw<-NULL
wgl.int<-make.null.sum(list(params$wgl$Int$eff$Const[simno],
	vs,vg,va,cw),n)
wgl.slp<-params$wgl$Slp$eff$Const[simno]}
else {
wgl<-getwgl(species)
if(nseas==4)wglseas<-c(mean(wgl$seas[1:3]),mean(wgl$seas[4:6]),mean(wgl$seas[7:9]),mean(wgl$seas[10:12]))
else wglseas<-wgl$seas
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
cost.sim.fish<-function(trip,n.land,p,cost.cell.lga,params,itrip,simno,length.list){
  page<-p$age.p[trip,]

haul.age<-zerovar1(length(page),params$age$Int$haul.tau[simno])
  haul.age<-haul.age-mean(haul.age)
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
    

l<-seq(length.list$minl,length.list$maxl,length.list$int)
lint1<-log(l-length.list$int/2)
lint2<-log(l+length.list$int/2)
lmid<-0.5*(lint1+lint2)
lnew<-seq(from=min(lint1),to=max(lint2),len=10000)
intno<-NULL
for(i in 1:length(lnew))intno<-c(intno,sum(lint2<=lnew[i]))
intno[10000]<-99
intno<-intno+1
  
lsd<-1/sqrt(params$lga$tau.obs[simno])
pla<-NULL
for(i in 1:length(plga))pla<-cbind(pla,page[i]*(pnorm(lint2,plga[i],lsd)-pnorm(lint1,plga[i],lsd)))
pla<-pla/sum(pla)
  pla<<-pla

#la<-matrix(rmultinom(1,n.land,pla),nrow=nrow(pla),ncol=ncol(pla))
#r<-3.4
#k<-1
#m<-10
  k<-params$COST.list$cens.mcmc$k[simno]
  m<-params$COST.list$cens.mcmc$m[simno]
  if(is.matrix(params$COST.list$cens.mcmc$r))rmean<-median(params$COST.list$cens.mcmc$r[,simno])
 else  r<-median(params$COST.list$cens.mcmc$r)
#pland<-k*pnorm(m*(lmid-r))
pnew<-k*pnorm(m*(lnew-r))
  pland<-tapply(pnew,intno,mean)
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
  agl.all.l<-agl.all.a<-wt<-agl<-NULL
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
if(!is.null(wglsd))wt.land[[i]]<-sum(exp(wglInt+wglSlp*log(lsamp)+rnorm(length(lsamp),0,wglsd)))
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
long.wt.sample<-wglInt+wglSlp*log(long.l.samp)+rnorm(length(long.l.samp),0,wglsd)
wt.samp[[i]]<-NA
if(!is.null(wglsd))wt.samp[[i]]<-sum(exp(long.wt.sample))
#agl.all.l<-c(agl.all.l,long.l.samp)
#agl.all.a<-c(agl.all.a,long.age.samp)
up.samp<-sample(1:length(long.l.samp),length(lsamp),replace=T)
agl.all.l<-c(agl.all.l,long.l.samp[up.samp])
agl.all.a<-c(agl.all.a,long.age.samp[up.samp])
la<-table(long.l.samp,long.age.samp)
   
l.freq.out[[i]]<-rowSums(la)
l.out[[i]]<-as.numeric(rownames(la))
}
la.table<-table(agl.all.l,agl.all.a)
agl<-make.agl(la.table,nasamp)
  wt<-exp(wglInt+wglSlp*log(agl$l.age)+rnorm(length(agl$l.age),0,wglsd))
}
list(l.age=agl$l.age,age=agl$age,wt=wt,l=l.out,lfreq=l.freq.out,
	wt.land=wt.land,wt.samp=wt.samp,nsize=nsize)
}
###################################################################
make.agl<-function(la.table,nasamp){
lamat2<-matrix(la.table,ncol=ncol(la.table))
l.freq<-rowSums(la.table)
l<-as.numeric(rownames(la.table))
anames<-as.integer(colnames(la.table))
l.age<-age<-NULL
for(j in 1:length(l)){
  asamp<-rep(anames,lamat2[j,])
	if(l.freq[j]<nasamp){
		l.age<-c(l.age,rep(l[j],l.freq[j]))
				age<-c(age,anames)}
	else
		{l.age<-c(l.age,rep(l[j],nasamp))                 
		if(length(asamp)>1)age<-c(age,sample(asamp,nasamp))
		if(length(asamp)==1)age<-c(age,rep(asamp,nasamp))
               }
				}
list(age=age,l.age=l.age)
}
###################################################################


make.disc.sample<-function(lamat,l.in,wglSlp,wglInt,wglsd,
	a.in,nlsamp,nasamp){
  nfish<-sum(lamat)
agl.all.l<-agl.all.a<-wt<-agl<-NULL
l.age.out<-age.out<-l.out<-l.freq.out<-nhaul<-wt.land<-wt.samp<-NULL
div<-c(4,9)
if(nfish>2){
is<-1
for(i in 1:length(div))is<-is+(log(nfish)>div[i])
  if(is==1)nhaul<-sample(c(2,5),1,prob=c(0.5,0.5))
  if(is==2)nhaul<-sample(c(2,3,4,5,6,7,10,11,12,13,14,15,17,18,19,21),1,
             prob=c(0.08,0.12,0.04,0.04,0.04,0.04,0.08,0.08,0.08,0.04,0.04,0.12,0.08,0.04,0.04,0.04))
  if(is==3)nhaul<-sample(c(4,5,8,9,11,12,14,16,20,21,23,26,36,37),1,
         prob=c(0.05,0.09,0.09,0.05,0.09,0.05,0.05,0.14,0.09,0.05,0.14,0.05,0.05,0.05))
     
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
wt.land[[i]]<-sum(exp(wglInt+wglSlp*log(lsamp)+rnorm(length(lsamp),0,wglsd)))
wt.samp[[i]]<-sum(exp(wglInt+wglSlp*log(long.l.samp)+rnorm(length(long.l.samp),0,wglsd)))
#agl.all.l<-c(agl.all.l,long.l.samp)
#agl.all.a<-c(agl.all.a,long.age.samp)
up.samp<-sample(1:length(long.l.samp),length(lsamp),replace=T)
agl.all.l<-c(agl.all.l,long.l.samp[up.samp])
agl.all.a<-c(agl.all.a,long.age.samp[up.samp])
}

la<-table(long.l.samp,long.age.samp)
lamat2<-matrix(la,ncol=ncol(la))
l.freq.out[[i]]<-rowSums(lamat2)
l.out[[i]]<-as.numeric(rownames(la))
ltot<-ltot[-sampuse]
atot<-atot[-sampuse]
}
la.table<-table(agl.all.l,agl.all.a)
  agl<-make.agl(la.table,nasamp)
  wt<-exp(wglInt+wglSlp*log(agl$l.age)+rnorm(length(agl$l.age),0,wglsd))
}
list(l.age=agl$l.age,age=agl$age,wt=wt,l=l.out,lfreq=l.freq.out,nhaul=nhaul,
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
	disc.commCat<-disc.catchCat<-land.trpCode<-disc.staNum<-disc.trpCode<-NULL
for(i in obs.first){
	land.lenCls<-c(land.lenCls,unlist(sample.land[[i]]$l))
	land.lenNum<-c(land.lenNum,unlist(sample.land[[i]]$lfreq))
	land.trpCode<-c(land.trpCode,rep(trpCode[i],
		length(unlist(sample.land[[i]]$l))))
	for(j in 1:length(sample.land[[i]]$l)){
		land.commCat<-c(land.commCat,rep(j,length(sample.land[[i]]$l[[j]])))}

if(type=="OBS"){
	disc.lenCls<-c(disc.lenCls,unlist(sample.disc[[i]]$l))
	disc.lenNum<-c(disc.lenNum,unlist(sample.disc[[i]]$lfreq))
  nh<-sample.disc[[i]]$nhaul
	if(!is.null(nh)){
	for(j in 1:nh)disc.staNum<-c(disc.staNum,
		rep(j,length(sample.disc[[i]]$l[[j]])))}

	disc.trpCode<-c(disc.trpCode,rep(trpCode[i],
		length(unlist(sample.disc[[i]]$l))))
}}
if(type=="OBS"){
disc.commCat<-rep(NA,length(disc.lenNum))
disc.catchCat<-rep("DIS",length(disc.trpCode))}
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
	disc.commCat<-disc.catchCat<-land.commCat<-disc.staNum<-NULL
for(i in obs.first){
	land.wt<-c(land.wt,unlist(sample.land[[i]]$wt.land))
	land.subSampWt<-c(land.subSampWt,unlist(sample.land[[i]]$wt.samp))
        if(length(sample.land[[i]]$wt.land)>0){
	for(j in 1:length(sample.land[[i]]$wt.land))land.commCat<-c(land.commCat,j)}
if(type=="OBS"){
	disc.wt<-c(disc.wt,unlist(sample.disc[[i]]$wt.land))
	disc.subSampWt<-c(disc.subSampWt,unlist(sample.disc[[i]]$wt.samp))
	nh<-sample.disc[[i]]$nhaul
	if(is.null(nh))nh<-0
	disc.trpCode<-c(disc.trpCode,rep(trpCode[i],nh))
	if(nh>0){for(j in 1:nh)disc.staNum<-c(disc.staNum,j)}
}}
land.staNum<-rep(NA,length(land.wt))
land.trpCode<-rep(trpCode,structure$n.size[obs.first])
if(type=="OBS"){
disc.commCat<-rep(NA,length(disc.wt))
disc.catchCat<-rep("DIS",length(disc.wt))}
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
disc.lenCls<-disc.age<-land.lenCls<-land.age<-disc.wt<-land.wt<-land.trpCode<-disc.trpCode<-NULL

for(i in obs.first){
if(type=="OBS"){
	disc.lenCls<-c(disc.lenCls,sample.disc[[i]]$l.age)
	disc.age<-c(disc.age,sample.disc[[i]]$age)
	disc.wt<-c(disc.wt,sample.disc[[i]]$wt)
	disc.trpCode<-c(disc.trpCode,rep(trpCode[i],length(sample.disc[[i]]$age)))
      }
	land.lenCls<-c(land.lenCls,sample.land[[i]]$l.age)
	land.age<-c(land.age,sample.land[[i]]$age)
	land.wt<-c(land.wt,sample.land[[i]]$wt)
	land.trpCode<-c(land.trpCode,rep(trpCode[i],length(sample.land[[i]]$age)))
}
catchCat<-c(rep("LAN",length(land.age)),rep("DIS",length(disc.age)))

ca.lenCls<-c(land.lenCls,disc.lenCls)
ca.age<-c(land.age,disc.age)
ca.wt<-c(land.wt,disc.wt)
ca.trpCode<-c(land.trpCode,disc.trpCode)
ca.spp<-rep(species,length(catchCat))
ca<-list(lenCls=ca.lenCls,age=ca.age,wt=ca.wt,trpCode=ca.trpCode,
	catchCat=catchCat,spp=ca.spp)

#-------------------------------------------------------------------
#hh
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
make.null.sum<-function(list,n){
x<-rep(0,n)
  for(i in 1:length(list)){
    if(!is.null(list[[i]]))x<-x+list[[i]]}
x
}
###################################################################
get.unobs.cell.params<-function(params,ncat,narea,simno){
model<-params$model$Int
unobs.cell.params<-NULL
if(model$cell){
ncov<-model$seas+model$area+model$gear
refmat<-unique(params$data$mcov$Int$cov[,2:(ncov+2)])
nobscell<-nrow(refmat)
year<-rep(1,nobscell)
if(model$seas)seas<-refmat[,1]
else seas<-year
if(model$gear)gear<-refmat[,model$seas+model$gear]
else gear<-year
if(model$area)area<-refmat[,model$seas+model$gear+model$area]
else area<-year

acov = list(year=year,seas=seas,gear=gear,area=area)
foo = make.cell.new(acov,model,ncat=ncat,nArea=narea,cond.u=TRUE) 
ncell = dim(foo$V.u)[1]
V.u.eig = eigen(foo$V.u,symmetric=TRUE)
r = sum(V.u.eig$val> (dim(foo$V.u)[1]*V.u.eig$val[1]*1e-10))
m<-sqrt(V.u.eig$val[1:r])
if(r>1){m<-sqrt(V.u.eig$val[1:r])
  C.u = V.u.eig$vec[,1:r]%*%diag(m)
unobs.cell.params = C.u%*%rnorm(r)}
if(r==1) {C.u<- V.u.eig$vec[,1]*m
unobs.cell.params<-C.u*rnorm(1)}
if(r==0)unobs.cell.params<-rep(0,ncell)

d<-dim(params$Int$eff$cell)
if(is.null(d))obs.cell.params<-params$Int$eff$cell
if(length(d)==2)obs.cell.params<-params$Int$eff$cell[,simno]
if(length(d)==3)obs.cell.params<-as.vector(t(params$Int$eff$cell[,,simno]))
Eu<-foo$E.u%*%obs.cell.params
unobs.cell.params = matrix(Eu+unobs.cell.params,ncol=ncat)
}
unobs.cell.params
}
###################################################################
