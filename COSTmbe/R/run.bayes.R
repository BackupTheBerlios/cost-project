run.bayes<-function(COST.data,fit=NULL,do.predict=T,species,timeStrata="quarter",burnin=2000,
                    thin=1,nmcmc=1000,ageMin=0,ageMax=20,
                    usewglrel=T,cov.list=NULL,pred.cov.real=NULL,arealist=NULL,l.int){
  obsdata<-mldata<-NULL
if(sum(COST.data@cs@tr$sampType=="M")>0)mldata<-subset(COST.data@cs,subset=(sampType=="M"))
if(sum(COST.data@cs@tr$sampType=="S")>0)obsdata<-subset(COST.data@cs,subset=(sampType=="S")) 
if(is.null(mldata)&(is.null(obsdata))){print("no data of type M or S");return}
input.data<-read.cost.data(obsdata,mldata,species,usewglrel)
nseas<-0
if(timeStrata=="quarter")nseas<-4
if(timeStrata=="month")nseas<-12
if(nseas==0){print("invalid time strata");return}
if(nseas==4){input.data$season_obs<-1+floor((input.data$season_obs-1)/3)    
input.data$season_mland<-1+floor((input.data$season_mland-1)/3)
input.data$seas<-1+floor((input.data$seas-1)/3) }
  
input.data$sampsize_disc[is.na(input.data$sampsize_disc)]<-0
input.data$haulsize_disc[is.na(input.data$haulsize_disc)]<-0
input.data$cens.mu=c(1,30.0,3.4)
input.data$cens.tau=c(1.0,1.0,1.0)
input.data$cens.pri=c(3.4,0.001,0.0001,0.001)
  seas.code<-sim.codecov(input.data$seas)
  gear.code<-sim.codecov(input.data$gear)
  area.code<-sim.codecov(input.data$area,arealist)
  adj<-get.adj(area.code$covlist)
aobs = NULL
  acov = list(year=input.data$year,seas=seas.code$cov,
	gear=gear.code$cov,area=area.code$cov)

  neigh = list(num=adj$num,adj=adj$adj) 

  simobj=list(aobs=aobs,acov=acov,neigh=neigh)
  class(simobj) <- "caa.data"
 if(is.null(fit)){
   if(is.null(cov.list)){print("cov.list empty");return}
models<-make.reg.models(cov.list)
  input.data<<-input.data

fit = cost.fit(simobj,input.data,burnin=burnin,numit.inner=thin,numit.outer=nmcmc,constr=1,seed=3421,
               ageMin=ageMin,ageMax=ageMax,nSeason=nseas,
  agemodel=models$agemodel,
               lgamodel=models$lgamodel,lgarel="log-linear",model1=T,model2=!usewglrel)
nHaul = input.data$n_trip_obs+input.data$n_trip_mland
fit$aobs=list(data=list(nBoats=nHaul))
if(usewglrel){
wglparams <- getwgl(species)
if(nseas==4)wglseas<-c(mean(wglparams$seas[1:3]),mean(wglparams$seas[4:6]),mean(wglparams$seas[7:9]),mean(wglparams$seas[10:12]))
else wglseas<-wglparams$seas
Int <- list(eff = list(const = wglparams$const, year = NULL, seas = wglseas, gear = NULL, area = NULL, cell = NULL),
            tau.area = NULL, ar = NULL, tau.cell = NULL, tau.haul = NULL)
Slp <- list(eff = list(const = wglparams$slope, year = NULL, seas = NULL, gear = NULL, area = NULL, cell = NULL),
            tau.area = NULL, ar = NULL, tau.cell = NULL, tau.haul = NULL)
fit<-insert.wgl.param(fit,Int,Slp,tau.obs=1000000) 
}
}
predict<-dbeObject.list<-NULL
if(do.predict){
seas.pred<-sim.codecov(pred.cov.real[,1],seas.code$covlist)
gear.pred<-sim.codecov(pred.cov.real[,2],gear.code$covlist)
area.pred<-sim.codecov(pred.cov.real[,3],area.code$covlist)
year.pred <-rep(1,length(seas.pred$cov))
  landings<-NULL
  if(timeStrata=="quarter")land.seas<-COST.data@cl@cl$quarter
  if(timeStrata=="month")land.seas<-COST.data@cl@cl$month
if(timeStrata!="quarter"&timeStrata!="month"){print("invalid time strata");return} 
  for(i in 1:length(seas.pred$cov)){ind<-(land.seas==seas.pred$cov[i])&
                 (COST.data@cl@cl$foCatEu5==gear.pred$cov[i])&
                 (COST.data@cl@cl$area==area.pred$cov[i])
    landings<-c(landings,COST.data@cl@cl$landWt[ind])}

landings<-1000*landings
predict <-  predict.fit.COST(fit,fit$COST.list,year.pred,seas.pred$cov,gear.pred$cov,
          area.pred$cov,landings,seas.pred$cov,
          t2.year=NULL,t2.seas=NULL,t2.gear=NULL,t2.area=NULL,
          burnin=0,nMC=100,l.int=l.int,
          par.haulsize=NULL)
dbeObject.list<-mbe2dbe(predict,species)
#if(nrow(pred.cov.real)==1) dbeObject.list<-mbe2dbe(predict,species)
#else dbeObject.list<-NULL
}
list(fit=fit,predict=predict,dbeObject.list=dbeObject.list)
}

#########################################################################
make.reg.models<-function(cov.list){
agemodel<-lgamodel<-wglmodel<-vector("list")
agemodel$Int<-lgamodel$Int<-wglmodel$Int<-lgamodel$Slp<-wglmodel$Slp<-vector("list")
agemodel$Int$year<-F
agemodel$Int$seas<-cov.list$ageseas
agemodel$Int$gear<-cov.list$agegear
agemodel$Int$area<-cov.list$agearea
agemodel$Int$cell<-F
agemodel$Int$haul<-T
agemodel$Hsz<-NULL
lgamodel$Int$year<-F
lgamodel$Int$seas<-cov.list$lgaseas
lgamodel$Int$gear<-cov.list$lgagear
lgamodel$Int$area<-cov.list$lgaarea
lgamodel$Int$cell<-F
lgamodel$Int$haul<-F
lgamodel$Hsz<-NULL
lgamodel$Slp$year<-F
lgamodel$Slp$seas<-F
lgamodel$Slp$gear<-F
lgamodel$Slp$area<-F
lgamodel$Slp$cell<-F
lgamodel$Slp$haul<-F
wglmodel$Int$year<-F
wglmodel$Int$seas<-cov.list$wglseas
wglmodel$Int$gear<-cov.list$wglgear
wglmodel$Int$area<-cov.list$wglarea
wglmodel$Int$cell<-F
wglmodel$Int$haul<-F
wglmodel$Hsz<-NULL
wglmodel$Slp$year<-F
wglmodel$Slp$seas<-F
wglmodel$Slp$gear<-F
wglmodel$Slp$area<-F
wglmodel$Slp$cell<-F
wglmodel$Slp$haul<-F
#if((cov.list$ageseas+cov.list$agegear+cov.list$agearea)>1)agemodel$Int$cell<-T
#if((cov.list$lgaseas+cov.list$lgagear+cov.list$lgaarea)>1)lgamodel$Int$cell<-T
#if((cov.list$wglseas+cov.list$wglgear+cov.list$wglarea)>1)wglmodel$Int$cell<-T
list(agemodel=agemodel,lgamodel=lgamodel,wglmodel=wglmodel)
}
#########################################################################

get.adj<-function(areanames)
{
  newnames<-roman.sub(unique(areanames))
  adj<-num<-NULL
  areas.ok<-unique(unlist(adjacent.areas))
  if(sum(is.element(areanames,areas.ok))==length(areanames)){
  for(a in newnames){
    adj1<-adjacent.areas[[a]]
    adj<-c(adj,intersect(adj1,newnames))
    num<-c(num,length(intersect(adj1,newnames)))}
  newadj<-rep(NA,length(adj))
  for(i in 1:length(unique(areanames)))newadj[adj==unique(newnames)[i]]<-i
adj<-newadj
  if(sum(num==0)>0)adj<-num<-NULL}
  else {
n<-length(areanames)
    num<-rep(2,n)
adj<-c(n,as.vector(matrix(c(2:n,1:(n-1)),nrow=2,byrow=T)),1)
  }
list(adj=adj,num=num,areanames=unique(areanames))
}
#########################################################################

roman.sub<-function(name){
  name<-sub("XIV","27.14.",name)
  name<-sub("VIIII","27.9.",name)
  name<-sub("VIII","27.8.",name)
  name<-sub("VII","27.7.",name)
  name<-sub("VI","27.6.",name)
  name<-sub("XIIII","27.14.",name)
  name<-sub("XIII","27.13.",name)
  name<-sub("XII","27.12.",name)
  name<-sub("XI","27.11.",name)
  name<-sub("IX","27.9.",name)
  name<-sub("IV","27.4.",name)
  name<-sub("V","27.5.",name)
  name<-sub("b1","b.1",name)
  name
}
    
#########################################################################
sim.codecov<-function(cov,covlist=NULL){
  newcov<-rep(NA,length(cov))
  if(is.null(covlist))covlist<-unique(cov)
  for (i in 1:length(covlist))newcov[cov==covlist[i]]<-i
  list(cov=newcov,covlist=covlist)
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
if(!is.null(dim(params$tau))){
  area.tau<-params$tau[1,]
  cell.tau<-params$tau[2,]
  haul.tau<-params$tau[3,]}
else {
  area.tau<-params$tau[1]
  cell.tau<-params$tau[2]
  haul.tau<-params$tau[3]}
}
if(sum(taulist)==2){
  id<-cumsum(taulist)
if(!is.null(dim(params$tau))){  
 if(model$area)area.tau<-params$tau[1,]
 if(model$cell) cell.tau<-params$tau[id[2],]
 if(model$haul)haul.tau<-params$tau[id[3],]}
  else {
 if(model$area)area.tau<-params$tau[1]
 if(model$cell) cell.tau<-params$tau[id[2]]
 if(model$haul)haul.tau<-params$tau[id[3]]}
}
params$area.tau<-area.tau
params$cell.tau<-cell.tau
params$haul.tau<-haul.tau
params
}
