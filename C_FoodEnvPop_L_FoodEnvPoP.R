##======================
rm(list=ls());
library(gtools); ## for permutations function
##======================

## calculate power set
powerset= function(x) 
{
  sets = lapply(1:(length(x)), function(i) combn(x, i, simplify = F))
  unlist(sets, recursive = F)
}
###calculate means for partitions
meanPART=function(PART)   ## PART:PARTITION
{
n=length(PART);
mean=vector(,n);
mean=unlist(lapply(PART,mean));
return(mean);
}


##############======================================#################
##############           DATA   AANLYSIS            #################
##############======================================#################

###-----------         data preparasion   -----------------------####         

FOOD0=read.csv(file.choose()); ##choose the csv file GFSI
ENV0=read.csv(file.choose());  ##choose the csv file EPI
POP0=read.csv(file.choose());  ##choose the csv file POP
FOOD=data.matrix(FOOD0);       ##convert csv table into a matrix format
ENV=data.matrix(ENV0);
POP=data.matrix(POP0);

COUNTRY=FOOD0[,1:2];           ##95 country names and codes
FOOD=FOOD[,3:11];            ##only numeric columns are fetched
ENV=ENV[,3:11];
POP=POP[,3:11];

##--calculate yearly net changes and growth rates of food security and envrionment
n1=dim(FOOD)[1];n2=dim(FOOD)[2];
diffFOOD=matrix(nrow=n1,ncol=n2-1);      ## yearly net changes
diffENV=matrix(nrow=n1,ncol=n2-1);
for(i in 1:n1)
{
   diffFOOD[i,]=diff(FOOD[i,]);
   diffENV[i,]=diff(ENV[i,]);
}
gr_FOOD=diffFOOD/FOOD[,1:8];              ##gr: yearly growtn rates
gr_ENV=diffENV/ENV[,1:8];

##--annual growth rate of population
difPOP=matrix(nrow=dim(POP)[1],ncol=dim(POP)[2]-1);
for(j in 1:(dim(POP)[2]-1))
{
difPOP[,j]=POP[,j+1]/POP[,j]-1;        ## yearly growth rates of pop
}

##--calculate the range based on growth rates of Food and Environment
min_FOOD=apply(gr_FOOD,2,min);
max_FOOD=apply(gr_FOOD,2,max);
min_ENV=apply(gr_ENV,2,min);
max_ENV=apply(gr_ENV,2,max);

###-----  partiton the real data to obtain parameters based on even scales ------###
                               
Part_FE=function(L,T) ## L:the number(or levels) of compartments for a partition
{
PARTITIONS_FOOD=as.list(numeric(L*T));dim(PARTITIONS_FOOD)=c(L,T);   ## partition the countries based on L-levelled food index
PARTITIONS_ENV=as.list(numeric(L*T));dim(PARTITIONS_ENV)=c(L,T);     ## partition the countries based on L-levelled env. index
for(t in 1:T)
{  
int_FOOD=seq(min_FOOD[t],max_FOOD[t],by=(max_FOOD[t]-min_FOOD[t])/(L-1)); ##int:intervals
cat_FOOD=findInterval(gr_FOOD[,t],int_FOOD);                       ##cat:categories
int_ENV=seq(min_ENV[t],max_ENV[t],by=(max_ENV[t]-min_ENV[t])/(L-1)); 
cat_ENV=findInterval(gr_ENV[,t],int_ENV);
   for(i in 1:L)                                                    ##partitioning
   {
   PARTITIONS_FOOD[[i,t]]=which(cat_FOOD==i);
   PARTITIONS_ENV[[i,t]]=which(cat_ENV==i);
   ##print(unlist(PARTITIONS_ENV[[i,t]]));
   }
}
return(list(PARTITIONS_FOOD,PARTITIONS_ENV));
}
#for example, a=Part_FE(L,T);a;


##--- the intersected country codes for the partitions based on food and env. ----###
## (some of them may be empty, or "NaN" or integer(0) in R)
INTERPartIndex=function(L,T) 
{
PARTITIONS_FOOD=Part_FE(L,T)[[1]];
PARTITIONS_ENV=Part_FE(L,T)[[2]];
InterPartIndex=as.list(numeric(T*L*L));dim(InterPartIndex)=c(T,L,L);
for(t in 1:T)
{
  for(i in 1:L)
  {
      for(j in 1:L)
      {
      InterPartIndex[[t,i,j]]=intersect(PARTITIONS_FOOD[[i,t]],PARTITIONS_ENV[[j,t]]);
      }
  }
}
return(InterPartIndex);
}
#for example, b=INTERPartIndex(L,T);b;


###------Size of compartments for each year's partition
SeedCompartmentSizes=function(L,T)
{
seedCompartmentSizes=list();
iNTERPartIndex=INTERPartIndex(L,T);
   for(t in 1:T)
   {
   seedCompartmentSizes[[t]]=matrix(unlist(lapply(iNTERPartIndex[t,,],length)),L,L);
   }
return(seedCompartmentSizes);
}
#for example, SeedCompartmentSizes(5,8);

##----------------(p,q) permutated Compartment sizes
PerCompSizes_t=function(t,pVec,qVvec,L,T) ##pVec, qVec are perms of {1,2,...,L}
{
seedCompartmentSizes_t=SeedCompartmentSizes(L,T)[[t]];
perCompSizes_t=seedCompartmentSizes_t[pVec,qVec];
return(perCompSizes_t);
}
#for example, pVec=c(5,3,1,2,4);qVec=c(2,5,1,4,3);L=5;T=8;t=1;
#PerCompSizes_t(t,pVec,qVvec,L,T);


###-------------Initialized Ordering wrt each year
IntialOrdering=function(L,T)
{
intialOrderingMatrix=matrix(0,23,T);
iNTERPartIndex=INTERPartIndex(L,T);
  for(t in 1:T)
  {
  intialOrderingMatrix[,t]=unlist(iNTERPartIndex[t,,]);
  }
return(intialOrderingMatrix);
}
#for example, IntialOrdering(5,8)

###----------------Permuted Compartment (countries)----------
PermCompCoun_t=function(t,pVec,qVvec,L,T)
{
permCompCoun_t=as.list(numeric(L*L));dim(permCompCoun_t)=c(L,L);
perCompSizes_t=PerCompSizes_t(t,pVec,qVvec,L,T);
index_t=cumsum(t(perCompSizes_t));
intialOrdering_t=IntialOrdering(L,T)[,t];
if(index_t[1]==0){permCompCoun_t[[1,1]]=0}else{permCompCoun_t[[1,1]]=intialOrdering_t[index_t[1]]};
  for(k in 2:(length(index_t)))
  {
  i=ceiling(k/L);
  j=k-L*(i-1);
  if(index_t[k]-index_t[k-1]==0){permCompCoun_t[[i,j]]=0}else{
  ele=intialOrdering_t[(index_t[k-1]+1):(index_t[k])];permCompCoun_t[[i,j]]=ele};
  } 
return(permCompCoun_t);
}
#for example,abc=PermCompCoun_t(t,pVec,qVvec,L,T);abc;


###-------Permuted Compartment (growth rate of population)--------
PermCompGrPop_t=function(t,pVec,qVvec,L,T)
{
permCompCoun_t=PermCompCoun_t(t,pVec,qVvec,L,T);
permCompGrPop_t=as.list(numeric(L*L));dim(permCompGrPop_t)=c(L,L);

  for(i in 1:L)
  {
  
       for(j in 1:L)
 {
       countries=permCompCoun_t[[i,j]];
       permCompCoun_t[[i,j]]=difPOP[countries,t];
       }
  }
return(permCompCoun_t);
}
#for example, bbb=PermCompGrPop_t(t,pVec,qVvec,L,T);bbb;


### --- means for permutated compartment (growth rate of population
MeanPermCompGrPop_t=function(t,pVec,qVvec,L,T)
{
PermCompGrPop_t=PermCompGrPop_t(t,pVec,qVvec,L,T);
meanPermCompGrPop_t=matrix(0,L,L);
  for(i in 1:L)
  {
     for(j in 1:L)
     {
     meanPermCompGrPop_t[i,j]=mean(PermCompGrPop_t[[i,j]]);
     }
  }
meanPermCompGrPop_t[is.nan(meanPermCompGrPop_t)]=0;
return(meanPermCompGrPop_t);
}
#for example, ccc=MeanPermCompGrPop_t(t,pVec,qVvec,L,T);ccc;

###----Food factor mean
FoodfactorMean_t=function(t,pVec,qVvec,L,T)
{
PermCompGrPop_t=PermCompGrPop_t(t,pVec,qVvec,L,T);
foodGroupSum_t=rep(0,L);
perCompSizes_t=rowSums(PerCompSizes_t(t,pVec,qVvec,L,T));
  for(i in 1:L)
  {
  foodGroupSum_t[i]=sum(unlist(PermCompGrPop_t[i,]))
  }
foodfactorMean_t=foodGroupSum_t/perCompSizes_t;
foodfactorMean_t[is.nan(foodfactorMean_t)]=0;
return(foodfactorMean_t);
}
#for example, g=FoodfactorMean_t(t,pVec,qVvec,L,T);g;

###----Environment factor mean
EnvfactorMean_t=function(t,pVec,qVvec,L,T)
{
PermCompGrPop_t=PermCompGrPop_t(t,pVec,qVvec,L,T);
envGroupSum_t=rep(0,L);
perCompSizes_t=colSums(PerCompSizes_t(t,pVec,qVvec,L,T));
  for(i in 1:L)
  {
  envGroupSum_t[i]=sum(unlist(PermCompGrPop_t[,i]))
  }
envfactorMean_t=envGroupSum_t/perCompSizes_t;
envfactorMean_t[is.nan(envfactorMean_t)]=0;
return(envfactorMean_t)
}
#for example, h=EnvfactorMean_t(t,pVec,qVvec,L,T);h;


## ---- Grand mean for a partition 
GrandMean_t=function(t,pVec,qVvec,L,T)
{
permCompGrPop_t=PermCompGrPop_t(t,pVec,qVvec,L,T);
grandMean_t=mean(unlist(permCompGrPop_t));
return(grandMean_t)
}
#for example, d=GrandMean_t(t,pVec,qVvec,L,T);d;


##------calculate AAE_P_gr_t-----------
AAE_P_gr_t=function(t,pVec,qVvec,L,T)
{
aAE_P_gr_t=matrix(0,L,L);
compartMean_t=MeanPermCompGrPop_t(t,pVec,qVvec,L,T);
permCompGrPop_t=PermCompGrPop_t(t,pVec,qVvec,L,T);
    for(i in 1:L)
    {
       for(j in 1:L)
       {
         if(permCompGrPop_t[i,j]=="numeric(0)"){aAE_P_gr_t[i,j]=0}else
         {aAE_P_gr_t[i,j]=sum(abs(permCompGrPop_t[[i,j]]-compartMean_t[i,j]))};      
       }
    }
return(sum(aAE_P_gr_t)/23);
}
#for example, e=AAE_P_gr_t(t,pVec,qVvec,L,T);e;


##-------calculate AAfV_P_gr_t-------------
AAfV_P_gr_t=function(t,pVec,qVvec,L,T)
{
perCompSizes_t=rowSums(PerCompSizes_t(t,pVec,qVvec,L,T));
foodfactorMean_t=FoodfactorMean_t(t,pVec,qVvec,L,T);
grandMean_t=GrandMean_t(t,pVec,qVvec,L,T);
aAfV_P_gr_t=abs(foodfactorMean_t-grandMean_t);
Foodsubsum=sum(perCompSizes_t*aAfV_P_gr_t);
return(Foodsubsum/L);
}
#for example, k=AAfV_P_gr_t(t,pVec,qVvec,L,T);k;

##-------calculate AAeV_P_gr_t-------------
AAeV_P_gr_t=function(t,pVec,qVvec,L,T)
{
perCompSizes_t=rowSums(PerCompSizes_t(t,pVec,qVvec,L,T));
envfactorMean_t=EnvfactorMean_t(t,pVec,qVvec,L,T);
grandMean_t=GrandMean_t(t,pVec,qVvec,L,T);
aAeV_P_gr_t=abs(envfactorMean_t-grandMean_t);
Envsubsum=sum(perCompSizes_t*aAeV_P_gr_t);
return(Envsubsum/L);
}
#for example, h=AAeV_P_gr_t(t,pVec,qVvec,L,T);h;


##-------calculate AAfeV_P_gr_t-------------
AAfeV_P_gr_t=function(t,pVec,qVvec,L,T)
{
perCompSizes_t=PerCompSizes_t(t,pVec,qVvec,L,T);
compartMean_t=MeanPermCompGrPop_t(t,pVec,qVvec,L,T);
foodfactorMean_t=FoodfactorMean_t(t,pVec,qVvec,L,T);
envfactorMean_t=EnvfactorMean_t(t,pVec,qVvec,L,T);
grandMean_t=GrandMean_t(t,pVec,qVvec,L,T);

a1=sweep(compartMean_t,1,foodfactorMean_t);a1;
a2=sweep(a1,2,envfactorMean_t);a2;
a3=a2+grandMean_t;
aAfeV_P_gr_t=sum(perCompSizes_t*abs(a3),na.rm=T);

return(aAfeV_P_gr_t/23);
}
#for example, m=AAfeV_P_gr_t(t,pVec,qVvec,L,T);m;


##---------calculate VfE_P_gr_t--------------
VfE_P_gr_t=function(t,pVec,qVvec,L,T)
{
aAE_P_gr_t=AAE_P_gr_t(t,pVec,qVvec,L,T);
aAfV_P_gr_t=AAfV_P_gr_t(t,pVec,qVvec,L,T);
vfE_P_gr_t=aAfV_P_gr_t/aAE_P_gr_t;
return(vfE_P_gr_t);
}
#for example, n=VfE_P_gr_t(t,pVec,qVvec,L,T);n;

##---------calculate VeE_P_gr_t--------------
VeE_P_gr_t=function(t,pVec,qVvec,L,T)
{
aAE_P_gr_t=AAE_P_gr_t(t,pVec,qVvec,L,T);
aAeV_P_gr_t=AAeV_P_gr_t(t,pVec,qVvec,L,T);
veE_P_gr_t=aAeV_P_gr_t/aAE_P_gr_t;
return(veE_P_gr_t);
}
#for example, q=VeE_P_gr_t(t,pVec,qVvec,L,T);q;


##---------calculate VfeE_P_gr_t--------------
VfeE_P_gr_t=function(t,pVec,qVvec,L,T)
{
aAE_P_gr_t=AAE_P_gr_t(t,pVec,qVvec,L,T);
aAfeV_P_gr_t=AAfeV_P_gr_t(t,pVec,qVvec,L,T);
vfeE_P_gr_t=aAfeV_P_gr_t/aAE_P_gr_t;
return(vfeE_P_gr_t);
}
#for example, r=VfeE_P_gr_t(t,pVec,qVvec,L,T);r;


##---------calculate VfVe_P_gr_t---------------
VfVe_P_gr_t=function(t,pVec,qVvec,L,T)
{
aAfV_P_gr_t=AAfV_P_gr_t(t,pVec,qVvec,L,T);
aAeV_P_gr_t=AAeV_P_gr_t(t,pVec,qVvec,L,T);
vfVe_P_gr_t=aAfV_P_gr_t/aAeV_P_gr_t;
return(vfVe_P_gr_t);
}
#for example, s=VfVe_P_gr_t(t,pVec,qVvec,L,T);s;

###################################################################
############-------values for parameters---------------

L=4;T=8;
vfE_P_gr_2014_Perm=veE_P_gr_2014_Perm=
vfeE_P_gr_2014_Perm=vfVe_P_gr_2014_Perm=
vfE_P_gr_2015_Perm=veE_P_gr_2015_Perm=
vfeE_P_gr_2015_Perm=vfVe_P_gr_2015_Perm=
vfE_P_gr_2016_Perm=veE_P_gr_2016_Perm=
vfeE_P_gr_2016_Perm=vfVe_P_gr_2016_Perm=
vfE_P_gr_2017_Perm=veE_P_gr_2017_Perm=
vfeE_P_gr_2017_Perm=vfVe_P_gr_2017_Perm=
vfE_P_gr_2018_Perm=veE_P_gr_2018_Perm=
vfeE_P_gr_2018_Perm=vfVe_P_gr_2018_Perm=
vfE_P_gr_2019_Perm=veE_P_gr_2019_Perm=
vfeE_P_gr_2019_Perm=vfVe_P_gr_2019_Perm=
vfE_P_gr_2020_Perm=veE_P_gr_2020_Perm=
vfeE_P_gr_2020_Perm=vfVe_P_gr_2020_Perm=
vfE_P_gr_2021_Perm=veE_P_gr_2021_Perm=
vfeE_P_gr_2021_Perm=vfVe_P_gr_2021_Perm=
rep(0,factorial(L)^2);


pm=expand.grid(1:factorial(L),1:factorial(L));
perm=permutations(L,L,1:L);perm;
for(k in 1:factorial(L)^2)
{
pVec=perm[pm[k,1],];
qVec=perm[pm[k,2],];
vfE_P_gr_2014_Perm[k]=VfE_P_gr_t(1,pVec,qVvec,L,T);
veE_P_gr_2014_Perm[k]=VeE_P_gr_t(1,pVec,qVvec,L,T);
vfeE_P_gr_2014_Perm[k]=VfeE_P_gr_t(1,pVec,qVvec,L,T);
vfVe_P_gr_2014_Perm[k]=VfVe_P_gr_t(1,pVec,qVvec,L,T);

vfE_P_gr_2015_Perm[k]=VfE_P_gr_t(2,pVec,qVvec,L,T);
veE_P_gr_2015_Perm[k]=VeE_P_gr_t(2,pVec,qVvec,L,T);
vfeE_P_gr_2015_Perm[k]=VfeE_P_gr_t(2,pVec,qVvec,L,T);
vfVe_P_gr_2015_Perm[k]=VfVe_P_gr_t(2,pVec,qVvec,L,T);

vfE_P_gr_2016_Perm[k]=VfE_P_gr_t(3,pVec,qVvec,L,T);
veE_P_gr_2016_Perm[k]=VeE_P_gr_t(3,pVec,qVvec,L,T);
vfeE_P_gr_2016_Perm[k]=VfeE_P_gr_t(3,pVec,qVvec,L,T);
vfVe_P_gr_2016_Perm[k]=VfVe_P_gr_t(3,pVec,qVvec,L,T);

vfE_P_gr_2017_Perm[k]=VfE_P_gr_t(4,pVec,qVvec,L,T);
veE_P_gr_2017_Perm[k]=VeE_P_gr_t(4,pVec,qVvec,L,T);
vfeE_P_gr_2017_Perm[k]=VfeE_P_gr_t(4,pVec,qVvec,L,T);
vfVe_P_gr_2017_Perm[k]=VfVe_P_gr_t(4,pVec,qVvec,L,T);

vfE_P_gr_2018_Perm[k]=VfE_P_gr_t(5,pVec,qVvec,L,T);
veE_P_gr_2018_Perm[k]=VeE_P_gr_t(5,pVec,qVvec,L,T);
vfeE_P_gr_2018_Perm[k]=VfeE_P_gr_t(5,pVec,qVvec,L,T);
vfVe_P_gr_2018_Perm[k]=VfVe_P_gr_t(5,pVec,qVvec,L,T);

vfE_P_gr_2019_Perm[k]=VfE_P_gr_t(6,pVec,qVvec,L,T);
veE_P_gr_2019_Perm[k]=VeE_P_gr_t(6,pVec,qVvec,L,T);
vfeE_P_gr_2019_Perm[k]=VfeE_P_gr_t(6,pVec,qVvec,L,T);
vfVe_P_gr_2019_Perm[k]=VfVe_P_gr_t(6,pVec,qVvec,L,T);

vfE_P_gr_2020_Perm[k]=VfE_P_gr_t(7,pVec,qVvec,L,T);
veE_P_gr_2020_Perm[k]=VeE_P_gr_t(7,pVec,qVvec,L,T);
vfeE_P_gr_2020_Perm[k]=VfeE_P_gr_t(7,pVec,qVvec,L,T);
vfVe_P_gr_2020_Perm[k]=VfVe_P_gr_t(7,pVec,qVvec,L,T);

vfE_P_gr_2021_Perm[k]=VfE_P_gr_t(8,pVec,qVvec,L,T);
veE_P_gr_2021_Perm[k]=VeE_P_gr_t(8,pVec,qVvec,L,T);
vfeE_P_gr_2021_Perm[k]=VfeE_P_gr_t(8,pVec,qVvec,L,T);
vfVe_P_gr_2021_Perm[k]=VfVe_P_gr_t(8,pVec,qVvec,L,T);
}



############################################
### presentations ####
###########################################
##--population growth rate~ food
par(mfrow=c(4,2));
plot(sort(vfE_P_gr_2014_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food, 2014");
plot(sort(vfE_P_gr_2015_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food, 2015");
plot(sort(vfE_P_gr_2016_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food, 2016");
plot(sort(vfE_P_gr_2017_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food, 2017");
plot(sort(vfE_P_gr_2018_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food, 2018");
plot(sort(vfE_P_gr_2019_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food, 2019");
plot(sort(vfE_P_gr_2020_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food, 2020");
plot(sort(vfE_P_gr_2021_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food, 2021");


quantile(vfE_P_gr_2014_Perm,0.95);
quantile(vfE_P_gr_2015_Perm,0.95);
quantile(vfE_P_gr_2016_Perm,0.95);
quantile(vfE_P_gr_2017_Perm,0.95);
quantile(vfE_P_gr_2018_Perm,0.95);
quantile(vfE_P_gr_2019_Perm,0.95);
quantile(vfE_P_gr_2020_Perm,0.95);
quantile(vfE_P_gr_2021_Perm,0.95);


quantile(vfE_P_gr_2014_Perm,0.9);
quantile(vfE_P_gr_2015_Perm,0.9);
quantile(vfE_P_gr_2016_Perm,0.9);
quantile(vfE_P_gr_2017_Perm,0.9);
quantile(vfE_P_gr_2018_Perm,0.9);
quantile(vfE_P_gr_2019_Perm,0.9);
quantile(vfE_P_gr_2020_Perm,0.9);
quantile(vfE_P_gr_2021_Perm,0.9);


vfE_P_gr_2014_Perm[1];
vfE_P_gr_2015_Perm[1];
vfE_P_gr_2016_Perm[1];
vfE_P_gr_2017_Perm[1];
vfE_P_gr_2018_Perm[1];
vfE_P_gr_2019_Perm[1];
vfE_P_gr_2020_Perm[1];
vfE_P_gr_2021_Perm[1];

##--------------
##--population growth rate~ environment

par(mfrow=c(4,2));
plot(sort(veE_P_gr_2014_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~env, 2014");
plot(sort(veE_P_gr_2015_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~env, 2015");
plot(sort(veE_P_gr_2016_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~env, 2016");
plot(sort(veE_P_gr_2017_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~env, 2017");
plot(sort(veE_P_gr_2018_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~env, 2018");
plot(sort(veE_P_gr_2019_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~env, 2019");
plot(sort(veE_P_gr_2020_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~env, 2020");
plot(sort(veE_P_gr_2021_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~env, 2021");


quantile(veE_P_gr_2014_Perm,0.95);
quantile(veE_P_gr_2015_Perm,0.95);
quantile(veE_P_gr_2016_Perm,0.95);
quantile(veE_P_gr_2017_Perm,0.95);
quantile(veE_P_gr_2018_Perm,0.95);
quantile(veE_P_gr_2019_Perm,0.95);
quantile(veE_P_gr_2020_Perm,0.95);
quantile(veE_P_gr_2021_Perm,0.95);


quantile(veE_P_gr_2014_Perm,0.9);
quantile(veE_P_gr_2015_Perm,0.9);
quantile(veE_P_gr_2016_Perm,0.9);
quantile(veE_P_gr_2017_Perm,0.9);
quantile(veE_P_gr_2018_Perm,0.9);
quantile(veE_P_gr_2019_Perm,0.9);
quantile(veE_P_gr_2020_Perm,0.9);
quantile(veE_P_gr_2021_Perm,0.9);


veE_P_gr_2014_Perm[1];
veE_P_gr_2015_Perm[1];
veE_P_gr_2016_Perm[1];
veE_P_gr_2017_Perm[1];
veE_P_gr_2018_Perm[1];
veE_P_gr_2019_Perm[1];
veE_P_gr_2020_Perm[1];
veE_P_gr_2021_Perm[1];

##----------------------
##--population growth rate~ interaction food and environment

par(mfrow=c(4,2));
plot(sort(vfeE_P_gr_2014_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food&env, 2014");
plot(sort(vfeE_P_gr_2015_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food&env, 2015");
plot(sort(vfeE_P_gr_2016_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food&env, 2016");
plot(sort(vfeE_P_gr_2017_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food&env, 2017");
plot(sort(vfeE_P_gr_2018_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food&env, 2018");
plot(sort(vfeE_P_gr_2019_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food&env, 2019");
plot(sort(vfeE_P_gr_2020_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food&env, 2020");
plot(sort(vfeE_P_gr_2021_Perm),xlab="sample sizes,L=3,T=8",ylab="s/n ratio",main="Pop. growth rate~food&env, 2021");


quantile(vfeE_P_gr_2014_Perm,0.95);
quantile(vfeE_P_gr_2015_Perm,0.95);
quantile(vfeE_P_gr_2016_Perm,0.95);
quantile(vfeE_P_gr_2017_Perm,0.95);
quantile(vfeE_P_gr_2018_Perm,0.95);
quantile(vfeE_P_gr_2019_Perm,0.95);
quantile(vfeE_P_gr_2020_Perm,0.95);
quantile(vfeE_P_gr_2021_Perm,0.95);


quantile(vfeE_P_gr_2014_Perm,0.9);
quantile(vfeE_P_gr_2015_Perm,0.9);
quantile(vfeE_P_gr_2016_Perm,0.9);
quantile(vfeE_P_gr_2017_Perm,0.9);
quantile(vfeE_P_gr_2018_Perm,0.9);
quantile(vfeE_P_gr_2019_Perm,0.9);
quantile(vfeE_P_gr_2020_Perm,0.9);
quantile(vfeE_P_gr_2021_Perm,0.9);


vfeE_P_gr_2014_Perm[1];
vfeE_P_gr_2015_Perm[1];
vfeE_P_gr_2016_Perm[1];
vfeE_P_gr_2017_Perm[1];
vfeE_P_gr_2018_Perm[1];
vfeE_P_gr_2019_Perm[1];
vfeE_P_gr_2020_Perm[1];
vfeE_P_gr_2021_Perm[1];

##-----------------------------

##-- food ~ environment
par(mfrow=c(4,2));
plot(sort(vfVe_P_gr_2014_Perm),xlab="sample sizes,L=3,T=8",ylab="relative ratio",main="food/env, 2014");
plot(sort(vfVe_P_gr_2015_Perm),xlab="sample sizes,L=3,T=8",ylab="relative ratio",main="food/env, 2015");
plot(sort(vfVe_P_gr_2016_Perm),xlab="sample sizes,L=3,T=8",ylab="relative ratio",main="food/env, 2016");
plot(sort(vfVe_P_gr_2017_Perm),xlab="sample sizes,L=3,T=8",ylab="relative ratio",main="food/env, 2017");
plot(sort(vfVe_P_gr_2018_Perm),xlab="sample sizes,L=3,T=8",ylab="relative ratio",main="food/env, 2018");
plot(sort(vfVe_P_gr_2019_Perm),xlab="sample sizes,L=3,T=8",ylab="relative ratio",main="food/env, 2019");
plot(sort(vfVe_P_gr_2020_Perm),xlab="sample sizes,L=3,T=8",ylab="relative ratio",main="food/env, 2020");
plot(sort(vfVe_P_gr_2021_Perm),xlab="sample sizes,L=3,T=8",ylab="relative ratio",main="food/env, 2021");


quantile(vfVe_P_gr_2014_Perm,0.95);
quantile(vfVe_P_gr_2015_Perm,0.95);
quantile(vfVe_P_gr_2016_Perm,0.95);
quantile(vfVe_P_gr_2017_Perm,0.95);
quantile(vfVe_P_gr_2018_Perm,0.95);
quantile(vfVe_P_gr_2019_Perm,0.95);
quantile(vfVe_P_gr_2020_Perm,0.95);
quantile(vfVe_P_gr_2021_Perm,0.95);


quantile(vfVe_P_gr_2014_Perm,0.9);
quantile(vfVe_P_gr_2015_Perm,0.9);
quantile(vfVe_P_gr_2016_Perm,0.9);
quantile(vfVe_P_gr_2017_Perm,0.9);
quantile(vfVe_P_gr_2018_Perm,0.9);
quantile(vfVe_P_gr_2019_Perm,0.9);
quantile(vfVe_P_gr_2020_Perm,0.9);
quantile(vfVe_P_gr_2021_Perm,0.9);


vfVe_P_gr_2014_Perm[1];
vfVe_P_gr_2015_Perm[1];
vfVe_P_gr_2016_Perm[1];
vfVe_P_gr_2017_Perm[1];
vfVe_P_gr_2018_Perm[1];
vfVe_P_gr_2019_Perm[1];
vfVe_P_gr_2020_Perm[1];
vfVe_P_gr_2021_Perm[1];


































## compute real AE, AAE 
InterPart=as.list(numeric(t*L*L));dim(InterPart)=c(t,L,L);
SSE_Part=as.list(numeric(t*L*L));dim(SSE_Part)=c(t,L,L);
mean_8=as.list(numeric(t*L*L));dim(mean_8)=c(t,L,L);
size_8=as.list(numeric(t*L*L));dim(size_8)=c(t,L,L);
sum_8=as.list(numeric(t*L*L));dim(sum_8)=c(t,L,L);
AB_8=as.list(numeric(t*L*L));dim(AB_8)=c(t,L,L);
AFa=as.list(numeric(t*L));dim(AFa)=c(t,L);      ##Abosolute value for factor A
AFb=as.list(numeric(t*L));dim(AFb)=c(t,L);
mean_i_8=as.list(numeric(t*L));dim(mean_i_8)=c(t,L);   ##mean for each level of factor A
mean_j_8=as.list(numeric(t*L));dim(mean_j_8)=c(t,L);   ##mean for each level of factor B


for(t in 1:8)
{
  for(i in 1:L)
  {
      for(j in 1:L)
      {
      ind= InterPartIndex[[t,i,j]];
          if(length(ind)!=0)
           { InterPart[[t,i,j]]=difPOP[ind]
             mean_8[[t,i,j]]=mean(InterPart[[t,i,j]]);
             size_8[[t,i,j]]=length(InterPart[[t,i,j]]);
             er=InterPart[[t,i,j]]-mean_8[[t,i,j]];
             SSE_Part[[t,i,j]]=sqrt(er*er);
             sum_8[[t,i,j]]=sum(unlist(SSE_Part[[t,i,j]]));
             exp1=mean_8[[t,i,j]]-mean(difPOP);
             AB_8[[t,i,j]]=exp1*size_8[[t,i,j]];
           }
          else
           {InterPart[[t,i,j]]=c(NaN);size_8[[t,i,j]]=c(NaN);
             SSE_Part[[t,i,j]]=c(NaN);
             sum_8[[t,i,j]]=c(NaN);
           }
      }
   AFa[[t,i]]=abs(sum(unlist(AB_8[t,i,]),na.rm=TRUE));
   AFb[[t,i]]=abs(sum(unlist(AB_8[t,,i]),na.rm=TRUE));
   mean_i_8[[t,i]]=sum(unlist(AB_8[t,i,]),na.rm=TRUE)/L;
   mean_j_8[[t,i]]=sum(unlist(AB_8[t,,i]),na.rm=TRUE)/L;
  }
}


##for SSab
cross=as.list(numeric(t*L*L));dim(cross)=c(t,L,L);
AFab=as.list(numeric(t*L*L));dim(AFab)=c(t,L,L);
for(t in 1:8)
{
  for(i in 1:L)
  {
      for(j in 1:L)
      {
      cross[[t,i,j]]=abs(mean_8[[t,i,j]]+mean(difPOP)-mean_i_8[[t,i]]-mean_j_8[[t,j]]);
      AFab[[t,i,j]]=cross[[t,i,j]]*size_8[[t,i,j]];
      }
   }
}



########## compute AAE,AAFa,AAFb,AAFab and s/n ratios


AE=vector(,8);
AAE=vector(,8);
AAFa=vector(,8);
AAFb=vector(,8);
AAFab=vector(,8);

for(t in 1:8)
{
AE[t]=sum(unlist(sum_8[t,,]),na.rm=TRUE);       ##AE:absolute error
AAE[t]=AE[t]/95;                               ##AAE:averaged absolute error
AAFa[t]=sum(unlist(as.matrix(AFa[t,])))/L;    ##AAFa:averaged absolute factor A
AAFb[t]=sum(unlist(as.matrix(AFb[t,])))/L;    ##AAFa:averaged absolute factor A
AAFab[t]=sum(unlist(as.matrix(AFab[t,,])),na.rm=TRUE)/23;
}

SaN=AAFa/AAE;    ##SaN:signal noise ration for Factor A
SbN=AAFb/AAE;    ##SaN:signal noise ration for Factor B
SabN=AAFab/AAE;  ##SaN:signal noise ration for Factor AB



###presentations of results
plot(SaN,type="b",ylim=c(0,20),xlab="Year(2014=1)",ylab="realized s/n",main="real SaN(black),SbN(red),SabN(blue)");
points(SbN,type="l", col="red");
points(SabN,type="l", col="blue");
####===============================================================#####
########            nonparametric distribution            ##############
####===============================================================#####

##------ set up the seed content matrix and seed size matrix ----------##

seed_23_8_mt=matrix(nrow=23,ncol=8);   
for(t in 1:8)
{
seedVec=vector();
for(k in 2:(L+1))
{
   for(i in 1:(k-1))
   {
      j=k-i;
      seedVec=c(seedVec,unlist(InterPartIndex[[t,i,j]]));
   }
}
for(k in (L+2):(2*L))
{
    for(i in (k-L):L)
    {
     j=k-i;
     seedVec=c(seedVec,unlist(InterPartIndex[[t,i,j]]));
    }
} 
seed_23_8_mt[,t]=seedVec;        ##  seed content matrix
}
seedSize=size_8;                 ##  seed size matrix

##------  sizes for all permutated seed size matrices, 8xfacLxfacL  -------------##

PERM_facL_L_mt=permutations(n=L,r=L,1:L,repeats.allowed=F);
Size_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(Size_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
  for(i in 1:factorial(L))
  {
       for(j in 1:factorial(L))
       {
       Size_8_facL_facL_ls[[t,i,j]]=seedSize[t,,][PERM_facL_L_mt[i,],PERM_facL_L_mt[j,]];
       }
  }
}

L=5;

##-- contents/partitioned countries for all permutated seed size matrices, 8xfacLxfacL -----------##

L_L_ls=as.list(numeric(L^2));dim(L_L_ls)=c(L,L);
Sample_ls=as.list(numeric(L^2));dim(Sample_ls)=c(L,L);

con=function(L_L_mat,sdt)   ##sdt:seed matrix for Year t
{
   for(k in 2:(L+1))
   {
        for(u in 1:(k-1))
        {v=k-u;va=unlist(L_L_mat[u,v]);
           if(va!="NaN"){L_L_ls[[u,v]]=sdt[1:va];
                      Sample_ls[[u,v]]=difPOP[,t][unlist(L_L_ls[[u,v]])];
                      sdt=setdiff(sdt,sdt[1:va]);}
           else{L_L_ls[[u,v]]=NaN;Sample_ls[[u,v]]=NaN;}
        }
    };
   for(k in (L+2):(2*L))
   {
        for(u in (k-L):L)
        {v=k-u;va=unlist(L_L_mat[u,v]);
       if(va!="NaN"){L_L_ls[[u,v]]=sdt[1:va];
                      Sample_ls[[u,v]]=difPOP[,t][unlist(L_L_ls[[u,v]])];
                      sdt=setdiff(sdt,sdt[1:va]);}
        else{L_L_ls[[u,v]]=NaN;Sample_ls[[u,v]]=NaN;}
         }
    } 
return(list(L_L_ls,Sample_ls));
}





Cont_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(Cont_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
sdt=seed_23_8_mt[,t];
     for(i in 1:factorial(L))
     {
           for(j in 1:factorial(L))
           {
           Cont_8_facL_facL_ls[[t,i,j]]=con(Size_8_facL_facL_ls[[t,i,j]],sdt)[2];
           }
     }

}


##----------compute in-group means, cross-group means, grand means-----------##

means=function(L,part_L_L_ls)          ##L: the number of compartments for each partition
{
mean_L_L_mx=matrix(nrow=L+1,ncol=L+1);
i=j=1;
  for(i in 1:L)
  {rowsum=n=0;
      for(j in 1:L)
      {
      if(part_L_L_ls[[i,j]]!="NaN"){mean_L_L_mx[i,j]=mean(part_L_L_ls[[i,j]]);
      rowsum=rowsum+sum(part_L_L_ls[[i,j]]);n=n+length(part_L_L_ls[[i,j]])}
      else{mean_L_L_mx[i,j]=NaN};
      }
   mean_L_L_mx[i,L+1]=rowsum/n;
  }

 for(j in 1:L)
{
colsum=n=0;
      for(i in 1:L)
      {
      if(part_L_L_ls[[i,j]]!="NaN"){colsum=colsum+sum(part_L_L_ls[[i,j]]);
      n=n+length(part_L_L_ls[[i,j]])};
      }
mean_L_L_mx[L+1,j]=colsum/n;
}
mean_L_L_mx[L+1,L+1]=sum(unlist(part_L_L_ls),na.rm=T)/23;
return(mean_L_L_mx);
}


mean_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(mean_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
     for(i in 1:factorial(L))
     {
           for(j in 1:factorial(L))
           {
           mean_8_facL_facL_ls[[t,i,j]]=means(L,Cont_8_facL_facL_ls[[t,i,j]][[1]]);
           }
     }
}

##----------------------- compute variances ---------------------------------##

abs_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));  ##absolute difference between/across groups
dim(abs_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
i=j=1;
     for(i in 1:factorial(L))
     {
           for(j in 1:factorial(L))
           {
           A=Cont_8_facL_facL_ls[[t,i,j]][[1]];
           B=mean_8_facL_facL_ls[[t,i,j]];
           C=as.list(numeric((L+1)^2));dim(C)=c(L+1,L+1);
                for(p in 1:L)
                {
                      for(q in 1:L)
                      {
                      if(A[[p,q]]!="NaN"){C[[p,q]]=abs(A[[p,q]]-B[1:L,1:L][p,q])}
                      else{C[[p,q]]=NaN};
                      }
                }
               
                for(u in 1:L)
                {
                if(B[u,L+1]!="NaN"){C[[u,L+1]]=abs(B[u,L+1]-B[L+1,L+1])}
                else{C[[u,L+1]]=NaN};
                }
                for(v in 1:L)
                {
                if(B[L+1,v]!="NaN"){C[[L+1,v]]=abs(B[L+1,v]-B[L+1,L+1])}
                else{C[[L+1,v]]=NaN};
                }

           abs_8_facL_facL_ls[[t,i,j]]=C;
           }
     }
}

##---------------  for product factor food x Environment  ----------------##

cross_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(cross_8_facL_facL_ls)=c(8,factorial(L),factorial(L));

for(t in 1:8)
{
    for(i in 1:factorial(L))
    {
          for(j in 1:factorial(L))
          {
          B=mean_8_facL_facL_ls[[t,i,j]];
          D=matrix(nrow=L,ncol=L);
              for(u in 1:L)
              {
                  for(v in 1:L)
                  {
                  if(B[u,v]!="NaN"){D[u,v]=abs(B[u,v]+B[L+1,L+1]-B[u,L+1]-B[L+1,v])}
                  else{D[u,v]=NaN};
                  }
              }
          cross_8_facL_facL_ls[[t,i,j]]=D;
          }
    }
}



### AAE: AVERAGE ABSOLUTE ERROR;
AAE_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(AAE_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
    for(i in 1:factorial(L))
    {
        for(j in 1:factorial(L))
        {
        AAE_8_facL_facL_ls[[t,i,j]]=sum(unlist(abs_8_facL_facL_ls[[t,i,j]]),na.rm=T)/23;
        }
    }
}


### size expanded
ExSize_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L))); ##Ex:expanded
dim(ExSize_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
    for(i in 1:factorial(L))
    {
          for(j in 1:factorial(L))
          {
          size=as.list(numeric((L+1)^2));dim(size)=c(L+1,L+1);
          size[1:L,1:L]=Size_8_facL_facL_ls[[t,i,j]];

              for(u in 1:L)
              {if(sum(unlist(size[u,1:L])=="NaN")!=L){size[[u,L+1]]=sum(unlist(size[u,1:L]),na.rm=T)}
               else{size[[u,L+1]]=NaN};
              }

              for(v in 1:L)
              {if(sum(unlist(size[1:L,v])=="NaN")!=L){size[[L+1,v]]=sum(unlist(size[1:L,v]),na.rm=T)}
               else{size[[u,L+1]]=NaN};
              }
          size[[L+1,L+1]]=sum(unlist(size[L+1,]));
          ExSize_8_facL_facL_ls[[t,i,j]]=size;
          }
    }
}




### AAVf: Average Absolute Variance for food 
AAVf_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(AAVf_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
    for(i in 1:factorial(L))
    {
          for(j in 1:factorial(L))
          {
           vec1=unlist(abs_8_facL_facL_ls[[t,i,j]][,L+1]);
           if("NaN"%in%vec1){vec1=vec1[-which(vec1=="NaN")]}else{vec1=vec1};
           vec2=unlist(ExSize_8_facL_facL_ls[[t,i,j]][,L+1]);
           if("NaN"%in%vec2){vec2=vec2[-which(vec2=="NaN")]}else{vec2=vec2};
           AAVf_8_facL_facL_ls[[t,i,j]]=vec1%*%vec2/L;
          }
     }
}

### AAVe: Average Absolute Variance for environment 
AAVe_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(AAVe_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
    for(i in 1:factorial(L))
    {
          for(j in 1:factorial(L))
          {
           vec1=unlist(abs_8_facL_facL_ls[[t,i,j]][L+1,]);
           if("NaN"%in%vec1){vec1=vec1[-which(vec1=="NaN")]}else{vec1=vec1};
           vec2=unlist(ExSize_8_facL_facL_ls[[t,i,j]][L+1,]);
           if("NaN"%in%vec2){vec2=vec2[-which(vec2=="NaN")]}else{vec2=vec2};
           AAVe_8_facL_facL_ls[[t,i,j]]=vec1%*%vec2/L;
          }
     }
}


### AAVfe: Average Absolute Variance for food x environment 
AAVfe_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(AAVfe_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
    for(i in 1:factorial(L))
    {
          for(j in 1:factorial(L))
          {
          sum=0;
          MAT=cross_8_facL_facL_ls[[t,i,j]];
          siz=Size_8_facL_facL_ls[[t,i,j]];
                for(u in 1:L)
                {
                      for(v in 1:L)
                      {
                      if(MAT[u,v]!="NaN"){sum=sum+MAT[u,v]*unlist(siz[u,v])}else{sum=sum};
                      }
                }
          AAVfe_8_facL_facL_ls[[t,i,j]]=sum/23;
          }
     }
}



### ratios: AAVf/AAE=VfE
VfE_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(VfE_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
    for(i in 1:factorial(L))
    {
          for(j in 1:factorial(L))
          {
           numerator=AAVf_8_facL_facL_ls[[t,i,j]];
           denominator= AAE_8_facL_facL_ls[[t,i,j]];
           VfE_8_facL_facL_ls[[t,i,j]]=numerator/denominator;
          }
    }
}

### ratios: AAVe/AAE=VeE
VeE_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(VeE_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
    for(i in 1:factorial(L))
    {
          for(j in 1:factorial(L))
          {
           numerator=AAVe_8_facL_facL_ls[[t,i,j]];
           denominator= AAE_8_facL_facL_ls[[t,i,j]];
           VeE_8_facL_facL_ls[[t,i,j]]=numerator/denominator;
          }
    }
}


### ratios: AAVfe/AAE=VfeE
VfeE_8_facL_facL_ls=as.list(numeric(8*factorial(L)*factorial(L)));
dim(VfeE_8_facL_facL_ls)=c(8,factorial(L),factorial(L));
for(t in 1:8)
{
    for(i in 1:factorial(L))
    {
          for(j in 1:factorial(L))
          {
           numerator=AAVfe_8_facL_facL_ls[[t,i,j]];
           denominator=AAE_8_facL_facL_ls[[t,i,j]];
           VfeE_8_facL_facL_ls[[t,i,j]]=numerator/denominator;
          }
    }
}


###----------------------- distributions ------------------------------###
food=sort(unlist(VfE_8_facL_facL_ls[1,,]));
m=max(food);
cdf=food/m;
plot(food,cdf);


env=sort(unlist(VeE_8_facL_facL_ls[1,,]));
m=max(env);
cdf=env/m;
plot(env,cdf);


food_env=sort(unlist(VeE_8_facL_facL_ls[1,,]));
m=max(food_env);
cdf=food_env/m;
plot(food_env,cdf);


###presentations of results





plot(unlist(VfE_8_facL_facL_ls[1,,]),xlab="samples"
,ylab="s/n: for food");

plot(unlist(VeE_8_facL_facL_ls[1,,]),xlab="samples"
,ylab="s/n: for environment");

plot(unlist(VfeE_8_facL_facL_ls[1,,]),xlab="samples",
ylab="s/n: the interaction between food and environment");





