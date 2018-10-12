#=====================================================================
library(spatstat)
library(tidyverse)
library(readxl)
library(coda)
# Import functions
source("SMR_functions.r")

#------------------------------------------------------------------
# Gudgenby Data
#-----------------------------------------------------------------
hist<- read_excel("data/Camera Detections modified.xlsx",sheet="Cats",skip=2)
IDS<- read_excel("data/Camera Detections modified.xlsx",sheet="catID")
cams<- read_excel("data/gudgenby_camsUTM_revised.xlsx")


#------------------------------------
history<- hist[,-1]

IDS<- IDS %>% filter(!(Age %in% c("Pup")))

IDk<- IDS %>% filter(Status=="known") %>% dplyr::select(ID)
IDu<- IDS %>% filter(Status=="unknown") %>% dplyr::select(ID)
IDk<- IDk$ID
IDu<- IDu$ID

h1<- build.capt(history, IDk, 1)
Yu<- build.sight(history, IDu, 1)

Yk<- capt.array(h1, dim(history)[2], dim(history)[1])


#------
win.graph(10,12)
par(mfrow=c(3,4), mar=c(0,0,0,0))
for(i in 1:12) {
  n<- rowSums(Yk[i,,])
  plot(foxlocs,main="")
  plot(foxlocs[n>0],pch=16,add=T)
  
}

#---------------------------------------------

camlocs<- cams %>% dplyr::select(x,y)
camlocs<- camlocs/1000 # km
cam.region<- ripras(camlocs$x,camlocs$y) #convex polygon around locations
cam.region<- dilation(cam.region, 5)  # add 5km buffer
camlocs<- as.ppp(camlocs,W=cam.region)
locs<- coords(camlocs)

M<-50
mmax<- 50
delta=c(0.1, 0.1, 2)
sigma.prior<- list("uniform",0, 5)
lam0.prior<- list("beta",2, 2)
xlim<- range(locs[,1])
ylim<- range(locs[,2])
#--------
#Unknown number of marks
inits<- function(){list(S=cbind(runif(M+mmax, xlim[1], xlim[2]), 
                        runif(M+mmax, ylim[1], ylim[2])), lam0=runif(1, 0.01, 0.5),
                        sigma=runif(1, 0.4, 2), psi=runif(1, 0.4, 0.6),psim=runif(1, 0.4, 0.6))}

mod<- SMR_dens(n=Yu, X=locs, y=Yk, M=M, mmax=mmax, obsmod = "bern", 
               niters=10000, region=cam.region, sigma.prior=sigma.prior,lam0.prio=lam0.prior,
               inits=inits(), delta=delta)

inds<- 2000:10000
out<- mod[inds,]

N<- out[,7]
A<- area.owin(cam.region)
D<- N/A	

mean(N)
median(N)  # Posterior mode of N
sd(N)
quantile(N, c(0.025,0.975))

mean(D)
median(D)  # Posterior mode of N
sd(D)
quantile(D, c(0.025,0.975))

mean(out[,2])
median(out[,2])  # Posterior mode of N
sd(out[,2])
quantile(out[,2], c(0.025,0.975))

mean(out[,1])
median(out[,1])  # Posterior mode of N
sd(out[,1])
quantile(out[,1], c(0.025,0.975))


win.graph(10,10)
par(mfrow=c(2,2),pty='s')
hist(N,nclass=50,main="",xlim=c(0,150),xlab="Parameter value",probability=T,las=0)
mtext(expression(paste("Abundance (",italic(hat(N)),")")))

hist(D,nclass=50,main="",xlim=c(0,0.25),xlab="Parameter value",probability=T,las=0)
mtext(expression(paste("Density (",italic(hat(D)),")")))

hist(out[,2],nclass=30,main="",xlim=c(0,0.1),xlab="Parameter value",probability=T)
mtext(expression(paste(g(0))))

hist(out[,1],nclass=30,main="",xlim=c(0,10),xlab="Parameter value",probability=T)
plot(function(x) dunif(x, 0, 10), 0, 10, xlim=c(0, 10),add=T,col="grey80",lty=2)
mtext(expression(paste(sigma)))
#=========================================
#
# Run parallel chains
#
#=========================================
library(parallel)

ni<- 20000
nb<- 10000
nc<- 3

ncores<- nc
logfile<- "mylog.txt"
cl<- makePSOCKcluster(ncores,outfile=logfile)
clusterSetRNGStream(cl)
clusterEvalQ(cl, {library(spatstat)})
clusterExport(cl, c("SMR_dens","calc.dist","prior.density","sample.prior"))

mod<- clusterCall(cl, SMR_dens, n=Yu, X=locs, y=Yk, M=M, mmax=mmax, obsmod = "bern", 
                  niters=ni, region=cam.region, sigma.prior=sigma.prior,lam0.prior=lam0.prior,
                  inits=inits(), delta=delta)

stopCluster(cl)

for(i in 1:nc) {
  mod[[i]]<- mod[[i]][-c(1:nb),]
  mod[[i]]<- as.mcmc(mod[[i]])
}

mod<- mcmc.list(mod)

Rhat<- gelman.diag(mod, transform = T, multivariate = F)
Summ<- summary(mod)
Rhat
Summ

library(ggmcmc)
library(gridExtra)
S<- ggs(mod)

win.graph(12,12)
p1<- ggs_density(S, family="N")
p2<- ggs_density(S, family="psi")
p3<- ggs_density(S, family="lam0")
p4<- ggs_density(S, family="sigma")

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)

