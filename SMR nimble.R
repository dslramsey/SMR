# JAGS model
library(readxl)
library(spatstat)
library(nimble)
library(dplyr)

#================================================================================
#
# Gudgenby Dingo density estimation
#
#================================================================================
source("Nimble Functions.r") # required for nimble model
source("SMR_functions.r")

hist<- read_excel("data/Camera Detections modified.xlsx",sheet="Wild dogs",skip=2)
IDS<- read_excel("data/Camera Detections modified.xlsx",sheet="dogID")
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

#----------------------------------------------
M<-50
mmax<- 50
buffer<- 5

locs<- cams %>% dplyr::select(x,y)
locs<- locs/1000 # km

xlim<- range(locs[,1])
ylim<- range(locs[,2])
xlim[1]<- xlim[1] - buffer
xlim[2]<- xlim[2] + buffer
ylim[1]<- ylim[1] - buffer
ylim[2]<- ylim[2] + buffer

A<- (xlim[2]-xlim[1]) * (ylim[2]-ylim[1])

ym<- apply(Yk, c(1,2), sum) 

n<- rowSums(Yu)
get.k <- function(x) length(x[!is.na(x)])
K<- apply(Yu, 1, get.k)
X<- as.matrix(locs)
m<- nrow(ym) + mmax
J<- nrow(Yu)

yaug<- matrix(0,mmax,J)
ym<- rbind(ym,yaug)

#-----------------------------------
code<- nimbleCode({
    
    for(i in 1:m) {
      wm[i] ~ dbern(psim)
      Sm[i,1] ~ dunif(xlim[1],xlim[2])
      Sm[i,2] ~ dunif(ylim[1],ylim[2])
      d2m[i,1:J]<- (Sm[i,1] - X[1:J,1])^2 + (Sm[i,2] - X[1:J,2])^2
      probm[i,1:J]<- g0 * exp(-d2m[i,1:J]/2/sigma^2) * wm[i]
      ym[i,1:J] ~ dbin_by_row(probm[i,1:J], K[1:J])
      
    }
    
    for(i in 1:M) {
      w[i] ~ dbern(psi)
      S[i,1] ~ dunif(xlim[1],xlim[2])
      S[i,2] ~ dunif(ylim[1],ylim[2])
    
      d2[i,1:J]<- (S[i,1] - X[1:J,1])^2 + (S[i,2] - X[1:J,2])^2
      prob[i,1:J]<- g0 * exp(-d2[i,1:J]/2/sigma^2) * w[i]
    }
    
    for(j in 1:J) {
      Ptrap[j]<- 1-prod(1-prob[1:M,j])
    }
    n[1:J] ~ dbin_by_row(Ptrap[1:J],K[1:J]) 
  
    sigma ~ dinvgamma(2.27, 2.14)
    #sigma ~ dunif(0,10)
    g0 ~ dbeta(2, 2)
    psi ~ dbeta(2,2)
    psim  ~ dbeta(2,2)
    Nu<- sum(w[1:M])
    Nm<- sum(wm[1:m])
    N<- Nu + Nm
    D<- N/A	
})

constants <- list(M=M,m=m,J=J)

data<- list(ym=ym,n=n,X=X,K=K,xlim=xlim,ylim=ylim,A=A)

sx.init<- runif(M, xlim[1], xlim[2])
sy.init<- runif(M, ylim[1], ylim[2])
S.init<- cbind(sx.init,sy.init)

smx.init<- runif(m, xlim[1], xlim[2])
smy.init<- runif(m, ylim[1], ylim[2])
Sm.init<- cbind(smx.init,smy.init)

inits<- list(sigma=1,g0=0.1,psi=0.5,psim=0.5,wm=c(rep(1,dim(Yk)[1]),rep(0,mmax)),
              w=rbinom(M,1,0.5),Sm=Sm.init,S=S.init)


## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check = TRUE)
Rmcmc<- compileNimble(Rmodel)
#------------------------------------------------------------------------------------------
mcmcspec <- configureMCMC(Rmodel, monitors=c("D","N","Nm","Nu","g0","sigma","w","wm","Sm","S"),
                          onlySlice = T)


mcmcspec$removeSamplers("w", print=FALSE)
mcmcspec$removeSamplers('S', print=FALSE)
mcmcspec$removeSamplers("wm", print=FALSE)
mcmcspec$removeSamplers('Sm', print=FALSE)

wNodes <- Rmodel$expandNodeNames("w")
sNodePairs <- split(matrix(Rmodel$expandNodeNames('S'), ncol = 2), 1:M)
for(i in seq_along(wNodes)) mcmcspec$addSampler(target = wNodes[i], type = custom_zs_sampler, 
                                                control = list(sNodes = sNodePairs[[i]], 
                                                birthProb = 0.5, deathProb = 0.5 ))

for(i in seq_along(sNodePairs)) mcmcspec$addSampler(target = sNodePairs[[i]], 
                                                    type = custom_s_with_indicator_sampler, 
                                                    control = list(scaleFactor = 0.5, 
                                                    scaleNode = 'sigma', indicatorNode = wNodes[i]))

wNodes <- Rmodel$expandNodeNames("wm")
sNodePairs <- split(matrix(Rmodel$expandNodeNames('Sm'), ncol = 2), 1:m)
for(i in seq_along(wNodes)) mcmcspec$addSampler(target = wNodes[i], type = custom_zs_sampler, 
                                                control = list(sNodes = sNodePairs[[i]], 
                                                 birthProb = 0.5, deathProb = 0.5 ))

for(i in seq_along(sNodePairs)) mcmcspec$addSampler(target = sNodePairs[[i]], 
                                                    type = custom_s_with_indicator_sampler, 
                                                    control = list(scaleFactor = 0.5, 
                                                    scaleNode = 'sigma', indicatorNode = wNodes[i]))


Cmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Cmcmc, project = Rmodel, resetFunctions = TRUE)

#========================
# code for binary sampler and block sampling for S
mcmcspec <- configureMCMC(Rmodel, monitors=c("N","D","g0","sigma","psi"), onlySlice = T)

mcmcspec$removeSamplers("w", print = FALSE)
Nodes <- Rmodel$expandNodeNames("w")
for(Node in Nodes) mcmcspec$addSampler(target = Node, type = "binary", print=FALSE)
## remove the default samplers for s
mcmcspec$removeSamplers('S', print = FALSE)
sNodePairs <- split( matrix(Rmodel$expandNodeNames('S'), ncol = 2), 1:M )
for(i in seq_along(sNodePairs)) mcmcspec$addSampler(target = sNodePairs[[i]], type = 'RW_block', 
                                                    control = list(adaptScaleOnly = TRUE), print=FALSE)
Cmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Cmcmc, project = Rmodel, resetFunctions = TRUE)

#----------------------
ni<- 110000
nb<- 10000
nt<- 10
nc<- 3

inits <- function(){list(sigma=1,g0=0.1,psi=0.5,psim=0.5,wm=c(rep(1,dim(Yk)[1]),rep(0,mmax)),
                                 w=rbinom(M,1,0.5),Sm=Sm.init,S=S.init)}

samp<- runMCMC(Cmodel, niter = ni, nburnin = nb, nchains = nc, thin=nt, inits = inits,
               samplesAsCodaMCMC = T)

dogs<- samp
saveRDS(dogs, "dogs.rds")

library(MCMCvis)
library(ggmcmc)
library(gridExtra)

MCMCsummary(dogs, params=c("D","N","Nm","Nu","g0","sigma"), digits=3, n.eff=TRUE,
            func=mymode, func_name = "Mode")

win.graph(12,12)
p1<- ggs_density(S, family="N") + xlim(0, 50)
p2<- ggs_density(S, family="D") + xlim(0, 0.5)
p3<- ggs_density(S, family="g0") + xlim(0, 1)
p4<- ggs_density(S, family="sigma") + xlim(0, 10)

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)


win.graph(12,12)
p1<- ggs_traceplot(S, family="N")
p2<- ggs_traceplot(S, family="D")
p3<- ggs_traceplot(S, family="g0")
p4<- ggs_traceplot(S, family="sigma")

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)

win.graph(12,7)
p1<- ggs_autocorrelation(S, family="g0")
p2<- ggs_autocorrelation(S, family="sigma")

grid.arrange(p1, p2, ncol=2)

#============================================================================
#
# Gudgenby cats density estimation
#
#============================================================================

hist<- read_excel("data/Camera Detections modified.xlsx",sheet="Cats",skip=2)
IDS<- read_excel("data/Camera Detections modified.xlsx",sheet="catID")
cams<- read_excel("data/gudgenby_camsUTM_revised.xlsx")

#------------------------------------
history<- hist[,-1]

IDk<- IDS %>% filter(Status=="known") %>% dplyr::select(ID)
IDu<- IDS %>% filter(Status=="unknown") %>% dplyr::select(ID)
IDk<- IDk$ID
IDu<- IDu$ID

h1<- build.capt(history, IDk, 1)
Yu<- build.sight(history, IDu, 1)

Yk<- capt.array(h1, dim(history)[2], dim(history)[1])

#----------------------------------------------
M<-100
mmax<- 100
buffer<- 5

xlim<- range(locs[,1])
ylim<- range(locs[,2])
xlim[1]<- xlim[1] - buffer
xlim[2]<- xlim[2] + buffer
ylim[1]<- ylim[1] - buffer
ylim[2]<- ylim[2] + buffer

A<- (xlim[2]-xlim[1]) * (ylim[2]-ylim[1])

ym<- apply(Yk, c(1,2), sum) 

n<- rowSums(Yu)
get.k <- function(x) length(x[!is.na(x)])
K<- apply(Yu, 1, get.k)
X<- as.matrix(locs)
m<- nrow(ym) + mmax
J<- nrow(Yu)

yaug<- matrix(0,mmax,J)
ym<- rbind(ym,yaug)

#----------------------------
constants <- list(M=M,m=m,J=J)

data<- list(ym=ym,n=n,X=X,K=K,xlim=xlim,ylim=ylim,A=A)

sx.init<- runif(M, xlim[1], xlim[2])
sy.init<- runif(M, ylim[1], ylim[2])
S.init<- cbind(sx.init,sy.init)

smx.init<- runif(m, xlim[1], xlim[2])
smy.init<- runif(m, ylim[1], ylim[2])
Sm.init<- cbind(smx.init,smy.init)

inits<- list(sigma=1,g0=0.1,psi=0.5,psim=0.5,wm=c(rep(1,dim(Yk)[1]),rep(0,mmax)),
             w=rbinom(M,1,0.5),Sm=Sm.init,S=S.init)

#-----------------------------
## 

constants <- list(M=M,m=m,J=J)

data<- list(ym=ym,n=n,X=X,K=K,xlim=xlim,ylim=ylim,A=A)

sx.init<- runif(M, xlim[1], xlim[2])
sy.init<- runif(M, ylim[1], ylim[2])
S.init<- cbind(sx.init,sy.init)

smx.init<- runif(m, xlim[1], xlim[2])
smy.init<- runif(m, ylim[1], ylim[2])
Sm.init<- cbind(smx.init,smy.init)

inits<- list(sigma=1,g0=0.1,psi=0.5,psim=0.5,wm=c(rep(1,dim(Yk)[1]),rep(0,mmax)),
             w=rbinom(M,1,0.5),Sm=Sm.init,S=S.init)


## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check = TRUE)
Rmcmc<- compileNimble(Rmodel)
#------------------------------------------------------------------------------------------
mcmcspec <- configureMCMC(Rmodel, monitors=c("D","N","Nm","Nu","g0","sigma","w","wm","Sm","S"),
                          onlySlice = T)


mcmcspec$removeSamplers("w", print=FALSE)
mcmcspec$removeSamplers('S', print=FALSE)
mcmcspec$removeSamplers("wm", print=FALSE)
mcmcspec$removeSamplers('Sm', print=FALSE)

wNodes <- Rmodel$expandNodeNames("w")
sNodePairs <- split(matrix(Rmodel$expandNodeNames('S'), ncol = 2), 1:M)
for(i in seq_along(wNodes)) mcmcspec$addSampler(target = wNodes[i], type = custom_zs_sampler, 
                                                control = list(sNodes = sNodePairs[[i]], 
                                                               birthProb = 0.5, deathProb = 0.5 ))

for(i in seq_along(sNodePairs)) mcmcspec$addSampler(target = sNodePairs[[i]], 
                                                    type = custom_s_with_indicator_sampler, 
                                                    control = list(scaleFactor = 0.5, 
                                                    scaleNode = 'sigma', indicatorNode = wNodes[i]))

wNodes <- Rmodel$expandNodeNames("wm")
sNodePairs <- split(matrix(Rmodel$expandNodeNames('Sm'), ncol = 2), 1:m)
for(i in seq_along(wNodes)) mcmcspec$addSampler(target = wNodes[i], type = custom_zs_sampler, 
                                                control = list(sNodes = sNodePairs[[i]], 
                                                birthProb = 0.5, deathProb = 0.5 ))

for(i in seq_along(sNodePairs)) mcmcspec$addSampler(target = sNodePairs[[i]], 
                                                    type = custom_s_with_indicator_sampler, 
                                                    control = list(scaleFactor = 0.5, 
                                                    scaleNode = 'sigma', indicatorNode = wNodes[i]))


Cmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Cmcmc, project = Rmodel, resetFunctions = TRUE)

#========================
# code for binary sampler and block sampling for S
mcmcspec <- configureMCMC(Rmodel, monitors=c("N","D","g0","sigma","psi"), onlySlice = T)

mcmcspec$removeSamplers("w", print = FALSE)
Nodes <- Rmodel$expandNodeNames("w")
for(Node in Nodes) mcmcspec$addSampler(target = Node, type = "binary", print=FALSE)
## remove the default samplers for s
mcmcspec$removeSamplers('S', print = FALSE)
sNodePairs <- split( matrix(Rmodel$expandNodeNames('S'), ncol = 2), 1:M )
for(i in seq_along(sNodePairs)) mcmcspec$addSampler(target = sNodePairs[[i]], type = 'RW_block', 
                                                    control = list(adaptScaleOnly = TRUE), print=FALSE)
Cmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Cmcmc, project = Rmodel, resetFunctions = TRUE)

#----------------------
ni<- 110000
nb<- 10000
nt<- 10
nc<- 3

inits <- function(){list(sigma=1,g0=0.1,psi=0.5,psim=0.5,wm=c(rep(1,dim(Yk)[1]),rep(0,mmax)),
                         w=rbinom(M,1,0.5),Sm=Sm.init,S=S.init)}

samp<- runMCMC(Cmodel, niter = ni, nburnin = nb, nchains = nc, thin=nt, inits = inits,
               samplesAsCodaMCMC = T)

cats<- samp
saveRDS(cats, "cats.rds")

library(MCMCvis)
library(ggmcmc)
library(gridExtra)

MCMCsummary(cats, params=c("D","N","Nm","Nu","g0","sigma"), digits=3, n.eff=TRUE,
            func=mymode, func_name = "Mode")

S<- ggs(cats)

win.graph(12,12)
p1<- ggs_density(S, family="N") + xlim(0, 100)
p2<- ggs_density(S, family="D") + xlim(0, 0.5)
p3<- ggs_density(S, family="g0") + xlim(0, 0.1)
p4<- ggs_density(S, family="sigma") + xlim(0, 10)

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)


win.graph(12,12)
p1<- ggs_traceplot(S, family="N")
p2<- ggs_traceplot(S, family="D")
p3<- ggs_traceplot(S, family="g0")
p4<- ggs_traceplot(S, family="sigma")

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)

win.graph(12,7)
p1<- ggs_autocorrelation(S, family="g0")
p2<- ggs_autocorrelation(S, family="sigma")

grid.arrange(p1, p2, ncol=2)
