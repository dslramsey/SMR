library(nimble)
library(readxl)
library(spatstat)

source("Nimble Functions.r") # required for nimble model

foxobs<- read.csv("data/gudgenby foxes.csv")
cams<- read_excel("data/gudgenby_camsUTM_revised.xlsx")

foxobs<- as.matrix(foxobs[,-1])

locs<- cams %>% dplyr::select(x,y)
locs<- locs/1000 # km

xlim<- range(locs[,1])
ylim<- range(locs[,2])
buffer<- 2
xlim[1]<- xlim[1] - buffer
xlim[2]<- xlim[2] + buffer
ylim[1]<- ylim[1] - buffer
ylim[2]<- ylim[2] + buffer

A<- (xlim[2]-xlim[1]) * (ylim[2]-ylim[1])

n<- rowSums(foxobs,na.rm=T)

X<- as.matrix(locs)
M=300
J<- nrow(locs)
get.k <- function(x) length(x[!is.na(x)])
K<- apply(foxobs, 1, get.k)
#=================================================
## define the model
code <- nimbleCode({
 
    for(i in 1:M) {
    w[i] ~ dbern(psi)
    S[i,1] ~ dunif(xlim[1],xlim[2])
    S[i,2] ~ dunif(ylim[1],ylim[2])
    
   
    d2[i,1:J]<- (S[i,1] - X[1:J,1])^2 + (S[i,2] - X[1:J,2])^2
    prob[i,1:J]<- exp(-d2[i,1:J]/2/sigma^2) * w[i]
    
    }
    
  for(j in 1:J) {
    Ptrap[j]<- g0*(1-prod(1-prob[1:M,j]))
  }
  n[1:J] ~ dbin_by_row(Ptrap[1:J],K[1:J]) 
  
  
  sigma ~ dinvgamma(2.27, 2.14)
  #sigma ~ dexp(1)
    g0 ~ dbeta(2, 2)
    psi ~ dbeta(2 ,2)
    N<- sum(w[1:M])
    D<- N/A
})
#-------------------------------------------------------------------------

constants <- list(M=M,J=J)

data <- list(n=n,X=X,xlim=xlim,ylim=ylim,A=A,K=K)

sx.init<- runif(M, xlim[1], xlim[2])
sy.init<- runif(M, ylim[1], ylim[2])
S.init<- cbind(sx.init,sy.init)

# Initial values
inits <- list(sigma=1,g0=0.1,psi=0.5,w=rbinom(M,1,0.5),S=S.init)

## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits, check = TRUE)
Rmcmc<- compileNimble(Rmodel)
#------------------------------------------------------------------------------------------
mcmcspec <- configureMCMC(Rmodel, monitors=c("N","D","g0","sigma","psi","S","w"), onlySlice = T)


mcmcspec$removeSamplers("w", print=FALSE)
mcmcspec$removeSamplers('S', print=FALSE)

wNodes <- Rmodel$expandNodeNames("w")
sNodePairs <- split( matrix(Rmodel$expandNodeNames('S'), ncol = 2), 1:M)
for(i in seq_along(wNodes)) mcmcspec$addSampler(target = wNodes[i], type = custom_zs_sampler, 
                            control = list(sNodes = sNodePairs[[i]], birthProb = 0.5, deathProb = 0.5 ))
for(i in seq_along(sNodePairs)) mcmcspec$addSampler(target = sNodePairs[[i]], 
                                type = custom_s_with_indicator_sampler, 
                      control = list(scaleFactor = 0.5, scaleNode = 'sigma', indicatorNode = wNodes[i]))

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

inits <- function(){list(sigma=1,g0=0.1,psi=0.5,w=rbinom(M,1,0.5),S=S.init)}

samp<- runMCMC(Cmodel, niter = ni, nburnin = nb, nchains = nc, thin=nt, inits = inits,
               samplesAsCodaMCMC = T)

foxes<- samp
saveRDS(foxes, "results/foxes.rds")

library(MCMCvis)
library(ggmcmc)
library(gridExtra)

MCMCsummary(foxes, params=c("D","N","g0","sigma"), digits=3, n.eff=TRUE,
            func=mymode, func_name = "Mode")



S<- ggs(samp)

win.graph(12,12)
p1<- ggs_density(S, family="N") + xlim(0, 300)
p2<- ggs_density(S, family="D") + xlim(0, 5)
p3<- ggs_density(S, family="g0") + xlim(0, 1)
p4<- ggs_density(S, family="sigma") + xlim(0, 2)

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

#-----------------------------------------------------  

load("foxHR")
load("cats")
load("dogs")

samp<- mcmc.list(sims[[1]],sims[[2]],sims[[3]])
parnames<- varnames(samp)
samp<- as.matrix(samp)

# expected foliar biomass
N<- samp[,grep("N",parnames)]
g0<- samp[,grep("g0",parnames)]
sigma<- samp[, grep("sigma",parnames)]
w<- samp[, grep("w",parnames)]
S<- samp[, grep("S",parnames)]

A<- (xlim[2]-xlim[1]) * (ylim[2]-ylim[1])
D<- N/A	

win.graph(10,10)
par(mfrow=c(2,2),pty='s')
hist(N,nclass=50,main="",xlim=c(0,300),xlab="Parameter value",ylab="Probability density",
     probability=T,las=0)
mtext(expression(paste("Abundance (",italic(hat(N)),")")))

plot(density(D),main="",xlim=c(0,2),xlab="Parameter value",ylab="Probability density",
     las=0,lwd=2,bty="n")
mtext(expression(paste("Density (",italic(hat(D)),")")))

plot(density(g0),main="",xlim=c(0,0.5),xlab="Parameter value",ylab="Probability density",
     las=0,lwd=2,bty="n")
mtext(expression(paste(g(0))))

plot(density(sigma),main="",xlim=c(0,5),xlab="Parameter value",ylab="Probability density",
     lwd=2,bty="n", ylim=c(0,4))
#plot(function(x) dunif(x, 0, 5), 0, 5, xlim=c(0, 3),add=T,col="grey80",lty=2)
plot(function(x) dgamma(x, 18, 30), 0, 5, xlim=c(0, 3),add=T,col="grey80",lty=2)
mtext(expression(paste(sigma)))


mymode(N)  # Posterior mode of N
mean(N)
median(N)
sd(N)
quantile(N,c(0.025,0.975))

# Density
mymode2(D) 
mean(D)
median(D)
sd(D)
quantile(D,c(0.025,0.975))

mean(g0)
median(g0)
sd(g0)
quantile(g0,c(0.025,0.975))

mean(sigma)
median(sigma)
sd(sigma)
quantile(sigma,c(0.025,0.975))


n<- rowSums(foxobs,na.rm=T)

win.graph(8,8)
obj<- list(Sx=S[,1:300],Sy=S[,301:600],z=w)
SCRdensity(obj,Xl=xlim[1],Xu=xlim[2],Yl=ylim[1],Yu=ylim[2],scalein=1,scaleout=1)
points(locs[,1],locs[,2])
points(locs[n>0,1],locs[n>0,2],pch=16)
mtext("Density",side=4,las=1,line=1,at=c(686,6047.5))
mtext("Easting",side=1,line=2)
mtext("Northing",side=2,line=2)
#------------------------
sims<- list(samp1=samp1,samp2=samp2,samp3=samp3)

save(sims, file="foxIP") #uninformative prior



