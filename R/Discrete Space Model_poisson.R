
library(readxl)
library(tidyverse)
library(sf)
library(nimble)
library(MCMCvis)

source("r/Functions.r") # required for nimble model

#-----------------------------------
map<- read_sf("Data/yanakie.geojson")

hist<- read_excel("Data/Wilsons_Prom_edited.xlsx",sheet="WilsonsPromData",na=" ")
IDS<- read_excel("Data/Wilsons_Prom_edited.xlsx",sheet="ListOfCatID")
cams<- read_excel("Data/Wilsons_Prom_edited.xlsx",sheet="CameraLocations")


#------------------------------------
hist<- hist %>% arrange(Cam)
cams<- cams %>% arrange(Cam)

history<- dplyr::select(hist,-Cam)

IDk<- IDS %>% filter(Status=="known") %>% dplyr::select(ID)
IDu<- IDS %>% filter(Status=="unknown") %>% dplyr::select(ID)
IDk<- IDk$ID
IDu<- IDu$ID


h1<- build.capt(history, IDk, 1)
Yu<- build.sight(history, IDu, 1)

Yk<- capt.array(h1, dim(history)[2], dim(history)[1], binary=FALSE)
#---------------------------------------------
# 
# Assemble data

M<-50  # data augmentation - unmarked
mmax<- 50 # data augmentation - marked
eps<- 500 # pixel size of state space

ym<- apply(Yk, c(1,2), sum, na.rm=T) 

n<- rowSums(Yu, na.rm=T)
get.k <- function(x) length(x[!is.na(x)])
K<- apply(Yu, 1, get.k)
m<- nrow(ym) + mmax
J<- nrow(Yu)

yaug<- matrix(0,mmax,J)
ym<- rbind(ym,yaug)


camlocs<- cams %>% dplyr::select(x=Easting,y=Northing)

cam.region<- make_grid(map, cell_diameter=eps, overlap="centre",xy=camlocs)
hex.cent<- st_centroid(cam.region)
X<- as.matrix(camlocs)/1000
grid<- st_coordinates(hex.cent)/1000
npix<- nrow(grid)
pix_area<- as.numeric(st_area(cam.region)[1]/1e6)
Area<- as.numeric(st_area(map)/1e6)

constants <- list(M=M,m=m,J=J,npix=npix,pix_area=pix_area, A=Area)

data <- list(ym=ym, n=n, grid=grid, X=X, K=K)


#================================================

code <- nimbleCode({
  
  # Priors
  log(lam0) <- alpha
  alpha ~ dnorm(0, sd=2)
  log(sigma) <- alpha1
  alpha1 ~ dnorm(0, sd=2)
  beta0 ~ dnorm(0, sd=2)
  
  for(j in 1:npix) {
    log(mum[j])<- beta0 + log(pix_area)
    log(muu[j])<- beta0 + log(pix_area)
    ppm[j]<- mum[j]/ENm
    ppu[j]<- muu[j]/ENu
  }
  ENm<- sum(mum[1:npix])
  ENu<- sum(muu[1:npix])
  psim <- ENm/m
  psiu <- ENu/M
  
  # Marked part
  for(i in 1:m) {
    z[i] ~ dbern(psim)
    s[i] ~ dcat(ppm[1:npix])
    g[i,1]<- grid[s[i],1]
    g[i,2]<- grid[s[i],2]
    prob[i,1:J]<- CalcDetection(g[i,1:2], lam0, sigma, X[1:J,1:2], J, z[i])
    ym[i,1:J] ~ dpois_by_row(prob[i,1:J],K[1:J])
  }
  
  # Unmarked part
  for(i in (m+1):(m+M)) {
    z[i] ~ dbern(psiu)
    s[i] ~ dcat(ppu[1:npix])
    g[i,1]<- grid[s[i], 1]
    g[i,2]<- grid[s[i], 2]
    prob[i,1:J]<- CalcDetection(g[i,1:2], lam0, sigma, X[1:J,1:2], J, z[i])
  }
  
  for(j in 1:J) {
    LamTot[j]<- sum(prob[(m+1):(m+M),j])
  }
  n[1:J] ~ dpois_by_row(LamTot[1:J],K[1:J])

  Nm<- sum(z[1:m])
  Nu<- sum(z[(m+1):(m+M)])
  N<- Nu + Nm
  D<- N/A	
})

#------------------------------------------------
# Initial values
sst <- sample(1:npix, (m+M), replace = TRUE)
zmst <- c(rep(1,dim(Yk)[1]),rep(0,mmax))
zust<- rbinom(M,1,0.9)
zst<- c(zmst,zust)

inits1 <- list(alpha= -4, alpha1=0, s=sst, z=zst, beta0 = -1)
             

## create the model object
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits1, check = FALSE)
Rmcmc<- compileNimble(Rmodel, showCompilerOutput = F)

conf <- configureMCMC(Rmodel)

## add a binary state sampler for each wm and w node
conf$removeSamplers(c("z"), print = FALSE)
Nodes <- Rmodel$expandNodeNames("z")
for(Node in Nodes) conf$addSampler(target = Node, type = "binary", print=FALSE)

conf$removeSamplers(c('alpha','alpha1','beta0'), print=FALSE)
conf$addSampler(target=c('alpha','alpha1','beta0'), type='AF_slice')

conf$resetMonitors()
conf$addMonitors(c("lam0","sigma","beta0","D","Nm","Nu","N"))

Cmodel <- buildMCMC(conf)

Cmcmc <- compileNimble(Cmodel, project = Rmodel, resetFunctions = T)

ni<- 2000  # for testing - run more iterations to improve convergence
nb<- 1000  # same
nc<- 3
nt<- 1
#--------------------


inits = function(){inits1}

samp<- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, thin = nt, inits = inits,  
               samplesAsCodaMCMC = TRUE)

MCMCsummary(samp, c("lam0","sigma","beta0","Nm","Nu","N","D"))


