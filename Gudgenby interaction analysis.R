library(readxl)
library(spatstat)
library(MCMCvis)
library(tidyverse)
# fit nonstationary marked Poisson process

cams<- read_excel("data/gudgenby_camsUTM_revised.xlsx")

cams<- cams[,c(2,3)]/1000 #km
xlim<- range(cams[,1])
ylim<- range(cams[,2])
buffer<- 2
xlim[1]<- xlim[1] - buffer
xlim[2]<- xlim[2] + buffer
ylim[1]<- ylim[1] - buffer
ylim[2]<- ylim[2] + buffer

site.region<- owin(xrange=xlim,yrange=ylim)

camlocs<- as.ppp(cams,W=site.region)

#---------------------
# Fox data
#

nsim<- 500
r=seq(0,2.5,0.1)

FD.res<- matrix(NA, nrow=length(r),ncol=nsim)
FC.res<- matrix(NA, nrow=length(r),ncol=nsim)
DC.res<- matrix(NA, nrow=length(r),ncol=nsim)

foxes<- readRDS("foxes.rds") # Use SPA_nimble.r to generate posteriors
dogs<- readRDS("dogs.rds") # use SMR_nimble.r to generate posteriors
cats<- readRDS("cats.rds") # use SMR_nimble.r to generate posteriors

for(i in 1:nsim) {
cat("Doing simulation ",i,"\n")
  
out<- MCMCpstr(foxes, c("S","w"), type="chains")

nsamp<- sample(1:30000, size=100, replace=F)

Sx<- out$S[,1,nsamp]
Sy<- out$S[,2,nsamp]
z<- out$w[,nsamp]

fox.x<- Sx[z==1]
fox.y<- Sy[z==1]

ok<- inside.owin(fox.x, fox.y, site.region)
fox.x<- fox.x[ok]
fox.y<- fox.y[ok]
#--------------------------
# Dog data

out<- MCMCpstr(dogs, c("Sm","wm","S","w"), type="chains")

Sx<- rbind(out$Sm[,1,nsamp],out$S[,1,nsamp])
Sy<- rbind(out$Sm[,2,nsamp],out$S[,2,nsamp])
z<- rbind(out$wm[,nsamp],out$w[,nsamp])

dog.x<- Sx[z==1]
dog.y<- Sy[z==1]

ok<- inside.owin(dog.x, dog.y, site.region)
dog.x<- dog.x[ok]
dog.y<- dog.y[ok]
#----------------------------
# cat data

out<- MCMCpstr(cats, c("Sm","wm","S","w"), type="chains")

Sx<- rbind(out$Sm[,1,nsamp],out$S[,1,nsamp])
Sy<- rbind(out$Sm[,2,nsamp],out$S[,2,nsamp])
z<- rbind(out$wm[,nsamp],out$w[,nsamp])

cat.x<- Sx[z==1]
cat.y<- Sy[z==1]

ok<- inside.owin(cat.x, cat.y, site.region)
cat.x<- cat.x[ok]
cat.y<- cat.y[ok]

X<- c(fox.x,dog.x,cat.x)
Y<- c(fox.y,dog.y,cat.y)
M<- factor(c(rep("F",length(fox.x)),rep("D",length(dog.x)),rep("C",length(cat.x))))

hr.ppp<- ppp(X, Y, marks=M,window=site.region)
hr.ppp<- as.ppp(hr.ppp)
hr.ppp<- hr.ppp[!duplicated(hr.ppp),]

FD<- markconnect(hr.ppp, "F","D", normalise=TRUE, r=r)
FC<- markconnect(hr.ppp, "F","C", normalise=TRUE, r=r)
DC<- markconnect(hr.ppp, "D","C", normalise=TRUE, r=r)

FD.res[,i]<- FD$iso
FC.res[,i]<- FC$iso
DC.res[,i]<- DC$iso

}

FD<- apply(FD.res,1,mean)
FD.cl<- apply(FD.res,1,quantile, c(0.025,0.975))
FC<- apply(FC.res,1,mean)
FC.cl<- apply(FC.res,1,quantile, c(0.025,0.975))
DC<- apply(DC.res,1,mean)
DC.cl<- apply(DC.res,1,quantile, c(0.025,0.975))

fd.results<- data.frame(r=r,pr=FD,lcl=FD.cl[1,],ucl=FD.cl[2,],Pair="Fox-Dog")
fc.results<- data.frame(r=r,pr=FC,lcl=FC.cl[1,],ucl=FC.cl[2,],Pair="Fox-Cat")
dc.results<- data.frame(r=r,pr=DC,lcl=DC.cl[1,],ucl=DC.cl[2,],Pair="Dog-Cat")
mc.results<- bind_rows(fd.results,fc.results,dc.results)


win.graph(12,10)
mc.results %>% ggplot(aes(r, pr, colour=Pair)) +
  geom_line() +
  geom_ribbon(aes(ymin=lcl,ymax=ucl, colour=Pair),alpha=0.25, fill="grey40",linetype="blank") +
  ylab(expression(paste(italic(p[i][j](r))))) +
  xlab(expression(paste(italic(r), (m)))) +
  geom_hline(aes(yintercept = 1),linetype=2, size=1) +
  theme_bw()


#------------------------------------------------------------------------------
#
# Calculate spatial probability of use by each species
#
#----------------------------------------------------------------------------
nsamp<- sample(1:30000, size=1000, replace=F) # 1000 posterior samples

out<- MCMCpstr(foxes, c("S","w"), type="chains")

nsamp<- sample(1:30000, size=100, replace=F)

Sx<- out$S[,1,nsamp]
Sy<- out$S[,2,nsamp]
z<- out$w[,nsamp]

fox.x<- Sx[z==1]
fox.y<- Sy[z==1]

ok<- inside.owin(fox.x, fox.y, site.region)
fox.x<- fox.x[ok]
fox.y<- fox.y[ok]
#--------------------------
# Dog data

out<- MCMCpstr(dogs, c("Sm","wm","S","w"), type="chains")

Sx<- rbind(out$Sm[,1,nsamp],out$S[,1,nsamp])
Sy<- rbind(out$Sm[,2,nsamp],out$S[,2,nsamp])
z<- rbind(out$wm[,nsamp],out$w[,nsamp])

dog.x<- Sx[z==1]
dog.y<- Sy[z==1]

ok<- inside.owin(dog.x, dog.y, site.region)
dog.x<- dog.x[ok]
dog.y<- dog.y[ok]
#----------------------------
# cat data

out<- MCMCpstr(cats, c("Sm","wm","S","w"), type="chains")

Sx<- rbind(out$Sm[,1,nsamp],out$S[,1,nsamp])
Sy<- rbind(out$Sm[,2,nsamp],out$S[,2,nsamp])
z<- rbind(out$wm[,nsamp],out$w[,nsamp])

cat.x<- Sx[z==1]
cat.y<- Sy[z==1]

ok<- inside.owin(cat.x, cat.y, site.region)
cat.x<- cat.x[ok]
cat.y<- cat.y[ok]

X<- c(fox.x,dog.x,cat.x)
Y<- c(fox.y,dog.y,cat.y)
M<- factor(c(rep("F",length(fox.x)),rep("D",length(dog.x)),rep("C",length(cat.x))))

hr.ppp<- ppp(X, Y, marks=M,window=site.region)
hr.ppp<- as.ppp(hr.ppp)
hr.ppp<- hr.ppp[!duplicated(hr.ppp),]


pred<- relrisk(hr.ppp)
#------------------------------
# Now some plots
#------------------------------

library(tidyverse)

pp<- lapply(pred, as.data.frame)
pp[[1]]$Species<- "Feral cat"
pp[[2]]$Species<- "Dingo"
pp[[3]]$Species<- "Fox"
pp<- do.call('rbind', pp)
names(pp)[3]<- "Prob"
pp$Species<- factor(pp$Species,levels=c("Dingo","Fox","Feral cat"))

win.graph(10,4)
pp %>% ggplot(aes(x, y)) +
  geom_raster(aes(fill=Prob)) +
  geom_point(aes(x, y), data=locs, shape=1, size=1.5) +
  scale_fill_distiller(palette = "Spectral",name="Probability") +
  scale_x_continuous(limits=c(xlim[1],xlim[2]),breaks=c(675,677.5,680)) +
  scale_y_continuous(limits=c(ylim[1],ylim[2])) +
  facet_wrap(~Species, nrow=1) +
  labs(x="Easting (km)",y="Northing (km)") +
  theme_bw()

pp<- as.data.frame(hr.ppp)
names(pp)[3]<- "Species"
pp$Species<- recode(pp$Species, "D"="Dingo","F"="Fox","C"="Feral cat")
pp$Species<- factor(pp$Species,levels=c("Dingo","Fox","Feral cat"))


win.graph(10,4)
pp %>% ggplot(aes(x, y)) +
  stat_bin_2d(aes(fill=..count../1000/0.5^2),geom = "tile", binwidth=c(0.5, 0.5),drop=F) +
  geom_point(aes(x, y), data=locs, shape=1, size=1.5) +
  facet_wrap(~Species,nrow=1) +
  scale_fill_distiller(palette = "Spectral", name="Density") +
  scale_x_continuous(limits=c(xlim[1],xlim[2]),breaks=c(675,677.5,680)) +
  scale_y_continuous(limits=c(ylim[1],ylim[2])) +
  labs(x="Easting (km)",y="Northing (km)") +
  theme_bw()

win.graph(10,4)
pp %>% ggplot(aes(x, y)) +
  stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE) +
  geom_point(aes(x, y), data=locs, shape=1, size=1.5) +
  facet_wrap(~Species) +
  scale_fill_distiller(palette = "Spectral") +
  labs(x="Easting",y="Northing") +
  theme_bw()
