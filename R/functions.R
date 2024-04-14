#----------------------------------------------------
#
#  functions
#
#----------------------------------------------------

# Binomial distribution for for vectors
dbin_by_row <- nimbleFunction(
  run = function(x = double(1), P = double(1), K = double(1), log = integer(0, default = 0)) {
    ans <- sum(dbinom(x, size=K, prob=P, log=TRUE))
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rbin_by_row  <- nimbleFunction(
  run = function(n = integer(0), P = double(1), K = double(1)) {
    ans <- rbinom(length(K), size=K, prob=P)
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dbin_by_row = list(
    BUGSdist = "dbin_by_row(P, K)",
    Rdist = "dbin_by_row(P, K)",
    range = c(0, Inf),
    types = c('value = double(1)', 'P = double(1)', 'K = double(1)'))
))
#-------------------------------
dpois_by_row <- nimbleFunction(
  run = function(x = double(1), lam = double(1), K = double(1), log = integer(0, default = 0)) {
    ans <- sum(dpois(x, lam*K, log=TRUE))
    returnType(double(0))
    if(log) return(ans)
    else return(exp(ans))
  })

rpois_by_row  <- nimbleFunction(
  run = function(n = integer(), lam = double(1), K = double(1)) {
    ans <- rpois(length(lam), lam*K)
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dpois_by_row = list(
    BUGSdist = "dpois_by_row(lam, K)",
    Rdist = "dpois_by_row(lam, K)",
    range = c(0, Inf),
    types = c('value = double(1)', 'lam = double(1)', 'K = double(1)'))
))

CalcDetection <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

#------------------------------
# Utility functions

#------------------------------------------------------------
# Build a capthist file from a matrix of detections

build.capt<- function(history, ID, sessid) {
  # cells of history with multiple ID must be separated by ","
  n<- length(ID)
  J<- dim(history)[1]
  K<- dim(history)[2]
  na.ind<- which(is.na(history), arr.ind=TRUE)
  chist<- list()
  ind<- 1
  for(i in 1:n) {
    id<- ID[i]
    for(j in 1:J) {
      for(k in 1:K) {
        ss<- strsplit(as.character(history[j,k]),",")[[1]]
        if(any(grepl(paste0("^",id,"$"),ss))) {
          chist[[ind]]<- data.frame(sessid,id,k,j)
          ind<- ind + 1
        }
      }
    }
  }
  chist.df<- do.call('rbind', chist)
  names(chist.df)<- c("Session","ID","Occasion","Detector")
  list(chist=chist.df,na.ind=na.ind)
}

#------------------------------------------------------------
# Build a unmarked sighting file from a matrix of detections

build.sight<- function(history, ID, sessid) {
  # cells of history with multiple ID must be separated by ","
  n<- length(ID)
  J<- dim(history)[1]
  K<- dim(history)[2]
  Yu<- matrix(0, J, K)
  na.mat<- is.na(history)
  Yu[na.mat]<- NA
  for(i in 1:n){
    id<- ID[i]
    for(j in 1:J) {
      for(k in 1:K) {
        ss<- strsplit(as.character(history[j,k]),",")[[1]]
        if(any(grepl(paste0("^",id,"$"),ss))) Yu[j,k]<- Yu[j,k] + 1
      }
    }
  }
  Yu
}

#----------------------------------------------
# Spread TrapID dataframe into detection array

capt.array<- function(ch, nocc, ntraps, binary=FALSE) {
  # ch is output from build.capt()
  trapid<- ch$chist
  id<- sort(unique(trapid$ID))
  n<- length(id)
  Yk<- array(0, c(n, ntraps, nocc))
  for(j in 1:nrow(ch$na.ind))
    Yk[, ch$na.ind[j,1], ch$na.ind[j,2]]<- NA # Fill in missing cells
  for(i in 1:n) {
    tmpid<- trapid[trapid$ID == id[i],]
    K<- nrow(tmpid)
    for(k in 1:K) {
      if(binary)
        Yk[i, tmpid$Detector[k], tmpid$Occasion[k]]<- 1
      else
        Yk[i, tmpid$Detector[k], tmpid$Occasion[k]]<- Yk[i, tmpid$Detector[k], tmpid$Occasion[k]] + 1
    }}
  Yk
}
#--------------------------------------------------------------------------------
calc.dist<- function(A, D) {
  #calculates distance matrix between devices and HR centres
  dx<- outer(A[,1], D[,1], FUN="-")
  dy<- outer(A[,2], D[,2], FUN="-")
  return(dx^2 + dy^2)
}

#---------------------------------------------------
# for discrete state space. input is an sf object
make_grid <- function(shp, cell_diameter, cell_area, square= FALSE, offset=c(0,0), 
                      xy=NULL, overlap = c("centre","any","all"), plot=TRUE) {
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  overlap<- match.arg(overlap)
  # generate array of hexagon centers
  g <- st_make_grid(shp, cellsize = cell_diameter, what="polygons", square=square, offset = offset)
  
  # clip to boundary of study area
  if(overlap == "centre")
    inside <- apply(st_within(st_centroid(g), shp), 1, any)
  else if(overlap == "any")
    inside <- apply(st_intersects(g, shp), 1, any)
  else if(overlap == "all")
    inside <- apply(st_within(g, shp), 1, any)
  g<- g[inside]
  g<- st_sf(HexID = 1:length(g), geometry=g)
  if(!is.null(xy)) {
    xy<- st_as_sf(data.frame(xy), coords=1:2, crs = st_crs(shp))
    iin<- apply(st_intersects(xy, g), 1, any)
    if(length(which(!iin)) > 0)
      cat("these traps are outside the mask area:",which(!iin))
  }
  if(plot) {
    p<- shp |> ggplot() +
      geom_sf(fill="grey70", color=NA, data=g) +
      geom_sf(fill=NA)   
    if(!is.null(xy) & (sum(iin) == nrow(xy)))
      p<- p + geom_sf(data=xy[iin,], shape=4, color="green")
    if(!is.null(xy) & (sum(iin) != nrow(xy))) {
      p<- p + geom_sf(data=xy[iin,], shape=4, color="green")
      p<- p +geom_sf(data=xy[!iin,], shape=4, color="red")
    }
    print(p)
  }
  return(g)
}


