
#-----------------------------------------------------------
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

capt.array<- function(ch, nocc, ntraps) {
  # ch is output from build.capt()
  trapid<- ch$chist
  id<- sort(unique(trapid$ID))
  n<- length(id)
  Yk<- array(0, c(n, ntraps, nocc))
  Yk[, ch$na.ind[,1], ch$na.ind[,2]]<- NA # Fill in missing cells
  for(i in 1:n) {
    tmpid<- trapid[trapid$ID == id[i],]
    K<- nrow(tmpid)
    for(k in 1:K) {
      Yk[i, tmpid$Detector[k], tmpid$Occasion[k]]<- 1
    }}
  Yk
}
