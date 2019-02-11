#----------------------------------------------------
#
#  Nimble functions
#
#----------------------------------------------------

# Binomial distribution for rom vectors
dbin_by_row <- nimbleFunction(
  run = function(x = double(1), P = double(1), K = double(1), log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dbinom(x[j], K[j], P[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rbin_by_row  <- nimbleFunction(
  run = function(n = integer(), P = double(1), K = double(1)) {
    declare(ans, double(1))
    J <- length(P)
    setSize(ans, J)
    for(j in 1:J)
      ans[j] <- rbinom(1, K[j], P[j])
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
  run = function(x = double(1), bigLambda = double(1), K = double(1), log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dpois(x[j], bigLambda[j]*K[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rpois_by_row  <- nimbleFunction(
  run = function(n = integer(), bigLambda = double(1), K = double(1)) {
    declare(ans, double(1))
    J<- length(bigLambda)
    setSize(ans, J)
    for(j in 1:J)
      ans[j] <- rpois(1, bigLambda[j]*K[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dpois_by_row = list(
    BUGSdist = "dpois_by_row(bigLambda, K)",
    Rdist = "dpois_by_row(bigLambda, K)",
    range = c(0, Inf),
    types = c('value = double(1)', 'bigLambda = double(1)', 'K = double(1)'))
))
#-------------------------------
TProb<- nimbleFunction(
  run=function(x = double(2)) {
    nr<- dim(x)[1]
    nc<- dim(x)[2]
    ans<- numeric(nc)
    for(i in 1:nc)
      ans[i]<- 1-prod(1-x[1:nr,i])
    returnType(double(1))
    return(ans)
  })
#------------------------
Tsum<- nimbleFunction(
  run=function(x = double(2), lam=double()) {
    nr<- dim(x)[1]
    nc<- dim(x)[2]
    ans<- numeric(nc)
    for(i in 1:nc)
      ans[i]<- lam * sum(x[1:nr,i])
    returnType(double(1))
    return(ans)
  })
#--------------------------------------------------------------------------------
# The following custom MCMC samplers are based on ideas for tuning the Nimble MCMC 
# algorithm for SCR models given at
# https://nature.berkeley.edu/~pdevalpine/SCR_NIMBLE_ideas/SCR_NIMBLE_ideas.html
#--------------------------------------------------------------------------------

custom_s_with_indicator_sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    scaleFactor <- control$scaleFactor
    scaleNode <- control$scaleNode
    calcNodes <- model$getDependencies(target)
    target1 <- target[1]
    target2 <- target[2]
    indicatorNode <- control$indicatorNode
  },
  run = function() {
    if(model[[indicatorNode]]==0) return()
    propSD <- model[[scaleNode]]*scaleFactor
    model[[target1]] <<- rnorm(1, mean = model[[target1]], sd = propSD)
    model[[target2]] <<- rnorm(1, mean = model[[target2]], sd = propSD)
    logMHR <- calculateDiff(model, calcNodes)
    jump <- decide(logMHR)
    if(jump)
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    else
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(
    reset = function () {}
  )
)

#--------------------------------------

custom_zs_sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    zNode <- target
    sNodes <- control$sNodes
    birthProb <- control$birthProb
    deathProb <- control$deathProb
    calcNodes <- model$getDependencies(c(zNode, sNodes))
  },
  run = function() {
    currentIndicatorValue <- model[[zNode]]
    currentLogProb <- getLogProb(model, calcNodes)
    if(currentIndicatorValue == 0) {
      if(runif(1,0,1) > birthProb) return()
      ## propose birth
      logProbReverseProposalValues <- getLogProb(model, sNodes)
      simulate(model, sNodes)
      logProbProposalValues <- calculate(model, sNodes)
      model[[zNode]] <<- 1
      proposalLogProb <- calculate(model, calcNodes)
      logMHR <- proposalLogProb - currentLogProb + log(deathProb) + logProbReverseProposalValues - log(birthProb) - logProbProposalValues
    } else {
      if(runif(1,0,1) > deathProb) return()
      ## propose death
      logProbReverseProposalValues <- getLogProb(model, sNodes)
      simulate(model, sNodes)
      logProbProposalValues <- calculate(model, sNodes)
      model[[zNode]] <<- 0
      proposalLogProb <- calculate(model, calcNodes)
      logMHR <- proposalLogProb - currentLogProb + log(birthProb) + logProbReverseProposalValues - log(deathProb) - logProbProposalValues
    }
    
    jump <- decide(logMHR)
    if(jump)
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    else
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(
    reset = function () {}
  )
)

#----------------------------

mymode<- function (x) 
{ # for continuous x only
  d <- density(x)
  mx<- which.max(d$y)
  d$x[mx]
}
