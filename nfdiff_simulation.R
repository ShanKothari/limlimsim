## For Tilman/Clark Ecology Theory and Concepts class F2015
## Species Packing With Niche and Fitness Differences:
## The Limits of Limiting Similarity

#######################################################
## uses traits and niche overlap to parametrize a simulation
## of Lotka-Volterra competition each step (for a certain amount of time)

## The distribution is based on extinction and population-dependent speciation
## There is no directionality to evolution (not sure how to add that)

## initial number of species
init.nsp<-10

## number of traits
n.traits<-3

## probability of speciation (multiply by number of species remaining)
pspec<-0.05

## trait change picked from Gaussian dist w/ sd defined by tcsd
tcsd<-0.01

## niche breadth, for calculating overlap
n.breadth<-0.05

## generations to run
gens<-100

## fitness calculation function
fitness<-function(traits){
  ## geometric mean
  ## trait value of 1 yields optimal fitness
  fit<-apply(traits,1,function(x) exp(mean(log(x))))
  return(fit)
}

## niche overlap calculation function
niche.overlap<-function(traits,n.breadth){
  ## overlap of each species with each other species
  ## assumes all species are Gaussian with equal niche breadth
  ## niche breadth here refers to st. dev. of Gaussian
  dist.mat<-as.matrix(dist(traits))
  overlap<-exp(-dist.mat^2/(n.breadth^2*2))
  return(overlap)
}

## population dynamics simulator
pop.size<-function(fitness,overlap,init.pop,runs){
  ## imposes a max carrying capacity of 1, since a_ii=1 for all i
  ## at fitness=0, population naturally trends toward 0
  lambda<-fitness+1
  pop<-rep(init.pop,length(fitness))
  run=1
  while(run<=runs){
    F<-lambda/(1+(pop %*% overlap))
    pop<-F*pop
    run=run+1
  }
  return(pop)
}

## determines which species speciate
speciation<-function(pop,pspec,spec.fit){
  rnums<-runif(length(pop),0,1)
  if(spec.fit){
    spec<-which(rnums<pspec*pop)
  }
  else{
    spec<-which(rnums<pspec)
  }
  return(spec)
}

run.sim<-function(init.nsp,n.traits,gens,pspec,tcsd,n.breadth,fit.diff=T,niche.diff=T,spec.fit=T){
  traits<-data.frame(sapply(1:n.traits,function(x) runif(n=init.nsp,min=0,max=1)))
  gen=1
  while(gen<=gens){
    ptm0<-proc.time()
    
    ## population simulation
    if(fit.diff){
      fit<-fitness(traits)
    } else {
      fit<-rep(1,nrow(traits))
    }
    if(niche.diff){
      overlap.mat<-niche.overlap(traits,n.breadth)
    } else {
      ## does this make sense? need to think
      overlap.mat<-matrix(rep(1,nrow(traits)^2),ncol=nrow(traits))
    }
    pop<-pop.size(fit,overlap.mat,init.pop=0.1,runs=1000)
    
    ## extinction
    extinct<-which(pop<0.01)
    if(length(extinct)>0){
      traits<-data.frame(traits[-extinct,])
      pop<-pop[-extinct]
    }
      
    if(gen<gens){
      ## speciation
      spec<-speciation(pop,pspec,spec.fit)
      traits<-rbind(traits,traits[spec,])
    
      ## trait evolution
      change<-matrix(rnorm(nrow(traits)*ncol(traits),0,tcsd),ncol=ncol(traits))
      traits<-traits+change
    
      ## make sure no trait values exceed 1
      overshoot<-which(traits>1,arr.ind=T)
      traits[overshoot]<-1
      undershoot<-which(traits<0,arr.ind=T)
      traits[undershoot]<-0.001
    }
    
    ## tracking time per iterate
    ptm1 <- proc.time() - ptm0
    jnk<-as.numeric(ptm1[3])
    cat('\n','It took ', jnk, "seconds to do iteration", gen, ", ",nrow(traits),"species remaining")
    gen<-gen+1
  }
  return(traits)
}

## nearest-neighbor distance calculator
nn.calc<-function(trait.mat){
  dist.mat<-as.matrix(dist(trait.out))
  diag(dist.mat)<-NA
  nn.dist<-apply(dist.mat,1,function(x) min(x,na.rm=T))
  return(nn.dist)
}