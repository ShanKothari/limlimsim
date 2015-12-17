#####################################################
## Iterative numerical solution of differential equations in one dimension
## These don't find how closely speacies can possibly pack; instead where
## success species invasions would result if they automatically evolved to
## trait optimum

n.breadth<-0.05
lower.b<-0
upper.b<-1

LV.eq<-function(trait.val,pt.df,n.breadth){
  niche.diffs<-exp(-(trait.val-pt.df$trait)^2/(n.breadth^2*2))
  ## this is the equilibrial LV population density
  pop<- -1*(trait.val-sum(niche.diffs*pt.df$pop))
  return(pop)
}

pt.df<-data.frame(trait=1,pop=1)
for(j in 1:15){
  opt<-optim(0.01,LV.eq,pt.df=pt.df,n.breadth=n.breadth,method="L-BFGS-B",lower=lower.b,upper=upper.b)
  pt.df[j+1,]<-c(opt$par,-1*opt$value)
  for(i in 1:50){
    for(i in 1:nrow(pt.df)){
      opt<-optim(pt.df$trait[i],LV.eq,pt.df=pt.df[-i,],n.breadth=n.breadth,method="L-BFGS-B",lower=lower.b,upper=upper.b)
      pt.df[i,]<-c(opt$par,-1*opt$value)
    }
  }
}

## using exponential decay fitness function

n.breadth<-0.1
lower.b<-0
upper.b<-Inf

LV.eq<-function(trait.val,pt.df,n.breadth){
  niche.diffs<-exp(-(trait.val-pt.df$trait)^2/(n.breadth^2*2))
  ## this is the equilibrial LV population density
  pop<- -1*(exp(-trait.val)-sum(niche.diffs*pt.df$pop))
  return(pop)
}

pt.df<-data.frame(trait=0,pop=1)
for(j in 1:5){
  opt<-optim(3,LV.eq,pt.df=pt.df,n.breadth=n.breadth,method="L-BFGS-B",lower=lower.b,upper=upper.b)
  pt.df[j+1,]<-c(opt$par,-1*opt$value)
  for(i in 1:50){
    for(i in 1:nrow(pt.df)){
      opt<-optim(pt.df$trait[i],LV.eq,pt.df=pt.df[-i,],n.breadth=n.breadth,method="L-BFGS-B",lower=lower.b,upper=upper.b)
      pt.df[i,]<-c(opt$par,-1*opt$value)
    }
  }
}

## alternately use this, together with min and which.min
## fit<-sapply(seq(0,1,0.001),function(x) LV.eq(x,pt.df,n.breadth=0.05))
## similarly, this can be used to make "fitness landscapes" for each sp given the traits of the rest
## subbing any sp for "1"
## fit<-sapply(seq(0,1,0.001),function(x) -1*LV.eq(x,pt.df[-1,],n.breadth))

n.points<-5
plot(pt.df[1:n.points,1],rep(0,n.points),ylim=c(0,1),xlim=c(0,1),main="Optimal trait value for each species in its own absence",xlab="trait",ylab="population size")
for(i in 1:n.points){points(seq(0,1,0.001),fit<-sapply(seq(0,1,0.001),function(x) -1*LV.eq(x,pt.df[-i,],n.breadth)),type="l")}
for(i in 1:n.points){abline(v=pt.df[i,1],col="red")}