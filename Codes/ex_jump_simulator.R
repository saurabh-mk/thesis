# jump simulation
ex.jumpsimulator=function(phy, alpha=0, sigmasq.brown=0.01, sigmasq.jump=1, jumps){
  
  .jumpsim=function(phy, alpha=0, sigmasq.brown=0.01, sigma.jump=0.2, lambda.jump=0.2){
    phy=reorder(phy)
    cedges=cumsum(phy$edge.length)
    tmax=sum(phy$edge.length)
    nn=phy$edge[,2]
    cs=c()
    tt=c()
    jumps=0
    if(lambda.jump>0){
      while(1){
        dt=rexp(1,lambda.jump)
        tt=c(tt,dt)
        cs=cumsum(tt)
        if(cs[length(cs)]>tmax){
          cs=cs[-length(cs)]
          jumps=length(cs)
          break()
        }
      }
    }
    
    jump.edges=rep(0, length(nn))
    if(jumps>0) {
      for(j in 1:length(cs)){
        tmp=min(which(cs[j]<cedges))
        jump.edges[tmp]=jump.edges[tmp]+1
      }
    }
    
    hist=as.data.frame(matrix(cbind(phy$edge, phy$edge.length, NA, jump.edges, NA, NA, NA), ncol=8))
    names(hist)=c("ancestor","descendant","edge","phenotype","jumps","time", "effect_Brown", "effect_Jump")
    
    root=Ntip(phy)+1
    for(i in 1:nrow(hist)){
      start=ifelse(hist$ancestor[i]==root, alpha, hist$phenotype[which(hist$descendant==hist$ancestor[i])])
      stime=ifelse(hist$ancestor[i]==root, 0, hist$time[which(hist$descendant==hist$ancestor[i])])
      t=hist$edge[i]
      hist$phenotype[i]=start+(Beffect<-rnorm(1, mean=0, sd=sqrt(sigmasq.brown*t)))
      hist$effect_Brown[i]=abs(Beffect^2)/hist$edge[i]
      hist$time[i]=stime+hist$edge[i]
      if((jmp<-hist$jumps[i])>0){
        hist$phenotype[i]=hist$phenotype[i]+(Jeffect<-sum(rnorm(jmp, mean=0, sd=sigma.jump)))
      } else {
        Jeffect=0
      }
      hist$effect_Jump[i]=(Jeffect/abs(Jeffect))*abs(Jeffect^2)/hist$edge[i]
    }
    
    return(hist)
  }
  
  lambda=jumps/sum(phy$edge.length)
  hist=.jumpsim(phy, alpha=alpha, sigmasq.brown=sigmasq.brown, sigma.jump=sqrt(sigmasq.jump), lambda.jump=lambda)
  dat=hist$phenotype[hist$descendant<=Ntip(phy)]
  names(dat)=phy$tip.label[hist$descendant[hist$descendant<=Ntip(phy)]]
  edges=hist$descendant[hist$jumps>0]
  hist$scl=hist$effect_Jump
  
  hist$cex=(hist$scl-min(hist$scl))/(max(hist$scl)-min(hist$scl))
  hist$cex=4*asin(sqrt(hist$cex))
  
  return(list(hist=hist, dat=dat, phy=phy))
}