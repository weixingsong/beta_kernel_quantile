# R-Codes for 
# Computing Bootstrap Bandwidths for Beta and Normal Kernel
# Quntile Estimators

# Input:
#   x: Data
#   pv: A number from (0,1) indicating which quantile to be found
#   B: Bootstrap sample size
#   hleft,hright: Interval where the bandwidth is sought

# Output:
#   Lb1: Bootstrap mse for Type-1 Beta Kernel 
#   Lb2: Bootstrap mse for Type-2 Beta Kernel 
#   Ln:  Bootstrap mse for Normal Kernel 
#   hb1: Bootstrap bandwidht for Type-1 Beta Kernel 
#   hb2: Bootstrap bandwidth for Type-2 Beta Kernel 
#   hn:  Bootstrap bandwidth for Type-1 Beta Kernel 


hboot=function(x,pv,B=200,hleft=0.0001,hright=0.5)
 {
  #set.seed(1234)
  Btemp=replicate(B,sort(sample(x,replace=TRUE)),simplify=FALSE)
  Bsample=sapply(Btemp,identity)
  nsize=length(x)
  iseq = seq(nsize)/nsize
  jseq = (seq(nsize)-1)/nsize         

  hseq=seq(hleft,hright,length=1000)
  Lb1=Lb2=Ln=Le=rep(0,length(hseq))
  k=1
  for(hv in hseq)
  {    
    B1v=pbeta(iseq,pv/hv+1,(1-pv)/hv+1)-
        pbeta(jseq,pv/hv+1,(1-pv)/hv+1)
    
    B2v=pbeta(iseq,pv/hv,(1-pv)/hv)-
        pbeta(jseq,pv/hv,(1-pv)/hv)
    
    Nv=pnorm(iseq,pv,hv) - pnorm(jseq,pv,hv)
    
    Ev=rep(0,nsize)
    for(kk in seq(nsize))
    {
      a=pv-hv
      b=pv+hv
      c=(kk-1)/nsize
      d=kk/nsize
      
      lend=max(a,c)
      rend=min(b,d)
      
      
      if(lend>=rend)
      {Ev[kk]=0}
      if(lend<rend)
      {
        bb=(rend-pv)/hv
        aa=(lend-pv)/hv
        Ev[kk]=((bb-aa)-(bb^3-aa^3)/3)*0.75
      }
    }

    Qb1t=t(Bsample)%*%B1v
    Qb2t=t(Bsample)%*%B2v
    QEpt=t(Bsample)%*%Ev
    Qnt=t(Bsample)%*%Nv
    Qet=quantile(x,pv)
    
    Lb1[k]=mean((Qb1t-Qet)^2)
    Lb2[k]=mean((Qb2t-Qet)^2)
    Le[k]=mean((QEpt-Qet)^2)
    Ln[k]=mean((Qnt-Qet)^2)
    k=k+1
  }
  
  hb1=hseq[Lb1==min(Lb1)]
  hb2=hseq[Lb2==min(Lb2)]
  he=hseq[Le==min(Le)]
  hn=hseq[Ln==min(Ln)]
  list(Lb1=Lb1,Lb2=Lb2,Ln=Ln,hb1=hb1,hb2=hb2,he=he,hn=hn)
}

