# Simulation Study for the paper

# Enhancing Quantile Function Estimation with Beta-Kernel Smoothing

# Danish Fire Data with transformation

# Table 3

library(QRM)
data(danish)
x=danish
x1=x[x>1]
x2=x[x>10]
result=matrix(0,nrow=6,ncol=5)
k=1
for (dat in c(1,2))
{
  for (pv in c(0.99,0.995,0.999))
  {
    if (dat==1)
    {x=x1}
    else
    {x=x2}  
    n=length(x)
    
    hb=0.01*n^(-2/3)
    hn=0.01*n^(-1/3)
    hv=0.01*n^(-2/3)
    he=0.01*n^(-1/3)
    
     booth=hboot(x,pv,B=500,hleft=0.00001,hright=0.5)
     hb=booth$hb1
     hv=booth$hb2
     he=booth$he
     hn=booth$hn
     cat(hb,hv,he,hn,"\n")
     
    xord = sort(x)
    iseq = seq(n)/n
    jseq = (seq(n)-1)/n  
    
    
    # Champernowne Transformation
    
    M = median(x)
    
    Lfun=function(a)
    {
      n*(log(a[1])+log((M+a[2])^a[1]-a[2]^a[1]))+
        (a[1]-1)*sum(log(x+a[2]))-
        2*sum(log((x+a[2])^a[1]+(M+a[2])^a[1]-2*a[2]^a[1]))
    }
    
    av = optim(par = c(0.03, 0.01),fn=Lfun)$par

    a1 = av[1]
    a2 = av[2]
    
    y=((x+a2)^a1-a2^a1)/((x+a2)^a1+(M+a2)^a1-2*a2^a1)
    

    yord = sort(y)
    

    Bv=pbeta(iseq,pv/hb+1,(1-pv)/hb+1)-
        pbeta(jseq,pv/hb+1,(1-pv)/hb+1)
    Nv=pnorm(iseq,pv,hn) - pnorm(jseq,pv,hn)
    Hv = pbeta(iseq,pv/hv,(1-pv)/hv)-
      pbeta(jseq,pv/hv,(1-pv)/hv)
    
    Ev=rep(0,n)
    for(kk in seq(n))
    {
      a=pv-he
      b=pv+he
      c=(kk-1)/n
      d=kk/n
      
      lend=max(a,c)
      rend=min(b,d)
      
      
      if(lend>=rend)
       {Ev[kk]=0}
      if(lend<rend)
       {
        bb=(rend-pv)/he
        aa=(lend-pv)/he
        Ev[kk]=((bb-aa)-(bb^3-aa^3)/3)*0.75
       }
    }
    
    Qbv = sum(yord*Bv)
    bfun=function(u)
      {
        ((u+a2)^a1-a2^a1)/((u+a2)^a1+(M+a2)^a1-2*a2^a1)-Qbv
      }
    Qbt=uniroot(bfun,interval=c(10,400))$root
      
    Qnv = sum(yord*Nv)
    nfun=function(u)
      {
        ((u+a2)^a1-a2^a1)/((u+a2)^a1+(M+a2)^a1-2*a2^a1)-Qnv
      }
    Qnt=uniroot(nfun,interval=c(0,400))$root
    
    Qhv = sum(yord*Hv)
    hfun=function(u)
    {
      ((u+a2)^a1-a2^a1)/((u+a2)^a1+(M+a2)^a1-2*a2^a1)-Qhv
    }
    Qht=uniroot(hfun,interval=c(10,400))$root
    
    QEpv = sum(yord*Ev)
    efun=function(u)
    {
      ((u+a2)^a1-a2^a1)/((u+a2)^a1+(M+a2)^a1-2*a2^a1)-QEpv
    }
    QEpt=uniroot(efun,interval=c(0,400))$root
    
      
    Qet= quantile(x,probs=c(pv))
    result[k,]=c(Qbt,Qht,QEpt,Qnt,Qet)
    k=k+1
  }
}

rowT1=paste("Threshold=1,","0.990:")
rowT2=paste("            ","0.995:")
rowT3=paste("            ","0.999:")
rowT4=paste("Threshold=2,","0.990:")
rowT5=paste("            ","0.995:")
rowT6=paste("            ","0.999:")
rowT=c(rowT1,rowT2,rowT3,rowT4,rowT5,rowT6)
colT=c("B1","B2","K","N","E")
dimnames(result)=list(rowT,colT)
round(result,3)

  