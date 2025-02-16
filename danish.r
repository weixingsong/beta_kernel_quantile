# Simulation Study for the paper

# Enhancing Quantile Function Estimation with Beta-Kernel Smoothing

# Danish Fire Data without transformation

# Table 2


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
     
     #hb=0.01*n^(-2/3)
     #hn=0.01*n^(-1/3)
     #hv=0.01*n^(-2/3)
     #he=0.01*n^(-1/3)
     
      booth=hboot(x,pv,B=500,hleft=0.00001,hright=0.5)
      hb=booth$hb1
      hv=booth$hb2
      hn=booth$hn
      he=booth$he
      cat(hb,hv,he,hn,"\n")
     
     xord = sort(x)
     iseq = seq(n)/n
     jseq = (seq(n)-1)/n         
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

     Qbt = sum(xord*Bv)
     Qnt = sum(xord*Nv)
     Qht = sum(xord*Hv)
     QEpt = sum(xord*Ev)
     Qet = quantile(x,probs=c(pv))
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
