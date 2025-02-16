# Simulation Study for the paper

# Estimating the Quantile Function Using Beta-Kernel Smoother

# Data are generated from a Beta-Distribution

set.seed(6789)
total=500

cv=0.02

pvseq=seq(0.01,0.99,length=20)


result = matrix(0,nrow=32,ncol=20)
par(mfrow=c(2,2))
rn = 1
for (n in c(50, 100, 200,500))
{ 
  for (beta in c(1.5, 2))
  {
    cn = 1
    for (pv in pvseq)
    {
      alpha = 1.5
      Qp = qbeta(pv,alpha,beta)
      Qnt = Qb1t = Qb2t = Qet = QEpt = rep(total)
      for (k in seq(total))
      {
        x = rbeta(n,alpha,beta)    
        
        hb1=hb2=cv*n^(-2/3)
        hn=cv*n^(-1/3)

        xord = sort(x)
        iseq = seq(n)/n
        jseq = (seq(n)-1)/n         
        
        B1v=pbeta(iseq,pv/hb1+1,(1-pv)/hb1+1)-
          pbeta(jseq,pv/hb1+1,(1-pv)/hb1+1)
        B2v=pbeta(iseq,pv/hb2,(1-pv)/hb2)-
          pbeta(jseq,pv/hb2,(1-pv)/hb2)
        Nv=pnorm(iseq,pv,hn) - pnorm(jseq,pv,hn)
        
        Ev=rep(0,n)
        for(kk in seq(n))
        {
          a=pv-hn
          b=pv+hn
          c=(kk-1)/n
          d=kk/n
          
          lend=max(a,c)
          rend=min(b,d)
          
          
          if(lend>=rend)
           {Ev[kk]=0}
          if(lend<rend)
           {
            funu=function(u)
             {
              0.75*(1-((u-pv)/hn)^2)/hn
             }
            Ev[kk]=integrate(funu,lend,rend)$value
           }
        }
       
        
        Qb1t[k] = sum(xord*B1v)
        Qb2t[k] = sum(xord*B2v)
        QEpt[k] = sum(xord*Ev)
        Qnt[k] = sum(xord*Nv)
        Qet[k] = quantile(x,probs=c(pv))
      }
      
      mseb1 = mean((Qb1t-Qp)^2)
      mseb2 = mean((Qb2t-Qp)^2)
      mseep = mean((QEpt-Qp)^2)
      msen = mean((Qnt-Qp)^2)
      
      msee = mean((Qet-Qp)^2)
      MSEr = c(mseb1,mseb2,mseep,msen)/msee
      result[rn:(rn+3),cn]=MSEr
      cn = cn+1
      
      # if((n==500)&((pv==0.05)|(pv==0.5)))
      # {
      #   pvb=paste("p=",toString(pv),",",sep="")
      #   pvc=paste("beta=",toString(beta),",",sep="")
      #   pva=paste(pvb,pvc)
      #   boxplot(Qb1t,Qb2t,Qnt,Qet,names=c('B1k','B2k','Nk','E'),xlab=pva)
      #   abline(h=Qp,col="red",lwd=2)
      # } 
    }
    rn=rn+4 
  }
}

rowT1 =paste("n=50, ","Beta(1.5,1.5)"," B1/E")
rowT2 =paste("      ","             "," B2/E")
rowT3 =paste("      ","             "," EP/E")
rowT4 =paste("      ","             "," N/E")
rowT5 =paste("n=50, ","Beta(1.5,2.0)"," B1/E")
rowT6 =paste("      ","             "," B2/E")
rowT7 =paste("      ","             "," EP/E")
rowT8 =paste("      ","             "," N/E")
rowT9 =paste("n=100,","Beta(1.5,1.5)"," B1/E")
rowT10 =paste("      ","             "," B2/E")
rowT11 =paste("      ","             "," EP/E")
rowT12 =paste("      ","             "," N/E")
rowT13=paste("n=100,","Beta(1.5,2.0)"," B1/E")
rowT14=paste("      ","             "," B2/E")
rowT15 =paste("      ","             "," EP/E")
rowT16=paste("      ","             "," N/E")
rowT17=paste("n=200,","Beta(1.5,1.5)"," B1/E")
rowT18=paste("      ","             "," B2/E")
rowT19 =paste("      ","             "," EP/E")
rowT20=paste("      ","             "," N/E")
rowT21=paste("n=200,","Beta(1.5,2.0)"," B1/E")
rowT22=paste("      ","             "," B2/E")
rowT23 =paste("      ","             "," EP/E")
rowT24=paste("      ","             "," N/E")
rowT25=paste("n=500,","Beta(1.5,1.5)"," B1/E")
rowT26=paste("      ","             "," B2/E")
rowT27 =paste("      ","             "," EP/E")
rowT28=paste("      ","             "," N/E")
rowT29=paste("n=500,","Beta(1.5,2.0)"," B1/E")
rowT30=paste("      ","             "," B2/E")
rowT31 =paste("      ","             "," EP/E")
rowT32=paste("      ","             "," N/E")
colT=pvseq

dimnames(result)=list(
  c(rowT1,rowT2,rowT3,rowT4,
    rowT5,rowT6,rowT7,rowT8,
    rowT9,rowT10,rowT11,rowT12,
    rowT13,rowT14,rowT15,rowT16,
    rowT17,rowT18,rowT19,rowT20,
    rowT21,rowT22,rowT23,rowT24,
    rowT25,rowT26,rowT27,rowT28,
    rowT29,rowT30,rowT31,rowT32
    ),
  colT)

cat("c=",cv,"\n")  
print(round(result,3))

nsize=c(50,50,100,100,200,200,500,500)
al=as.character(1.5)
beta=c(1.5,2,1.5,2,1.5,2,1.5,2)
pvseq=seq(0.01,0.99,length=20)

#datatext=write.table(result, file="C:\\Users\\weixing\\Dropbox\\Research\\Research(2022)\\Smoothing Quantile\\Simulation\\matrix_data.txt", row.names=FALSE, col.names=FALSE)

par(mfrow=c(2,3))
for(j in c(1,3,5,2,4,6))
{
    ns=as.character(nsize[j])
    bt=as.character(beta[j])
    dm =t(result[((j-1)*4+1):((j-1)*4+4),])

    # Correct xlab using bquote()
    matplot(pvseq, dm, lty=c(1,2,1,2), 
            type="l", lwd=1, bty = "l", col=c("red","blue","black","darkgrey"),
            xlab=bquote("n=" ~ .(ns) ~ ", " ~ alpha ~ "=" ~ .(al) ~ ", " ~ beta ~ "=" ~ .(bt)),
            ylab="Relative Efficiency", ylim=c(0.5,1.1))

    # Legend without border
    #legend("bottom", legend=c("B1/E","B2/E","N/E"), lty=1:3, bty="n", col=c("black","black","black"))
}