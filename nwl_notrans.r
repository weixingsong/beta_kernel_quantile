# Simulation Study for the paper

# Enhancing Quantile Function Estimation with Beta-Kernel Smoothing

# Normal-Webull-Lognormal: Relative Efficiency Ratio Plot 
# No Transformation
# Figures 9-12

set.seed(6789)
cv = 0.03
n = 200

pvseq=seq(0.01,0.99,length=20)

total=500
par(mfrow=c(3,3))
result = matrix(0,nrow=12,ncol=20)
cn = 1
for(pv in pvseq)
{ 
  #if(abs(pv-0.7)<=0.05)
  #{
  #  x11()
  #  par(mfrow=c(3,3))
  #}
  rn = 1
  for (model in c(1,2,3))
  {
    Qpn = qnorm(pv,5,1)
    Qpw = qweibull(pv,1.5,1)
    Qpl = qlnorm(pv,0,0.5)
    Qp = Qpn*(model==1)+Qpw*(model==2)+Qpl*(model==3)
    
    Qnt = Qbt = Qet = Qht =  QEpt = rep(total)
    for (k in seq(total))
    {
      xn = rnorm(n,5,1)
      xw = rweibull(n,1.5,1)
      xl = rlnorm(n,0,0.5)
      x = xn*(model==1)+xw*(model==2)+xl*(model==3)
    
      hb = cv*n^(-2/3)
      hn = cv*n^(-1/3)
      hv = cv*n^(-2/3)
      he = cv*n^(-1/3)
      
      xord = sort(x)
      iseq = seq(n)/n
      jseq = (seq(n)-1)/n         
      
      Bv = pbeta(iseq,pv/hb+1,(1-pv)/hb+1)-
        pbeta(jseq,pv/hb+1,(1-pv)/hb+1)
      Nv = pnorm(iseq,pv,hn) - pnorm(jseq,pv,hn)
      Hv = pbeta(iseq,pv/hv,(1-pv)/hv)-
        pbeta(jseq,pv/hv,(1-pv)/hv)
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
          bb=(rend-pv)/he
          aa=(lend-pv)/he
          Ev[kk]=((bb-aa)-(bb^3-aa^3)/3)*0.75
        }
      }
      
      Qbt[k] = sum(xord*Bv)
     
      Qnt[k] = sum(xord*Nv)
      
      Qht[k] = sum(xord*Hv)
      
      QEpt[k] = sum(xord*Ev)
      
      Qet[k] = quantile(x,probs=c(pv))
    }
    
    
    
    mseb = mean((Qbt-Qp)^2)    
    mseh = mean((Qht-Qp)^2)
    mseep = mean((QEpt-Qp)^2)
    msen = mean((Qnt-Qp)^2)
    msee = mean((Qet-Qp)^2)
    
    MSEr = c(mseb,mseh,mseep,msen)/msee
    
    # if ((pv == 0.2)|(pv == 0.3)|(pv == 0.5)|(pv == 0.7)|(pv == 0.8)|(pv == 0.95))
    # {
    #   if(model==1)
    #   {
    #     distri=' Normal'
    #   }
    #   else if (model==2)
    #   {
    #     distri=" Weibull" 
    #   }  
    #   else
    #   {
    #     distri=" Lognormal" 
    #   }
    #   pva=paste("p=",toString(pv),",",sep="")
    #   pva=paste(pva,distri)
    #   boxplot(Qbt,Qnt,Qet,names=c('Bk','Nk','E'),xlab=pva)
    #   abline(h=Qp,lwd=2,lty=2)  
    # }
    
    result[rn:(rn+3),cn]=MSEr
    rn = rn+4
  }  
  cn=cn+1
  cat(cn,"\n")
}  


result=round(result,3)
pvseq= c(0.05,seq(0.1,0.9,by=0.1),0.95)
#result1=result[1:3,]
#dimnames(result1)=list(c('Normal: MSEb1/MSEe','        MSEb2/MSEe','         MSEn/MSEe'),pvseq)
#result2=result[4:6,]
#dimnames(result2)=list(c('Weibull: MSEb1/MSEe','         MSEb2/MSEe','          MSEn/MSEe'),pvseq)
#result3=result[7:9,]
#dimnames(result3)=list(c('Lognormal: MSEb1/MSEe','           MSEb2/MSEe','            MSEn/MSEe'),pvseq)
#rbind(result1,result2,result3)


al=as.character(1.5)
beta=c(1.5,2,1.5,2,1.5,2,1.5,2)
pvseq=seq(0.01,0.99,length=20)
cv=as.character(cv)
ns=as.character(n)

#datatext=write.table(result, file="C:\\Users\\weixing\\Dropbox\\Research\\Research(2022)\\Smoothing Quantile\\Simulation\\matrix_data.txt", row.names=FALSE, col.names=FALSE)

par(mfrow=c(1,3))
for(j in c(1,2,3))
{
    dist = c("N(5,1)", "Weibull(1,1.5)", "LN(0,0.5)")
    dm = t(result[((j-1)*4+1):((j-1)*4+4),])

    # Correct xlab using bquote()
    matplot(pvseq, dm, lty=c(1,2,1,2), type="l", lwd=1,  bty = "l",
            col=c("red","blue","black","darkgrey"),
            xlab=bquote("n=" ~ .(ns) ~ ", "~ "c=" ~ .(cv) ~ ", " ~ .(dist[j])),  
            ylab="Relative Efficiency", ylim=c(0.6,1.1)) 

    # Legend without border
    #legend("bottom", legend=c("B1/E","B2/E","N/E"), lty=1:3, bty="n", col=c("black","black","black"))
}


