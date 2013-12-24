plot.pcovr <-
  function(x, ...){
    modsel <- x$modsel
    alpha <- x$alpha
    R <- x$R
    vec <- x$Rvalues
    a <- x$Alphavalues
    
    if (modsel=="sim"){
      Qy2 <- x$Qy2
      plot(a,Qy2[,1],xlab="Weighting Parameter",ylab="Qy2",col=1,type="l",ylim=c(0,1))
      for (r in 1:length(vec)){
        points(a,Qy2[,r],col=r,type="l",lty=r)
      }
      points(alpha,Qy2[which(a==alpha),which(vec==R)],pch=2,col=vec[vec==R])
      text(alpha,Qy2[which(a==alpha),which(vec==R)]+.05,"optimal",col=vec[vec==R])
      legend("topleft",as.character(vec),title="Number of Components",col=1:length(vec),lty=1:length(vec))
      
    } else if (modsel=="seqRcv"){
      Qy2 <- x$Qy2
      plot(vec,Qy2,xlab="Number of Components",ylab="Qy2",col=1,type="l",pch=1,ylim=c(0,1))
      points(R,Qy2[vec==R],pch=2,col=1)
      text(R,Qy2[vec==R]-.05,"optimal")
      legend("topleft",c("Qy2"),col=c(1),lty=c(1))
      
    } else if (modsel=="seq"){
      VAF <- x$VAFsum
      plot(vec,VAF,xlab="Number of Components",ylab="VAFsum",col=1,type="l",pch=1,ylim=c(0,1))
      points(R,VAF[vec==R],pch=2,col=1)
      text(R,VAF[vec==R]-.05,"optimal")
      legend("topleft",c("VAFsum"),col=c(1),lty=c(1))
      
    } else if (modsel=="seqAcv"){
      par(mfrow = c(2,1))
      VAF <- x$VAFsum
      Qy2 <- data.matrix(x$Qy2)
      plot(vec,VAF,xlab="Number of Components",ylab="VAFsum",col=1,type="l",pch=1,ylim=c(0,1))
      points(R,VAF[vec==R],pch=2,col=1)
      text(R,VAF[vec==R]-.05,"optimal")
      legend("topleft",c("VAFsum"),pch=c(1),col=c(1),lty=c(1))
      plot(a,Qy2,xlab="Weighting parameter",ylab="Qy2",col=1,type="l",pch=1,ylim=c(0,1))
      points(alpha,t(Qy2)[a==alpha],pch=2,col=1)
      text(alpha,Qy2[a==alpha]-.05,"optimal")
      legend("topleft",c("Qy2"),col=c(1),lty=c(1))
    } 
  }