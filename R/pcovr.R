pcovr <- 
  function(X,Y,modsel="seq",Rmin=1,Rmax=NULL,weight=NULL,rot="varimax", target=NULL, prep="stand", ratio=1, fold="LeaveOneOut"){
    library(GPArotation)
    library(ThreeWay)
    J <- ncol(X)
    N <- nrow(X)
    K <- ncol(Y)
    
    if (is.null(weight)){
      a <- seq(from=.01, to=1, by=.01)
    } else {
      a <- weight
    }
    
    if (is.null(Rmax)){
      Rmax <- J/3
    }
    vec <- Rmin:Rmax
    
    
    # PREPROCESSING
    Z <- switch(prep, 
                stand=nrm2(scale(cbind(X,Y), center = TRUE, scale = FALSE))*N^(1/2), 
                cent=scale(cbind(X,Y), center = TRUE, scale = FALSE))
    X <- Z[,1:J]
    Y <- array(Z[,(J+1):(J+K)],c(N,K))
    
    
    # MODEL SELECTION
    AlphaMaxLik <- (SUM(X)$ssq/ (SUM(X)$ssq + SUM(Y)$ssq*ratio))
    alpha <- a[which.min(abs(a - AlphaMaxLik))]
    
    if (modsel=="sim"){
      CV <- array(NA,c(length(a),length(vec)))
      for (l in 1:length(vec)){
        for (w in 1:length(a)){
          para <- pcovr_est(X,Y, r = vec[l], a = a[w], cross = TRUE,fold)
          CV[w,l] <- para$Qy2
        }
      }
      alpha <- a[arrayInd(which.max(CV),dim(CV))[1]]
      R <- vec[arrayInd(which.max(CV),dim(CV))[2]]
      
      plot(a,CV[,1],xlab="Weighting Parameter",ylab="Qy2",col=1,type="b",pch=1,ylim=c(0,1))
      for (r in 1:length(vec)){
        points(a,CV[,r],col=r,type="b",lty=r,pch=1)
      }
      points(alpha,CV[which(a==alpha),which(vec==R)],pch=2,col=vec[vec==R])
      text(alpha,CV[which(a==alpha),which(vec==R)]+.05,"optimal",col=vec[vec==R])
      legend("topleft",as.character(vec),title="Number of Components",col=1:length(vec),lty=1:length(vec),pch=1:length(vec))
      
    } else if (modsel=="seqRcv"){
      CV <- array(NA,c(1,length(vec)))
      for (l in 1:length(vec)){
        para <- pcovr_est(X,Y,vec[l],alpha,TRUE,fold)
        CV[1,l] <- para$Qy2
      }
      R <- vec[which.max(CV)]
      
      plot(vec,CV,xlab="Number of Components",ylab="Qy2",col=1,type="b",pch=1,ylim=c(0,1))
      points(R,CV[vec==R],pch=2,col=1)
      text(R,CV[vec==R]-.05,"optimal")
      legend("topleft",c("Qy2"),pch=c(1),col=c(1),lty=c(1))
      
    } else {
      VAF <- array(NA,c(1,length(vec)))
      scr <- array(NA,c(1,length(vec)))
      for (l in 1:length(vec)){
        para <- pcovr_est(X,Y,vec[l], a = alpha, cross = FALSE,fold)
        VAF[1,l] <- alpha*para$Rx2 + (1-alpha)*para$Ry2
      }
      if (length(vec)>2){
        for (u in 2:(length(vec)-1)){
          scr[,u]=(VAF[u]-VAF[u-1])/(VAF[u+1]-VAF[u])
        }
      } else{
        scr <- VAF
      }
      R <- vec[which.max(scr)]
      
      plot(vec,VAF,xlab="Number of Components",ylab="VAFsum",col=1,type="b",pch=1,ylim=c(0,1))
      points(R,VAF[vec==R],pch=2,col=1)
      text(R,VAF[vec==R]-.05,"optimal")
      legend("topleft",c("VAFsum"),pch=c(1),col=c(1),lty=c(1))
    }
    
    if (modsel=="seqAcv"){
      CV <- array(NA,c(length(a),1))
      for (w in 1:length(a)){
        para <- pcovr_est(X,Y,R,a[w],TRUE,fold)
        CV[w,1] <- para$Qy2
      }
      alpha <- a[which.max(CV)]
      
      plot(a,CV,xlab="Weighting parameter",ylab="Qy2",col=1,type="b",pch=1,ylim=c(0,1))
      points(alpha,CV[a==alpha],pch=2,col=1)
      text(alpha,CV[a==alpha]-.05,"optimal")
      legend("topleft",c("Qy2"),pch=c(1),col=c(1),lty=c(1))
    }
    
    
    # ANALYSES
    results <- list()
    para <- pcovr_est(X,Y,r = R, a = alpha, cross = FALSE,fold)
    Te <- para$Te * N^(1/2)
    W <- para$W * N^(1/2)
    Px <- para$Px * N^(-1/2)
    Py <- para$Py * N^(-1/2)
    
    #rotate parameters
    if ((R>1) & (rot!="none")){
      Bpx <- t(Px)
      if (rot=="quartimin"){
        rotation <- GPFoblq(Bpx, Tmat=diag(ncol(Bpx)))
      }
      if (rot=="varimax"){
        rotation <- GPForth(Bpx, Tmat=diag(ncol(Bpx)))
      }
      if (rot=="target"){
        rotation <- targetT(Bpx, Tmat=diag(ncol(Bpx)), target)
      }
      rotE <- t(ginv(rotation$Th))
      Px <- t(rotation$loadings)
      Te <- Te %*% t(ginv(rotE))
      W <- W %*% t(ginv(rotE))
      Py <- t(rotE) %*% Py
    }
    
    #results
    Rx2 <- para$Rx2
    Ry2 <- para$Ry2
    
    if (modsel == "sim"){
      data.frame(CV)
      rownames(CV) <- a
      colnames(CV) <- vec
      results <- list(Px=Px,Py=Py,Te=Te,W=W,Rx2=Rx2,Ry2=Ry2,Qy2=CV,alpha=alpha,R=R)
    } else if (modsel=="seqRcv"){
      data.frame(CV)
      colnames(CV) <- vec
      results <- list(Px=Px,Py=Py,Te=Te,W=W,Rx2=Rx2,Ry2=Ry2,Qy2=CV,alpha=alpha,R=R)
    } else if (modsel=="seqAcv"){
      data.frame(CV)
      rownames(CV) <- a
      data.frame(VAF)
      colnames(VAF) <- vec
      results <- list(Px=Px,Py=Py,Te=Te,W=W,Rx2=Rx2,Ry2=Ry2,Qy2=CV,VAFsum=VAF,alpha=alpha,R=R)
    } else {
      data.frame(VAF)
      colnames(VAF) <- vec
      results <- list(Px=Px,Py=Py,Te=Te,W=W,Rx2=Rx2,Ry2=Ry2,VAFsum=VAF,alpha=alpha,R=R)
    }
    return(results)
  }

