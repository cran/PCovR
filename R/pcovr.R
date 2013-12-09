pcovr <-
function(X,Y,modsel="seq",Rmin=1,Rmax=ncol(X)/3,R=NULL,weight=NULL,rot="varimax", target=NULL, prepX="stand",prepY="stand", ratio=1, fold="LeaveOneOut"){
  UseMethod("pcovr")
}
