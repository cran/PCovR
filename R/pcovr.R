pcovr <-
function(X,Y,modsel="seq",Rmin=1,Rmax=NULL,weight=NULL,rot="varimax", target=NULL, prep="stand", ratio=1, fold="LeaveOneOut"){
  UseMethod("pcovr")
}
