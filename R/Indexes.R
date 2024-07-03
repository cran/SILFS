
vectornorm2 <- function(x){
  temp <- 0
  for (i in 1:length(x)){
    temp <- temp + x[i]^2
  }
  v <- sqrt(temp)
  return(v)
}



CM_return <- function(x,y){
  if((x+y)==0){
    return(0)
  } else{
    return(x/(x+y))
  }
}




CM <- function(xhat,x,eps=0.05){
  x=abs(x)
  xhat=abs(xhat)
  result=matrix(0,1,6)
  colnames(result)=c("TPR","TNR","FPR","FNR","FDR","Precision")
  nz_loc=which(x>=eps)

  TP=sum(xhat[nz_loc]>=eps)
  TN=sum(xhat[-nz_loc]<eps)
  FP=sum(xhat[-nz_loc]>=eps)
  FN=sum(xhat[nz_loc]<eps)

  result[1,1]=TP/(TP+FN)
  result[1,2]=TN/(TN+FP)
  result[1,3]=FP/(FP+TN)
  result[1,4]=FN/(TP+FN)
  result[1,5]=CM_return(FP,TP)    #FP/(FP+TP)
  result[1,6]=CM_return(TP,FP)
  return(result)
}



EVMSE<-function(tv,est){
  return(vectornorm2(tv-est)/sqrt(length(tv)))
}
