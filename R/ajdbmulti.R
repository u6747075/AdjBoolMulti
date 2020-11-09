
#'Adjacency Matrix and Boolean Matrix Multiplication
#'
#'This function calculates adjacency matrix or boolean matrix
#'with desired power, "n". This function can have 2 matrices
#'to multiply rather than 1 matrix.
#'
#' @param mtr squared matrix or proper matrix
#' @param mtr2 proper matrix where dim(mtr)[2] == dim(mtr2)[1]
#' @param n power
#' @param type calculate adjacency matrix if 1, and boolean matrix if 2
#' @return returns multiplied (or powered) adjacency matrix if type chosen was 1,
#' and boolean matrix if type chosen was 2
#' @author Mahiro Y (new R user)
#'
#' @note although the speed of the algorithm of adjacency matrix calculation
#' is 1000 times slower than normal calculation done by'expm',
#' there is no problem about the speed as long as if you don't run huge iteration or
#' recursion.
#' @export

mtrmulti = function(mtr,mtr2=mtr,n=2,type=1){
  isSquare = dim(mtr)[2]==dim(mtr2)[1]
  isNull = any(c(any(is.na(mtr)), any(is.na(mtr2))))
  isNnumeric = is.numeric(n)
  if (!isSquare) {
    stop("Matrix dimension erro")
  }
  if (isNull) {
    stop("one (or more) matrix contains NA")
  }
  if (!isNnumeric) {
    stop("n is not numeric")
  }

  if (type==1) {
    if (any(mtr != mtr2)) {
      mtr = calcMM.W(mtr,mtr2)
    }

    return(MM.W(mtr,n=n))
  }
  else if(type==2){
    if (any(mtr != mtr2)) {
      mtr = mtr%*%mtr2
    }

    return(MM.N(mtr,n))
  }
  else{stop("null type")}
}

calcMM.W=function(mtr,mtr2=mtr){
  rtn = NULL
  for (i in 1:dim(mtr)[1]) {
    for (j in 1:dim(mtr2)[2]) {
      rtn=c(rtn,min(mtr[i,]+mtr2[,j]))
    }
  }
  newM=matrix(rtn,nrow = dim(mtr)[2],byrow = T)
  return(newM)
}
MM.W = function(mtr,mtr2=mtr,n){
  if (n==0) {
    return(mtr)
  }
  nums = NULL
  while (n>1) {
    dd = floor(log(n,base=2))
    n= n-2^dd
    nums = c(nums,dd)
  }

  rtn = NULL
  for(i in nums){
    newM = mtr
    for (j in 1:i) {
      newM=calcMM.W(newM)
    }
    if(any(is.null(rtn))){rtn=newM}
    else{rtn=calcMM.W(rtn,newM)}

  }
  return(rtn)
}

MM.N = function(mtr,n){

  rtn=mtr%^%n
  rtn = (rtn>=1)
  return(1*rtn)
}
