#################################################################################
## functions for simulation
## Nov 2018
#################################################################################


#' Simulate from the type one state-space Model on Stiefel manifold.
#'
#' This function simulates from the type one model on Stiefel manifold.
#' See Details part below.
#'
#' The type one model on Stiefel manifold takes the form:
#' \deqn{\boldsymbol{y}_t \quad = \quad \boldsymbol{\alpha}_t \boldsymbol{\beta} ' \boldsymbol{x}_t + \boldsymbol{B} \boldsymbol{z}_t + \boldsymbol{\varepsilon}_t}
#' \deqn{\boldsymbol{\alpha}_{t+1} | \boldsymbol{\alpha}_{t} \quad \sim \quad ML (p, r, \boldsymbol{\alpha}_{t} \boldsymbol{D})}
#' where \eqn{\boldsymbol{y}_t} is a \eqn{p}-vector of the dependent variable,
#' \eqn{\boldsymbol{x}_t} and \eqn{\boldsymbol{z}_t} are explanatory variables wit dimension \eqn{q_1} and \eqn{q_2},
#' \eqn{\boldsymbol{x}_t} and \eqn{\boldsymbol{z}_t} have no overlap,
#' matrix \eqn{\boldsymbol{B}} is the coefficients for \eqn{\boldsymbol{z}_t},
#' \eqn{\boldsymbol{\varepsilon}_t} is the error vector.
#'
#' The matrices \eqn{\boldsymbol{\alpha}_t} and \eqn{\boldsymbol{\beta}} have dimensions \eqn{p \times r} and \eqn{q_1 \times r}, respectively.
#' Note that \eqn{r} is strictly smaller than both \eqn{p} and \eqn{q_1}.
#' \eqn{\boldsymbol{\alpha}_t} and \eqn{\boldsymbol{\beta}} are both non-singular matrices.
#' \eqn{\boldsymbol{\alpha}_t} is time-varying while \eqn{\boldsymbol{\beta}} is time-invariant.
#'
#' Furthermore, \eqn{\boldsymbol{\alpha}_t} fulfills the condition \eqn{\boldsymbol{\alpha}_t' \boldsymbol{\alpha}_t = \boldsymbol{I}_r},
#' and therefor it evolves on the Stiefel manifold.
#'
#' \eqn{ML (p, r, \boldsymbol{\alpha}_{t} \boldsymbol{D})} denotes the Matrix Langevin distribution or matrix von Mises-Fisher distribution on the Stiefel manifold.
#' Its density function takes the form
#' \deqn{f(\boldsymbol{\alpha_{t+1}} ) = \frac{ \mathrm{etr} \left\{ \boldsymbol{D} \boldsymbol{\alpha}_{t}' \boldsymbol{\alpha_{t+1}} \right\} }{ _{0}F_1 (\frac{p}{2}; \frac{1}{4}\boldsymbol{D}^2 ) }}
#' where \eqn{\mathrm{etr}} denotes \eqn{\mathrm{exp}(\mathrm{tr}())},
#' and \eqn{_{0}F_1 (\frac{p}{2}; \frac{1}{4}\boldsymbol{D}^2 )} is the (0,1)-type hypergeometric function for matrix.
#'
#' Note that the function does not add intercept automatically.
#'
#'
#' @param iT the sample size.
#' @param mX the matrix containing X_t with dimension \eqn{T \times q_1}.
#' @param mZ the matrix containing Z_t with dimension \eqn{T \times q_2}.
#' @param mY initial values of the dependent variable for \code{ik-1} up to 0. If \code{mY = NULL}, then no lagged dependent variables in regressors.
#' @param alpha_0 the initial alpha, \eqn{p \times r}.
#' @param beta the \eqn{\beta} matrix, iqx+ip*ik, y_{1,t-1},y_{1,t-2},...,y_{2,t-1},y_{2,t-2},...
#' @param mB the coefficient matrix \eqn{\boldsymbol{B}} before \code{mZ} with dimension \eqn{p \times q_2}.
#' @param Omega covariance matrix of the errors.
#' @param vD vector of the diagonals of \eqn{D}.
#' @param burnin burn-in sample size (matrix Langevin).
#'
#' @return A list containing the sampled data and the dynamics of alpha.
#'
#' The object is a list containing the following components:
#' \item{dData}{a data.frame of the sampled data}
#' \item{aAlpha}{an array of the \eqn{\boldsymbol{\alpha}_{t}} with the dimension \eqn{T \times p \times r}}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @keywords simulation
#'
#' @examples
#'
#' iT = 50 # sample size
#' ip = 2 # dimension of the dependent variable
#' ir = 1 # rank number
#' iqx=2 # number of variables in X
#' iqz=2 # number of variables in Z
#' ik = 1 # lag length
#'
#' if(iqx==0) mX=NULL else mX = matrix(rnorm(iT*iqx),iT, iqx)
#' if(iqz==0) mZ=NULL else mZ = matrix(rnorm(iT*iqz),iT, iqz)
#' if(ik==0) mY=NULL else mY = matrix(0, ik, ip)
#'
#' alpha_0 = matrix(c(runif_sm(num=1,ip=ip,ir=ir)), ip, ir)
#' beta = matrix(c(runif_sm(num=1,ip=ip*ik+iqx,ir=ir)), ip*ik+iqx, ir)
#' if(ip*ik+iqz==0) mB=NULL else mB = matrix(c(runif_sm(num=1,ip=(ip*ik+iqz)*ip,ir=1)), ip, ip*ik+iqz)
#' vD = 50
#'
#' ret = SimModel1(iT=iT, mX=mX, mZ=mZ, mY=mY, alpha_0=alpha_0, beta=beta, mB=mB, vD=vD)
#'
#' @export
SimModel1 <- function(iT, mX=NULL, mZ=NULL, mY=NULL, alpha_0, beta, mB=NULL, Omega=NULL, vD, burnin=100)
{
  # the dimensions
  ip = dim(alpha_0)[1]
  ir = dim(beta)[2]
  if(ip<=ir) stop(simpleError("p is smaller or equal to r!"))

  iqx = 0; iqz = 0
  if(!is.null(mX)) iqx = dim(mX)[2]
  if(!is.null(mZ)) iqz = dim(mZ)[2]

  if(!is.null(mY)){
    ik = dim(mY)[1]; idx = (0):(-ik+1)
  } else{
    ik = 0; idx = NULL
  }

  # the lagged dependent variable will be put into X
  iqxx = iqx+ik*ip
  if(iqxx<=ir) stop(simpleError("q_1 is smaller or equal to r!"))

  # the lagged dependent variable will be put into Z
  iqzz = iqz+ik*ip

  # initialize the data.frame
  tmp = rbind(mY,matrix(0,iT,ip))
  colnames(tmp) = paste('Y',1:ip,sep='')
  rownames(tmp) = c(rev(idx),1:iT)
  dData = data.frame(tmp)

  if(!is.null(mX)){
    tmp = rbind(matrix(0,ik,iqx), mX)
    if(is.null(colnames(mX))) colnames(tmp) = paste('X',1:iqx,sep='')
    else  colnames(tmp) = colnames(mX)
    dData = cbind(dData,tmp)
  }

  if(!is.null(mZ)){
    tmp = rbind(matrix(0,ik,iqz), mZ)
    if(is.null(colnames(mZ))) colnames(tmp) = paste('Z',1:iqz,sep='')
    else  colnames(tmp) = colnames(mZ)
    dData = cbind(dData,tmp)
  }

  if(is.null(Omega)) Omega = diag(ip)
  mE = rbind(matrix(0,ik,ip),matrix(rnorm(iT*ip),iT,ip))
  mE = mE %*% chol(Omega)
  colnames(mE) = paste('E',1:ip,sep='')
  dData = cbind(dData,mE)

  # the simulated alpha always starts from 2
  aAlpha = array(0,dim=c(iT+1,ip,ir))
  aAlpha[1,,] = alpha_0

  mD = diag(vD,ir,ir)

  # start simulation
  if(is.null(mB)){
    for(iter in 1:iT){
      tmp1 = matrix(aAlpha[iter,,],ip,ir)
      mF = tmp1%*%mD
      tmp1 = rmLB_sm(num=burnin,mJ=matrix(0,ip,ip),mH=matrix(0,ir,ir),mC=mF,mX=tmp1,ir=ir)
      tmp1 = matrix(tmp1[burnin,,],ip,ir)
      aAlpha[iter+1,,] = tmp1
      tmp2 = c(as.matrix(dData[idx+iter+ik-1,1:ip]))
      # lag first, then X or Z
      dData[iter+ik,1:ip] = tmp1%*%t(beta)%*%c(tmp2,mX[iter,]) + c(mE[iter+ik,])
    }
  }else{
    for(iter in 1:iT){
      tmp1 = matrix(aAlpha[iter,,],ip,ir)
      mF = tmp1%*%mD
      tmp1 = rmLB_sm(num=burnin,mJ=matrix(0,ip,ip),mH=matrix(0,ir,ir),mC=mF,mX=tmp1,ir=ir)
      tmp1 = matrix(tmp1[burnin,,],ip,ir)
      aAlpha[iter+1,,] = tmp1
      tmp2 = c(as.matrix(dData[idx+iter+ik-1,1:ip]))
      # lag first, then X or Z
      dData[iter+ik,1:ip] = tmp1%*%t(beta)%*%c(tmp2,mX[iter,]) + mB%*%c(tmp2,mZ[iter,]) + c(mE[iter+ik,])
    }
  }

  return(list(dData=dData,aAlpha=array(aAlpha[2:(iT+1),,],dim=c(iT,ip,ir))))
}


#' Simulate from the type two state-space Model on Stiefel manifold.
#'
#' This function simulates from the type two model on Stiefel manifold.
#' See Details part below.
#'
#' The type two model on Stiefel manifold takes the form:
#' \deqn{\boldsymbol{y}_t \quad = \quad \boldsymbol{\alpha} \boldsymbol{\beta}_t ' \boldsymbol{x}_t + \boldsymbol{B}' \boldsymbol{z}_t + \boldsymbol{\varepsilon}_t}
#' \deqn{\boldsymbol{\beta}_{t+1} | \boldsymbol{\beta}_{t} \quad \sim \quad ML (q_1, r, \boldsymbol{\beta}_{t} \boldsymbol{D})}
#' where \eqn{\boldsymbol{y}_t} is a \eqn{p}-vector of the dependent variable,
#' \eqn{\boldsymbol{x}_t} and \eqn{\boldsymbol{z}_t} are explanatory variables wit dimension \eqn{q_1} and \eqn{q_2},
#' \eqn{\boldsymbol{x}_t} and \eqn{\boldsymbol{z}_t} have no overlap,
#' matrix \eqn{\boldsymbol{B}} is the coefficients for \eqn{\boldsymbol{z}_t},
#' \eqn{\boldsymbol{\varepsilon}_t} is the error vector.
#'
#' The matrices \eqn{\boldsymbol{\alpha}} and \eqn{\boldsymbol{\beta}_t} have dimensions \eqn{p \times r} and \eqn{q_1 \times r}, respectively.
#' Note that \eqn{r} is strictly smaller than both \eqn{p} and \eqn{q_1}.
#' \eqn{\boldsymbol{\alpha}} and \eqn{\boldsymbol{\beta}_t} are both non-singular matrices.
#' \eqn{\boldsymbol{\beta}_t} is time-varying while \eqn{\boldsymbol{\alpha}} is time-invariant.
#'
#' Furthermore, \eqn{\boldsymbol{\beta}_t} fulfills the condition \eqn{\boldsymbol{\beta}_t' \boldsymbol{\beta}_t = \boldsymbol{I}_r},
#' and therefor it evolves on the Stiefel manifold.
#'
#' \eqn{ML (p, r, \boldsymbol{\beta}_t \boldsymbol{D})} denotes the Matrix Langevin distribution or matrix von Mises-Fisher distribution on the Stiefel manifold.
#' Its density function takes the form
#' \deqn{f(\boldsymbol{\beta_{t+1}} ) = \frac{ \mathrm{etr} \left\{ \boldsymbol{D} \boldsymbol{\beta}_{t}' \boldsymbol{\beta_{t+1}} \right\} }{ _{0}F_1 (\frac{p}{2}; \frac{1}{4}\boldsymbol{D}^2 ) }}
#' where \eqn{\mathrm{etr}} denotes \eqn{\mathrm{exp}(\mathrm{tr}())},
#' and \eqn{_{0}F_1 (\frac{p}{2}; \frac{1}{4}\boldsymbol{D}^2 )} is the (0,1)-type hypergeometric function for matrix.
#'
#' Note that the function does not add intercept automatically.
#'
#'
#' @param iT the sample size.
#' @param mX the matrix containing X_t with dimension \eqn{T \times q_1}.
#' @param mZ the matrix containing Z_t with dimension \eqn{T \times q_2}.
#' @param mY initial values of the dependent variable for \code{ik-1} up to 0. If \code{mY = NULL}, then no lagged dependent variables in regressors.
#' @param beta_0 the initial beta, iqx+ip*ik, y_{1,t-1},y_{1,t-2},...,y_{2,t-1},y_{2,t-2},....
#' @param alpha the \eqn{\alpha} matrix, \eqn{p \times r}.
#' @param mB the coefficient matrix \eqn{\boldsymbol{B}} before \code{mZ} with dimension \eqn{p \times q_2}.
#' @param Omega covariance matrix of the errors.
#' @param vD vector of the diagonals of \eqn{D}.
#' @param burnin burn-in sample size (matrix Langevin).
#'
#' @return A list containing the sampled data and the dynamics of beta.
#'
#' The object is a list containing the following components:
#' \item{dData}{a data.frame of the sampled data}
#' \item{aBeta}{an array of the \eqn{\boldsymbol{\beta}_t} with the dimension \eqn{T \times q_1 \times r}}
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @keywords simulation
#'
#' @examples
#'
#' iT = 50
#' ip = 2
#' ir = 1
#' iqx =3
#' iqz=2
#' ik = 1
#'
#' if(iqx==0) mX=NULL else mX = matrix(rnorm(iT*iqx),iT, iqx)
#' if(iqz==0) mZ=NULL else mZ = matrix(rnorm(iT*iqz),iT, iqz)
#' if(ik==0) mY=NULL else mY = matrix(0, ik, ip)
#'
#' alpha = matrix(c(runif_sm(num=1,ip=ip,ir=ir)), ip, ir)
#' beta_0 = matrix(c(runif_sm(num=1,ip=ip*ik+iqx,ir=ir)), ip*ik+iqx, ir)
#' if(ip*ik+iqz==0) mB=NULL else mB = matrix(c(runif_sm(num=1,ip=(ip*ik+iqz)*ip,ir=1)), ip, ip*ik+iqz)
#' vD = 50
#'
#' ret = SimModel2(iT=iT, mX=mX, mZ=mZ, mY=mY, alpha=alpha, beta_0=beta_0, mB=mB, vD=vD)
#'
#' @export
SimModel2 <- function(iT, mX=NULL, mZ=NULL, mY=NULL, beta_0, alpha, mB=NULL, Omega=NULL, vD, burnin=100)
{
  if(!is.null(mY)){
    ik = dim(mY)[1]; idx = (0):(-ik+1)
  } else{
    ik = 0; idx = NULL
  }

  # the dimensions
  ip = dim(alpha)[1]
  ir = dim(alpha)[2]
  if(ip<=ir) stop(simpleError("p is smaller or equal to r!"))

  iqx = 0; iqz = 0
  if(!is.null(mX)) iqx = dim(mX)[2]
  if(!is.null(mZ)) iqz = dim(mZ)[2]

  if(!is.null(mY)){
    ik = dim(mY)[1]; idx = (0):(-ik+1)
  } else{
    ik = 0; idx = NULL
  }

  # the lagged dependent variable will be put into X
  iqxx = iqx+ik*ip
  if(iqxx<=ir) stop(simpleError("q_1 is smaller or equal to r!"))

  # the lagged dependent variable will be put into Z
  iqzz = iqz+ik*ip

  # initialize the data.frame
  tmp = rbind(mY,matrix(0,iT,ip))
  colnames(tmp) = paste('Y',1:ip,sep='')
  rownames(tmp) = c(rev(idx),1:iT)
  dData = data.frame(tmp)

  if(!is.null(mX)){
    tmp = rbind(matrix(0,ik,iqx), mX)
    if(is.null(colnames(mX))) colnames(tmp) = paste('X',1:iqx,sep='')
    else  colnames(tmp) = colnames(mX)
    dData = cbind(dData,tmp)
  }

  if(!is.null(mZ)){
    tmp = rbind(matrix(0,ik,iqz), mZ)
    if(is.null(colnames(mZ))) colnames(tmp) = paste('Z',1:iqz,sep='')
    else  colnames(tmp) = colnames(mZ)
    dData = cbind(dData,tmp)
  }

  if(is.null(Omega)) Omega = diag(ip)
  mE = rbind(matrix(0,ik,ip),matrix(rnorm(iT*ip),iT,ip))
  mE = mE %*% chol(Omega)
  colnames(mE) = paste('E',1:ip,sep='')
  dData = cbind(dData,mE)

  # the simulated beta always starts from 2
  aBeta = array(0,dim=c(iT+1,iqxx,ir))
  aBeta[1,,] = beta_0

  mD = diag(vD,ir,ir)

  # start simulation
  if(is.null(mB)){
    for(iter in 1:iT){
      tmp1 = matrix(aBeta[iter,,],iqxx,ir)
      mF = tmp1%*%mD
      tmp1 = rmLB_sm(num=burnin,mJ=matrix(0,iqxx,iqxx),mH=matrix(0,ir,ir),mC=mF,mX=tmp1,ir=ir)
      tmp1 = matrix(tmp1[burnin,,],iqxx,ir)
      aBeta[iter+1,,] = tmp1
      tmp2 = c(as.matrix(dData[idx+iter+ik-1,1:ip]))
      # lag first, then X or Z
      dData[iter+ik,1:ip] = alpha%*%t(tmp1)%*%c(tmp2,mX[iter,]) + c(mE[iter+ik,])
    }
  }else{
    for(iter in 1:iT){
      tmp1 = matrix(aBeta[iter,,],iqxx,ir)
      mF = tmp1%*%mD
      tmp1 = rmLB_sm(num=burnin,mJ=matrix(0,iqxx,iqxx),mH=matrix(0,ir,ir),mC=mF,mX=tmp1,ir=ir)
      tmp1 = matrix(tmp1[burnin,,],iqxx,ir)
      aBeta[iter+1,,] = tmp1
      tmp2 = c(as.matrix(dData[idx+iter+ik-1,1:ip]))
      # lag first, then X or Z
      dData[iter+ik,1:ip] = alpha%*%t(tmp1)%*%c(tmp2,mX[iter,]) + mB%*%c(tmp2,mZ[iter,]) + c(mE[iter+ik,])
    }
  }

  return(list(dData=dData,aBeta=array(aBeta[2:(iT+1),,],dim=c(iT,iqxx,ir))))
}
