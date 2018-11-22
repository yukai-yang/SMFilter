#################################################################################
## filtering algorithms
## Nov 2018
#################################################################################


#' Filtering algorithm for the type one model.
#'
#' This function implements the filtering algorithm for the type one model.
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
#'
#' @param mY the matrix containing Y_t with dimension \eqn{T \times p}.
#' @param mX the matrix containing X_t with dimension \eqn{T \times q_1}.
#' @param mZ the matrix containing Z_t with dimension \eqn{T \times q_2}.
#' @param beta the \eqn{\beta} matrix.
#' @param mB the coefficient matrix \eqn{\boldsymbol{B}} before \code{mZ} with dimension \eqn{p \times q_2}.
#' @param Omega covariance matrix of the errors.
#' @param vD vector of the diagonals of \eqn{D}.
#' @param U0 initial value of the alpha sequence.
#' @param method a string representing the optimization method from c('max_1','max_2','max_3','min_1','min_2').
#'
#' @return an array \code{aAlpha} containing the modal orientations of alpha in the prediction step.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @keywords filtering
#'
#' @examples
#'
#' iT = 50
#' ip = 2
#' ir = 1
#' iqx = 4
#' iqz=0
#' ik = 0
#' Omega = diag(ip)*.1
#'
#' if(iqx==0) mX=NULL else mX = matrix(rnorm(iT*iqx),iT, iqx)
#' if(iqz==0) mZ=NULL else mZ = matrix(rnorm(iT*iqz),iT, iqz)
#' if(ik==0) mY=NULL else mY = matrix(0, ik, ip)
#'
#' alpha_0 = matrix(c(runif_sm(num=1,ip=ip,ir=ir)), ip, ir)
#' beta = matrix(c(runif_sm(num=1,ip=ip*ik+iqx,ir=ir)), ip*ik+iqx, ir)
#' mB=NULL
#' vD = 100
#'
#' ret = SimModel1(iT=iT, mX=mX, mZ=mZ, mY=mY, alpha_0=alpha_0, beta=beta, mB=mB, vD=vD, Omega=Omega)
#' mYY=as.matrix(ret$dData[,1:ip])
#' fil = FilterModel1(mY=mYY, mX=mX, mZ=mZ, beta=beta, mB=mB, Omega=Omega, vD=vD, U0=alpha_0)
#'
#' @export
FilterModel1 <- function(mY, mX, mZ, beta, mB=NULL, Omega, vD, U0, method='max_1')
{
  iT = nrow(mY)
  ip = ncol(mY)
  ir = ncol(beta)

  aAlpha = array(0,dim=c(iT+1,ip,ir))
  mJ = chol2inv(chol(Omega))
  mD = diag(vD,ir,ir)

  if(is.null(mB)){
    mBZ = matrix(0,iT,ip)
  }else{
    mBZ = NULL
    for(iter in 1:iT) mBZ = rbind(mBZ, mB%*%mZ[iter,])
  }

  if(length(method) > 1) stop(simpleError("Don't choose multiple methods!"))
  opt = switch(match(method,c('max_1','max_2','max_3','min_1','min_2')),
               MaxFx1, MaxFx2, MaxFx3, MinLine, MinLineOpt)
  if(is.null(opt)) stop(simpleError("The method does not match from c('max_1','max_2','max_3','min_1','min_2')!"))

  mU = U0
  aAlpha[1,,] = mU
  for(iter in 1:iT){
    mH = -.5 * t(beta)%*%c(mX[iter,])%*%t(mX[iter,])%*%beta
    mC = mU%*%mD + mJ%*% c(mY[iter,]-mBZ[iter,])%*%t(mX[iter,])%*%beta
    tmp = svd(mC)
    tmp = tmp$u%*%t(tmp$v)

    # maximize and update mU
    mU = opt(mX=tmp,mA=mH,mB=mJ,mC=mC)$mX

    aAlpha[iter+1,,] = mU
  }
  return(aAlpha)
}


#' Filtering algorithm for the type two model.
#'
#' This function implements the filtering algorithm for the type two model.
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
#'
#' @param mY the matrix containing Y_t with dimension \eqn{T \times p}.
#' @param mX the matrix containing X_t with dimension \eqn{T \times q_1}.
#' @param mZ the matrix containing Z_t with dimension \eqn{T \times q_2}.
#' @param alpha the \eqn{\alpha} matrix.
#' @param mB the coefficient matrix \eqn{\boldsymbol{B}} before \code{mZ} with dimension \eqn{p \times q_2}.
#' @param Omega covariance matrix of the errors.
#' @param vD vector of the diagonals of \eqn{D}.
#' @param U0 initial value of the alpha sequence.
#' @param method a string representing the optimization method from c('max_1','max_2','max_3','min_1','min_2').
#'
#' @return an array \code{aAlpha} containing the modal orientations of alpha in the prediction step.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @keywords filtering
#'
#' @examples
#'
#' iT = 50
#' ip = 2
#' ir = 1
#' iqx = 4
#' iqz=0
#' ik = 0
#' Omega = diag(ip)*.1
#'
#' if(iqx==0) mX=NULL else mX = matrix(rnorm(iT*iqx),iT, iqx)
#' if(iqz==0) mZ=NULL else mZ = matrix(rnorm(iT*iqz),iT, iqz)
#' if(ik==0) mY=NULL else mY = matrix(0, ik, ip)
#'
#' alpha = matrix(c(runif_sm(num=1,ip=ip,ir=ir)), ip, ir)
#' beta_0 = matrix(c(runif_sm(num=1,ip=ip*ik+iqx,ir=ir)), ip*ik+iqx, ir)
#' mB=NULL
#' vD = 100
#'
#' ret = SimModel2(iT=iT, mX=mX, mZ=mZ, mY=mY, alpha=alpha, beta_0=beta_0, mB=mB, vD=vD)
#' mYY=as.matrix(ret$dData[,1:ip])
#' fil = FilterModel2(mY=mYY, mX=mX, mZ=mZ, alpha=alpha, mB=mB, Omega=Omega, vD=vD, U0=beta_0)
#'
#' @export
FilterModel2 <- function(mY, mX, mZ, alpha, mB=NULL, Omega, vD, U0, method='max_1')
{
  iT = nrow(mY)
  ip = ncol(mY)
  iqxx = ncol(mX)
  ir = ncol(alpha)

  aBeta = array(0,dim=c(iT+1,iqxx,ir))
  iOmg = chol2inv(chol(Omega))
  mH = -.5*t(alpha)%*%iOmg%*%alpha
  mU = U0
  aBeta[1,,] = mU
  mD = diag(vD,ir,ir)

  if(is.null(mB)){
    mBZ = matrix(0,iT,ip)
  }else{
    mBZ = NULL
    for(iter in 1:iT) mBZ = rbind(mBZ, mB%*%mZ[iter,])
  }

  if(length(method) > 1) stop(simpleError("Don't choose multiple methods!"))
  opt = switch(match(method,c('max_1','max_2','max_3','min_1','min_2')),
               MaxFx1, MaxFx2, MaxFx3, MinLine, MinLineOpt)
  if(is.null(opt)) stop(simpleError("The method does not match from c('max_1','max_2','max_3','min_1','min_2')!"))

  for(iter in 1:iT){
    mJ = mX[iter,] %*% t(mX[iter,])
    mC = mU%*%mD + mX[iter,]%*%t(mY[iter,]-mBZ[iter,])%*%iOmg%*%alpha

    # maximize and update mU
    mU = opt(mX=mU,mA=mH,mB=mJ,mC=mC)$mX

    aBeta[iter+1,,] = mU
  }
  return(aBeta)
}
