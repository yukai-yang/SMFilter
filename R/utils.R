#################################################################################
## utility functions
#################################################################################

vnum = "1.0.0"
location = "GitHub"

# simple cat
cat0 <- function(...)
{
  words = list(...)
  for(tmp in words) cat(tmp)
  cat("\n")
}

#' Show the version number of some information.
#'
#' This function shows the version number and some information of the package.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @keywords utils
#'
#' @export
version <- function(){
  cat0("SMFilter version ", vnum, " (Red Filter) from ",location)
}

# exponential trace
etr <- function(mx)
{
  return(exp(sum(diag(mx))))
}

# compute X (X'X)^{-1/2}
morth <- function(mx)
{
  ir = dim(mx)[2]
  tmp = svd(t(mx)%*%mx)
  return(mx%*%(tmp$u%*%diag(1/sqrt(tmp$d),ir)%*%t(tmp$u)))
}

# former rUStiefel

#' Sample from the uniform distribution on the Stiefel manifold.
#'
#' This function draws a sample from the uniform distribution on the Stiefel manifold.
#'
#' The Stiefel manifold with dimension \eqn{p} and \eqn{r} (\eqn{p \geq r}) is a space whose points are \eqn{r}-frames in \eqn{R^p}.
#' A set of \eqn{r} orthonormal vectors in \eqn{R^p} is called an \eqn{r}-frame in \eqn{R^p}.
#' The Stiefel manifold is a collection of \eqn{p \times r} full rank matrices \eqn{X} such that \eqn{X'X = I_r}.
#'
#' @param num number of observations or sample size.
#' @param ip the first dimension \eqn{p} of the matrix.
#' @param ir the second dimension \eqn{r} of the matrix.
#'
#' @return an array with dimension \code{num}, \code{ip} and \code{ir} containing a sample of draws from the uniform distribution on the Stiefel manifold.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#' @keywords utils
#'
#' @examples
#'
#' runif_sm(10,4,2)
#'
#' @export
runif_sm <- function(num,ip,ir)
{
  mtmp = array(rnorm(ip*ir*num),dim=c(num,ip,ir))
  ret = t(apply(mtmp,1,morth))
  dim(ret) = c(num,ip,ir)
  return(ret)
}

# compute the integral: int etr(F'X) [dX]
# which is equal to 0F1(p/2, F'F/4)
# where mF p by r matrix
# mF = X0%*%diag(vD)
int0F1 <- function(vD, ip, ir, X0, aUS)
{
  mF = X0%*%diag(vD,ir)
  ftmp <- function(mX){
    return(etr(matrix(t(mF)%*%mX,ir,ir)))
  }
  tmp = apply(aUS,1,ftmp)
  return(mean(tmp))
}

# QR decomposition for mX
# then return the complement part of Q
cmpx <- function(mX)
{
  tmp = qr(mX)
  return(qr.Q(tmp, complete = TRUE)[,-(1L:tmp$rank)])
}

# compute the modal orientation given a sample aX from Stiefel manifold
mod_orientation <- function(aX,iN,ip,ir)
{
  ftmp <- function(nter) return(tcrossprod(aX[nter,,]))
    # return(aX[nter,,]%*%t(aX[nter,,]))

  tmp = sapply(1:iN,ftmp)
  dim(tmp) = c(ip,ip,iN)
  tmp = apply(tmp,c(1,2),mean)
  return(svd(tmp,ir)$u)
}

# target density is of the form
#   dbeta(t,1/2,k)*exp(a*t) * exp(b*sqrt(t))*(exp(c*sqrt(t))+exp(-c*sqrt(t)))
#  envelope density is of the form
#   dbeta(t,1/2,g)   where g<=k
rtheta_lb <- function(dk, da, db, dc)
{
  du=2; dg=dk; lrth=0; lrmx=0

  if(da>0){
    dg = max(1/(1+log(2+da)), dk-da)
    if((dk-dg)>0) lrmx = da-dk+dg + (dk-dg)*log((dk-dg)/da)
  }
  if(db<=0) lrmx = lrmx+dc+log(1+exp(-2*dc))
  if(db>0){
    ct = dc+log(.5*(1+exp(-2*dc)))
    lrmx = lrmx+sqrt(db*db+ct*ct) + log(2)
  }

  while(log(du)>lrth-lrmx){
    du = runif(1)
    dth = rbeta(1,.5,dg)
    lrth = da*dth+(dk-dg)*log(1-dth)+sqrt(1-dth)*db+sqrt(dth)*dc+log(1+exp(-2*sqrt(dth)*dc))
  }

  return(dth)
}

# steps (a) to (d) in Hoff
ry_lb <- function(vy, vl, vd)
{
  iN = length(vy)
  dk= .5*(iN-1)

  vy2 = vy**2
  for(iter in 1:iN){ # random order?
    omyi = 1/(1-vy2[iter])
    smyi = sqrt(omyi)
    ta = vl[iter] + (vl[iter]*vy2[iter] - sum(vl*vy2))*omyi
    tb = (sum(vy*vd) -vy[iter]*vd[iter])*smyi

    theta = rtheta_lb(dk,ta,tb,abs(vd[iter])) # rtheta_lb, a scalor

    vy = vy*sqrt(1-theta)*smyi
    vy2 = vy**2

    pp = 1/(1+exp(2*sqrt(theta)*vd[iter]))
    vy[iter] = sqrt(theta) * sample(c(-1,1),size=1,prob=c(pp,1-pp))
    vy2[iter] = theta
  }

  return(vy)
}

# former RVBL

#' Sample from the vector Langevin-Bingham on the Stiefel manifold.
#'
#' This function draws a sample from the vector Langevin-Bingham on the Stiefel manifold.
#'
#' The vector Langevin-Bingham distribution on the Stiefel manifold has the density kernel:
#' \deqn{f(X) \propto \mathrm{etr}\{ x' A x + c' x \}}
#' where \eqn{x} satisfies \eqn{x'x = 1}, and \eqn{A} is a symmetric matrix.
#'
#' @param num number of observations or sample size.
#' @param mA the matrix A which is symmetric ip*ip matrix.
#' @param vc the vector c with dimension ip.
#' @param vx the vector x, the initial value.
#'
#' @return an array containing a sample of draws from the vector Langevin-Bingham on the Stiefel manifold.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#'
#' @section References:
#' Hoff, P. D. (2009) "Simulation of the Matrix Bingham—von Mises—Fisher Distribution, With Applications to Multivariate and Relational Data", Journal of Computational and Graphical Statistics, Vol. 18, pp. 438-456.
#'
#' @keywords utils
#'
#' @export
rvlb_sm <- function(num, mA, vc, vx)
{
  tmp = eigen(mA,symmetric=T)
  mE = tmp$vectors
  vl = tmp$values
  #tmp = svd(mA)
  #mE = tmp$u
  #vl = tmp$d

  aret = array(0,dim=c(num,length(vx),1))

  for(iter in 1:num){
    vy = t(mE) %*% vx
    vd = t(mE) %*% vc
    vx = mE %*% ry_lb(vy,vl,vd) # keep this C interface function
    vx = vx/sqrt(sum(vx^2))
    aret[iter,,] = vx
  }

  return(aret)
}


# fomer RMBL

#' Sample from the matrix Langevin-Bingham on the Stiefel manifold.
#'
#' This function draws a sample from the matrix Langevin-Bingham on the Stiefel manifold.
#'
#' The matrix Langevin-Bingham distribution on the Stiefel manifold has the density kernel:
#' \deqn{f(X) \propto \mathrm{etr}\{ H X' J X + C' X \}}
#' where \eqn{X} satisfies \eqn{X'X = I_r}, and \eqn{H} and \eqn{J} are symmetric matrices.
#'
#' @param num number of observations or sample size.
#' @param mJ symmetric ip*ip matrix
#' @param mH symmetric ir*ir matrix
#' @param mC ip*ir matrix
#' @param mX ip*ir matrix, the initial value
#' @param ir ir
#'
#' @return an array containing a sample of draws from the matrix Langevin-Bingham on the Stiefel manifold.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#'
#' #' @section References:
#' Hoff, P. D. (2009) "Simulation of the Matrix Bingham—von Mises—Fisher Distribution, With Applications to Multivariate and Relational Data", Journal of Computational and Graphical Statistics, Vol. 18, pp. 438-456.
#'
#' @keywords utils
#'
#' @export
rmLB_sm <- function(num, mJ, mH, mC, mX, ir)
{
  ip = dim(mX)[1]
  # ip must be strictly bigger than ir

  aret = array(0,dim=c(num,ip,ir))

  if(ir==1){
    aret = rvlb_sm(num=num, mA=mJ*c(mH), vc=c(mC), vx=c(mX))
  }else{
    for(iter in 1:num){
      for (r in sample(1:ir)) {
        mN = cmpx(mX[,-r]) # null space
        An = mH[r,r] * t(mN) %*% (mJ) %*% mN
        cn = t(mN) %*% mC[,r]
        mX[,r] = mN %*% rvlb_sm(1, An, cn, t(mN)%*%mX[,r])
      }
      aret[iter,,] = mX
    }
  }

  return(aret)
}


# return the thetas/pi between two vectors
thetas <- function(mY, vX)
{
  tmp = c(mY %*% vX)
  tmax <- function(xx){
    if(xx>1) return(1)
    if(xx<(-1)) return(-1)
    return(xx)
  }
  tmp = sapply(tmp,tmax)
  return(acos(tmp)/pi)
}


FNorm2 <- function(mX)
{
  return(sum(diag(t(mX)%*%mX)))
}


#' Compute the squared Frobenius distance between two matrices.
#'
#' This function Compute the squared Frobenius distance between two matrices.
#'
#' The Frobenius distance between two matrices is defined to be
#' \deqn{d(X, Y) = \sqrt{ \mathrm{tr} \{ A' A \} }}
#' where \eqn{A = X - Y}.
#'
#' The Frobenius distance is a possible measure of the distance between two points on the Stiefel manifold.
#'
#' @param mX a \eqn{p \times r} matrix where \eqn{p \geq r}.
#' @param mY another \eqn{p \times r} matrix where \eqn{p \geq r}.
#'
#' @return the Frobenius distance.
#'
#' @author Yukai Yang, \email{yukai.yang@@statistik.uu.se}
#'
#' @keywords utils
#'
#' @examples
#'
#' FDist2(runif_sm(1,4,2)[1,,], runif_sm(1,4,2)[1,,])
#'
#' @export
FDist2 <- function(mX,mY)
{
  tmp = mX - mY
  return(FNorm2(tmp))
}
