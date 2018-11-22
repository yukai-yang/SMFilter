#################################################################################
## Functions for optimization
## Nov 2018
#################################################################################

# QR decompositon
qrQ <- function(mX)
{
  tmp = qr(mX)
  return(qr.Q(tmp)%*%diag(sign(diag(qr.R(tmp))),dim(mX)[2]))
}


GetYY <- function(tau, mG, mX)
{
  mU = cbind(mG, mX)
  mV = cbind(mX, -mG)
  mW = mU%*%t(mV)
  vd = dim(mU)
  if(vd[1]>=vd[2]){
    tmp = solve(t(mV)%*%mU*tau*.5 + diag(1,vd[2]))
    mY = mX - tau*mU%*% tmp %*%t(mV)%*%mX
    gY = - ( diag(1,vd[1]) - tau*.5*mU%*%tmp%*%t(mV) ) %*% mW %*% ( mX + mY )*.5
  }else{
    tmp = solve(diag(1,vd[1])+mW*tau*.5)
    mY = tmp %*% (diag(1,vd[1])-mW*tau*.5)%*%mX
    gY = -tmp %*% mW %*% ( mX + mY )*.5
  }
  return(list(mY=mY, gY=gY))
}


# maximize the function
# tr{A X' B X + C' X} w.r.t. X
# max_1
MaxFx1 <- function(mX, mA, mB, mC, tol=1e-8, maxiter=50)
{
  numit = 1
  fn1 = sum(diag(mA%*%t(mX)%*%mB%*%mX+t(mC)%*%mX))

  fn = NULL; lX = list()
  while(numit < maxiter){
    tmp = svd(mA%*%t(mX)%*%mB+t(mC))
    mY = tmp$v%*%t(tmp$u)
    fn2 = sum(diag(mA%*%t(mY)%*%mB%*%mY+t(mC)%*%mY))

    if(abs(fn2-fn1)<tol) {fn = c(fn,fn1); lX[[length(lX)+1]] = mX; break}
    if(fn2<fn1) {fn = c(fn,fn1); lX[[length(lX)+1]] = mX}
    else{ fn = c(fn,fn2); lX[[length(lX)+1]] = mY }

    numit = numit+1
    mX = mY
    fn1 = fn2
  }
  mY = lX[[which.max(fn)]]

  return(list(mX=mY, val=max(fn), numit=numit))
}


# maximize the function
# tr{A X' B X + C' X} w.r.t. X
# max_2
MaxFx2 <- function(mX, mA, mB, mC, method='CG', tol=1e-10)
{
  ip = dim(mB)[1]
  ir = dim(mA)[1]
  ir2 = ir+ir

  fnFx <- function(mX){
    # compute tr{ A X' B X + C' X}
    return(sum(diag(mA%*%t(mX)%*%mB%*%mX + t(mC)%*%mX)))
  }

  if(ip>=(ir2)){
    futil1 <- function(tau,mU,mV){
      tmp = solve(t(mV)%*%mU*tau*.5 + diag(1,ir2,ir2))
      return(mX - tau*mU%*% tmp %*%t(mV)%*%mX)
    }
    futil2 <- function(tau,mU,mV){
      tmp = solve(t(mV)%*%mU*tau*.5 + diag(1,ir2,ir2))
      mY = mX - tau*mU%*% tmp %*%t(mV)%*%mX
      return(-( diag(1,ip) - tau*.5*mU%*%tmp%*%t(mV) ) %*% mU%*%t(mV) %*% ( mX + mY )*.5)
    }
  }else{
    futil1 <- function(tau,mU,mV){
      mW = mU%*%t(mV)
      tmp = solve(diag(1,ip)+mW*tau*.5)
      return(tmp %*% (diag(1,ip)-mW*tau*.5)%*%mX)
    }
    futil2 <- function(tau,mU,mV){
      mW = mU%*%t(mV)
      tmp = solve(diag(1,ip)+mW*tau*.5)
      mY = tmp %*% (diag(1,ip)-mW*tau*.5)%*%mX
      return(-tmp %*% mW %*% ( mX + mY )*.5)
    }
  }

  ftarget <- function(tau){
    return(fnFx(mX=futil1(tau=tau,mU=mU,mV=mV)))
  }

  fgr <- function(tau){
    return( sum(diag(t(mG)%*%futil2(tau=tau,mU=mU,mV=mV))) )
  }

  tfn = fnFx(mX)/ip
  niter = 1
  while(T){
    mG = mB%*%mX%*%mA*2 + mC
    mU = cbind(mG, mX)
    mV = cbind(mX, -mG)

    grad = mG - mX%*%t(mG)%*%mX
    #if(FNorm2(grad)<tol) break

    ret = optim(par=0, fn=ftarget, gr=fgr, method=method,control=list(fnscale=-1))
    mX = GetYY(tau=ret$par,mG=mG,mX=mX)$mY

    if(fnFx(mX)/ip - tfn < tol) break
    tfn = fnFx(mX)/ip

    niter = niter+1
  }

  return(list(mX=mX,grad=grad,numit=niter))
}


# maximize the function
# tr{A X' B X + C' X} w.r.t. X
# max_3
MaxFx3 <- function(mX, mA, mB, mC, method='Brent', tol=1e-10)
{
  ip = dim(mB)[1]
  ir = dim(mA)[1]

  fnFx <- function(mX){
    # compute tr{ A X' B X + C' X}
    return(sum(diag(mA%*%t(mX)%*%mB%*%mX + t(mC)%*%mX)))
  }

  ftarget <- function(tau){
    tmp = qrQ(mX+grad*tau)
    return(fnFx(mX=tmp))
  }

  tfn = fnFx(mX)/ip
  niter = 1
  while(T){
    mG = mB%*%mX%*%mA*2 + mC
    #grad = mG - mX%*%t(mG)%*%mX
    grad = mG - mX%*%(t(mG)%*%mX + t(mX)%*%mG)*.5

    ret = optim(par=0,fn=ftarget,method=method,
                lower=-10,upper=10,control=list(fnscale=-1))
    mX = qrQ(mX+grad*ret$par)

    if(fnFx(mX)/ip - tfn < tol) break
    tfn = fnFx(mX)/ip

    niter = niter+1
  }

  return(list(mX=mX,grad=grad,numit=niter))
}


# apply the Armijo step line search algorithm
# note that it is a minimization
# min tr{ AX'BX + C'X }
# min_2
MinLine <- function(mX, mA, mB, mC, tol=1e-8, maxiter=100, alpha=1, beta=.5, sigma=.1)
{
  mA = -mA; mC = -mC

  ip = dim(mB)[1]
  ir = dim(mA)[1]

  fnFx <- function(mX){
    # compute tr{ A X' B X + C' X}
    return(sum(diag(mA%*%t(mX)%*%mB%*%mX + t(mC)%*%mX))/ip)
  }

  ftarget <- function(tau){
    tmp = qrQ(mX-grad*tau); return(fnFx(tmp))
  }

  tfn = fnFx(mX)
  niter = 1

  while(T){
    if(niter == maxiter) break

    mG = (mB%*%mX%*%mA*2 + mC)/ip
    #grad = mG - mX%*%t(mG)%*%mX
    grad = mG - mX%*%(t(mG)%*%mX + t(mX)%*%mG)*.5

    tau = alpha
    counter = 0
    while( ftarget(tau) - tfn >= 0 && counter <20 ){
      tau = tau*beta; counter=counter+1
      #print(tau)
    }
    if(counter==20) tau = 0

    #if(ftarget(tau == tfn)) break

    tmp = mX-grad*tau
    mX = qrQ(tmp)

    if(abs(fnFx(mX) - tfn) < tol) break
    tfn = fnFx(mX)

    niter = niter+1
  }

  return(list(mX=mX,val=fnFx(mX),grad=grad,numit=niter))
}


# apply the Armijo step line search algorithm
# note that it is a minimization
# min tr{ AX'BX + C'X }
# min_3
MinLineOpt <- function(mX, mA, mB, mC, tol=1e-8, maxiter=100)
{
  mA = -mA; mC = -mC

  ip = dim(mB)[1]
  ir = dim(mA)[1]

  fnFx <- function(mX){
    # compute tr{ A X' B X + C' X}
    return(sum(diag(mA%*%t(mX)%*%mB%*%mX + t(mC)%*%mX))/ip)
  }

  ftarget <- function(tau) return(fnFx(mX=retr(tt=tau)))

  retr <- function(tt) return(qrQ(mX+grad*tt))

  tfn = fnFx(mX)
  niter = 1

  while(T){
    if(niter == maxiter) break

    mG = (mB%*%mX%*%mA*2 + mC)/ip
    #grad = mG - mX%*%t(mG)%*%mX
    grad = mG - mX%*%(t(mG)%*%mX + t(mX)%*%mG)*.5

    ret=optim(par=0,fn=ftarget,method="CG")

    mX = retr(tt=ret$par)

    if(abs(fnFx(mX) - tfn) < tol) break
    tfn = fnFx(mX)

    niter = niter+1
  }

  return(list(mX=mX,val=fnFx(mX)*ip,grad=grad*ip,numit=niter))
}

