#' @export

mean.transpose = function(x,tol = 1e-3)
{
  n = nrow(x)
  p = ncol(x)
  if(sum(is.na(x))==0)
    {
      nu = colMeans(x)
      xnew = scale(x,center=TRUE,scale=FALSE)
      mu = rowMeans(xnew)
      xnew = t(scale(t(xnew),center=TRUE,scale=FALSE))
    }
  else
    {
      ind = 1
      xnew = x
      mut = nut = 0
      while(abs(ind)>tol)
        {
          mu = rowMeans(xnew,na.rm=TRUE)
          xnew = xnew - mu
          nu = colMeans(xnew,na.rm=TRUE)
          xnew = t(t(xnew) - nu)
          ind = sum(mu) + sum(nu)
          mut = mut + mu
          nut = nut + nu
        }
      xnew[is.na(x)] = 0
      mu = mut
      nu = nut
    }
  M = t(t((matrix(0,n,p) + mu)) + nu)
  return(list(x=x,xcen=xnew,mu=mu,nu=nu,M=M))
}

cov.transpose = function(xc,rhor,rhoc,qr=2,qc=2,seed=1,thr=1e-3,maxit=1e3,trace=FALSE,thr.glasso=1e-3,maxit.glasso=1e3,pen.diag=TRUE)
{
  set.seed(seed)
  n = nrow(xc)
  p = ncol(xc)
  if(qr==2 & qc==2)
    {
      svdx = svd(xc,nu=n,LINPACK=TRUE)
      V = svdx$v
      d = svdx$d
      U = svdx$u
      r = length(which(d>1e-4))
      d = d[1:r]
      c1 = -4*rhoc*p^2
      c2 = 32*rhor*rhoc*p + d^4*(n-p)
      c3 = 4*rhor*(d^4 - 16*rhor*rhoc)
      beta = sqrt((-c2 - sqrt(c2^2 - 4*c1*c3))/(2*c1))
      theta = (d^2*beta)/(p*beta^2 - 4*rhor)
      if(sum(theta)>1e10)
        {
          ind = which(theta>1e6)
          r = min(ind) - 1
        }
      if(r<p){theta[(r+1):p] = rep(2*sqrt(rhoc/n),p-r) }
      beta[(r+1):n] = rep(2*sqrt(rhor/p),n-r)
      Sigmahat = U%*%diag(beta)%*%t(U)
      sigi = U%*%diag(1/beta)%*%t(U)
      Deltahat = V%*%as.matrix(diag(theta))%*%t(V)
      delti = V%*%diag(1/theta)%*%t(V)
    }
  if(qr==1 & qc==1)
    {
      delti = glasso(t(xc)%*%xc/n,rhoc,maxit=maxit.glasso,thr=thr.glasso,penalize.diag=pen.diag)$wi
      sigi = glasso(xc%*%t(xc)/p,rhor,maxit=maxit.glasso,thr=thr.glasso,penalize.diag=pen.diag)$wi
      ind = 1; iter = 0;
      sr = xc%*%delti%*%t(xc)/n
      while(ind>thr & iter<maxit)
        {
          iter = iter + 1
          oldS = sigi
          oldD = delti
          gr = glasso(sr,2*rhor/p,maxit=maxit.glasso,thr=thr.glasso,penalize.diag=pen.diag)
          sigi = gr$wi
          sc = t(xc)%*%sigi%*%xc/n
          gc = glasso(sc,2*rhoc/n,maxit=maxit.glasso,thr=thr.glasso,penalize.diag=pen.diag)
          delti = gc$wi
          sr = xc%*%delti%*%t(xc)/p
          if(trace)
            {
              Sigmahat = gr$w
              Deltahat = gc$w
              cat(MNloglike(xc,M=matrix(0,n,p),Sig=Sigmahat,Delt=Deltahat,rhor=rhor,rhoc=rhoc,qr=qr,qc=qc),fill=TRUE)
            }
          ind = (sum((oldD - delti)^2) + sum((oldS - sigi)^2))/(sum((oldS)^2) + sum((oldD)^2))
        }
      if(!trace)
        {
          Sigmahat = gr$w
          Deltahat = gc$w
        }
    }
  if(qr==1 & qc==2)
    {
      svdx = svd(xc,LINPACK=TRUE)
      V = svdx$v
      theta = runif(p)
      ind = 1
      iter = 0
      solve.quad=function(a,b,cc){ return((-b+sqrt(b^2-4*a*cc))/(2*a))}
      delti = V%*%diag(1/theta)%*%t(V)
      sigi = matrix(0,n,n)
      while(ind>thr & iter<maxit)
        {
          oldD = delti
          oldS = sigi
          sr = xc%*%delti%*%t(xc)/p
          gr = glasso(sr,2*rhor/p,maxit=maxit.glasso,thr=thr.glasso,penalize.diag=pen.diag)
          sigi = gr$wi
          sc = t(xc)%*%sigi%*%xc
          esc = eigen(sc)
          theta = (esc$values + sqrt(esc$values^2 + 16*n*rhoc))/(2*n)
          delti = esc$vectors%*%diag(1/theta)%*%t(esc$vectors)
          if(trace)
            {
              Sigmahat = gr$w
              Deltahat = esc$vectors%*%diag(theta)%*%t(esc$vectors)
              cat(MNloglike(xc,M=matrix(0,n,p),Sig=Sigmahat,Delt=Deltahat,rhor=rhor,rhoc=rhoc,qr=qr,qc=qc),fill=TRUE)
            }
          ind = (sum((oldD - delti)^2) + sum((oldS - sigi)^2))/(sum((oldS)^2) + sum((oldD)^2))          
        }
      if(!trace)
        {
          Sigmahat = gr$w
          Deltahat = esc$vectors%*%diag(theta)%*%t(esc$vectors)
        }
    }
  if(qr==2 & qc==1)
    {
      svdx = svd(xc,nu=n,LINPACK=TRUE)
      U = svdx$u
      beta = runif(n)
      ind = 1
      iter = 0
      solve.quad=function(a,b,cc){ return((-b+sqrt(b^2-4*a*cc))/(2*a))}
      sigi = U%*%diag(1/beta)%*%t(U)
      delti = matrix(0,p,p)
      while(ind>thr & iter<maxit)
        {
          oldD = delti
          oldS = sigi
          sc = t(xc)%*%sigi%*%xc/n
          gc = glasso(sc,2*rhoc/n,maxit=maxit.glasso,thr=thr.glasso)
          delti = gc$wi
          sr = xc%*%delti%*%t(xc)
          svr = svd(sr,nu=n,LINPACK=TRUE)
          vals = svr$d
          U = svr$u
          beta = (vals + sqrt(vals^2 + 16*p*rhor))/(2*p)
          sigi = U%*%diag(1/beta)%*%t(U)
          if(trace)
            {
              Deltahat = gc$w
              Sigmahat = U%*%diag(beta)%*%t(U)
              cat(MNloglike(xc,M=matrix(0,n,p),Sig=Sigmahat,Delt=Deltahat,rhor=rhor,rhoc=rhoc,qr=qr,qc=qc),fill=TRUE)
            }
          ind = (sum((oldD - delti)^2) + sum((oldS - sigi)^2))/(sum((oldS)^2) + sum((oldD)^2))          
        }
      if(!trace)
        {
          Deltahat = gc$w
          Sigmahat = U%*%diag(beta)%*%t(U)
        }

    }
  return(list(Sigmahat=Sigmahat,Deltahat=Deltahat,Sigmaihat=sigi,Deltaihat=delti))
}


MNloglike = function(x,M,Sig,Delt,rhor,rhoc,qr=2,qc=2,Sigi=NULL,Delti=NULL)
{
  n = nrow(x)
  p = ncol(x)
  if(length(Sigi)==0){  Sigi = solve(Sig) }
  if(length(Delti)==0){   Delti = solve(Delt) }
  if(qr==2){ tr = sum((rhor*Sigi)^2)} else { tr = sum(abs(rhor*Sigi))}
  if(qc==2){ tc = sum((rhoc*Delti)^2)} else { tc = sum(abs(rhoc*Delti))}
  if(det(Sigi)==0){ t1 = log(1e-300)
                  }else{ t1 = log(det(Sigi))}
  if(sum(is.na(x))==0)
    {
      val = (p/2)*t1 + (n/2)*log(det(Delti)) - (1/2)*sum(diag(Sigi%*%(x-M)%*%Delti%*%t(x-M))) - tr - tc
    }
  else
    {
      xs = x
      xc = x - M
      tlds = 0
      for(j in 1:p)
        {
          oj = !is.na(x[,j])
          sigioj = solve(Sig[oj,oj])
          tlds = tlds + log(det(sigioj))
          xs[oj,j] = chol(sigioj)%*%xc[oj,j]
        }
      tldd = trt = 0
      for(i in 1:n)
        {
          oi = !is.na(x[i,])
          deltioi = solve(Delt[oi,oi])
          tldd = tldd + log(det(deltioi))
          trt = trt + sum(diag(xs[i,oi]%*%t(xs[i,oi])%*%deltioi))
        }
      val = (1/2)*(tlds + tldd - trt) - tr - tc
    }
  return(val)
}


MVNloglike = function(x,mu,Sig,rho,q)
{
  n = nrow(x)
  Sigi = ginv(Sig)
  if(q==1){ tl = rho*sum(abs(Sigi)) }
  else { tl = rho*sum(Sigi^2) }
  if(sum(is.na(x))==0)
    {
      xc = t(t(x) - mu)
      val = (n/2)*log(det(Sigi)) - (1/2)*sum(diag(xc%*%Sigi%*%t(xc))) - tl
    }
  else
    {
      tld = tm = 0
      for(i in 1:n)
        {
          oi = !is.na(x[i,])
          Sigioi = ginv(Sig[oi,oi])
          tld = tld + log(det(Sigioi))
	  if (!is.null(x[i,oi]-mu[oi]))
	  {
		if(ncol(as.matrix(x[i,oi] - mu[oi]))==nrow(Sigioi))
		{
          		temp = as.matrix((x[i,oi] - mu[oi]))%*%Sigioi
          		tm = tm + temp%*%as.matrix(t(x[i,oi] - mu[oi]))
		}
		else
		{
			temp = as.matrix(t(x[i,oi] - mu[oi]))%*%Sigioi
          		tm = tm + temp%*%as.matrix((x[i,oi] - mu[oi]))
		}
	  }
        }
      val = (1/2)*(tld -tm) - tl
    }
  return(val)
}

TREC = function(x,sig,delt,sigi,delti,M,thr=1e-3,maxit=1e3)
{
  n = nrow(x)
  p = ncol(x)
  xhat = x
  xhat[is.na(x)] = M[is.na(x)]
  ind = 1
  iter = 1
  rmi = (1:n)[apply(is.na(x),1,sum)>0]
  cmj = (1:p)[apply(is.na(x),2,sum)>0]
  nrmi = apply(is.na(x),1,sum)
  ncmj = apply(is.na(x),2,sum)
  while(ind>thr & iter<maxit)
    {
      oldx = xhat
      for(i in rmi)
        {
          ey = M[i,] + (-sigi[i,-i]/sigi[i,i])%*%as.matrix((xhat[-i,] - M[-i,]))
          covy = delt/sigi[i,i]
          mj = (1:p)[is.na(x[i,])]
          if(nrmi[i]>=round(p/2))
            {
              swpz = matinv(covy,(1:p)[-mj])
              xhat[i,mj] = ey[mj] + swpz[mj,-mj,drop=FALSE]*(xhat[i,-mj] - ey[-mj])
			  xhat[i,mj] = round(xhat[i,mj],2)
            }
          else
            {
              swpz = matinv(delti,(1:p)[mj])
              xhat[i,mj] = ey[mj] - swpz[mj,-mj,drop=FALSE]*(xhat[i,-mj] - ey[-mj])
            }
        }
      for(j in cmj)
        {
          ey = M[,j] + as.matrix((xhat[,-j] - M[,-j]))%*%(-delti[-j,j]/delti[j,j])
          covy = sig/delti[j,j]
          mi = (1:n)[is.na(x[,j])]
          if(ncmj[j]>=round(n/2))
            {
              swpz = matinv(covy,(1:n)[-mi])
              xhat[mi,j] = ey[mi] + swpz[mi,-mi,drop=FALSE]*(xhat[-mi,j] - ey[-mi])
			  xhat[i,mj] = round(xhat[i,mj],2)
            }
          else
            {
              swpz = matinv(sigi*delti[j,j],(1:n)[mi])
              xhat[mi,j] = ey[mi] - swpz[mi,-mi,drop=FALSE]*(xhat[-mi,j] - ey[-mi])
            }
        }
      ind = sum((oldx - xhat)^2)/sum(oldx^2)
      iter = iter + 1
    }
  return(list(xhat=xhat,iter=iter))
}

REC = function(x,trace=TRUE,maxit=10,thr.em=1e-3)
{
  n = nrow(x)
  p = ncol(x)
  missi = (1:n)[apply(is.na(x),1,sum)>0]
  if(n>1){muhat = colMeans(x,na.rm=TRUE)} else{muhat = rep(0,p)}
  xc = t(t(x) - muhat)
  xc[is.na(x)] = 0
  Sighat = t(xc)%*%xc/n
  tlike = 0
  iter = 0
  if(trace)
    {
      tlike[iter + 1] = MVNloglike(x,muhat,Sighat,0,2)
      cat(tlike[iter+1],fill=T)
    }
  ind = 1
  xhat = t(t(xc) + muhat)
  while(ind>thr.em & iter<maxit)
    {
      oldxh = xhat
      iter = iter + 1
      covc = matrix(0,p,p)
      for(i in missi)
        {
          mi = is.na(x[i,])
          swp = matinv(Sighat,(1:p)[!mi])
          xhat[i,mi] = muhat[mi] + swp[mi,!mi,drop=FALSE]%*%t(xhat[i,!mi,drop=FALSE] - muhat[!mi])
          covc[mi,mi] = covc[mi,mi] + Sighat[mi,mi] - swp[mi,!mi,drop=FALSE]%*%Sighat[!mi,mi,drop=FALSE]
        }
      if(n>1){muhat = colMeans(xhat)}
      xhc = t(t(xhat) - muhat)
      Sighat = t(xhc)%*%xhc/n + covc/n
      if(trace)
        {
          xh = xhat
          xh[is.na(x)] = NA
          tlike[iter + 1] = MVNloglike(xh,muhat,Sighat,0,2)
          cat(tlike[iter+1],fill=T)
        }
      ind = sum((xhat - oldxh)^2)/sum(x^2,na.rm=TRUE)
    }
  if(!trace){tlike=NULL}
  return(list(xhat=xhat,Sig=Sighat,mu=muhat,loglike=tlike))

}


