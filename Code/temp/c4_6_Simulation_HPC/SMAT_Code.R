SMAT <-
function(y,x,s,w=NULL,status=NULL,prev=NULL,working="unstr",conv=1e-4,maxiter=5e2) {
  d=w # weighting
  w=s # SNP effect
  
	if(!is.matrix(y)) { stop("y needs to be a matrix.") }
	if(!is.matrix(x)) { stop("x needs to be a matrix.") }

  cor.y = cor(y,use="pairwise.complete.obs")
  if(sum(cor.y<=0)>0) {
    warning("At least one pair of the phenotypes is not positively correlated; consider transformation.  As an alternative, consider M-DF test implemented within library geepack.")
  }

  if(is.vector(w)) w=as.matrix(w)
	#if(!is.matrix(w)) { stop("s needs to be a matrix.") }

	# check input types and dimensions
	n = nrow(y)
	M = ncol(y)
	p = ncol(x)
	q = 1 # set to 1 for the moment
	if(nrow(x)!=n | nrow(w)!=n) { stop("Row numbers of y, x, s do not match.") }
	
	if(!is.null(d)) {
		if(length(d)!=n) { stop("Row numbers of w and y do not match.") }
		is.wt = "user"; 
	} else if(!is.null(status)){ 
                status=as.vector(status)
                if(length(status)!=n) { stop("Row numbers of status and y do not match.") }
                if(is.null(prev)) { stop("prev should be provided if status is provided.") }
                if(length(setdiff(status,c(NA,1,0)))>0) { stop("status should only contain 0, 1, or NA.") }
		d = (status==1)*prev[1] + (status==0)*(1-prev[1]);  # prev should be a scalar
                is.wt = "casecontrol";
	} else {
            if(!is.null(prev)) {
                warning("Equal weight is used as status is not provided.")
            }
            d = rep(1,n)
            is.wt = "equal"
        } 

  # use only non-missing data
  index.observed=which( rowSums(is.na(y))==0 & 
                        rowSums(is.na(x))==0 & 
                        is.na(d)==FALSE & 
                        rowSums(is.na(w))==FALSE )
  if(length(index.observed)==0) stop("No complete observations are found!")
  y=y[index.observed,,drop=F]
  x=x[index.observed,,drop=F]
  d=d[index.observed]
  w=as.matrix(w[index.observed,,drop=F])

  if(is.wt=="casecontrol") {
    status=status[index.observed]
    pn=mean(status)
    d = d / ( (status==1)*pn + (status==0)*(1-pn) )
  }

  # use the new dimensions
  n = nrow(y)
  M = ncol(y)
  p = ncol(x)
    
  sum.d = sum(d) # udpate

  # convergence indicator
  is.conv = TRUE

	# precomputation
	if(is.wt!="equal") {
		y = diag(sqrt(d))%*%y
		x = diag(sqrt(d))%*%x
		w = diag(sqrt(d))%*%w	
	}
	yy=crossprod(y,y)
	yx=crossprod(y,x)
	yw=crossprod(y,w)
	xx=crossprod(x,x)
	xw=crossprod(x,w)
	xy=t(yx)
	ww=crossprod(w,w)
	wx=t(xw)
	wy=t(yw)

	# initial values of parameters
	gamhat = rep(0,p*M+q)
	sigma2 = apply(y,2,sd)^2 # sample variance
	psi = diag(sigma2)
	eyeM = diag(M)
	lM = matrix(1,ncol=1,nrow=M)
	if(working=="ind") { 
		R=eyeM 
	} else if(working=="exch") { 
		theta=mean(lowerTriangle(cor(y)))
		R=matrix(theta,M,M)
		diag(R)=1 
	} else { 
		R=cor(y) 
	} # unstructured, sample correlation

	lower = rep(1,M)

	# big while loop for estimation ---------------------------
	max.diff = 1
	iter = 1
	while(max.diff>conv & iter<maxiter) {

		invR = solve(R)

		#------------------
		#  update gamhat 
		#------------------
		gamhat.old = gamhat

		phi = diag(1/sqrt(c(sigma2))) # not psi
		invR.phi = invR %*% phi
		l.invR.phi = t(lM) %*% invR.phi
		l.invR = t(lM) %*% invR
		invR.l = t(l.invR)
		l.invR.l = l.invR %*% lM
		XinvRX = rbind(cbind(kronecker(invR,xx),kronecker(invR.l,xw)),
                       cbind(kronecker(l.invR,wx),kronecker(ww,l.invR.l)))
		XinvRy = rbind(matrix(t(invR.phi%*%yx),ncol=1),
                       matrix(l.invR.phi%*%yw,ncol=1))

		gamhat=solve(XinvRX,XinvRy)
		beta = matrix(gamhat[1:(M*p)],nrow=p)
		alpha = matrix(gamhat[(M*p+1):(M*p+q)],nrow=q)

		#------------------
		#  update sigma2
		#------------------
		sigma2.old = sigma2

		for(j in 1:M) {
 			old2=sigma2[j]
			incre=old2
			while(abs(incre)/old2>=conv) {
				old2  = sigma2[j]
				old   = sqrt(sigma2[j])
				bj    = beta[,j]
				upper = yy[j,j]-old*yx[j,]%*%bj-old*yw[j]%*%alpha-sum.d*old2
				lower[j] = sum.d+t(bj)%*%xx%*%bj/2+t(alpha)%*%wx%*%bj+t(alpha)%*%ww%*%alpha/2
				incre = upper/lower[j]
				sigma2[j] = old2 + incre
			}
		}
		
		#------------------
		#  update R
		#------------------
		if(working=="exch" | working=="unstr") {
			phi = diag(1/sqrt(c(sigma2))) # not psi, updated by new sigma2
  			a = alpha%*%t(lM)
			rr = phi%*%yy%*%phi - phi%*%yx%*%beta - t(beta)%*%t(yx)%*%phi - phi%*%yw%*%a - t(a)%*%t(yw)%*%phi + t(beta)%*%xx%*%beta + t(a)%*%ww%*%a + t(beta)%*%xw%*%a + t(a)%*%wx%*%beta
			sum.rr = sum(lowerTriangle(rr))
			if(working=="exch") {
				theta = sum.rr / (n*M*(M-1)/2 - p*M - q - M) # double check this
				R = matrix(theta,M,M)
				diag(R) = 1
			} else {
				rr = rr/n
				R = rr
				diag(R) = 1
			}
		}
		
		# compute differences
		max.diff = max(abs(gamhat-gamhat.old)/abs(gamhat.old),abs(sigma2-sigma2.old)/sigma2.old)
		iter = iter+1
		
	} # end of big while loop
	# big while loop for estimation completed -----------------
  
  # check if maximum #iter is reched
  if(iter>=maxiter) { is.conv = FALSE; warning("Maximum number of iterations reached. Please increase maxiter and rerun the function.") }
  
	# invR
	invR = solve(R)
	
	# y.star
	y.star = y%*%diag(1/sqrt(sigma2))
	# residual matrix at estimated values (nxM)
	resid = y.star - x%*%beta - w%*%alpha%*%t(lM)
	
	# score * t(score) in summation
	U.sig = (y.star * resid - d %*% t(lM))/n # (nxM)  # NOTE THIS 1/n
	U.beta = matrix(0,nrow=n,ncol=p*M)
	U.alpha = matrix(0,nrow=n,ncol=q)
	resid.invR = resid %*% invR # (nxM)
	UU = matrix(0,nrow=p*M+q+M,ncol=p*M+q+M)
  U = matrix(0,nrow=p*M+q+M,ncol=n) # useful in testing for common exposure
	for(i in 1:n) {
		Ui.sig = U.sig[i,]
		Ui.beta = as.vector( x[i,] %*% t(resid.invR[i,]) )
		U.beta[i,] = Ui.beta
		Ui.alpha = w[i,] %*% t(resid.invR[i,]) %*% lM
		U.alpha[i,] = Ui.alpha
		Ui = c(Ui.sig, Ui.beta, Ui.alpha)
		UU = UU + tcrossprod(Ui)
    U[,i] = Ui
	}
	
	# hessian matrix on 1-1 corner
	H.sig.sig = diag(lower/n/sigma2)
  # hessian matrix on 2-2 corner
	H.gam.gam = XinvRX
	# hessian matrix on 1-2 corner
	H.sig.gam.left = matrix(0,nrow=M,ncol=p) # useful matrix
	H.sig.gam.right = matrix(0,nrow=M,ncol=q) # useful matrix
	for(i in 1:M) {
		H.sig.gam.left[i,] = t(beta[,i])%*%xx + t(alpha)%*%wx
		H.sig.gam.right[i,] = t(beta[,i])%*%xw + t(alpha)%*%ww 
	}
  H.sig.gam.left = vdiag(H.sig.gam.left)
	H.sig.gam = cbind(H.sig.gam.left,H.sig.gam.right) / n
  # hessian matrix on 2-1 corner
  H.gam.sig.upper = t(H.sig.gam.left) %*% diag(1/sigma2) %*% invR / 2 # pM x M
  H.gam.sig.lower = t(vdiag(H.sig.gam.right)) %*% diag(1/sigma2) %*% invR %*% lM / 2 # Mq x 1
  H.gam.sig = matrix(0,nrow=p*M+q,ncol=M)
  for(i in 1:M) {
    H.gam.sig[1:(p*M),i] = matrix(H.gam.sig.upper[(1:p)+(i-1)*p,],nrow=p*M)
    H.gam.sig[(p*M+1):(p*M+q),i] = matrix(H.gam.sig.lower[(1:q)+(i-1)*q,],nrow=q)
  }
	# hessian matrix
	H = rbind(cbind(H.sig.sig,H.sig.gam),cbind(H.gam.sig,H.gam.gam))
	
	# information matrix
	info = t(H)%*%solve(UU)%*%H
	
	# covariance matrix
	cov.mat = solve(info)
	
	# s.e.
	se.vec = sqrt(diag(cov.mat))
	se.sig = se.vec[1:M]
	se.gam = se.vec[(M+1):(p*M+q+M)]
  se.beta = se.gam[1:(p*M)]
  se.alpha = se.gam[(p*M+1):(p*M+q)]
  
  # p-value
  pval.alpha = 2*(1-pnorm(abs(alpha/se.alpha)))

  # score test for common exposure
  Lambda = diag(M)
  Lambda = Lambda[-1,]
  Lam.invR = Lambda%*%invR
  U2 = matrix(0,nrow=M-1,ncol=n)
  UU2 = matrix(0,nrow=M+p*M+q,ncol=M-1) # G12, Note that G11=UU
  #U2U = t(UU2) # G21
  U2U2 = matrix(0,nrow=M-1,ncol=M-1) # G22
  for(i in 1:n) {
    U2[,i] = c(w[i,]) * Lam.invR %*% resid[i,]
    UU2 = UU2 + tcrossprod(U[,i],U2[,i])
    #U2U = U2U + tcrossprod(U2[,i],U[,i])
    U2U2 = U2U2 + tcrossprod(U2[,i])
  }
  U2U = t(UU2)
  sum.U2 = rowSums(U2)
  # Note: big.U = rbind(U,U2) (not in summation), dimension = (p*M+q+M+M-1) x n

  # A matrix 
  A.left = matrix(0,nrow=M-1,ncol=M)
  A.right = matrix(0,nrow=M-1,ncol=p*M+q)
  for(i in 1:n) {
    A.left = A.left + 0.5 * c(w[i,]) * Lam.invR %*% diag(1/sigma2) %*% diag(c( crossprod(x[i,],beta) + c(w[i,])*c(alpha) ))
    A.right = A.right + c(w[i,]) * Lam.invR %*% cbind(kronecker(eyeM,t(x[i,])), c(w[i,])*lM)
  }
  A = cbind(A.left,A.right) # Note: same as the original SAS

  # SIGMA matrix
  AinvH = A%*%solve(H)
  tAinvH = t(AinvH)
  SIG = U2U2 - AinvH%*%UU2 - U2U%*%tAinvH + AinvH%*%UU%*%tAinvH

  # score statistic & p-value
  score.test = t(sum.U2) %*% solve(SIG) %*% sum.U2
  pval.common = 1-pchisq(score.test,M-1)

	# output a list
  out = list(beta=beta,se.beta=matrix(se.beta,nrow=p),
             alpha=c(alpha),se.alpha=c(se.alpha),
             sigma2=sigma2,se.sigma2=se.sig,
             pvalue=c(pval.alpha),pvalue.common=c(pval.common),is.conv=is.conv)
  
  if(!is.null(colnames(x))) {
    rownames(out$beta)=rownames(out$se.beta)=colnames(x) # dim p
  } else {
    rownames(out$beta)=rownames(out$se.beta)=paste("predictor",1:p,sep="") # dim p
  }
  
  if(!is.null(colnames(y))) {
    colnames(out$beta)=colnames(out$se.beta)=colnames(y) # dim M
  } else {
    colnames(out$beta)=colnames(out$se.beta)=paste("outcome",1:M,sep="") # dim M
  }
  
  return(out)
}

lowerTriangle <-
function(x, diag=F) {
	x[lower.tri(x,diag=diag)]
}

vdiag <-
function(v) {
	nr = nrow(v)
	nc = ncol(v)
	mat = matrix(0,nrow=nr,ncol=nc*nr)
	for(i in 1:nrow(v)) {
		mat[i,(1:nc)+(i-1)*nc] = v[i,]
	}
	return(mat)
}

SMAT_temp <-function(y,x,s,w=NULL,status=NULL,prev=NULL,working="unstr",conv=1e-4,maxiter=5e2) {
  d=w # weighting
  w=s # SNP effect
  
	if(!is.matrix(y)) { stop("y needs to be a matrix.") }
	if(!is.matrix(x)) { stop("x needs to be a matrix.") }

  cor.y = cor(y,use="pairwise.complete.obs")
  if(sum(cor.y<=0)>0) {
    warning("At least one pair of the phenotypes is not positively correlated; consider transformation.  As an alternative, consider M-DF test implemented within library geepack.")
  }

  if(is.vector(w)) w=as.matrix(w)
	#if(!is.matrix(w)) { stop("s needs to be a matrix.") }

	# check input types and dimensions
	n = nrow(y)
	M = ncol(y)
	p = ncol(x)
	q = 1 # set to 1 for the moment
	if(nrow(x)!=n | nrow(w)!=n) { stop("Row numbers of y, x, s do not match.") }
	
	if(!is.null(d)) {
		if(length(d)!=n) { stop("Row numbers of w and y do not match.") }
		is.wt = "user"; 
	} else if(!is.null(status)){ 
                status=as.vector(status)
                if(length(status)!=n) { stop("Row numbers of status and y do not match.") }
                if(is.null(prev)) { stop("prev should be provided if status is provided.") }
                if(length(setdiff(status,c(NA,1,0)))>0) { stop("status should only contain 0, 1, or NA.") }
		d = (status==1)*prev[1] + (status==0)*(1-prev[1]);  # prev should be a scalar
                is.wt = "casecontrol";
	} else {
            if(!is.null(prev)) {
                warning("Equal weight is used as status is not provided.")
            }
            d = rep(1,n)
            is.wt = "equal"
        } 

  # use only non-missing data
  index.observed=which( rowSums(is.na(y))==0 & 
                        rowSums(is.na(x))==0 & 
                        is.na(d)==FALSE & 
                        rowSums(is.na(w))==FALSE )
  if(length(index.observed)==0) stop("No complete observations are found!")
  y=y[index.observed,,drop=F]
  x=x[index.observed,,drop=F]
  d=d[index.observed]
  w=as.matrix(w[index.observed,,drop=F])

  if(is.wt=="casecontrol") {
    status=status[index.observed]
    pn=mean(status)
    d = d / ( (status==1)*pn + (status==0)*(1-pn) )
  }

  # use the new dimensions
  n = nrow(y)
  M = ncol(y)
  p = ncol(x)
    
  sum.d = sum(d) # udpate

  # convergence indicator
  is.conv = TRUE

	# precomputation
	if(is.wt!="equal") {
		y = diag(sqrt(d))%*%y
		x = diag(sqrt(d))%*%x
		w = diag(sqrt(d))%*%w	
	}
	yy=crossprod(y,y)
	yx=crossprod(y,x)
	yw=crossprod(y,w)
	xx=crossprod(x,x)
	xw=crossprod(x,w)
	xy=t(yx)
	ww=crossprod(w,w)
	wx=t(xw)
	wy=t(yw)

	# initial values of parameters
	gamhat = rep(0,p*M+q)
	sigma2 = apply(y,2,sd)^2 # sample variance
	psi = diag(sigma2)
	eyeM = diag(M)
	lM = matrix(1,ncol=1,nrow=M)
	if(working=="ind") { 
		R=eyeM 
	} else if(working=="exch") { 
		theta=mean(lowerTriangle(cor(y)))
		R=matrix(theta,M,M)
		diag(R)=1 
	} else { 
		R=cor(y) 
	} # unstructured, sample correlation

	lower = rep(1,M)

	# big while loop for estimation ---------------------------
	max.diff = 1
	iter = 1
	while(max.diff>conv & iter<maxiter) {

		invR = solve(R)

		#------------------
		#  update gamhat 
		#------------------
		gamhat.old = gamhat

		phi = diag(1/sqrt(c(sigma2))) # not psi
		invR.phi = invR %*% phi
		l.invR.phi = t(lM) %*% invR.phi
		l.invR = t(lM) %*% invR
		invR.l = t(l.invR)
		l.invR.l = l.invR %*% lM
		XinvRX = rbind(cbind(kronecker(invR,xx),kronecker(invR.l,xw)),
                       cbind(kronecker(l.invR,wx),kronecker(ww,l.invR.l)))
		XinvRy = rbind(matrix(t(invR.phi%*%yx),ncol=1),
                       matrix(l.invR.phi%*%yw,ncol=1))

		gamhat=solve(XinvRX,XinvRy)
		beta = matrix(gamhat[1:(M*p)],nrow=p)
		alpha = matrix(gamhat[(M*p+1):(M*p+q)],nrow=q)

		#------------------
		#  update sigma2
		#------------------
		sigma2.old = sigma2

		for(j in 1:M) {
 			old2=sigma2[j]
			incre=old2
			while(abs(incre)/old2>=conv) {
				old2  = sigma2[j]
				old   = sqrt(sigma2[j])
				bj    = beta[,j]
				upper = yy[j,j]-old*yx[j,]%*%bj-old*yw[j]%*%alpha-sum.d*old2
				lower[j] = sum.d+t(bj)%*%xx%*%bj/2+t(alpha)%*%wx%*%bj+t(alpha)%*%ww%*%alpha/2
				incre = upper/lower[j]
				sigma2[j] = old2 + incre
			}
		}
		
		#------------------
		#  update R
		#------------------
		if(working=="exch" | working=="unstr") {
			phi = diag(1/sqrt(c(sigma2))) # not psi, updated by new sigma2
  			a = alpha%*%t(lM)
			rr = phi%*%yy%*%phi - phi%*%yx%*%beta - t(beta)%*%t(yx)%*%phi - phi%*%yw%*%a - t(a)%*%t(yw)%*%phi + t(beta)%*%xx%*%beta + t(a)%*%ww%*%a + t(beta)%*%xw%*%a + t(a)%*%wx%*%beta
			sum.rr = sum(lowerTriangle(rr))
			if(working=="exch") {
				theta = sum.rr / (n*M*(M-1)/2 - p*M - q - M) # double check this
				R = matrix(theta,M,M)
				diag(R) = 1
			} else {
				rr = rr/n
				R = rr
				diag(R) = 1
			}
		}
		
		# compute differences
		max.diff = max(abs(gamhat-gamhat.old)/abs(gamhat.old),abs(sigma2-sigma2.old)/sigma2.old)
		iter = iter+1
		
	} # end of big while loop
	# big while loop for estimation completed -----------------
  
  # check if maximum #iter is reched
  if(iter>=maxiter) { is.conv = FALSE; warning("Maximum number of iterations reached. Please increase maxiter and rerun the function.") }
  
	# invR
	invR = solve(R)
	
	# y.star
	y.star = y%*%diag(1/sqrt(sigma2))
	# residual matrix at estimated values (nxM)
	resid = y.star - x%*%beta - w%*%alpha%*%t(lM)
	
	# score * t(score) in summation
	U.sig = (y.star * resid - d %*% t(lM))/n # (nxM)  # NOTE THIS 1/n
	U.beta = matrix(0,nrow=n,ncol=p*M)
	U.alpha = matrix(0,nrow=n,ncol=q)
	resid.invR = resid %*% invR # (nxM)
	UU = matrix(0,nrow=p*M+q+M,ncol=p*M+q+M)
  U = matrix(0,nrow=p*M+q+M,ncol=n) # useful in testing for common exposure
	for(i in 1:n) {
		Ui.sig = U.sig[i,]
		Ui.beta = as.vector( x[i,] %*% t(resid.invR[i,]) )
		U.beta[i,] = Ui.beta
		Ui.alpha = w[i,] %*% t(resid.invR[i,]) %*% lM
		U.alpha[i,] = Ui.alpha
		Ui = c(Ui.sig, Ui.beta, Ui.alpha)
		UU = UU + tcrossprod(Ui)
    U[,i] = Ui
	}
	
	# hessian matrix on 1-1 corner
	H.sig.sig = diag(lower/n/sigma2)
  # hessian matrix on 2-2 corner
	H.gam.gam = XinvRX
	# hessian matrix on 1-2 corner
	H.sig.gam.left = matrix(0,nrow=M,ncol=p) # useful matrix
	H.sig.gam.right = matrix(0,nrow=M,ncol=q) # useful matrix
	for(i in 1:M) {
		H.sig.gam.left[i,] = t(beta[,i])%*%xx + t(alpha)%*%wx
		H.sig.gam.right[i,] = t(beta[,i])%*%xw + t(alpha)%*%ww 
	}
  H.sig.gam.left = vdiag(H.sig.gam.left)
	H.sig.gam = cbind(H.sig.gam.left,H.sig.gam.right) / n
  # hessian matrix on 2-1 corner
  H.gam.sig.upper = t(H.sig.gam.left) %*% diag(1/sigma2) %*% invR / 2 # pM x M
  H.gam.sig.lower = t(vdiag(H.sig.gam.right)) %*% diag(1/sigma2) %*% invR %*% lM / 2 # Mq x 1
  H.gam.sig = matrix(0,nrow=p*M+q,ncol=M)
  for(i in 1:M) {
    H.gam.sig[1:(p*M),i] = matrix(H.gam.sig.upper[(1:p)+(i-1)*p,],nrow=p*M)
    H.gam.sig[(p*M+1):(p*M+q),i] = matrix(H.gam.sig.lower[(1:q)+(i-1)*q,],nrow=q)
  }
	# hessian matrix
	H = rbind(cbind(H.sig.sig,H.sig.gam),cbind(H.gam.sig,H.gam.gam))
	
	# information matrix
	info = t(H)%*%solve(UU)%*%H
	
	# covariance matrix
	cov.mat = solve(info)
	
	# s.e.
	se.vec = sqrt(diag(cov.mat))
	se.sig = se.vec[1:M]
	se.gam = se.vec[(M+1):(p*M+q+M)]
  se.beta = se.gam[1:(p*M)]
  se.alpha = se.gam[(p*M+1):(p*M+q)]
  
  # p-value
  pval.alpha = 2*(1-pnorm(abs(alpha/se.alpha)))

  # score test for common exposure
  Lambda = diag(M)
  Lambda = Lambda[-1,]
  Lam.invR = Lambda%*%invR
  U2 = matrix(0,nrow=M-1,ncol=n)
  UU2 = matrix(0,nrow=M+p*M+q,ncol=M-1) # G12, Note that G11=UU
  #U2U = t(UU2) # G21
  U2U2 = matrix(0,nrow=M-1,ncol=M-1) # G22
  for(i in 1:n) {
    U2[,i] = c(w[i,]) * Lam.invR %*% resid[i,]
    UU2 = UU2 + tcrossprod(U[,i],U2[,i])
    #U2U = U2U + tcrossprod(U2[,i],U[,i])
    U2U2 = U2U2 + tcrossprod(U2[,i])
  }
  U2U = t(UU2)
  sum.U2 = rowSums(U2)
  # Note: big.U = rbind(U,U2) (not in summation), dimension = (p*M+q+M+M-1) x n

  # A matrix 
  A.left = matrix(0,nrow=M-1,ncol=M)
  A.right = matrix(0,nrow=M-1,ncol=p*M+q)
  for(i in 1:n) {
    A.left = A.left + 0.5 * c(w[i,]) * Lam.invR %*% diag(1/sigma2) %*% diag(c( crossprod(x[i,],beta) + c(w[i,])*c(alpha) ))
    A.right = A.right + c(w[i,]) * Lam.invR %*% cbind(kronecker(eyeM,t(x[i,])), c(w[i,])*lM)
  }
  A = cbind(A.left,A.right) # Note: same as the original SAS

  # SIGMA matrix
  AinvH = A%*%solve(H)
  tAinvH = t(AinvH)
  SIG = U2U2 - AinvH%*%UU2 - U2U%*%tAinvH + AinvH%*%UU%*%tAinvH

  # score statistic & p-value
  score.test = t(sum.U2) %*% solve(SIG) %*% sum.U2
  pval.common = 1-pchisq(score.test,M-1)

	# output a list
  out = list(resid=resid,beta=beta,se.beta=matrix(se.beta,nrow=p),
             alpha=c(alpha),se.alpha=c(se.alpha),
             sigma2=sigma2,se.sigma2=se.sig,
             pvalue=c(pval.alpha),pvalue.common=c(pval.common),is.conv=is.conv)
  
  if(!is.null(colnames(x))) {
    rownames(out$beta)=rownames(out$se.beta)=colnames(x) # dim p
  } else {
    rownames(out$beta)=rownames(out$se.beta)=paste("predictor",1:p,sep="") # dim p
  }
  
  if(!is.null(colnames(y))) {
    colnames(out$beta)=colnames(out$se.beta)=colnames(y) # dim M
  } else {
    colnames(out$beta)=colnames(out$se.beta)=paste("outcome",1:M,sep="") # dim M
  }
  
  return(out)
}


SMAT_resid <-
function(y,x,s,w=NULL,status=NULL,prev=NULL,working="unstr",conv=1e-4,maxiter=5e2) {
  d=w # weighting
  w=s # SNP effect
  
	if(!is.matrix(y)) { stop("y needs to be a matrix.") }
	if(!is.matrix(x)) { stop("x needs to be a matrix.") }

  cor.y = cor(y,use="pairwise.complete.obs")
  if(sum(cor.y<=0)>0) {
    warning("At least one pair of the phenotypes is not positively correlated; consider transformation.  As an alternative, consider M-DF test implemented within library geepack.")
  }

  if(is.vector(w)) w=as.matrix(w)
	#if(!is.matrix(w)) { stop("s needs to be a matrix.") }

	# check input types and dimensions
	n = nrow(y)
	M = ncol(y)
	p = ncol(x)
	q = 1 # set to 1 for the moment
	if(nrow(x)!=n | nrow(w)!=n) { stop("Row numbers of y, x, s do not match.") }
	
	if(!is.null(d)) {
		if(length(d)!=n) { stop("Row numbers of w and y do not match.") }
		is.wt = "user"; 
	} else if(!is.null(status)){ 
                status=as.vector(status)
                if(length(status)!=n) { stop("Row numbers of status and y do not match.") }
                if(is.null(prev)) { stop("prev should be provided if status is provided.") }
                if(length(setdiff(status,c(NA,1,0)))>0) { stop("status should only contain 0, 1, or NA.") }
		d = (status==1)*prev[1] + (status==0)*(1-prev[1]);  # prev should be a scalar
                is.wt = "casecontrol";
	} else {
            if(!is.null(prev)) {
                warning("Equal weight is used as status is not provided.")
            }
            d = rep(1,n)
            is.wt = "equal"
        } 

  # use only non-missing data
  index.observed=which( rowSums(is.na(y))==0 & 
                        rowSums(is.na(x))==0 & 
                        is.na(d)==FALSE & 
                        rowSums(is.na(w))==FALSE )
  if(length(index.observed)==0) stop("No complete observations are found!")
  y=y[index.observed,,drop=F]
  x=x[index.observed,,drop=F]
  d=d[index.observed]
  w=as.matrix(w[index.observed,,drop=F])

  if(is.wt=="casecontrol") {
    status=status[index.observed]
    pn=mean(status)
    d = d / ( (status==1)*pn + (status==0)*(1-pn) )
  }

  # use the new dimensions
  n = nrow(y)
  M = ncol(y)
  p = ncol(x)
    
  sum.d = sum(d) # udpate

  # convergence indicator
  is.conv = TRUE

	# precomputation
	if(is.wt!="equal") {
		y = diag(sqrt(d))%*%y
		x = diag(sqrt(d))%*%x
		w = diag(sqrt(d))%*%w	
	}
	yy=crossprod(y,y)
	yx=crossprod(y,x)
	yw=crossprod(y,w)
	xx=crossprod(x,x)
	xw=crossprod(x,w)
	xy=t(yx)
	ww=crossprod(w,w)
	wx=t(xw)
	wy=t(yw)

	# initial values of parameters
	gamhat = rep(0,p*M+q)
	sigma2 = apply(y,2,sd)^2 # sample variance
	psi = diag(sigma2)
	eyeM = diag(M)
	lM = matrix(1,ncol=1,nrow=M)
	if(working=="ind") { 
		R=eyeM 
	} else if(working=="exch") { 
		theta=mean(lowerTriangle(cor(y)))
		R=matrix(theta,M,M)
		diag(R)=1 
	} else { 
		R=cor(y) 
	} # unstructured, sample correlation

	lower = rep(1,M)

	# big while loop for estimation ---------------------------
	max.diff = 1
	iter = 1
	while(max.diff>conv & iter<maxiter) {

		invR = solve(R)

		#------------------
		#  update gamhat 
		#------------------
		gamhat.old = gamhat

		phi = diag(1/sqrt(c(sigma2))) # not psi
		invR.phi = invR %*% phi
		l.invR.phi = t(lM) %*% invR.phi
		l.invR = t(lM) %*% invR
		invR.l = t(l.invR)
		l.invR.l = l.invR %*% lM
		XinvRX = rbind(cbind(kronecker(invR,xx),kronecker(invR.l,xw)),
                       cbind(kronecker(l.invR,wx),kronecker(ww,l.invR.l)))
		XinvRy = rbind(matrix(t(invR.phi%*%yx),ncol=1),
                       matrix(l.invR.phi%*%yw,ncol=1))

		gamhat=solve(XinvRX,XinvRy)
		beta = matrix(gamhat[1:(M*p)],nrow=p)
		alpha = matrix(gamhat[(M*p+1):(M*p+q)],nrow=q)

		#------------------
		#  update sigma2
		#------------------
		sigma2.old = sigma2

		for(j in 1:M) {
 			old2=sigma2[j]
			incre=old2
			while(abs(incre)/old2>=conv) {
				old2  = sigma2[j]
				old   = sqrt(sigma2[j])
				bj    = beta[,j]
				upper = yy[j,j]-old*yx[j,]%*%bj-old*yw[j]%*%alpha-sum.d*old2
				lower[j] = sum.d+t(bj)%*%xx%*%bj/2+t(alpha)%*%wx%*%bj+t(alpha)%*%ww%*%alpha/2
				incre = upper/lower[j]
				sigma2[j] = old2 + incre
			}
		}
		
		#------------------
		#  update R
		#------------------
		if(working=="exch" | working=="unstr") {
			phi = diag(1/sqrt(c(sigma2))) # not psi, updated by new sigma2
  			a = alpha%*%t(lM)
			rr = phi%*%yy%*%phi - phi%*%yx%*%beta - t(beta)%*%t(yx)%*%phi - phi%*%yw%*%a - t(a)%*%t(yw)%*%phi + t(beta)%*%xx%*%beta + t(a)%*%ww%*%a + t(beta)%*%xw%*%a + t(a)%*%wx%*%beta
			sum.rr = sum(lowerTriangle(rr))
			if(working=="exch") {
				theta = sum.rr / (n*M*(M-1)/2 - p*M - q - M) # double check this
				R = matrix(theta,M,M)
				diag(R) = 1
			} else {
				rr = rr/n
				R = rr
				diag(R) = 1
			}
		}
		
		# compute differences
		max.diff = max(abs(gamhat-gamhat.old)/abs(gamhat.old),abs(sigma2-sigma2.old)/sigma2.old)
		iter = iter+1
		
	} # end of big while loop
	# big while loop for estimation completed -----------------
  
  # check if maximum #iter is reched
  if(iter>=maxiter) { is.conv = FALSE; warning("Maximum number of iterations reached. Please increase maxiter and rerun the function.") }
  
	# invR
	invR = solve(R)
	
	# y.star
	y.star = y%*%diag(1/sqrt(sigma2))
	# residual matrix at estimated values (nxM)
	resid = y.star - x%*%beta - w%*%alpha%*%t(lM)
	
	# score * t(score) in summation
	U.sig = (y.star * resid - d %*% t(lM))/n # (nxM)  # NOTE THIS 1/n
	U.beta = matrix(0,nrow=n,ncol=p*M)
	U.alpha = matrix(0,nrow=n,ncol=q)
	resid.invR = resid %*% invR # (nxM)
	UU = matrix(0,nrow=p*M+q+M,ncol=p*M+q+M)
  U = matrix(0,nrow=p*M+q+M,ncol=n) # useful in testing for common exposure
	for(i in 1:n) {
		Ui.sig = U.sig[i,]
		Ui.beta = as.vector( x[i,] %*% t(resid.invR[i,]) )
		U.beta[i,] = Ui.beta
		Ui.alpha = w[i,] %*% t(resid.invR[i,]) %*% lM
		U.alpha[i,] = Ui.alpha
		Ui = c(Ui.sig, Ui.beta, Ui.alpha)
		UU = UU + tcrossprod(Ui)
    U[,i] = Ui
	}
	
	# hessian matrix on 1-1 corner
	H.sig.sig = diag(lower/n/sigma2)
  # hessian matrix on 2-2 corner
	H.gam.gam = XinvRX
	# hessian matrix on 1-2 corner
	H.sig.gam.left = matrix(0,nrow=M,ncol=p) # useful matrix
	H.sig.gam.right = matrix(0,nrow=M,ncol=q) # useful matrix
	for(i in 1:M) {
		H.sig.gam.left[i,] = t(beta[,i])%*%xx + t(alpha)%*%wx
		H.sig.gam.right[i,] = t(beta[,i])%*%xw + t(alpha)%*%ww 
	}
  H.sig.gam.left = vdiag(H.sig.gam.left)
	H.sig.gam = cbind(H.sig.gam.left,H.sig.gam.right) / n
  # hessian matrix on 2-1 corner
  H.gam.sig.upper = t(H.sig.gam.left) %*% diag(1/sigma2) %*% invR / 2 # pM x M
  H.gam.sig.lower = t(vdiag(H.sig.gam.right)) %*% diag(1/sigma2) %*% invR %*% lM / 2 # Mq x 1
  H.gam.sig = matrix(0,nrow=p*M+q,ncol=M)
  for(i in 1:M) {
    H.gam.sig[1:(p*M),i] = matrix(H.gam.sig.upper[(1:p)+(i-1)*p,],nrow=p*M)
    H.gam.sig[(p*M+1):(p*M+q),i] = matrix(H.gam.sig.lower[(1:q)+(i-1)*q,],nrow=q)
  }
	# hessian matrix
	H = rbind(cbind(H.sig.sig,H.sig.gam),cbind(H.gam.sig,H.gam.gam))
	
	# information matrix
	info = t(H)%*%solve(UU)%*%H
	
	# covariance matrix
	cov.mat = solve(info)
	
	# s.e.
	se.vec = sqrt(diag(cov.mat))
	se.sig = se.vec[1:M]
	se.gam = se.vec[(M+1):(p*M+q+M)]
  se.beta = se.gam[1:(p*M)]
  se.alpha = se.gam[(p*M+1):(p*M+q)]
  
  # p-value
  pval.alpha = 2*(1-pnorm(abs(alpha/se.alpha)))

  # score test for common exposure
  Lambda = diag(M)
  Lambda = Lambda[-1,]
  Lam.invR = Lambda%*%invR
  U2 = matrix(0,nrow=M-1,ncol=n)
  UU2 = matrix(0,nrow=M+p*M+q,ncol=M-1) # G12, Note that G11=UU
  #U2U = t(UU2) # G21
  U2U2 = matrix(0,nrow=M-1,ncol=M-1) # G22
  for(i in 1:n) {
    U2[,i] = c(w[i,]) * Lam.invR %*% resid[i,]
    UU2 = UU2 + tcrossprod(U[,i],U2[,i])
    #U2U = U2U + tcrossprod(U2[,i],U[,i])
    U2U2 = U2U2 + tcrossprod(U2[,i])
  }
  U2U = t(UU2)
  sum.U2 = rowSums(U2)
  # Note: big.U = rbind(U,U2) (not in summation), dimension = (p*M+q+M+M-1) x n

  # A matrix 
  A.left = matrix(0,nrow=M-1,ncol=M)
  A.right = matrix(0,nrow=M-1,ncol=p*M+q)
  for(i in 1:n) {
    A.left = A.left + 0.5 * c(w[i,]) * Lam.invR %*% diag(1/sigma2) %*% diag(c( crossprod(x[i,],beta) + c(w[i,])*c(alpha) ))
    A.right = A.right + c(w[i,]) * Lam.invR %*% cbind(kronecker(eyeM,t(x[i,])), c(w[i,])*lM)
  }
  A = cbind(A.left,A.right) # Note: same as the original SAS

  # SIGMA matrix
  AinvH = A%*%solve(H)
  tAinvH = t(AinvH)
  SIG = U2U2 - AinvH%*%UU2 - U2U%*%tAinvH + AinvH%*%UU%*%tAinvH

  # score statistic & p-value
  score.test = t(sum.U2) %*% solve(SIG) %*% sum.U2
  pval.common = 1-pchisq(score.test,M-1)

  resid2<-resid*d

  
	# output a list
  out = list(resid2=resid2,resid=resid,beta=beta,se.beta=matrix(se.beta,nrow=p),
             alpha=c(alpha),se.alpha=c(se.alpha),
             sigma2=sigma2,se.sigma2=se.sig,
             pvalue=c(pval.alpha),pvalue.common=c(pval.common),is.conv=is.conv)
  
  if(!is.null(colnames(x))) {
    rownames(out$beta)=rownames(out$se.beta)=colnames(x) # dim p
  } else {
    rownames(out$beta)=rownames(out$se.beta)=paste("predictor",1:p,sep="") # dim p
  }
  
  if(!is.null(colnames(y))) {
    colnames(out$beta)=colnames(out$se.beta)=colnames(y) # dim M
  } else {
    colnames(out$beta)=colnames(out$se.beta)=paste("outcome",1:M,sep="") # dim M
  }
  
  return(out)
}

