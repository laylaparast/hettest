delta.e.estimate = function(sone=NULL,szero=NULL, szerop, yzerop, extrapolate = TRUE,mat = NULL, n1=NULL, n0=NULL) {
	if(!is.null(mat)){
		sone = mat[1:n1]
		szero = mat[(n1+1):(n1+n0)]
	}
	h.paper4 = bw.nrd(szerop)*(length(szerop)^(-0.11))
	mu.s0 = sapply(szero,pred.smooth.2,kernel.use=szerop, bw=h.paper4, outcome=yzerop)
  	if(sum(is.na(mu.s0))>0 & extrapolate){
  		print(paste("Note: ", sum(is.na(mu.s0)), " values extrapolated."))
    	c.mat = cbind(szero, mu.s0)
    	for(o in 1:length(mu.s0)) {
    		if(is.na(mu.s0[o])){
    			distance = abs(szero - szero[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			mu.s0[o] = new.est[1]   #in case there are multiple matches
    		}
    	}
		}
	mu.s1 = sapply(sone,pred.smooth.2,kernel.use=szerop, bw=h.paper4, outcome=yzerop)
  	if(sum(is.na(mu.s1))>0 & extrapolate){
  		print(paste("Note: ", sum(is.na(mu.s1)), " values extrapolated."))
    	c.mat = cbind(sone, mu.s1)
    	for(o in 1:length(mu.s1)) {
    		if(is.na(mu.s1[o])){
    			distance = abs(sone - sone[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			mu.s1[o] = new.est[1]   #in case there are multiple matches
    		}
    	}
		}
	delta.e = mean(mu.s1) - mean(mu.s0)	
	sd.e = sqrt(var(mu.s1)/length(sone)+var(mu.s0)/length(szero))
	test.statistic.e = delta.e/sd.e
	p.value.e = 2*(1- pnorm(abs(test.statistic.e)))
	return(list("delta.e" = delta.e, "sd.e" = sd.e, "test.statistic.e" = test.statistic.e, "p.value.e" = p.value.e))
}

pred.smooth.2 <-function(kernel.use,kernel.apply, bw,outcome) { 	
    return(sum(Kern.FUN(kernel.use,kernel.apply,bw=bw)*outcome)/sum(Kern.FUN(kernel.use,kernel.apply,bw=bw)))
 
  }

Kern.FUN <- function(zz,zi,bw) 
  { 
    out = (VTM(zz,length(zi))- zi)/bw
	dnorm(out)/bw
           
  }
  
 VTM<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }

kf=function(x, h){return(dnorm(x/h)/h)}

 
delta.h.estimate = function(sone=NULL,szero=NULL,wone=NULL, wzero=NULL, szerop, wzerop, yzerop, extrapolate = TRUE,mat = NULL, n1=NULL, n0=NULL) {
	c.adj=1
	if(!is.null(mat)){
		sone = mat[1:n1]
		szero = mat[(n1+1):(n1+n0)]
		wone = mat[(n1+n0+1):(n1+n0+n1)]
		wzero = mat[(n1+n0+n1+1):(n1+n0+n1+n0)]
	}

	h.paper2 =c.adj*(sqrt(4)*bw.nrd(szerop)*length(szerop)^(1/5))*(length(szerop)^(-0.4))  #hy2s1
	h.paper3 =c.adj*(sqrt(4)*bw.nrd(wzerop)*length(wzerop)^(1/5))*(length(wzerop)^(-0.4)) #hy2w1
	
	#treated group, current study
	h.paper0 =c.adj*(bw.nrd(wone)*length(wone)^(1/5))*(length(wone)^(-0.4))   #hsw1
	temp.1 = sapply(wone, mu1swy.h,s.use=szerop, w.use=wzerop, y.use=yzerop, s.apply=sone, h.paper3=h.paper3, h.paper2=h.paper2)
	weight.11=Kern.FUN(zz=wone, zi=wone, bw=h.paper0)
	m1.1 = apply(temp.1*t(weight.11),2,sum)/apply(t(weight.11),2,sum)
	S.1 = diag(temp.1)
	
	temp.0 = sapply(wzero, mu1swy.h,s.use=szerop, w.use=wzerop, y.use=yzerop, s.apply=szero, h.paper3=h.paper3, h.paper2=h.paper2)
	weight.00=Kern.FUN(zz=wzero, zi=wzero, bw=h.paper0)
	m0.0 = apply(temp.0*t(weight.00),2,sum)/apply(t(weight.00),2,sum)	
	S.0 = diag(temp.0)
	if(sum(is.na(m0.0))>0) {print(paste("Number of extreme values excluded:",sum(is.na(m0.0))))}
		
	weight.10=Kern.FUN(zz=wone, zi=wzero, bw=h.paper0)
    m1.W0 = (weight.10 %*% S.1)/apply(weight.10,1,sum)
	weight.01=Kern.FUN(zz=wzero, zi=wone, bw=h.paper0)
    m0.W1 = (weight.01 %*% S.0)/apply(weight.01,1,sum)
	n1 = length(sone); n0 = length(szero)
	
	first.term = m1.W0 - m0.0
	second.term = m1.1 - m0.W1
	num.na = sum(is.na(first.term)) + sum(is.na(second.term))
	
	if(num.na>0) {print(paste("Number of extreme values excluded:",num.na))}
	
	delta.h = (1/(n1+n0-num.na))*(sum(first.term, na.rm = T) + sum(second.term, na.rm = T))
	
	#for variance
	mu1 = mean(m1.1,na.rm = T)
	mu0 = mean(m0.0, na.rm = T)
	pi.1 = n1/(n1+n0); pi.0 = n0/(n1+n0)
	var.first.term = 1/(n1^2)*sum((S.1 - pi.0*m1.1 - pi.1*m0.W1-pi.1*(mu1-mu0))^2)
	var.second.term = 1/(n0^2)*sum((S.0 - pi.0*m1.W0 - pi.1*m0.0-pi.1*(mu0-mu1))^2)
	var.term = var.first.term + var.second.term
	
	sd.h = sqrt(var.term)
	test.statistic.h = delta.h/sd.h
	p.value.h = 2*(1- pnorm(abs(test.statistic.h)))
	return(list("delta.h" = delta.h, "sd.h" = sd.h, "test.statistic.h" = test.statistic.h, "p.value.h" = p.value.h))
}


mu1swy.h=function(s.use,w.use,y.use,s.apply, w.grd, h.paper3, h.paper2)
{m=length(s.apply)
 res=rep(NA, m)
 for(i in 1:m)
   {weight=kf(s.apply[i]-s.use, h.paper2)*kf(w.grd-w.use, h.paper3)
    res[i]=sum(weight*y.use)/sum(weight)
    }
  if(sum(is.na(res))!=0) {
  	ind = which(is.na(res))
  	for(i in ind) 
   		{mat = cbind(abs(s.apply-s.apply[i]), res) 
   		 mmm = which(mat[,1] == min(mat[-ind,1]))[1]   
   		 res[i] = res[mmm]  
   		}}
return(res)
}
 
