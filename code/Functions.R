library(MASS)

##main function

##input 
##W the covariates; X the observed response; D the treatment assignment
##p dimension of brain regions; n number of subjects; m number of repeated observation
##q dimension of covariates; B bootstrap sampling size; aph significant level ; c augmentation proportion

##output 
##Tstat: a p by p matrix, Tstat(i,j) is the estimated treatment effect of the connection between region i and region j;
##result: a p by p matrix, result(i,j) is 0 or 1, result(i,j)=1 means treatment effect on connection between
##region i and region j is significant; result(i,j)=0 means the treatment effect is not significant

corr_estimate<-function(W,X,D,p,n,m,q,B,aph,c){
  
  #estimation for correlation matrix 
  
  Y.estimate_<-apply(X,3,cor,simplify = FALSE)
  Y.estimate<-array(unlist(Y.estimate_),dim=c(p,p,n))
  
  #logistic regression for propensity score
  
  logistReg = glm(D ~ W + 0, family = binomial)
  beta.estimate = logistReg$coefficient
  Prob.estimate = as.vector(exp(W %*% beta.estimate) / (1 + exp(W %*% beta.estimate)))
  
  #point estimation of the treatment effect
  
  weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)
  tau.estimate<-apply(Y.estimate* weight[slice.index(Y.estimate, 3)],c(1,2),mean)
  
  #estimation for eta, the influence function and variance
  Fisher.Information = 0
  for (i in 1 : n){
    Fisher.Information = Fisher.Information + Prob.estimate[i] * (1 - Prob.estimate[i]) * W[i, ] %*% t(W[i, ])
  }
  Fisher.Information = Fisher.Information / n
  Theta = solve(Fisher.Information)
  
  weight.H = D * (1 - Prob.estimate) / Prob.estimate + (1 - D) * Prob.estimate / (1 - Prob.estimate)
  H = array(0, c(p, p, q))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      H[j1, j2, ] = colMeans(weight.H * Y.estimate[j1, j2, ] * W)
    }
  }
  
  eta = array(0, c(p, p, n))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      for (i in 1 : n){
        eta[j1, j2, i] = Y.estimate[j1, j2, i] * weight[i] - H[j1, j2, ] %*% Theta %*% W[i, ] * (D[i] - Prob.estimate[i])
      }
    }
  }
  
  theta.var<-apply( (eta-tau.estimate[slice.index(eta,c(1,2))])^2,c(1,2),mean)
  
  
  #standardized treatment effect
  Tstat = sqrt(n) * abs(tau.estimate) * theta.var^(-0.5)
  #set variance term, the diagonal to be zero
  Index = matrix(1, p, p)-diag(1,p,p)
  Tstat0 = Tstat * Index
  

  #Proposed multiple testing procedure
  z = array(0, c(p, p, B))
  for(b in 1 : B){
    g = rnorm(n)
    temp =apply( (eta-tau.estimate[slice.index(eta,c(1,2))])*g[slice.index(eta, 3)],c(1,2) , sum)
    z[, , b] = theta.var^(-0.5) * temp / sqrt(n) * Index
  }
  
  result<-matrix(0,p,p)
  z.initial = z
  Tstat.initial = Tstat0
  repeat{
    Tstat.max = max(abs(Tstat.initial))
    index.temp = which(abs(Tstat0) == Tstat.max, arr.ind = TRUE)
    z.max = apply(abs(z.initial), 3, max)
    z.max.quan = quantile(z.max, 1-aph)
    if (Tstat.max < z.max.quan) break
    for (i in 1 : dim(index.temp)[1]){
      Tstat.initial[index.temp[i, 1], index.temp[i, 2]] = 0
      z.initial[index.temp[i, 1], index.temp[i, 2], ] = 0
      result[index.temp[i,1],index.temp[i,2]]=1
    }
  }
  
  #Augmentation
  size<-sum(result)/2
  num_add<-floor(c*size/(1-c))
  if(num_add>=1){
    test_replace<-sort(Tstat.initial,decreasing=TRUE)
    for(i in 1:(2*num_add)){
      index.temp = which(abs(Tstat0) == test_replace[i], arr.ind = TRUE)
      for (i in 1 : dim(index.temp)[1]){
        result[index.temp[i,1],index.temp[i,2]]=1
      }
    }}
  
  re<-list(Tstat,result)
  return(re)}



##functions to generate simulation data

##covariance matrix generating function
ARcov = function(p, rho){
  cov0 = matrix(0, p, p)
  for(i in 1 : p){
    for (j in 1 : p){
      cov0[i, j] = rho^(abs(i - j))
    }
  }
  return(cov0)
}


##function for data generating: block diagonal setting, treatment effect 
##the nonzero treatment effects appear in the blocks of size s0 along the main diagonal of the correlation matrice

##input
##detailed meaning of these parameters can be found in the article
##p is dimension of brain regions; Raw0 is the parameter for correlation structure; q is dimension for covariates
##Beta is the true logistic model parameter; gamma decides the relation between the treatment effect and covariates
##n subject number; m time series length; s0 is the size of `block'`;  raw is strength of time dependence

##output: a list of W,D,X
##W: simulated covariates; X: simulated response; D: simulated treatment assignment

Sample_1<-function(p,Raw0,q,Beta,gamma,n,m,s0,raw){
  Sigma0 = ARcov(p, Raw0)
  Sigmaw = ARcov(q, 0.7) * 0.5
  W = mvrnorm(n, rep(0, q), Sigmaw)
  Prob = as.vector(exp(W %*% Beta) / (1 + exp(W %*% Beta)))
  
  #--- Data simulation---
  D = c()
  for(i in 1 : n){
    D[i] = rbinom(1, 1, Prob[i])
  }
  
  epsilon = rnorm(n, 0, 0.1)
  delta = as.vector(exp(abs(W %*% gamma)) / (1 + exp(abs(W %*% gamma)))) - 0.5 + epsilon 
  # need to take absolute value here otherwise the mean of delta is close to zero, no average treatment effect.
  X = array(0, c(m, p, n))
  E = array(0, c(m, p, n))
  
  #block diagonal signal structure
  k=p/s0
  signa.vector=array(0,c(k,p))
  for(i in 1:k){
    signa.vector[i,(((i-1)*s0+1):(i*s0))]<-rep(1,s0)
  }
  signal<-t(signa.vector)%*%signa.vector
  signal=signal-diag(rep(1,p))
  for (i in 1 : n){
    Y0 = Sigma0
    Y1 = Sigma0 + signal* abs(delta[i])*Raw0*2
    Y = D[i] * Y1 + (1 - D[i]) * Y0
    
    # for time dependent data
    E[, , i] = mvrnorm(m, rep(0,p), Y)
    for(j in 1:m){
      if(j==1){X[j, , i]=E[j, , i]}
      if(j>=2){X[j, , i]=raw*X[j-1, , i]+sqrt(1-raw^2)*E[j, , i]}
    }
  }
  outcome<-list(W,D,X)
  return(outcome)}




##function for data generating: super diagonal setting
##the nonzero causal effects are only on the s0th to s1th super diagonals of the correlation matrice

##input
##detailed meaning of these parameters can be found in the article
##p is dimension of brain regions; Raw0 is the parameter for correlation structure; q is dimension for covariates
##Beta is the true logistic model parameter; gamma decides the relation between the treatment effect and covariates
##n subject number; m time series length; raw is strength of time dependence

##output: a list of W,D,X
##W: simulated covariates; X: simulated response; D: simulated treatment assignment


Sample_2<-function(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw){
  Sigma0 = ARcov(p, Raw0)
  Sigmaw = ARcov(q, 0.7) * 0.5
  W = mvrnorm(n, rep(0, q), Sigmaw)
  Prob = as.vector(exp(W %*% Beta) / (1 + exp(W %*% Beta)))
  
  #--- Data simulation---
  D = c()
  for(i in 1 : n){
    D[i] = rbinom(1, 1, Prob[i])
  }
  
  epsilon = rnorm(n, 0, 0.1)
  delta = as.vector(exp(abs(W %*% gamma)) / (1 + exp(abs(W %*% gamma)))) - 0.5 + epsilon 
  # need to take absolute value here otherwise the mean of delta is close to zero, no average treatment effect.
  X = array(0, c(m, p, n))
  E = array(0, c(m, p, n))
  
  #super diagonal
  signal=matrix(0,nrow=p,ncol=p)
  for(k in s0:s1){
    for(i in 1:(p-k+1)){
      signal[i,i+k-1]=1
      signal[i+k-1,i]=1}
  }
  
  for (i in 1 : n){
    Y0 = Sigma0 
    Y1 = Sigma0 + 2*signal* abs(delta[i])*Raw0
    Y = D[i] * Y1 + (1 - D[i]) * Y0
    
    # for time dependent data
    E[, , i] = mvrnorm(m, rep(0,p), Y)
    for(j in 1:m){
      if(j==1){X[j, , i]=E[j, , i]}
      if(j>=2){X[j, , i]=raw*X[j-1, , i]+sqrt(1-raw^2)*E[j, , i]}
    }
  }
  outcome<-list(W,D,X)
  return(outcome)}



#This function is used for data generating: Super diagonal setting
#When Raw0 is large, there is a prob that covariance matrix is not positive definite.
#`Sample_3' is similar to function `Sample_2'`
#This function will skip those non-positive-definite matrices and generate another one, which is the only difference from `Sample_2`


Sample_3<-function(p,Raw0,q,Beta,gamma,n,m,s0,s1,raw){
  Sigma0 = ARcov(p, Raw0)
  Sigmaw = ARcov(q, 0.7) * 0.5
  W = mvrnorm(n, rep(0, q), Sigmaw)
  Prob = as.vector(exp(W %*% Beta) / (1 + exp(W %*% Beta)))
  
  #--- Data simulation---
  D = c()
  for(i in 1 : n){
    D[i] = rbinom(1, 1, Prob[i])
  }
  
  epsilon = rnorm(n, 0, 0.1)
  delta = as.vector(exp(abs(W %*% gamma)) / (1 + exp(abs(W %*% gamma)))) - 0.5 + epsilon 
  # need to take absolute value here otherwise the mean of delta is close to zero, no average treatment effect.
  X = array(0, c(m, p, n))
  E = array(0, c(m, p, n))
  
  signal=matrix(0,nrow=p,ncol=p)
  for(k in s0:s1){
    for(i in 1:(p-k+1)){
      signal[i,i+k-1]=1
      signal[i+k-1,i]=1}
  }
  
  for (i in 1 : n){
    while(1){
      Y0 = Sigma0 
      Y1 = Sigma0 + 2*signal* abs(delta[i])*Raw0
      Y = D[i] * Y1 + (1 - D[i]) * Y0
      if(all(eigen(Y)$values>0)==1){break}
      epsilon_1 = rnorm(1, 0, 0.1)
      delta[i] = as.vector(exp(abs(W %*% gamma)) / (1 + exp(abs(W %*% gamma))))[i] - 0.5 + epsilon_1
      
    }
    # for time dependent data
    E[, , i] = mvrnorm(m, rep(0,p), Y)
    for(j in 1:m){
      if(j==1){X[j, , i]=E[j, , i]}
      if(j>=2){X[j, , i]=raw*X[j-1, , i]+sqrt(1-raw^2)*E[j, , i]}
    }
  }
  outcome<-list(W,D,X)
  return(outcome)}





## Functions below are extension of the main function `corr_estimate` for simulation studies
## After the estimation of treatment effect matrix and step-down multiple testing procedure, 
## they calculate the true positive rate, the false positive rate and fdp


## This function uses the proposed step-down procedure and calculate tpr, fpr and fdp under block diagonal setting

##input 
##W the covariates; X the observed response; D the treatment
##p dimension of brain regions; n number of subjects; m number of repeated observation
##q dimension of covariates; B bootstrap sampling size; aph significant level ; c augmentation proportion
##s0 is the block size of the block diagonal setting

##output 
##p_true_rate: true positive rate
##p_false_rate: false positive rate
##fdr: false positive proportion

corr_estimate_diag<-function(W,X,D,p,n,m,q,B,aph,c,s0){
 
  #causal effect estimation 
  
  Y.estimate_<-apply(X,3,cor,simplify = FALSE)
  Y.estimate<-array(unlist(Y.estimate_),dim=c(p,p,n))
  
  
  logistReg = glm(D ~ W + 0, family = binomial)
  beta.estimate = logistReg$coefficient
  Prob.estimate = as.vector(exp(W %*% beta.estimate) / (1 + exp(W %*% beta.estimate)))
  
  weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)
  
  tau.estimate<-apply(Y.estimate* weight[slice.index(Y.estimate, 3)],c(1,2),mean)
  
  
  
  #estimation for eta, the influence function and variance
  Fisher.Information = 0
  for (i in 1 : n){
    Fisher.Information = Fisher.Information + Prob.estimate[i] * (1 - Prob.estimate[i]) * W[i, ] %*% t(W[i, ])
  }
  Fisher.Information = Fisher.Information / n
  Theta = solve(Fisher.Information)
  
  weight.H = D * (1 - Prob.estimate) / Prob.estimate + (1 - D) * Prob.estimate / (1 - Prob.estimate)
  H = array(0, c(p, p, q))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      H[j1, j2, ] = colMeans(weight.H * Y.estimate[j1, j2, ] * W)
    }
  }
  
  eta = array(0, c(p, p, n))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      for (i in 1 : n){
        eta[j1, j2, i] = Y.estimate[j1, j2, i] * weight[i] - H[j1, j2, ] %*% Theta %*% W[i, ] * (D[i] - Prob.estimate[i])
      }
    }
  }
  
  theta.var<-apply( (eta-tau.estimate[slice.index(eta,c(1,2))])^2,c(1,2),mean)
  
  
  #standardized treatment effect
  Tstat = sqrt(n) * abs(tau.estimate) * theta.var^(-0.5)
  #set variance term, the diagonal to be zero
  
  Index = matrix(1, p, p)-diag(1,p,p)
  
  Tstat0 = Tstat * Index
  
  
  #Proposed multiple testing procedure
  z = array(0, c(p, p, B))
  for(b in 1 : B){
    g = rnorm(n)
    temp =apply( (eta-tau.estimate[slice.index(eta,c(1,2))])*g[slice.index(eta, 3)],c(1,2) , sum)
    z[, , b] = theta.var^(-0.5) * temp / sqrt(n) * Index
  }
  
  result<-matrix(0,p,p)
  z.initial = z
  Tstat.initial = Tstat0
  repeat{
    Tstat.max = max(abs(Tstat.initial))
    index.temp = which(abs(Tstat0) == Tstat.max, arr.ind = TRUE)
    z.max = apply(abs(z.initial), 3, max)
    z.max.quan = quantile(z.max, 1-aph)
    if (Tstat.max < z.max.quan) break
    for (i in 1 : dim(index.temp)[1]){
      Tstat.initial[index.temp[i, 1], index.temp[i, 2]] = 0
      z.initial[index.temp[i, 1], index.temp[i, 2], ] = 0
      result[index.temp[i,1],index.temp[i,2]]=1
    }
  }
  
  #Augment
  size<-sum(result)/2
  num_add<-floor(c*size/(1-c))
  if(num_add>=1){
    test_replace<-sort(Tstat.initial,decreasing=TRUE)
    for(i in 1:(2*num_add)){
      index.temp = which(abs(Tstat0) == test_replace[i], arr.ind = TRUE)
      for (i in 1 : dim(index.temp)[1]){
        result[index.temp[i,1],index.temp[i,2]]=1
      }
    }}
  
  
  
  #This part for calculation of simulation results
  p_true<-0
  p_false<-0
  fdr<-0
  
  k=p/s0
  signa.vector=array(0,c(k,p))
  for(i in 1:k){
    signa.vector[i,(((i-1)*s0+1):(i*s0))]<-rep(1,s0)
  }
  signal<-t(signa.vector)%*%signa.vector
  
  for(i in 1:p){
    for(j in 1:p){
      if(j>i&&signal[i,j]==1&&result[i,j]==1){p_true=p_true+1}
      if(j>i&&signal[i,j]==0&&result[i,j]==1){p_false=p_false+1}
    }}
  
  if((p_true+p_false)>0){fdr<-p_false/(p_false+p_true)}else{
    fdr<-0}
  p_true_rate<-p_true/((sum(signal)-p)/2)
  p_false_rate<-p_false/((p*p-sum(signal))/2)
  
  #print results
  #print(Tstat)
  #print(result)
  #print(p_true_rate)
  #print(p_false_rate)
  #print(fdr)
  re<-c(p_true_rate,p_false_rate,fdr)
  return(re)}




## This function uses the BH procedure and calculate tpr, fpr and fdp under block diagonal setting

##input 
##W the covariates; X the observed response; D the treatment
##p dimension of brain regions; n number of subjects; m number of repeated observation
##q dimension of covariates; B bootstrap sampling size; aph significant level ; c augmentation proportion
##s0 is the block size of the block diagonal setting

##output 
##p_true_rate: true positive rate
##p_false_rate: false positive rate
##fdr: false positive proportion


corr_estimate_bh<-function(W,X,D,p,n,m,q,B,aph,c,s0){
  
  #causal effect estimation 
  
  Y.estimate_<-apply(X,3,cor,simplify = FALSE)
  Y.estimate<-array(unlist(Y.estimate_),dim=c(p,p,n))
  
  logistReg = glm(D ~ W + 0, family = binomial)
  beta.estimate = logistReg$coefficient
  Prob.estimate = as.vector(exp(W %*% beta.estimate) / (1 + exp(W %*% beta.estimate)))
  
  weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)
  
  tau.estimate<-apply(Y.estimate* weight[slice.index(Y.estimate, 3)],c(1,2),mean)
  
  
  #estimation for eta, the influence function and variance
  Fisher.Information = 0
  for (i in 1 : n){
    Fisher.Information = Fisher.Information + Prob.estimate[i] * (1 - Prob.estimate[i]) * W[i, ] %*% t(W[i, ])
  }
  Fisher.Information = Fisher.Information / n
  Theta = solve(Fisher.Information)
  weight.H = D * (1 - Prob.estimate) / Prob.estimate + (1 - D) * Prob.estimate / (1 - Prob.estimate)
  H = array(0, c(p, p, q))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      H[j1, j2, ] = colMeans(weight.H * Y.estimate[j1, j2, ] * W)
    }
  }
  
  eta = array(0, c(p, p, n))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      for (i in 1 : n){
        eta[j1, j2, i] = Y.estimate[j1, j2, i] * weight[i] - H[j1, j2, ] %*% Theta %*% W[i, ] * (D[i] - Prob.estimate[i])
      }
    }
  }

  theta.var<-apply( (eta-tau.estimate[slice.index(eta,c(1,2))])^2,c(1,2),mean)
  
  
  #standardized treatment effect
  Tstat = sqrt(n) * abs(tau.estimate) * theta.var^(-0.5)
  
  #P-value
  P=2*pnorm(Tstat,0,1,lower.tail = FALSE)
  P1=P
  for(i in 1:p){
    for(j in 1:p){
      if(j<=i){P1[i,j]=1}
    }
  }
  #BH-procedure
  test=sort(P1)
  num=1
  while(test[num]<num*aph/(p*(p-1)/2)){
    num=num+1
  }
  num=num-1
  
  result<-matrix(0,p,p)
  if(num>1){
    for(i in 1:num){
      index.temp = which(P == test[i], arr.ind = TRUE)
      for (j in 1 : dim(index.temp)[1]){
        result[index.temp[j,1],index.temp[j,2]]=1
      }
    }
  }
  
  #calculation for simulation results
  p_true<-0
  p_false<-0
  fdr<-0
  
  k=p/s0
  signa.vector=array(0,c(k,p))
  for(i in 1:k){
    signa.vector[i,(((i-1)*s0+1):(i*s0))]<-rep(1,s0)
  }
  signal<-t(signa.vector)%*%signa.vector
  
  for(i in 1:p){
    for(j in 1:p){
      if(j>i&&signal[i,j]==1&&result[i,j]==1){p_true=p_true+1}
      if(j>i&&signal[i,j]==0&&result[i,j]==1){p_false=p_false+1}
    }}
  
  if((p_true+p_false)>0){fdr<-p_false/(p_false+p_true)}else{
    fdr<-0}
  p_true_rate<-p_true/((sum(signal)-p)/2)
  p_false_rate<-p_false/((p*p-sum(signal))/2)
  
  #print results
  #print(Tstat)
  #print(result)
  #print(p_true_rate)
  #print(p_false_rate)
  #print(fdr)
  re<-c(p_true_rate,p_false_rate,fdr)
  return(re)}





## This function uses the proposed step-down procedure and calculate tpr, fpr and fdp under super diagonal setting

##input 
##W the covariates; X the observed response; D the treatment
##p dimension of brain regions; n number of subjects; m number of repeated observation
##q dimension of covariates; B bootstrap sampling size; aph significant level ; c augmentation proportion
##s0,s1: the nonzero causal effects are only on the s0th to s1th super diagonals of the correlation matrice

##output 
##p_true_rate: true positive rate
##p_false_rate: false positive rate
##fdr: false positive proportion


corr_estimate_diag_2<-function(W,X,D,p,n,m,q,B,aph,c,s0,s1){
  
  #causal effect estimation 
  
  Y.estimate_<-apply(X,3,cor,simplify = FALSE)
  Y.estimate<-array(unlist(Y.estimate_),dim=c(p,p,n))
  
  
  logistReg = glm(D ~ W + 0, family = binomial)
  beta.estimate = logistReg$coefficient
  Prob.estimate = as.vector(exp(W %*% beta.estimate) / (1 + exp(W %*% beta.estimate)))
  
  weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)
  
  tau.estimate<-apply(Y.estimate* weight[slice.index(Y.estimate, 3)],c(1,2),mean)
  
  
  
  #estimation for eta, the influence function and variance
  Fisher.Information = 0
  for (i in 1 : n){
    Fisher.Information = Fisher.Information + Prob.estimate[i] * (1 - Prob.estimate[i]) * W[i, ] %*% t(W[i, ])
  }
  Fisher.Information = Fisher.Information / n
  Theta = solve(Fisher.Information)
  
  weight.H = D * (1 - Prob.estimate) / Prob.estimate + (1 - D) * Prob.estimate / (1 - Prob.estimate)
  H = array(0, c(p, p, q))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      H[j1, j2, ] = colMeans(weight.H * Y.estimate[j1, j2, ] * W)
    }
  }
  
  eta = array(0, c(p, p, n))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      for (i in 1 : n){
        eta[j1, j2, i] = Y.estimate[j1, j2, i] * weight[i] - H[j1, j2, ] %*% Theta %*% W[i, ] * (D[i] - Prob.estimate[i])
      }
    }
  }
  
  theta.var<-apply( (eta-tau.estimate[slice.index(eta,c(1,2))])^2,c(1,2),mean)
  
  
  #standardized treatment effect
  Tstat = sqrt(n) * abs(tau.estimate) * theta.var^(-0.5)
  #set variance term, the diagonal to be zero
  
  Index = matrix(1, p, p)-diag(1,p,p)
  
  Tstat0 = Tstat * Index
  
  
  #Proposed multiple testing procedure
  z = array(0, c(p, p, B))
  for(b in 1 : B){
    g = rnorm(n)
    temp =apply( (eta-tau.estimate[slice.index(eta,c(1,2))])*g[slice.index(eta, 3)],c(1,2) , sum)
    z[, , b] = theta.var^(-0.5) * temp / sqrt(n) * Index
  }
  
  result<-matrix(0,p,p)
  z.initial = z
  Tstat.initial = Tstat0
  repeat{
    Tstat.max = max(abs(Tstat.initial))
    index.temp = which(abs(Tstat0) == Tstat.max, arr.ind = TRUE)
    z.max = apply(abs(z.initial), 3, max)
    z.max.quan = quantile(z.max, 1-aph)
    if (Tstat.max < z.max.quan) break
    for (i in 1 : dim(index.temp)[1]){
      Tstat.initial[index.temp[i, 1], index.temp[i, 2]] = 0
      z.initial[index.temp[i, 1], index.temp[i, 2], ] = 0
      result[index.temp[i,1],index.temp[i,2]]=1
    }
  }
  
  
  size<-sum(result)/2
  num_add<-floor(c*size/(1-c))
  if(num_add>=1){
    test_replace<-sort(Tstat.initial,decreasing=TRUE)
    for(i in 1:(2*num_add)){
      index.temp = which(abs(Tstat0) == test_replace[i], arr.ind = TRUE)
      for (i in 1 : dim(index.temp)[1]){
        result[index.temp[i,1],index.temp[i,2]]=1
      }
    }}
  
  
  p_true<-0
  p_false<-0
  fdr<-0
  
  
  signal=matrix(0,nrow=p,ncol=p)
  for(j in s0:s1){
    for(i in 1:(p-j+1)){
      signal[i,i+j-1]=1
      signal[i+j-1,i]=1}
  }
  
  
  
  for(i in 1:p){
    for(j in 1:p){
      if(j>i&&signal[i,j]==1&&result[i,j]==1){p_true=p_true+1}
      if(j>i&&signal[i,j]==0&&result[i,j]==1){p_false=p_false+1}
    }}
  
  if((p_true+p_false)>0){fdr<-p_false/(p_false+p_true)}else{
    fdr<-0}
  p_true_rate<-p_true/((sum(signal))/2)
  p_false_rate<-p_false/((p*p-p-sum(signal))/2)
  
  
  #print(Tstat)
  #print(result)
  #print(p_true_rate)
  #print(p_false_rate)
  #print(fdr)
  re<-c(p_true_rate,p_false_rate,fdr)
  return(re)}




## This function uses the BH procedure and calculate tpr, fpr and fdp under super diagonal setting

##input 
##W the covariates; X the observed response; D the treatment
##p dimension of brain regions; n number of subjects; m number of repeated observation
##q dimension of covariates; B bootstrap sampling size; aph significant level ; c augmentation proportion
##s0,s1: the nonzero causal effects are only on the s0th to s1th super diagonals of the correlation matrice

##output 
##p_true_rate: true positive rate
##p_false_rate: false positive rate
##fdr: false positive proportion


corr_estimate_bh_2<-function(W,X,D,p,n,m,q,B,aph,c,s0,s1){
  
  #causal effect estimation 
  
  Y.estimate_<-apply(X,3,cor,simplify = FALSE)
  Y.estimate<-array(unlist(Y.estimate_),dim=c(p,p,n))
  
  logistReg = glm(D ~ W + 0, family = binomial)
  beta.estimate = logistReg$coefficient
  Prob.estimate = as.vector(exp(W %*% beta.estimate) / (1 + exp(W %*% beta.estimate)))
  
  weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)
  
  tau.estimate<-apply(Y.estimate* weight[slice.index(Y.estimate, 3)],c(1,2),mean)
  
  
  #estimation for eta, the influence function and variance
  Fisher.Information = 0
  for (i in 1 : n){
    Fisher.Information = Fisher.Information + Prob.estimate[i] * (1 - Prob.estimate[i]) * W[i, ] %*% t(W[i, ])
  }
  Fisher.Information = Fisher.Information / n
  Theta = solve(Fisher.Information)
  weight.H = D * (1 - Prob.estimate) / Prob.estimate + (1 - D) * Prob.estimate / (1 - Prob.estimate)
  H = array(0, c(p, p, q))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      H[j1, j2, ] = colMeans(weight.H * Y.estimate[j1, j2, ] * W)
    }
  }
  
  eta = array(0, c(p, p, n))
  for (j1 in 1 : p){
    for (j2 in 1 : p){
      for (i in 1 : n){
        eta[j1, j2, i] = Y.estimate[j1, j2, i] * weight[i] - H[j1, j2, ] %*% Theta %*% W[i, ] * (D[i] - Prob.estimate[i])
      }
    }
  }
  
  
  theta.var<-apply( (eta-tau.estimate[slice.index(eta,c(1,2))])^2,c(1,2),mean)
  
  Tstat = sqrt(n) * abs(tau.estimate) * theta.var^(-0.5)
  P=2*pnorm(Tstat,0,1,lower.tail = FALSE)
  
  P1=P
  for(i in 1:p){
    for(j in 1:p){
      if(j<=i){P1[i,j]=1}
    }
  }
  test=sort(P1)
  num=1
  while(test[num]<num*aph/(p*(p-1)/2)){
    num=num+1
  }
  num=num-1
  
  
  
  result<-matrix(0,p,p)
  if(num>1){
    for(i in 1:num){
      index.temp = which(P == test[i], arr.ind = TRUE)
      for (j in 1 : dim(index.temp)[1]){
        result[index.temp[j,1],index.temp[j,2]]=1
      }
    }
  }
  
  
  p_true<-0
  p_false<-0
  fdr<-0
  
  signal=matrix(0,nrow=p,ncol=p)
  for(j in s0:s1){
    for(i in 1:(p-j+1)){
      signal[i,i+j-1]=1
      signal[i+j-1,i]=1}
  }
  
  
  
  
  for(i in 1:p){
    for(j in 1:p){
      if(j>i&&signal[i,j]==1&&result[i,j]==1){p_true=p_true+1}
      if(j>i&&signal[i,j]==0&&result[i,j]==1){p_false=p_false+1}
    }}
  
  if((p_true+p_false)>0){fdr<-p_false/(p_false+p_true)}else{
    fdr<-0}
  p_true_rate<-p_true/((sum(signal))/2)
  p_false_rate<-p_false/((p*p-p-sum(signal))/2)
  
  
  #print(Tstat)
  #print(result)
  #print(p_true_rate)
  #print(p_false_rate)
  #print(fdr)
  re<-c(p_true_rate,p_false_rate,fdr)
  return(re)}






















