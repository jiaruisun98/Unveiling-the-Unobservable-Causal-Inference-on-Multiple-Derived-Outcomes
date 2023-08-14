library(openxlsx)
source("Functions.R")

##This code is for the real case study. 

#Connection_diff1 provides the Nerwork in `Table1' and `Figure S2'
#Tau_sig provides the `Estimated Effect` in `Table1`
#CI_low and CI_high provides the 95% CI in `Table1`
#Ave-trt and Ace-cl provide Ave-trt and Ace-cl in `Table1`

##read in brian region response data
data <- data.frame()
data <- array(data=0,dim=c(175,116,79))
#79 subjects,116ROIs,175 time points

for(i in 50952:50962){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}

for(i in 50964:51003){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951-1] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}


for(i in 51006:51021){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951-3] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}

for(i in 51023:51030){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951-4] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}


for(i in 51032:51035){
  path <- paste0("NYU_00",i,"_rois_aal",".1D")
  data[,,i-50951-5] <-as.matrix(read.table(file = path, header = TRUE))[1:175,]
}

X<-data

#read in covariance
Covariance<-read.csv("phenotypic_NYU.csv")
#Covariance<-read.xlsx("COV_NYU.xlsx")
Covariance<-Covariance[which(Covariance[,"DX_GROUP"]==1),]
D<-Covariance[,"CURRENT_MED_STATUS"]


#deleting subjects with missing data
W<-Covariance[,c(5,6,8,9)]
W<-as.matrix(W)
X<-X[,,which(W[,3]!=-9999)]
D<-as.vector(D)
D<-D[which(W[,3]!=-9999)]
W<-W[which(W[,3]!=-9999),]
X<-X[,,-c(66,76)]
W<-W[-c(66,76),]
D<-D[-c(66,76)]
ID1=which(D==1)
ID0=which(D==0)

##data:
##W the covariates; X the observed response; D the treatment

#parameters
##p dimension of brain regions; n number of subjects; m number of repeated observation
##q dimension of covariates; B bootstrap sampling size; aph significant level ; c augmentation proportion
p<-116
n<-76
m<-175
q<-4
B<-1000
aph<-0.05
c<-0.1

set.seed(12345678)

#############################
##Treatment effect estimation and multiple testing procedure of proposed method
#############################

##causal effect estimation 
#estimate correlation matrix
Y.estimate_<-apply(X,3,cor,simplify = FALSE)
Y.estimate<-array(unlist(Y.estimate_),dim=c(p,p,n))

#logistic regression for propensity score
logistReg = glm(D ~ W + 0, family = binomial)
beta.estimate = logistReg$coefficient
Prob.estimate = as.vector(exp(W %*% beta.estimate) / (1 + exp(W %*% beta.estimate)))
weight = D / Prob.estimate - (1 - D) / (1 - Prob.estimate)

#IPW, treatment effect estimation
tau.estimate<-apply(Y.estimate* weight[slice.index(Y.estimate, 3)],c(1,2),mean)

weight1<-weight[ID1]
tau.estimate1<-apply(Y.estimate[, , ID1]* weight1[slice.index(Y.estimate[, , ID1], 3)],c(1,2),sum)/n
weight0<-weight[ID0]
tau.estimate0<-apply(Y.estimate[, , ID0]* weight0[slice.index(Y.estimate[, , ID0], 3)],c(1,2),sum)/n

#estimation for eta, the influence function and variance
Fisher.Information = 0
for (i in 1 : n){
  Fisher.Information = Fisher.Information + Prob.estimate[i] * (1 - Prob.estimate[i]) * W[i, ] %*% t(W[i, ])}

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
result<-matrix(0,p,p)

z = array(0, c(p, p, B))
for(b in 1 : B){
  g = rnorm(n)
  temp =apply( (eta-tau.estimate[slice.index(eta,c(1,2))])*g[slice.index(eta, 3)],c(1,2) , sum)
  z[, , b] = theta.var^(-0.5) * temp / sqrt(n) * Index
}
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

#detected signals
Connection1<-result
Connection_diff1<-which(Connection1 ==1, arr.ind = TRUE)
print(Connection_diff1)

#CI for the causal effect for detected connections
Tau_sig=tau.estimate[Connection_diff1]
theta_sig=theta.var[Connection_diff1]
CI_low=Tau_sig-sqrt(theta_sig)/sqrt(n)*qnorm(0.975,0,1)
CI_high=Tau_sig+sqrt(theta_sig)/sqrt(n)*qnorm(0.975,0,1)

#Ave-trt+Ave-cl
Ave_trt<-tau.estimate1[Connection_diff1]
Ave_cl<-tau.estimate0[Connection_diff1]

###################
##Non-causal method: two-sample t-test directly on the subject level sample correlation
###################

#parameters
##p dimension of brain regions; n number of subjects; m number of repeated observation
##q dimension of covariates; B bootstrap sampling size; aph significant level ; c augmentation proportion
p<-116
n<-76
m<-175
q<-4
B<-1000
aph<-0.05
c<-0.1


##Correlation estimation
Y.estimate_<-apply(X,3,cor,simplify = FALSE)
Y.estimate<-array(unlist(Y.estimate_),dim=c(p,p,n))

Y.estimate.y = Y.estimate[,,which(D==1)]
Y.estimate.n = Y.estimate[,,which(D==0)]

m1 = sum(D==1)
m2 = sum(D==0)

##Difference between treatment and control
Difference = apply(Y.estimate.y,c(1,2),sum)*1/sqrt(m1)-
  apply(Y.estimate.n,c(1,2),sum)*sqrt(m1)/m2
Difference = abs(Difference)

Y.estimate.y.average = apply(Y.estimate.y,c(1,2),mean)
Y.estimate.n.average = apply(Y.estimate.n,c(1,2),mean)

Index = matrix(1, p, p)-diag(1,p,p)
Difference0 = Difference * Index

##proposed procedure 
z = array(0, c(p, p, B))
for(b in 1 : B){
  g1 = rnorm(m1)
  #temp = 0
  #for (i in 1 : m1){
  #  temp = temp + g1[i] * (Y.estimate.y[, , i] - Y.estimate.y.average)/sqrt(m1)
  #}
  g2 = rnorm(m2)
  #for (i in 1 : m2){
  #  temp = temp + g2[i] * (Y.estimate.n[, , i] - Y.estimate.n.average)*sqrt(m1)/m2
  #}
  temp =apply( (Y.estimate.y-Y.estimate.y.average[slice.index(Y.estimate.y,c(1,2))])*g1[slice.index(Y.estimate.y, 3)],c(1,2) , sum)/sqrt(m1)-
    apply( (Y.estimate.n-Y.estimate.n.average[slice.index(Y.estimate.n,c(1,2))])*g2[slice.index(Y.estimate.n, 3)],c(1,2) , sum)*sqrt(m1)/m2

  z[, , b] = temp * Index
}

result<-matrix(0,p,p)

Tstat.initial = Difference0
z.initial = z

repeat{
  Tstat.max = max(abs(Tstat.initial))
  index.temp = which(abs(Difference0) == Tstat.max, arr.ind = TRUE)
  z.max = apply(abs(z.initial), 3, max)
  z.max.quan = quantile(z.max, 1-aph)
  if (Tstat.max < z.max.quan) break
  # reject = rbind(reject, index.temp)
  for (i in 1 : dim(index.temp)[1]){
    Tstat.initial[index.temp[i, 1], index.temp[i, 2]] = 0
    z.initial[index.temp[i, 1], index.temp[i, 2], ] = 0
    result[index.temp[i,1],index.temp[i,2]]=1
  }
  #  cat("Max Stat = ", Tstat.max, "quantile = ", z.max.quan, "\n")
}

#Agumentation

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


result_diff2<-which(result ==1, arr.ind = TRUE)
print(result_diff2)
##There are no significant connections

#############################################
##Treatment effect estimation and multiple testing procedure of BH method
#############################################

#estimate correlation matrix
Y.estimate_<-apply(X,3,cor,simplify = FALSE)
Y.estimate<-array(unlist(Y.estimate_),dim=c(p,p,n))

#logistic regression for propensity score
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



##Use the BH procedure instead of the proposed procedure

Tstat = sqrt(n) * abs(tau.estimate) * theta.var^(-0.5)
P=2*pnorm(Tstat,0,1,lower.tail = FALSE)

P1=P
for(i in 1:p){
  for(j in 1:p){
    if(j<=i){P1[i,j]=1} }}

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


#detected signals
Connection_bh<-result
Connection_diff_bh<-which(Connection_bh ==1, arr.ind = TRUE)
print(Connection_diff_bh)

#CI for the causal effect for detected connections
#Tau_sig=tau.estimate[Connection_diff1]
#theta_sig=theta.var[Connection_diff1]
#CI_low=Tau_sig-sqrt(theta_sig)/sqrt(n)*qnorm(0.975,0,1)
#CI_high=Tau_sig+sqrt(theta_sig)/sqrt(n)*qnorm(0.975,0,1)

#Ave-trt+Ave-cl
#Ave_trt<-tau.estimate1[Connection_diff1]
#Ave_cl<-tau.estimate0[Connection_diff1]


















