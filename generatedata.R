## Simulation

generatedata <- function(n, frailty=0, normal=TRUE){
  
  X1 = rbinom(n, 1, 0.5)
  X2 = rbinom(n, 1, 0.5)
  X3 = rbinom(n, 1, 0.5)
  X = cbind(X1, X2, X3)
  if (!normal) {
    e = log(rgamma(n, 4, 4))
  } else {
    e = rnorm(n, 0, 1)
  }
  
  # treated
  a1 = 0.01; b1 = 0.04-0.02; beta1 = c(0.2,0.2,-0.2)
  a2 = 0.02; b2 = 0.02+0.01; beta2 = c(0.2,-0.2,-0.2)
  a3 = 0.01+0.01; b3 = 0.02; beta3 = c(0.1,0.1,0.1)
  d12 = exp(-0.3); d13 = exp(0.2); d21 = exp(-0.2); d23 = exp(0.5)
  Xb1 = exp(as.numeric(X%*%beta1)+frailty*e)
  Xb2 = exp(as.numeric(X%*%beta2)+frailty*e)
  Xb3 = exp(as.numeric(X%*%beta3)+frailty*e)
  
  # first jump
  E1 = rep(0,n)
  a = a1*Xb1 + a2*Xb2 + a3*Xb3; b = b1*Xb1 + b2*Xb2 + b3*Xb3
  T1 = (-b+sqrt(b^2-2*a*log(runif(n))))/a
  prob1 = Xb1*(a1*T1+b1); prob2 = Xb2*(a2*T1+b2); prob3 = Xb3*(a3*T1+b3)
  e3 = rbinom(n, 1, prob3/(prob1+prob2+prob3))
  E1[e3==1] = 3
  e1 = rbinom(n, 1, prob1/(prob1+prob2))
  E1[E1==0] = e1[E1==0]
  E1[E1==0] = 2
  
  # second jump from "1"
  T2 = T1
  E2 = rep(0,n)
  a = a2*Xb2*d12 + a3*Xb3*d13; b = b2*Xb2*d12 + b3*Xb3*d13
  c = a/2*T1^2 + b*T1
  t2 = (-b+sqrt(b^2-2*a*(log(runif(n))-c)))/a
  prob2 = Xb2*(a2*t2+b2)*d12; prob3 = Xb3*(a3*t2+b3)*d13
  e3 = rbinom(n, 1, prob3/(prob2+prob3))
  E2[E1==1&e3==1] = 3
  E2[E1==1&e3==0] = 2
  T2[E1==1] = t2[E1==1]
  # second jump from "2"
  a = a1*Xb1*d21 + a3*Xb3*d23; b = b1*Xb1*d21 + b3*Xb3*d23
  c = a/2*T1^2 + b*T1
  t2 = (-b+sqrt(b^2-2*a*(log(runif(n))-c)))/a
  prob1 = Xb1*(a1*t2+b1)*d21; prob3 = Xb3*(a3*t2+b3)*d23
  e3 = rbinom(n, 1, prob3/(prob1+prob3))
  E2[E1==2&e3==1] = 3
  E2[E1==2&e3==0] = 1
  T2[E1==2] = t2[E1==2]
  T2[T2<T1] = T1[T2<T1] + 0.0001
  
  # third jump
  T3 = T2
  E3 = rep(0,n)
  a = a3*Xb3*d23*d13; b = b3*Xb3*d23*d13
  c = a/2*T2^2 + b*T2
  t3 = (-b+sqrt(b^2-2*a*(log(runif(n))-c)))/a
  T3[E2==1|E2==2] = t3[E2==1|E2==2]
  T2[T2<T1] = T1[T2<T1] + 0.0001
  E3[E2==1|E2==2] = 3
  
  T_1 = T1*(E1==1) + T2*(E2==1)
  T_2 = T1*(E1==2) + T2*(E2==2)
  T_3 = T1*(E1==3) + T2*(E2==3) + T3*(E3==3)
  D_1 = (E1==1)|(E2==1)
  D_2 = (E1==2)|(E2==2)
  D_3 = (E1==3)|(E2==3)|(E3==3)
  T_1[D_1==0] = T_3[D_1==0]
  T_2[D_2==0] = T_3[D_2==0]
  
  # censor
  Xc = exp(-as.numeric(X%*%c(-0.1,0,0)))
  T_c = 5 + sqrt(-20*Xc*log(runif(n)))
  T_c[T_c>15] = 15
  D_1_1 = D_1 * (T_1<=T_c)
  D_2_1 = D_2 * (T_2<=T_c)
  D_3_1 = D_3 * (T_3<=T_c)
  T_3_1 = T_3 * D_3_1 + T_c * (1-D_3_1)
  T_1_1 = T_1 * D_1_1 + T_3_1 * (1-D_1_1)
  T_2_1 = T_2 * D_2_1 + T_3_1 * (1-D_2_1)
  
  # control
  a1 = 0.01; b1 = 0.04; beta1 = c(0.2,0.2,-0.2)
  a2 = 0.02; b2 = 0.02; beta2 = c(0.2,-0.2,-0.2)
  a3 = 0.01; b3 = 0.02; beta3 = c(0.1,0.1,0.1)
  d12 = exp(-0.3); d13 = exp(0.1); d21 = exp(-0.2); d23 = exp(0.5)
  Xb1 = exp(as.numeric(X%*%beta1)+frailty*e)
  Xb2 = exp(as.numeric(X%*%beta2)+frailty*e)
  Xb3 = exp(as.numeric(X%*%beta3)+frailty*e)
  
  # first jump
  E1 = rep(0,n)
  a = a1*Xb1 + a2*Xb2 + a3*Xb3; b = b1*Xb1 + b2*Xb2 + b3*Xb3
  T1 = (-b+sqrt(b^2-2*a*log(runif(n))))/a
  prob1 = Xb1*(a1*T1+b1); prob2 = Xb2*(a2*T1+b2); prob3 = Xb3*(a3*T1+b3)
  e3 = rbinom(n, 1, prob3/(prob1+prob2+prob3))
  E1[e3==1] = 3
  e1 = rbinom(n, 1, prob1/(prob1+prob2))
  E1[E1==0] = e1[E1==0]
  E1[E1==0] = 2
  
  # second jump from "1"
  T2 = T1
  E2 = rep(0,n)
  a = a2*Xb2*d12 + a3*Xb3*d13; b = b2*Xb2*d12 + b3*Xb3*d13
  c = a/2*T1^2 + b*T1
  t2 = (-b+sqrt(b^2-2*a*(log(runif(n))-c)))/a
  prob2 = Xb2*(a2*t2+b2)*d12; prob3 = Xb3*(a3*t2+b3)*d13
  e3 = rbinom(n, 1, prob3/(prob2+prob3))
  E2[E1==1&e3==1] = 3
  E2[E1==1&e3==0] = 2
  T2[E1==1] = t2[E1==1]
  # second jump from "2"
  a = a1*Xb1*d21 + a3*Xb3*d23; b = b1*Xb1*d21 + b3*Xb3*d23
  c = a/2*T1^2 + b*T1
  t2 = (-b+sqrt(b^2-2*a*(log(runif(n))-c)))/a
  prob1 = Xb1*(a1*t2+b1)*d21; prob3 = Xb3*(a3*t2+b3)*d23
  e3 = rbinom(n, 1, prob3/(prob1+prob3))
  E2[E1==2&e3==1] = 3
  E2[E1==2&e3==0] = 1
  T2[E1==2] = t2[E1==2]
  T2[T2<T1] = T1[T2<T1] + 0.0001
  
  # third jump
  T3 = T2
  E3 = rep(0,n)
  a = a3*Xb3*d23*d13; b = b3*Xb3*d23*d13
  c = a/2*T2^2 + b*T2
  t3 = (-b+sqrt(b^2-2*a*(log(runif(n))-c)))/a
  T3[E2==1|E2==2] = t3[E2==1|E2==2]
  T2[T2<T1] = T1[T2<T1] + 0.0001
  E3[E2==1|E2==2] = 3
  
  T_1 = T1*(E1==1) + T2*(E2==1)
  T_2 = T1*(E1==2) + T2*(E2==2)
  T_3 = T1*(E1==3) + T2*(E2==3) + T3*(E3==3)
  D_1 = (E1==1)|(E2==1)
  D_2 = (E1==2)|(E2==2)
  D_3 = (E1==3)|(E2==3)|(E3==3)
  T_1[D_1==0] = T_3[D_1==0]
  T_2[D_2==0] = T_3[D_2==0]
  
  # censor
  Xc = exp(-as.numeric(X%*%c(-0.1,-0.1,0)))
  T_c = 6 + sqrt(-25*Xc*log(runif(n)))
  T_c[T_c>15] = 15
  D_1_0 = D_1 * (T_1<=T_c)
  D_2_0 = D_2 * (T_2<=T_c)
  D_3_0 = D_3 * (T_3<=T_c)
  T_3_0 = T_3 * D_3_0 + T_c * (1-D_3_0)
  T_1_0 = T_1 * D_1_0 + T_3_0 * (1-D_1_0)
  T_2_0 = T_2 * D_2_0 + T_3_0 * (1-D_2_0)
  
  # treatment and observed data
  ps = 1 / (1+exp(-as.numeric(-0.5+X%*%c(-0.3, 0.4, 0.5))))
  A = rbinom(n, 1, ps)
  T_1 = T_1_1 * A + T_1_0 * (1-A)
  T_2 = T_2_1 * A + T_2_0 * (1-A)
  T_3 = T_3_1 * A + T_3_0 * (1-A)
  D_1 = D_1_1 * A + D_1_0 * (1-A)
  D_2 = D_2_1 * A + D_2_0 * (1-A)
  D_3 = D_3_1 * A + D_3_0 * (1-A)
  # Effect targets on 1; out of 2.
  
  return(list(T1=T_1,D1=D_1,T2=T_2,D2=D_2,T3=T_3,D3=D_3,A=A,X=X))
  
}

trueCIF <- function(N=10000, frailty=0, normal=TRUE){
  tseq = seq(0,15,0.01)
  X1 = rbinom(N, 1, 0.5)
  X2 = rbinom(N, 1, 0.5)
  X3 = rbinom(N, 1, 0.5)
  X = cbind(X1, X2, X3)
  if (!normal) {
    e = log(rgamma(N, 4, 4))
  } else {
    e = rnorm(N, 0, 1)
  }
  
  # treated
  a1 = 0.01; b1 = 0.04-0.02; beta1 = c(0.2,0.2,-0.2)
  a2 = 0.02; b2 = 0.02+0.01; beta2 = c(0.2,-0.2,-0.2)
  a3 = 0.01+0.01; b3 = 0.02; beta3 = c(0.1,0.1,0.1)
  d12 = exp(-0.3); d13 = exp(0.2); d21 = exp(-0.2); d23 = exp(0.5)
  Xb1 = exp(as.numeric(X%*%beta1)+frailty*e)
  Xb2 = exp(as.numeric(X%*%beta2)+frailty*e)
  Xb3 = exp(as.numeric(X%*%beta3)+frailty*e)
  # hazards
  lam_og1 = sapply(tseq, function(t) (a1*t+b1)*Xb1)
  lam_or1 = sapply(tseq, function(t) (a2*t+b2)*Xb2)
  lam_od1 = sapply(tseq, function(t) (a3*t+b3)*Xb3)
  lam_ogr1 = lam_or1 * d12
  lam_ogd1 = lam_od1 * d13
  lam_org1 = lam_og1 * d21
  lam_ord1 = lam_od1 * d23
  lam_orgd1 = lam_ogrd1 = lam_od1 * d23 * d13
  # cif
  Lam_o1 = t(apply(lam_od1+lam_og1+lam_or1, 1, cumsum))*0.01
  Lam_og1 = t(apply(lam_ogd1+lam_ogr1, 1, cumsum))*0.01
  Lam_or1 = t(apply(lam_ord1+lam_org1, 1, cumsum))*0.01
  Lam_ogr1 = t(apply(lam_ogrd1, 1, cumsum))*0.01
  Lam_org1 = t(apply(lam_orgd1, 1, cumsum))*0.01
  dF_or1 = exp(-Lam_o1)*lam_or1
  dF_og1 = exp(-Lam_o1)*lam_og1
  dF_od1 = exp(-Lam_o1)*lam_od1
  F_or1 = t(apply(dF_or1, 1, cumsum))*0.01
  F_og1 = t(apply(dF_og1, 1, cumsum))*0.01
  F_od1 = t(apply(dF_od1, 1, cumsum))*0.01
  dF_og_1 = t(apply(dF_og1*exp(Lam_og1), 1, cumsum))*exp(-Lam_og1)*0.01
  dF_or_1 = t(apply(dF_or1*exp(Lam_or1), 1, cumsum))*exp(-Lam_or1)*0.01
  dF_ogr1 = dF_og_1*lam_ogr1
  F_ogr1 = t(apply(dF_ogr1, 1, cumsum))*0.01
  dF_org1 = dF_or_1*lam_org1
  F_org1 = t(apply(dF_org1, 1, cumsum))*0.01
  dF_ord1 = dF_or_1*lam_ord1
  F_ord1 = t(apply(dF_ord1, 1, cumsum))*0.01
  dF_ogd1 = dF_og_1*lam_ogd1
  F_ogd1 = t(apply(dF_ogd1, 1, cumsum))*0.01
  dF_org_1 = t(apply(dF_org1*exp(Lam_org1), 1, cumsum))*exp(-Lam_org1)*0.01
  dF_ogr_1 = t(apply(dF_ogr1*exp(Lam_ogr1), 1, cumsum))*exp(-Lam_ogr1)*0.01
  dF_orgd1 = dF_org_1*lam_orgd1
  F_orgd1 = t(apply(dF_orgd1, 1, cumsum))*0.01
  dF_ogrd1 = dF_ogr_1*lam_ogrd1
  F_ogrd1 = t(apply(dF_ogrd1, 1, cumsum))*0.01
  F_1 = F_od1 + F_ogd1 + F_ord1 + F_ogrd1 + F_orgd1
  
  # control
  a1 = 0.01; b1 = 0.04; beta1 = c(0.2,0.2,-0.2)
  a2 = 0.02; b2 = 0.02; beta2 = c(0.2,-0.2,-0.2)
  a3 = 0.01; b3 = 0.02; beta3 = c(0.1,0.1,0.1)
  d12 = exp(-0.3); d13 = exp(0.1); d21 = exp(-0.2); d23 = exp(0.5)
  Xb1 = exp(as.numeric(X%*%beta1)+frailty*e)
  Xb2 = exp(as.numeric(X%*%beta2)+frailty*e)
  Xb3 = exp(as.numeric(X%*%beta3)+frailty*e)
  # hazards
  lam_og0 = sapply(tseq, function(t) (a1*t+b1)*Xb1)
  lam_or0 = sapply(tseq, function(t) (a2*t+b2)*Xb2)
  lam_od0 = sapply(tseq, function(t) (a3*t+b3)*Xb3)
  lam_ogr0 = lam_or0 * d12
  lam_ogd0 = lam_od0 * d13
  lam_org0 = lam_og0 * d21
  lam_ord0 = lam_od0 * d23
  lam_orgd0 = lam_ogrd0 = lam_od0 * d23 * d13
  # cif
  Lam_o0 = t(apply(lam_od0+lam_og0+lam_or0, 1, cumsum))*0.01
  Lam_og0 = t(apply(lam_ogd0+lam_ogr0, 1, cumsum))*0.01
  Lam_or0 = t(apply(lam_ord0+lam_org0, 1, cumsum))*0.01
  Lam_ogr0 = t(apply(lam_ogrd0, 1, cumsum))*0.01
  Lam_org0 = t(apply(lam_orgd0, 1, cumsum))*0.01
  dF_or0 = exp(-Lam_o0)*lam_or0
  dF_og0 = exp(-Lam_o0)*lam_og0
  dF_od0 = exp(-Lam_o0)*lam_od0
  F_or0 = t(apply(dF_or0, 1, cumsum))*0.01
  F_og0 = t(apply(dF_og0, 1, cumsum))*0.01
  F_od0 = t(apply(dF_od0, 1, cumsum))*0.01
  dF_og_0 = t(apply(dF_og0*exp(Lam_og0), 1, cumsum))*exp(-Lam_og0)*0.01
  dF_or_0 = t(apply(dF_or0*exp(Lam_or0), 1, cumsum))*exp(-Lam_or0)*0.01
  dF_ogr0 = dF_og_0*lam_ogr0
  F_ogr0 = t(apply(dF_ogr0, 1, cumsum))*0.01
  dF_org0 = dF_or_0*lam_org0
  F_org0 = t(apply(dF_org0, 1, cumsum))*0.01
  dF_ord0 = dF_or_0*lam_ord0
  F_ord0 = t(apply(dF_ord0, 1, cumsum))*0.01
  dF_ogd0 = dF_og_0*lam_ogd0
  F_ogd0 = t(apply(dF_ogd0, 1, cumsum))*0.01
  dF_org_0 = t(apply(dF_org0*exp(Lam_org0), 1, cumsum))*exp(-Lam_org0)*0.01
  dF_ogr_0 = t(apply(dF_ogr0*exp(Lam_ogr0), 1, cumsum))*exp(-Lam_ogr0)*0.01
  dF_orgd0 = dF_org_0*lam_orgd0
  F_orgd0 = t(apply(dF_orgd0, 1, cumsum))*0.01
  dF_ogrd0 = dF_ogr_0*lam_ogrd0
  F_ogrd0 = t(apply(dF_ogrd0, 1, cumsum))*0.01
  F_0 = F_od0 + F_ogd0 + F_ord0 + F_ogrd0 + F_orgd0
  
  # conditional interventional
  lam_og_a = lam_og1
  lam_od_a = lam_od1
  lam_or_a = lam_or0
  lam_ogr_a = lam_ogr0
  lam_ogd_a = lam_ogd1
  lam_org_a = lam_org1
  lam_ord_a = lam_ord1
  lam_ogrd_a = lam_ogrd1
  lam_orgd_a = lam_orgd1
  Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))*0.01
  Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))*0.01
  Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))*0.01
  Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))*0.01
  Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))*0.01
  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_og_a = exp(-Lam_o_a)*lam_og_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_or_a = t(apply(dF_or_a, 1, cumsum))*0.01
  F_og_a = t(apply(dF_og_a, 1, cumsum))*0.01
  F_od_a = t(apply(dF_od_a, 1, cumsum))*0.01
  dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)*0.01
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)*0.01
  dF_ogr_a = dF_og__a*lam_ogr_a
  F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))*0.01
  dF_org_a = dF_or__a*lam_org_a
  F_org_a = t(apply(dF_org_a, 1, cumsum))*0.01
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))*0.01
  dF_ogd_a = dF_og__a*lam_ogd_a
  F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))*0.01
  dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)*0.01
  dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)*0.01
  dF_orgd_a = dF_org__a*lam_orgd_a
  F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))*0.01
  dF_ogrd_a = dF_ogr__a*lam_ogrd_a
  F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))*0.01
  F_cin = F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a
  
  # interventional
  prev_g1 = (F_og1-F_ogr1-F_ogd1)/(1-F_od1-F_or1-F_ogd1-F_ogr1)
  lam_or1 = lam_ogr1 = prev_g1*lam_ogr1+(1-prev_g1)*lam_or1
  prev_g0 = (F_og0-F_ogr0-F_ogd0)/(1-F_od0-F_or0-F_ogd0-F_ogr0)
  lam_or0 = lam_ogr0 = prev_g0*lam_ogr0+(1-prev_g0)*lam_or0
  
  Lam_o1 = t(apply(lam_od1+lam_og1+lam_or1, 1, cumsum))*0.01
  Lam_og1 = t(apply(lam_ogd1+lam_ogr1, 1, cumsum))*0.01
  Lam_or1 = t(apply(lam_ord1+lam_org1, 1, cumsum))*0.01
  Lam_ogr1 = t(apply(lam_ogrd1, 1, cumsum))*0.01
  Lam_org1 = t(apply(lam_orgd1, 1, cumsum))*0.01
  dF_or1 = exp(-Lam_o1)*lam_or1
  dF_og1 = exp(-Lam_o1)*lam_og1
  dF_od1 = exp(-Lam_o1)*lam_od1
  F_or1 = t(apply(dF_or1, 1, cumsum))*0.01
  F_og1 = t(apply(dF_og1, 1, cumsum))*0.01
  F_od1 = t(apply(dF_od1, 1, cumsum))*0.01
  dF_og_1 = t(apply(dF_og1*exp(Lam_og1), 1, cumsum))*exp(-Lam_og1)*0.01
  dF_or_1 = t(apply(dF_or1*exp(Lam_or1), 1, cumsum))*exp(-Lam_or1)*0.01
  dF_ogr1 = dF_og_1*lam_ogr1
  F_ogr1 = t(apply(dF_ogr1, 1, cumsum))*0.01
  dF_org1 = dF_or_1*lam_org1
  F_org1 = t(apply(dF_org1, 1, cumsum))*0.01
  dF_ord1 = dF_or_1*lam_ord1
  F_ord1 = t(apply(dF_ord1, 1, cumsum))*0.01
  dF_ogd1 = dF_og_1*lam_ogd1
  F_ogd1 = t(apply(dF_ogd1, 1, cumsum))*0.01
  dF_org_1 = t(apply(dF_org1*exp(Lam_org1), 1, cumsum))*exp(-Lam_org1)*0.01
  dF_ogr_1 = t(apply(dF_ogr1*exp(Lam_ogr1), 1, cumsum))*exp(-Lam_ogr1)*0.01
  dF_orgd1 = dF_org_1*lam_orgd1
  F_orgd1 = t(apply(dF_orgd1, 1, cumsum))*0.01
  dF_ogrd1 = dF_ogr_1*lam_ogrd1
  F_ogrd1 = t(apply(dF_ogrd1, 1, cumsum))*0.01
  F_1o = F_od1 + F_ogd1 + F_ord1 + F_ogrd1 + F_orgd1
  
  Lam_o0 = t(apply(lam_od0+lam_og0+lam_or0, 1, cumsum))*0.01
  Lam_og0 = t(apply(lam_ogd0+lam_ogr0, 1, cumsum))*0.01
  Lam_or0 = t(apply(lam_ord0+lam_org0, 1, cumsum))*0.01
  Lam_ogr0 = t(apply(lam_ogrd0, 1, cumsum))*0.01
  Lam_org0 = t(apply(lam_orgd0, 1, cumsum))*0.01
  dF_or0 = exp(-Lam_o0)*lam_or0
  dF_og0 = exp(-Lam_o0)*lam_og0
  dF_od0 = exp(-Lam_o0)*lam_od0
  F_or0 = t(apply(dF_or0, 1, cumsum))*0.01
  F_og0 = t(apply(dF_og0, 1, cumsum))*0.01
  F_od0 = t(apply(dF_od0, 1, cumsum))*0.01
  dF_og_0 = t(apply(dF_og0*exp(Lam_og0), 1, cumsum))*exp(-Lam_og0)*0.01
  dF_or_0 = t(apply(dF_or0*exp(Lam_or0), 1, cumsum))*exp(-Lam_or0)*0.01
  dF_ogr0 = dF_og_0*lam_ogr0
  F_ogr0 = t(apply(dF_ogr0, 1, cumsum))*0.01
  dF_org0 = dF_or_0*lam_org0
  F_org0 = t(apply(dF_org0, 1, cumsum))*0.01
  dF_ord0 = dF_or_0*lam_ord0
  F_ord0 = t(apply(dF_ord0, 1, cumsum))*0.01
  dF_ogd0 = dF_og_0*lam_ogd0
  F_ogd0 = t(apply(dF_ogd0, 1, cumsum))*0.01
  dF_org_0 = t(apply(dF_org0*exp(Lam_org0), 1, cumsum))*exp(-Lam_org0)*0.01
  dF_ogr_0 = t(apply(dF_ogr0*exp(Lam_ogr0), 1, cumsum))*exp(-Lam_ogr0)*0.01
  dF_orgd0 = dF_org_0*lam_orgd0
  F_orgd0 = t(apply(dF_orgd0, 1, cumsum))*0.01
  dF_ogrd0 = dF_ogr_0*lam_ogrd0
  F_ogrd0 = t(apply(dF_ogrd0, 1, cumsum))*0.01
  F_0o = F_od0 + F_ogd0 + F_ord0 + F_ogrd0 + F_orgd0
  
  lam_og_a = lam_og1
  lam_od_a = lam_od1
  lam_or_a = lam_or0
  lam_ogr_a = lam_ogr0
  lam_ogd_a = lam_ogd1
  lam_org_a = lam_org1
  lam_ord_a = lam_ord1
  lam_ogrd_a = lam_ogrd1
  lam_orgd_a = lam_orgd1
  Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))*0.01
  Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))*0.01
  Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))*0.01
  Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))*0.01
  Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))*0.01
  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_og_a = exp(-Lam_o_a)*lam_og_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_or_a = t(apply(dF_or_a, 1, cumsum))*0.01
  F_og_a = t(apply(dF_og_a, 1, cumsum))*0.01
  F_od_a = t(apply(dF_od_a, 1, cumsum))*0.01
  dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)*0.01
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)*0.01
  dF_ogr_a = dF_og__a*lam_ogr_a
  F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))*0.01
  dF_org_a = dF_or__a*lam_org_a
  F_org_a = t(apply(dF_org_a, 1, cumsum))*0.01
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))*0.01
  dF_ogd_a = dF_og__a*lam_ogd_a
  F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))*0.01
  dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)*0.01
  dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)*0.01
  dF_orgd_a = dF_org__a*lam_orgd_a
  F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))*0.01
  dF_ogrd_a = dF_ogr__a*lam_ogrd_a
  F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))*0.01
  F_int = F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a
  
  return(list(tt=tseq, F_1=colMeans(F_1), F_0=colMeans(F_0), F_cin=colMeans(F_cin),
              F_1o=colMeans(F_1o), F_0o=colMeans(F_0o), F_int=colMeans(F_int)))
}

matchy = function(x,y,newx){
  newy = sapply(newx, function(t) y[which(x==t)])
  newy = as.numeric(newy)
  newy[is.na(newy)] = 0
  return(newy)
}

matcht = function(x,y,newx){
  newy = sapply(newx, function(t) y[max(which(x<=t))])
  newy[is.na(newy)] = 0
  return(newy)
}