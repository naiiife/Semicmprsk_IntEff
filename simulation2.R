## Simulation

library(survival)
library(coxme)
source('generatedata.R')
source('fit_p.R')

lid = Sys.getenv("SLURM_ARRAY_TASK_ID")
if (lid=="") lid=0
set.seed(lid)

n = 600
dat = generatedata(n, frailty=1)
A = dat$A
X = dat$X
Tg = dat$T1
Tr = dat$T2
Td = dat$T3
Dg = dat$D1
Dr = dat$D2
Dd = dat$D3
xseq = 1:15

# npmle fit
fit = RandInt(Tg,Dg,Tr,Dr,Td,Dd,A,X)
res_np = fit$res_np
time_np = fit$time_np

# frailty fit
fit = RandInt_fr(Tg,Dg,Tr,Dr,Td,Dd,A,X)
res_fr = fit$res_fr
time_fr = fit$time_fr
sigma = fit$sigma

## bootstrap
B = 200
res_np.b = res_fr.b = array(NA, dim=c(10,15,B))
sigma.b = NULL
for (b in 1:B){
ss1 = sample(which(dat$A==1),replace=TRUE)
ss0 = sample(which(dat$A==0),replace=TRUE)
ss = c(ss1,ss0)
A = dat$A[ss]
X = as.matrix(dat$X[ss,])
Tg=dat$T1[ss];Dg=dat$D1[ss];Tr=dat$T2[ss]
Dr=dat$D2[ss];Td=dat$T3[ss];Dd=dat$D3[ss]
fit.b = try(RandInt(Tg,Dg,Tr,Dr,Td,Dd,A,X))
if (!'try-error'%in%class(fit.b)){
res_np.b[,,b] = fit.b$res_np
}
fit.b = try(RandInt_fr(Tg,Dg,Tr,Dr,Td,Dd,A,X))
if (!'try-error'%in%class(fit.b)){
res_fr.b[,,b] = fit.b$res_fr
sigma.b = rbind(sigma.b, fit$sigma)
}
}
se_np = apply(res_np.b, c(1,2), sd, na.rm=TRUE)
se_fr = apply(res_fr.b, c(1,2), sd, na.rm=TRUE)
se_sig = apply(sigma.b, 2, sd, na.rm=TRUE)

save(res_np, se_np, time_np, res_fr, sigma, se_fr, se_sig, time_fr, 
     file=paste0('cif',lid,'.Rdata'))
