#setwd('./')
source('generatedata.R')
library(xtable)
xseq=1:15
true = trueCIF(N=10000, frailty=1, normal=TRUE)
F_1.t = matcht(true$tt, true$F_1, xseq)
F_0.t = matcht(true$tt, true$F_0, xseq)
F_1o.t = matcht(true$tt, true$F_1o, xseq)
F_0o.t = matcht(true$tt, true$F_0o, xseq)
F_int.t = matcht(true$tt, true$F_int, xseq)
F_cin.t = matcht(true$tt, true$F_cin, xseq)
F_ide.t = F_int.t - F_0o.t
F_iie.t = F_1o.t - F_int.t
F_cde.t = F_cin.t - F_0.t
F_cie.t = F_1.t - F_cin.t
F_1.np = F_0.np = F_1o.np = F_0o.np = F_int.np = F_cin.np = NULL
F_1.fr = F_0.fr = F_1o.fr = F_0o.fr = F_int.fr = F_cin.fr = NULL
se_1.np = se_0.np = se_1o.np = se_0o.np = se_int.np = se_cin.np = NULL
se_1.fr = se_0.fr = se_1o.fr = se_0o.fr = se_int.fr = se_cin.fr = NULL
F_ide.np = F_iie.np = F_cde.np = F_cie.np = F_ide.fr = F_iie.fr = F_cde.fr = F_cie.fr = NULL
se_ide.np = se_iie.np = se_cde.np = se_cie.np = se_ide.fr = se_iie.fr = se_cde.fr = se_cie.fr = NULL
sig.fr = se_sig.fr = time.np = time.fr = NULL

for (lid in 1:1000){
  load(paste0('res2_600/cif',lid,'.Rdata'))
  F_1.np = rbind(F_1.np, res_np[1,])
  F_0.np = rbind(F_0.np, res_np[2,])
  F_1o.np = rbind(F_1o.np, res_np[3,])
  F_0o.np = rbind(F_0o.np, res_np[4,])
  F_int.np = rbind(F_int.np, res_np[5,])
  F_cin.np = rbind(F_cin.np, res_np[6,])
  F_ide.np = rbind(F_ide.np, res_np[7,])
  F_iie.np = rbind(F_iie.np, res_np[8,])
  F_cde.np = rbind(F_cde.np, res_np[9,])
  F_cie.np = rbind(F_cie.np, res_np[10,])
  se_1.np = rbind(se_1.np, se_np[1,])
  se_0.np = rbind(se_0.np, se_np[2,])
  se_1o.np = rbind(se_1o.np, se_np[3,])
  se_0o.np = rbind(se_0o.np, se_np[4,])
  se_int.np = rbind(se_int.np, se_np[5,])
  se_cin.np = rbind(se_cin.np, se_np[6,])
  se_ide.np = rbind(se_ide.np, se_np[7,])
  se_iie.np = rbind(se_iie.np, se_np[8,])
  se_cde.np = rbind(se_cde.np, se_np[9,])
  se_cie.np = rbind(se_cie.np, se_np[10,])
  time.np = append(time.np, as.numeric(time_np))
  time.fr = append(time.fr, as.numeric(time_fr))
  
  F_1.fr = rbind(F_1.fr, res_fr[1,])
  F_0.fr = rbind(F_0.fr, res_fr[2,])
  F_1o.fr = rbind(F_1o.fr, res_fr[3,])
  F_0o.fr = rbind(F_0o.fr, res_fr[4,])
  F_int.fr = rbind(F_int.fr, res_fr[5,])
  F_cin.fr = rbind(F_cin.fr, res_fr[6,])
  F_ide.fr = rbind(F_ide.fr, res_fr[7,])
  F_iie.fr = rbind(F_iie.fr, res_fr[8,])
  F_cde.fr = rbind(F_cde.fr, res_fr[9,])
  F_cie.fr = rbind(F_cie.fr, res_fr[10,])
  se_1.fr = rbind(se_1.fr, se_fr[1,])
  se_0.fr = rbind(se_0.fr, se_fr[2,])
  se_1o.fr = rbind(se_1o.fr, se_fr[3,])
  se_0o.fr = rbind(se_0o.fr, se_fr[4,])
  se_int.fr = rbind(se_int.fr, se_fr[5,])
  se_cin.fr = rbind(se_cin.fr, se_fr[6,])
  se_ide.fr = rbind(se_ide.fr, se_fr[7,])
  se_iie.fr = rbind(se_iie.fr, se_fr[8,])
  se_cde.fr = rbind(se_cde.fr, se_fr[9,])
  se_cie.fr = rbind(se_cie.fr, se_fr[10,])
  sig.fr = rbind(sig.fr, sigma)
  se_sig.fr = rbind(se_sig.fr, se_sig)
}

bias_F1.np = (colMeans(F_1.np, na.rm=TRUE) - F_1.t) *100
bias_F0.np = (colMeans(F_0.np, na.rm=TRUE) - F_0.t) *100
bias_F1o.np = (colMeans(F_1o.np, na.rm=TRUE) - F_1o.t) *100
bias_F0o.np = (colMeans(F_0o.np, na.rm=TRUE) - F_0o.t) *100
bias_Fint.np = (colMeans(F_int.np, na.rm=TRUE) - F_int.t) *100
bias_Fcin.np = (colMeans(F_cin.np, na.rm=TRUE) - F_cin.t) *100
bias_Fide.np = (colMeans(F_ide.np, na.rm=TRUE) - F_ide.t) *100
bias_Fiie.np = (colMeans(F_iie.np, na.rm=TRUE) - F_iie.t) *100
bias_Fcde.np = (colMeans(F_cde.np, na.rm=TRUE) - F_cde.t) *100
bias_Fcie.np = (colMeans(F_cie.np, na.rm=TRUE) - F_cie.t) *100
bias_F1.fr = (colMeans(F_1.fr, na.rm=TRUE) - F_1.t) *100
bias_F0.fr = (colMeans(F_0.fr, na.rm=TRUE) - F_0.t) *100
bias_F1o.fr = (colMeans(F_1o.fr, na.rm=TRUE) - F_1o.t) *100
bias_F0o.fr = (colMeans(F_0o.fr, na.rm=TRUE) - F_0o.t) *100
bias_Fint.fr = (colMeans(F_int.fr, na.rm=TRUE) - F_int.t) *100
bias_Fcin.fr = (colMeans(F_cin.fr, na.rm=TRUE) - F_cin.t) *100
bias_Fide.fr = (colMeans(F_ide.fr, na.rm=TRUE) - F_ide.t) *100
bias_Fiie.fr = (colMeans(F_iie.fr, na.rm=TRUE) - F_iie.t) *100
bias_Fcde.fr = (colMeans(F_cde.fr, na.rm=TRUE) - F_cde.t) *100
bias_Fcie.fr = (colMeans(F_cie.fr, na.rm=TRUE) - F_cie.t) *100

sd_F1.np = apply(F_1.np,2,sd, na.rm=TRUE) *100
sd_F0.np = apply(F_0.np,2,sd, na.rm=TRUE) *100
sd_F1o.np = apply(F_1o.np,2,sd, na.rm=TRUE) *100
sd_F0o.np = apply(F_0o.np,2,sd, na.rm=TRUE) *100
sd_Fint.np = apply(F_int.np,2,sd, na.rm=TRUE) *100
sd_Fcin.np = apply(F_cin.np,2,sd, na.rm=TRUE) *100
sd_Fide.np = apply(F_ide.np,2,sd, na.rm=TRUE) *100
sd_Fiie.np = apply(F_iie.np,2,sd, na.rm=TRUE) *100
sd_Fcde.np = apply(F_cde.np,2,sd, na.rm=TRUE) *100
sd_Fcie.np = apply(F_cie.np,2,sd, na.rm=TRUE) *100
sd_F1.fr = apply(F_1.fr,2,sd, na.rm=TRUE) *100
sd_F0.fr = apply(F_0.fr,2,sd, na.rm=TRUE) *100
sd_F1o.fr = apply(F_1o.fr,2,sd, na.rm=TRUE) *100
sd_F0o.fr = apply(F_0o.fr,2,sd, na.rm=TRUE) *100
sd_Fint.fr = apply(F_int.fr,2,sd, na.rm=TRUE) *100
sd_Fcin.fr = apply(F_cin.fr,2,sd, na.rm=TRUE) *100
sd_Fide.fr = apply(F_ide.fr,2,sd, na.rm=TRUE) *100
sd_Fiie.fr = apply(F_iie.fr,2,sd, na.rm=TRUE) *100
sd_Fcde.fr = apply(F_cde.fr,2,sd, na.rm=TRUE) *100
sd_Fcie.fr = apply(F_cie.fr,2,sd, na.rm=TRUE) *100

se_F1.np = apply(se_1.np,2,mean, na.rm=TRUE) *100
se_F0.np = apply(se_0.np,2,mean, na.rm=TRUE) *100
se_F1o.np = apply(se_1o.np,2,mean, na.rm=TRUE) *100
se_F0o.np = apply(se_0o.np,2,mean, na.rm=TRUE) *100
se_Fint.np = apply(se_int.np,2,mean, na.rm=TRUE) *100
se_Fcin.np = apply(se_cin.np,2,mean, na.rm=TRUE) *100
se_Fide.np = apply(se_ide.np,2,mean, na.rm=TRUE) *100
se_Fiie.np = apply(se_iie.np,2,mean, na.rm=TRUE) *100
se_Fcde.np = apply(se_cde.np,2,mean, na.rm=TRUE) *100
se_Fcie.np = apply(se_cie.np,2,mean, na.rm=TRUE) *100
se_F1.fr = apply(se_1.fr,2,mean, na.rm=TRUE) *100
se_F0.fr = apply(se_0.fr,2,mean, na.rm=TRUE) *100
se_F1o.fr = apply(se_1o.fr,2,mean, na.rm=TRUE) *100
se_F0o.fr = apply(se_0o.fr,2,mean, na.rm=TRUE) *100
se_Fint.fr = apply(se_int.fr,2,mean, na.rm=TRUE) *100
se_Fcin.fr = apply(se_cin.fr,2,mean, na.rm=TRUE) *100
se_Fide.fr = apply(se_ide.fr,2,mean, na.rm=TRUE) *100
se_Fiie.fr = apply(se_iie.fr,2,mean, na.rm=TRUE) *100
se_Fcde.fr = apply(se_cde.fr,2,mean, na.rm=TRUE) *100
se_Fcie.fr = apply(se_cie.fr,2,mean, na.rm=TRUE) *100

cv_F1.np = apply(t(F_1.np-1.96*se_1.np)<=F_1.t & t(F_1.np+1.96*se_1.np)>=F_1.t,1,mean, na.rm=TRUE)
cv_F0.np = apply(t(F_0.np-1.96*se_0.np)<=F_0.t & t(F_0.np+1.96*se_0.np)>=F_0.t,1,mean, na.rm=TRUE)
cv_F1o.np = apply(t(F_1o.np-1.96*se_1o.np)<=F_1o.t & t(F_1o.np+1.96*se_1o.np)>=F_1o.t,1,mean, na.rm=TRUE)
cv_F0o.np = apply(t(F_0o.np-1.96*se_0o.np)<=F_0o.t & t(F_0o.np+1.96*se_0o.np)>=F_0o.t,1,mean, na.rm=TRUE)
cv_Fint.np = apply(t(F_int.np-1.96*se_int.np)<=F_int.t & t(F_int.np+1.96*se_int.np)>=F_int.t,1,mean, na.rm=TRUE)
cv_Fcin.np = apply(t(F_cin.np-1.96*se_cin.np)<=F_cin.t & t(F_cin.np+1.96*se_cin.np)>=F_cin.t,1,mean, na.rm=TRUE)
cv_Fide.np = apply(t(F_ide.np-1.96*se_ide.np)<=F_ide.t & t(F_ide.np+1.96*se_ide.np)>=F_ide.t,1,mean, na.rm=TRUE)
cv_Fiie.np = apply(t(F_iie.np-1.96*se_iie.np)<=F_iie.t & t(F_iie.np+1.96*se_iie.np)>=F_iie.t,1,mean, na.rm=TRUE)
cv_Fcde.np = apply(t(F_cde.np-1.96*se_cde.np)<=F_cde.t & t(F_cde.np+1.96*se_cde.np)>=F_cde.t,1,mean, na.rm=TRUE)
cv_Fcie.np = apply(t(F_cie.np-1.96*se_cie.np)<=F_cie.t & t(F_cie.np+1.96*se_cie.np)>=F_cie.t,1,mean, na.rm=TRUE)
cv_F1.fr = apply(t(F_1.fr-1.96*se_1.fr)<=F_1.t & t(F_1.fr+1.96*se_1.fr)>=F_1.t,1,mean, na.rm=TRUE)
cv_F0.fr = apply(t(F_0.fr-1.96*se_0.fr)<=F_0.t & t(F_0.fr+1.96*se_0.fr)>=F_0.t,1,mean, na.rm=TRUE)
cv_F1o.fr = apply(t(F_1o.fr-1.96*se_1o.fr)<=F_1o.t & t(F_1o.fr+1.96*se_1o.fr)>=F_1o.t,1,mean, na.rm=TRUE)
cv_F0o.fr = apply(t(F_0o.fr-1.96*se_0o.fr)<=F_0o.t & t(F_0o.fr+1.96*se_0o.fr)>=F_0o.t,1,mean, na.rm=TRUE)
cv_Fint.fr = apply(t(F_int.fr-1.96*se_int.fr)<=F_int.t & t(F_int.fr+1.96*se_int.fr)>=F_int.t,1,mean, na.rm=TRUE)
cv_Fcin.fr = apply(t(F_cin.fr-1.96*se_cin.fr)<=F_cin.t & t(F_cin.fr+1.96*se_cin.fr)>=F_cin.t,1,mean, na.rm=TRUE)
cv_Fide.fr = apply(t(F_ide.fr-1.96*se_ide.fr)<=F_ide.t & t(F_ide.fr+1.96*se_ide.fr)>=F_ide.t,1,mean, na.rm=TRUE)
cv_Fiie.fr = apply(t(F_iie.fr-1.96*se_iie.fr)<=F_iie.t & t(F_iie.fr+1.96*se_iie.fr)>=F_iie.t,1,mean, na.rm=TRUE)
cv_Fcde.fr = apply(t(F_cde.fr-1.96*se_cde.fr)<=F_cde.t & t(F_cde.fr+1.96*se_cde.fr)>=F_cde.t,1,mean, na.rm=TRUE)
cv_Fcie.fr = apply(t(F_cie.fr-1.96*se_cie.fr)<=F_cie.t & t(F_cie.fr+1.96*se_cie.fr)>=F_cie.t,1,mean, na.rm=TRUE)

## Plot: interventional CIFs
plot(0:12, c(0,F_1o.t)[1:13], type='l', xlim=c(0,12), ylim=c(0,1), lwd=1.5,
     xlab='Time', ylab='Cumulative incidence', main='Setting 2: Estimated and true CIFs')
points(0:12, c(0,colMeans(F_1o.np))[1:13], type='l', lty=2, lwd=1.5,col='brown')
points(0:12, c(0,colMeans(F_1o.np))[1:13], pch=1, col='brown')
points(0:12, c(0,colMeans(F_1o.fr))[1:13], type='l', lty=3, lwd=1.5, col='darkcyan')
points(0:12, c(0,colMeans(F_1o.fr))[1:13], pch=3, col='darkcyan')
legend('topleft', cex=0.9, legend=c('True F(t;1,1)','No-frailty model','Frailty model'), lty=c(1,2,3),
       lwd=c(1.5,1.5,1.5),pch=c(NA,1,3),col=c('black','brown','darkcyan'))

plot(0:12, c(0,F_0o.t)[1:13], type='l', xlim=c(0,12), ylim=c(0,1), lwd=1.5,
     xlab='Time', ylab='Cumulative incidence', main='Setting 2: Estimated and true CIFs')
points(0:12, c(0,colMeans(F_0o.np))[1:13], type='l', lty=2, lwd=1.5,col='brown')
points(0:12, c(0,colMeans(F_0o.np))[1:13], pch=1, col='brown')
points(0:12, c(0,colMeans(F_0o.fr))[1:13], type='l', lty=3, lwd=1.5, col='darkcyan')
points(0:12, c(0,colMeans(F_0o.fr))[1:13], pch=3, col='darkcyan')
legend('topleft', cex=0.9, legend=c('True F(t;0,0)','No-frailty model','Frailty model'), lty=c(1,2,3),
       lwd=c(1.5,1.5,1.5),pch=c(NA,1,3),col=c('black','brown','darkcyan'))

plot(0:12, c(0,F_int.t)[1:13], type='l', xlim=c(0,12), ylim=c(0,1), lwd=1.5,
     xlab='Time', ylab='Cumulative incidence', main='Setting 2: Estimated and true CIFs')
points(0:12, c(0,colMeans(F_int.np))[1:13], type='l', lty=2, lwd=1.5,col='brown')
points(0:12, c(0,colMeans(F_int.np))[1:13], pch=1, col='brown')
points(0:12, c(0,colMeans(F_int.fr))[1:13], type='l', lty=3, lwd=1.5, col='darkcyan')
points(0:12, c(0,colMeans(F_int.fr))[1:13], pch=3, col='darkcyan')
legend('topleft', cex=0.9, legend=c('True F(t;0,1)','No-frailty model','Frailty model'), lty=c(1,2,3),
       lwd=c(1.5,1.5,1.5),pch=c(NA,1,3),col=c('black','brown','darkcyan'))

plot(0:12, c(0,F_1.t)[1:13], type='l', xlim=c(0,12), ylim=c(0,1), lwd=1.5,
     xlab='Time', ylab='Cumulative incidence', main='Setting 2: Estimated and true CIFs')
points(0:12, c(0,colMeans(F_1.np))[1:13], type='l', lty=2, lwd=1.5,col='brown')
points(0:12, c(0,colMeans(F_1.np))[1:13], pch=1, col='brown')
points(0:12, c(0,colMeans(F_1.fr))[1:13], type='l', lty=3, lwd=1.5, col='darkcyan')
points(0:12, c(0,colMeans(F_1.fr))[1:13], pch=3, col='darkcyan')
legend('topleft', cex=0.9, legend=c('True F(t;1)','No-frailty model','Frailty model'), lty=c(1,2,3),
       lwd=c(1.5,1.5,1.5),pch=c(NA,1,3),col=c('black','brown','darkcyan'))

plot(0:12, c(0,F_0.t)[1:13], type='l', xlim=c(0,12), ylim=c(0,1), lwd=1.5,
     xlab='Time', ylab='Cumulative incidence', main='Setting 2: Estimated and true CIFs')
points(0:12, c(0,colMeans(F_0.np))[1:13], type='l', lty=2, lwd=1.5,col='brown')
points(0:12, c(0,colMeans(F_0.np))[1:13], pch=1, col='brown')
points(0:12, c(0,colMeans(F_0.fr))[1:13], type='l', lty=3, lwd=1.5, col='darkcyan')
points(0:12, c(0,colMeans(F_0.fr))[1:13], pch=3, col='darkcyan')
legend('topleft', cex=0.9, legend=c('True F(t;0)','No-frailty model','Frailty model'), lty=c(1,2,3),
       lwd=c(1.5,1.5,1.5),pch=c(NA,1,3),col=c('black','brown','darkcyan'))

plot(0:12, c(0,F_cin.t)[1:13], type='l', xlim=c(0,12), ylim=c(0,1), lwd=1.5,
     xlab='Time', ylab='Cumulative incidence', main='Setting 2: Estimated and true CIFs')
points(0:12, c(0,colMeans(F_cin.np))[1:13], type='l', lty=2, lwd=1.5,col='brown')
points(0:12, c(0,colMeans(F_cin.np))[1:13], pch=1, col='brown')
points(0:12, c(0,colMeans(F_cin.fr))[1:13], type='l', lty=3, lwd=1.5, col='darkcyan')
points(0:12, c(0,colMeans(F_cin.fr))[1:13], pch=3, col='darkcyan')
legend('topleft', cex=0.9, legend=c(expression(True~F[CI]*'(t;0,1)'),'No-frailty model','Frailty model'), lty=c(1,2,3),
       lwd=c(1.5,1.5,1.5),pch=c(NA,1,3),col=c('black','brown','darkcyan'))


table.np = rbind(bias_F1.np, sd_F1.np, se_F1.np, cv_F1.np, 
                 bias_F0.np, sd_F0.np, se_F0.np, cv_F0.np, 
                 bias_F1o.np, sd_F1o.np, se_F1o.np, cv_F1o.np, 
                 bias_F0o.np, sd_F0o.np, se_F0o.np, cv_F0o.np, 
                 bias_Fint.np, sd_Fint.np, se_Fint.np, cv_Fint.np, 
                 bias_Fcin.np, sd_Fcin.np, se_Fcin.np, cv_Fcin.np,
                 bias_Fide.np, sd_Fide.np, se_Fide.np, cv_Fide.np, 
                 bias_Fiie.np, sd_Fiie.np, se_Fiie.np, cv_Fiie.np,
                 bias_Fcde.np, sd_Fcde.np, se_Fcde.np, cv_Fcde.np, 
                 bias_Fcie.np, sd_Fcie.np, se_Fcie.np, cv_Fcie.np
                 )
table.fr = rbind(bias_F1.fr, sd_F1.fr, se_F1.fr, cv_F1.fr, 
                 bias_F0.fr, sd_F0.fr, se_F0.fr, cv_F0.fr, 
                 bias_F1o.fr, sd_F1o.fr, se_F1o.fr, cv_F1o.fr, 
                 bias_F0o.fr, sd_F0o.fr, se_F0o.fr, cv_F0o.fr, 
                 bias_Fint.fr, sd_Fint.fr, se_Fint.fr, cv_Fint.fr, 
                 bias_Fcin.fr, sd_Fcin.fr, se_Fcin.fr, cv_Fcin.fr,
                 bias_Fide.fr, sd_Fide.fr, se_Fide.fr, cv_Fide.fr, 
                 bias_Fiie.fr, sd_Fiie.fr, se_Fiie.fr, cv_Fiie.fr,
                 bias_Fcde.fr, sd_Fcde.fr, se_Fcde.fr, cv_Fcde.fr, 
                 bias_Fcie.fr, sd_Fcie.fr, se_Fcie.fr, cv_Fcie.fr
                 )
table.np = round(table.np[,c(2,4,6,8,10)], 3)
table.fr = round(table.fr[,c(2,4,6,8,10)], 3)
xtable(cbind(table.np,table.fr), digits=3, caption='Bias, SD, SE, and CP')
