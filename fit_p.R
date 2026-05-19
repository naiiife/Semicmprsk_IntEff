RandInt <- function(Tg,Dg,Tr,Dr,Td,Dd,A,X) {
  
  xseq = 1:15
  time.start = Sys.time()
  tt = sort(unique(c(Td[Dd==1],Tr[Dr==1],Tg[Dg==1])))
  l = length(tt)
  X = as.matrix(X)
  p = ncol(X)
  n = nrow(X)
  
  df = data.frame(id=1:n, Td=Td, Dd=Dd, Tr=Tr, Dr=Dr, Tg=Tg, Dg=Dg, X=X, A=A)
  xvars = grep("^X", names(df), value=TRUE)
  mg_d = tmerge(data1=df, data2=df, id=id, death=event(Td,Dd))
  mg_d = tmerge(mg_d, df, id=id, Gevent=tdc(Tg), Revent=tdc(Tr))
  df = data.frame(id=1:n, Tr=Tr, Dr=Dr, Tg=Tg, Dg=Dg, X=X, A=A)
  df$Tr[df$Tr>df$Tg] = df$Tg[df$Tr>df$Tg]; df$Dr[df$Tr>df$Tg] = 0
  mg_g = tmerge(data1=df, data2=df, id=id, gvhd=event(df$Tg,df$Dg))
  mg_g = tmerge(mg_g, df, id=id, Revent=tdc(df$Tr))
  df = data.frame(id=1:n, Tr=Tr, Dr=Dr, Tg=Tg, Dg=Dg, X=X, A=A)
  df$Tg[df$Tg>df$Tr] = df$Tr[df$Tg>df$Tr]; df$Dg[df$Tg>df$Tr] = 0
  mg_r = tmerge(data1=df, data2=df, id=id, relapse=event(df$Tr,df$Dr))
  mg_r = tmerge(mg_r, df, id=id, Gevent=tdc(df$Tg))
  mg_d = data.frame(id=mg_d$id,X_g=mg_d[,xvars]*0,X_r=mg_d[,xvars]*0,X_d=mg_d[,xvars],A=mg_d$A,
                    tstart=mg_d$tstart,tstop=mg_d$tstop,
                    event=ifelse(mg_d$death==1,'death','censoring'),
                    Gevent_d=mg_d$Gevent,Revent_d=mg_d$Revent,
                    Gevent_r=0,Revent_g=0,terevent='death')
  mg_g = data.frame(id=mg_g$id,X_g=mg_g[,xvars],X_r=mg_g[,xvars]*0,X_d=mg_g[,xvars]*0,A=mg_g$A,
                    tstart=mg_g$tstart,tstop=mg_g$tstop,
                    event=ifelse(mg_g$gvhd==1,'gvhd','censoring'),
                    Gevent_d=0,Revent_d=0,
                    Gevent_r=0,Revent_g=mg_g$Revent,terevent='gvhd')
  mg_r = data.frame(id=mg_r$id,X_g=mg_r[,xvars]*0,X_r=mg_r[,xvars],X_d=mg_r[,xvars]*0,A=mg_r$A,
                    tstart=mg_r$tstart,tstop=mg_r$tstop,
                    event=ifelse(mg_r$relapse==1,'relapse','censoring'),
                    Gevent_d=0,Revent_d=0,
                    Gevent_r=mg_r$Gevent,Revent_g=0,terevent='relapse')
  mg = rbind(mg_d,mg_g,mg_r)
  mg = mg[order(mg$id),]
  xvars = grep("^X", names(mg), value=TRUE)
  cox_formula = reformulate(c("Gevent_r","Gevent_d","Revent_g","Revent_d",xvars,
                              "strata(terevent)"), 
                            response="Surv(tstart,tstop,event!='censoring')")
  fit1 = coxph(cox_formula, data=mg, subset=(A==1))
  fit0 = coxph(cox_formula, data=mg, subset=(A==0))
  coefg = rbind(cbind(fit1$coefficients,sqrt(diag(vcov(fit1))),
                      fit0$coefficients,sqrt(diag(vcov(fit0))))[c(4+1:p,3),],
                Gevent=NA)
  coefr = rbind(cbind(fit1$coefficients,sqrt(diag(vcov(fit1))),
                      fit0$coefficients,sqrt(diag(vcov(fit0))))[c(4+p+1:p,1),],
                Revent=NA)
  coefd = cbind(fit1$coefficients,sqrt(diag(vcov(fit1))),
                fit0$coefficients,sqrt(diag(vcov(fit0))))[c(4+2*p+1:p,4,2),]
  coefs = round(cbind(coefg,coefr[c(1:p,p+2,p+1),],coefd), 3)
  rownames(coefs) = c(grep("^X", names(df), value=TRUE),'Revent','Gevent')
  colnames(coefs) = paste(rep(c('G','R','D'),each=4), 
                          rep(c('est','se'),6), rep(c('Z=1','Z=0'),each=2))
  
  F_1 = F_0 = F_1o = F_0o = Fd_spe = Fd_int = Fd_cin = 0
  Xbd1 = as.numeric(X%*%coefd[1:p,1])
  delta_gd1 = coefd[p+2,1]; delta_rd1 = coefd[p+1,1]
  Xbg1 = as.numeric(X%*%coefg[1:p,1])
  delta_rg1 = coefg[p+1,1]
  Xbr1 = as.numeric(X%*%coefr[1:p,1])
  delta_gr1 = coefr[p+1,1]
  Xbd0 = as.numeric(X%*%coefd[1:p,3])
  delta_gd0 = coefd[p+2,3]; delta_rd0 = coefd[p+1,3]
  Xbg0 = as.numeric(X%*%coefg[1:p,3])
  delta_rg0 = coefg[p+1,3]
  Xbr0 = as.numeric(X%*%coefr[1:p,3])
  delta_gr0 = coefr[p+1,3]
  lam1 = basehaz(fit1)
  lamd1 = lam1[lam1$strata=='death',]
  lam_d1 = diff(c(0,matcht(lamd1$time, lamd1$hazard, tt)))
  lamg1 = lam1[lam1$strata=='gvhd',]
  lam_g1 = diff(c(0,matcht(lamg1$time, lamg1$hazard, tt)))
  lamr1 = lam1[lam1$strata=='relapse',]
  lam_r1 = diff(c(0,matcht(lamr1$time, lamr1$hazard, tt)))
  lam0 = basehaz(fit0)
  lamd0 = lam0[lam0$strata=='death',]
  lam_d0 = diff(c(0,matcht(lamd0$time, lamd0$hazard, tt)))
  lamg0 = lam0[lam0$strata=='gvhd',]
  lam_g0 = diff(c(0,matcht(lamg0$time, lamg0$hazard, tt)))
  lamr0 = lam0[lam0$strata=='relapse',]
  lam_r0 = diff(c(0,matcht(lamr0$time, lamr0$hazard, tt)))
  
  lam_od1 = sapply(1:l, function(t) lam_d1[t]*exp(Xbd1))
  lam_ogd1 = lam_od1 * exp(delta_gd1)
  lam_ord1 = lam_od1 * exp(delta_rd1)
  lam_ogrd1 = lam_orgd1 = lam_od1 * exp(delta_gd1+delta_rd1)
  lam_og1 = sapply(1:l, function(t) lam_g1[t]*exp(Xbg1))
  lam_org1 = lam_og1 *exp(delta_rg1)
  lam_or1 = sapply(1:l, function(t) lam_r1[t]*exp(Xbr1))
  lam_ogr1 = lam_or1 * exp(delta_gr1)
  lam_od0 = sapply(1:l, function(t) lam_d0[t]*exp(Xbd0))
  lam_ogd0 = lam_od0 * exp(delta_gd0)
  lam_ord0 = lam_od0 * exp(delta_rd0)
  lam_ogrd0 = lam_orgd0 = lam_od0 * exp(delta_gd0+delta_rd0)
  lam_og0 = sapply(1:l, function(t) lam_g0[t]*exp(Xbg0))
  lam_org0 = lam_og0 *exp(delta_rg0)
  lam_or0 = sapply(1:l, function(t) lam_r0[t]*exp(Xbr0))
  lam_ogr0 = lam_or0 * exp(delta_gr0)
  
  # potential incidence
  #prev_g = (F_og_1-F_ogr_1-F_ogd_1)/(1-F_od_1-F_or_1-F_ogd_1-F_ogr_1)
  #lam_or1 = lam_ogr1 = prev_g*lam_ogr1+(1-prev_g)*lam_or1
  Lam_o_1 = t(apply(lam_od1+lam_og1+lam_or1, 1, cumsum))
  Lam_og_1 = t(apply(lam_ogd1+lam_ogr1, 1, cumsum))
  Lam_or_1 = t(apply(lam_ord1+lam_org1, 1, cumsum))
  Lam_ogr_1 = t(apply(lam_ogrd1, 1, cumsum))
  Lam_org_1 = t(apply(lam_orgd1, 1, cumsum))
  dF_or_1 = exp(-Lam_o_1)*lam_or1
  dF_og_1 = exp(-Lam_o_1)*lam_og1
  dF_od_1 = exp(-Lam_o_1)*lam_od1
  F_or_1 = t(apply(dF_or_1, 1, cumsum))
  F_og_1 = t(apply(dF_og_1, 1, cumsum))
  F_od_1 = t(apply(dF_od_1, 1, cumsum))
  dF_og__1 = t(apply(dF_og_1*exp(Lam_og_1), 1, cumsum))*exp(-Lam_og_1)
  dF_or__1 = t(apply(dF_or_1*exp(Lam_or_1), 1, cumsum))*exp(-Lam_or_1)
  dF_ogr_1 = dF_og__1*lam_ogr1
  F_ogr_1 = t(apply(dF_ogr_1, 1, cumsum))
  dF_org_1 = dF_or__1*lam_org1
  F_org_1 = t(apply(dF_org_1, 1, cumsum))
  dF_ord_1 = dF_or__1*lam_ord1
  F_ord_1 = t(apply(dF_ord_1, 1, cumsum))
  dF_ogd_1 = dF_og__1*lam_ogd1
  F_ogd_1 = t(apply(dF_ogd_1, 1, cumsum))
  dF_org__1 = t(apply(dF_org_1*exp(Lam_org_1), 1, cumsum))*exp(-Lam_org_1)
  dF_ogr__1 = t(apply(dF_ogr_1*exp(Lam_ogr_1), 1, cumsum))*exp(-Lam_ogr_1)
  dF_orgd_1 = dF_org__1*lam_orgd1
  F_orgd_1 = t(apply(dF_orgd_1, 1, cumsum))
  dF_ogrd_1 = dF_ogr__1*lam_ogrd1
  F_ogrd_1 = t(apply(dF_ogrd_1, 1, cumsum))
  F_1 = colMeans(F_od_1+F_ord_1+F_ogd_1+F_orgd_1+F_ogrd_1)
  
  #prev_g = (F_og_0-F_ogr_0-F_ogd_0)/(1-F_od_0-F_or_0-F_ogd_0-F_ogr_0)
  #lam_or0 = lam_ogr0 = prev_g*lam_ogr0+(1-prev_g)*lam_or0
  Lam_o_0 = t(apply(lam_od0+lam_og0+lam_or0, 1, cumsum))
  Lam_og_0 = t(apply(lam_ogd0+lam_ogr0, 1, cumsum))
  Lam_or_0 = t(apply(lam_ord0+lam_org0, 1, cumsum))
  Lam_ogr_0 = t(apply(lam_ogrd0, 1, cumsum))
  Lam_org_0 = t(apply(lam_orgd0, 1, cumsum))
  dF_or_0 = exp(-Lam_o_0)*lam_or0
  dF_og_0 = exp(-Lam_o_0)*lam_og0
  dF_od_0 = exp(-Lam_o_0)*lam_od0
  F_or_0 = t(apply(dF_or_0, 1, cumsum))
  F_og_0 = t(apply(dF_og_0, 1, cumsum))
  F_od_0 = t(apply(dF_od_0, 1, cumsum))
  dF_og__0 = t(apply(dF_og_0*exp(Lam_og_0), 1, cumsum))*exp(-Lam_og_0)
  dF_or__0 = t(apply(dF_or_0*exp(Lam_or_0), 1, cumsum))*exp(-Lam_or_0)
  dF_ogr_0 = dF_og__0*lam_ogr0
  F_ogr_0 = t(apply(dF_ogr_0, 1, cumsum))
  dF_org_0 = dF_or__0*lam_org0
  F_org_0 = t(apply(dF_org_0, 1, cumsum))
  dF_ord_0 = dF_or__0*lam_ord0
  F_ord_0 = t(apply(dF_ord_0, 1, cumsum))
  dF_ogd_0 = dF_og__0*lam_ogd0
  F_ogd_0 = t(apply(dF_ogd_0, 1, cumsum))
  dF_org__0 = t(apply(dF_org_0*exp(Lam_org_0), 1, cumsum))*exp(-Lam_org_0)
  dF_ogr__0 = t(apply(dF_ogr_0*exp(Lam_ogr_0), 1, cumsum))*exp(-Lam_ogr_0)
  dF_orgd_0 = dF_org__0*lam_orgd0
  F_orgd_0 = t(apply(dF_orgd_0, 1, cumsum))
  dF_ogrd_0 = dF_ogr__0*lam_ogrd0
  F_ogrd_0 = t(apply(dF_ogrd_0, 1, cumsum))
  F_0 = colMeans(F_od_0+F_ord_0+F_ogd_0+F_orgd_0+F_ogrd_0)
  
  # counterfactual incidence, separable
  # a_or=a_gr=a_og=a_rg=0,a_od=a_rd=a_gd=1
  lam_og_a = lam_og0
  lam_od_a = lam_od1
  lam_or_a = lam_or0
  lam_ogr_a = lam_ogr0
  lam_ogd_a = lam_ogd1
  lam_org_a = lam_org0
  lam_ord_a = lam_ord1
  lam_ogrd_a = lam_ogrd1
  lam_orgd_a = lam_orgd1
  Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
  Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
  Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
  Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
  Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_og_a = exp(-Lam_o_a)*lam_og_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_or_a = t(apply(dF_or_a, 1, cumsum))
  F_og_a = t(apply(dF_og_a, 1, cumsum))
  F_od_a = t(apply(dF_od_a, 1, cumsum))
  dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
  dF_ogr_a = dF_og__a*lam_ogr_a
  F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
  dF_org_a = dF_or__a*lam_org_a
  F_org_a = t(apply(dF_org_a, 1, cumsum))
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))
  dF_ogd_a = dF_og__a*lam_ogd_a
  F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
  dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
  dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
  dF_orgd_a = dF_org__a*lam_orgd_a
  F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
  dF_ogrd_a = dF_ogr__a*lam_ogrd_a
  F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
  Fd_spe = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)
  
  # counterfactual incidence, conditionally interventional
  # a_or=a_gr=0,a_od=a_og=a_rd=a_gd=a_rg=1
  lam_og_a = lam_og1
  lam_od_a = lam_od1
  lam_or_a = lam_or0
  lam_ogr_a = lam_ogr0
  lam_ogd_a = lam_ogd1
  lam_org_a = lam_org1
  lam_ord_a = lam_ord1
  lam_ogrd_a = lam_ogrd1
  lam_orgd_a = lam_orgd1
  Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
  Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
  Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
  Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
  Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_og_a = exp(-Lam_o_a)*lam_og_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_or_a = t(apply(dF_or_a, 1, cumsum))
  F_og_a = t(apply(dF_og_a, 1, cumsum))
  F_od_a = t(apply(dF_od_a, 1, cumsum))
  dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
  dF_ogr_a = dF_og__a*lam_ogr_a
  F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
  dF_org_a = dF_or__a*lam_org_a
  F_org_a = t(apply(dF_org_a, 1, cumsum))
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))
  dF_ogd_a = dF_og__a*lam_ogd_a
  F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
  dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
  dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
  dF_orgd_a = dF_org__a*lam_orgd_a
  F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
  dF_ogrd_a = dF_ogr__a*lam_ogrd_a
  F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
  Fd_cin = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)
  
  # potential incidence, overall
  lam_og_a = lam_og1
  lam_od_a = lam_od1
  lam_ogd_a = lam_ogd1
  lam_org_a = lam_org1
  lam_ord_a = lam_ord1
  lam_ogrd_a = lam_ogrd1
  lam_orgd_a = lam_orgd1
  prev_g = (F_og_1-F_ogr_1-F_ogd_1)/(1-F_od_1-F_or_1-F_ogd_1-F_ogr_1)
  lam_or_a = lam_ogr_a = prev_g*lam_ogr1+(1-prev_g)*lam_or1
  Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
  Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
  Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
  Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
  Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_og_a = exp(-Lam_o_a)*lam_og_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_or_a = t(apply(dF_or_a, 1, cumsum))
  F_og_a = t(apply(dF_og_a, 1, cumsum))
  F_od_a = t(apply(dF_od_a, 1, cumsum))
  dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
  dF_ogr_a = dF_og__a*lam_ogr_a
  F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
  dF_org_a = dF_or__a*lam_org_a
  F_org_a = t(apply(dF_org_a, 1, cumsum))
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))
  dF_ogd_a = dF_og__a*lam_ogd_a
  F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
  dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
  dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
  dF_orgd_a = dF_org__a*lam_orgd_a
  F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
  dF_ogrd_a = dF_ogr__a*lam_ogrd_a
  F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
  F_1o = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)
  
  # lam_or = lam_ogr
  lam_og_a = lam_og0
  lam_od_a = lam_od0
  lam_ogd_a = lam_ogd0
  lam_org_a = lam_org0
  lam_ord_a = lam_ord0
  lam_ogrd_a = lam_ogrd0
  lam_orgd_a = lam_orgd0
  prev_g = (F_og_0-F_ogr_0-F_ogd_0)/(1-F_od_0-F_or_0-F_ogd_0-F_ogr_0)
  lam_or_a = lam_ogr_a = prev_g*lam_ogr0+(1-prev_g)*lam_or0
  Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
  Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
  Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
  Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
  Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_og_a = exp(-Lam_o_a)*lam_og_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_or_a = t(apply(dF_or_a, 1, cumsum))
  F_og_a = t(apply(dF_og_a, 1, cumsum))
  F_od_a = t(apply(dF_od_a, 1, cumsum))
  dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
  dF_ogr_a = dF_og__a*lam_ogr_a
  F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
  dF_org_a = dF_or__a*lam_org_a
  F_org_a = t(apply(dF_org_a, 1, cumsum))
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))
  dF_ogd_a = dF_og__a*lam_ogd_a
  F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
  dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
  dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
  dF_orgd_a = dF_org__a*lam_orgd_a
  F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
  dF_ogrd_a = dF_ogr__a*lam_ogrd_a
  F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
  F_0o = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)
  
  # counterfactual incidence, interventional
  lam_og_a = lam_og1
  lam_od_a = lam_od1
  lam_ogd_a = lam_ogd1
  lam_org_a = lam_org1
  lam_ord_a = lam_ord1
  lam_ogrd_a = lam_ogrd1
  lam_orgd_a = lam_orgd1
  prev_g = (F_og_0-F_ogr_0-F_ogd_0)/(1-F_od_0-F_or_0-F_ogd_0-F_ogr_0)
  lam_or_a = lam_ogr_a = prev_g*lam_ogr0+(1-prev_g)*lam_or0
  Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
  Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
  Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
  Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
  Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
  dF_or_a = exp(-Lam_o_a)*lam_or_a
  dF_og_a = exp(-Lam_o_a)*lam_og_a
  dF_od_a = exp(-Lam_o_a)*lam_od_a
  F_or_a = t(apply(dF_or_a, 1, cumsum))
  F_og_a = t(apply(dF_og_a, 1, cumsum))
  F_od_a = t(apply(dF_od_a, 1, cumsum))
  dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
  dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
  dF_ogr_a = dF_og__a*lam_ogr_a
  F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
  dF_org_a = dF_or__a*lam_org_a
  F_org_a = t(apply(dF_org_a, 1, cumsum))
  dF_ord_a = dF_or__a*lam_ord_a
  F_ord_a = t(apply(dF_ord_a, 1, cumsum))
  dF_ogd_a = dF_og__a*lam_ogd_a
  F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
  dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
  dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
  dF_orgd_a = dF_org__a*lam_orgd_a
  F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
  dF_ogrd_a = dF_ogr__a*lam_ogrd_a
  F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
  Fd_int = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)
  
  F_1 = matcht(tt, F_1, xseq)
  F_0 = matcht(tt, F_0, xseq)
  F_1o = matcht(tt, F_1o, xseq)
  F_0o = matcht(tt, F_0o, xseq)
  Fd_int = matcht(tt, Fd_int, xseq)
  Fd_cin = matcht(tt, Fd_cin, xseq)
  F_ide = Fd_int - F_0o
  F_iie = F_1o - Fd_int
  F_cde = Fd_cin - F_0
  F_cie = F_1 - Fd_cin
  res_np = rbind(F_1, F_0, F_1o, F_0o, Fd_int, Fd_cin, F_ide, F_iie, F_cde, F_cie)
  rownames(res_np) = c('Potential 1', 'Potential 0', 'Overall 1', 'Overall 0',
                       'Interv', 'Cond Interv',
                       'IDE','IIE','CDE','CIE')
  colnames(res_np) = xseq
  time.end = Sys.time()
  time_np = difftime(time.end, time.start, unit='secs')
  
  return(list(res_np=res_np, time_np=time_np))
}

library(statmod)
library(MASS)
library(coxme)
ngrid = 8
quad = gauss.quad(ngrid, "hermite")
b_list = quad$nodes
w_list = quad$weights

RandInt_fr <- function(Tg,Dg,Tr,Dr,Td,Dd,A,X) {
  
  xseq = 1:15
  time.start = Sys.time()
  tt = sort(unique(c(Td[Dd==1],Tr[Dr==1],Tg[Dg==1])))
  l = length(tt)
  X = as.matrix(X)
  p = ncol(X)
  n = nrow(X)
  
  df = data.frame(id=1:n, Td=Td, Dd=Dd, Tr=Tr, Dr=Dr, Tg=Tg, Dg=Dg, X=X, A=A)
  xvars = grep("^X", names(df), value=TRUE)
  mg_d = tmerge(data1=df, data2=df, id=id, death=event(Td,Dd))
  mg_d = tmerge(mg_d, df, id=id, Gevent=tdc(Tg), Revent=tdc(Tr))
  df = data.frame(id=1:n, Tr=Tr, Dr=Dr, Tg=Tg, Dg=Dg, X=X, A=A)
  df$Tr[df$Tr>df$Tg] = df$Tg[df$Tr>df$Tg]; df$Dr[df$Tr>df$Tg] = 0
  mg_g = tmerge(data1=df, data2=df, id=id, gvhd=event(df$Tg,df$Dg))
  mg_g = tmerge(mg_g, df, id=id, Revent=tdc(df$Tr))
  df = data.frame(id=1:n, Tr=Tr, Dr=Dr, Tg=Tg, Dg=Dg, X=X, A=A)
  df$Tg[df$Tg>df$Tr] = df$Tr[df$Tg>df$Tr]; df$Dg[df$Tg>df$Tr] = 0
  mg_r = tmerge(data1=df, data2=df, id=id, relapse=event(df$Tr,df$Dr))
  mg_r = tmerge(mg_r, df, id=id, Gevent=tdc(df$Tg))
  mg_d = data.frame(id=mg_d$id,X_g=mg_d[,xvars]*0,X_r=mg_d[,xvars]*0,X_d=mg_d[,xvars],A=mg_d$A,
                    tstart=mg_d$tstart,tstop=mg_d$tstop,
                    event=ifelse(mg_d$death==1,'death','censoring'),
                    Gevent_d=mg_d$Gevent,Revent_d=mg_d$Revent,
                    Gevent_r=0,Revent_g=0,terevent='death')
  mg_g = data.frame(id=mg_g$id,X_g=mg_g[,xvars],X_r=mg_g[,xvars]*0,X_d=mg_g[,xvars]*0,A=mg_g$A,
                    tstart=mg_g$tstart,tstop=mg_g$tstop,
                    event=ifelse(mg_g$gvhd==1,'gvhd','censoring'),
                    Gevent_d=0,Revent_d=0,
                    Gevent_r=0,Revent_g=mg_g$Revent,terevent='gvhd')
  mg_r = data.frame(id=mg_r$id,X_g=mg_r[,xvars]*0,X_r=mg_r[,xvars],X_d=mg_r[,xvars]*0,A=mg_r$A,
                    tstart=mg_r$tstart,tstop=mg_r$tstop,
                    event=ifelse(mg_r$relapse==1,'relapse','censoring'),
                    Gevent_d=0,Revent_d=0,
                    Gevent_r=mg_r$Gevent,Revent_g=0,terevent='relapse')
  mg = rbind(mg_d,mg_g,mg_r)
  mg = mg[order(mg$id),]
  xvars = grep("^X", names(mg), value=TRUE)
  cox_formula = reformulate(c("Gevent_r","Gevent_d","Revent_g","Revent_d",xvars,
                              "strata(terevent)","(1|id)"), 
                            response="Surv(tstart,tstop,event!='censoring')")
  fit1 = coxme(cox_formula, data=mg, subset=(A==1))
  fit0 = coxme(cox_formula, data=mg, subset=(A==0))
  s1 = sqrt(VarCorr(fit1)$id); s0 = sqrt(VarCorr(fit0)$id)
  cox_formula = reformulate(c("Gevent_r","Gevent_d","Revent_g","Revent_d",xvars,
                              "strata(terevent)","offset(frailty_offset)"),
                            response="Surv(tstart,tstop,event!='censoring')")
  b1 = ranef(fit1)$id
  frailty_offset = b1[as.character(mg[mg$A==1,]$id)]
  fit1_coxph = coxph(cox_formula, data=mg[mg$A==1,], init=fixef(fit1),
                     control=coxph.control(iter.max=0))
  lam1 = basehaz(fit1_coxph, centered=FALSE)
  b0 = ranef(fit0)$id
  frailty_offset = b0[as.character(mg[mg$A==0,]$id)]
  fit0_coxph = coxph(cox_formula, data=mg[mg$A==0,], init=fixef(fit0),
                     control=coxph.control(iter.max = 0))
  lam0 = basehaz(fit0_coxph, centered=FALSE)
  
  coefg = rbind(cbind(fit1$coefficients,sqrt(diag(vcov(fit1))),
                      fit0$coefficients,sqrt(diag(vcov(fit0))))[c(4+1:p,3),],
                Gevent=NA)
  coefr = rbind(cbind(fit1$coefficients,sqrt(diag(vcov(fit1))),
                      fit0$coefficients,sqrt(diag(vcov(fit0))))[c(4+p+1:p,1),],
                Revent=NA)
  coefd = cbind(fit1$coefficients,sqrt(diag(vcov(fit1))),
                fit0$coefficients,sqrt(diag(vcov(fit0))))[c(4+2*p+1:p,4,2),]
  
  Xbd1 = as.numeric(X%*%coefd[1:p,1])
  delta_gd1 = coefd[p+2,1]; delta_rd1 = coefd[p+1,1]
  Xbg1 = as.numeric(X%*%coefg[1:p,1])
  delta_rg1 = coefg[p+1,1]
  Xbr1 = as.numeric(X%*%coefr[1:p,1])
  delta_gr1 = coefr[p+1,1]
  Xbd0 = as.numeric(X%*%coefd[1:p,3])
  delta_gd0 = coefd[p+2,3]; delta_rd0 = coefd[p+1,3]
  Xbg0 = as.numeric(X%*%coefg[1:p,3])
  delta_rg0 = coefg[p+1,3]
  Xbr0 = as.numeric(X%*%coefr[1:p,3])
  delta_gr0 = coefr[p+1,3]
  lamd1 = lam1[lam1$strata=='death',]
  lam_d1 = diff(c(0,matcht(lamd1$time, lamd1$hazard, tt)))
  lamg1 = lam1[lam1$strata=='gvhd',]
  lam_g1 = diff(c(0,matcht(lamg1$time, lamg1$hazard, tt)))
  lamr1 = lam1[lam1$strata=='relapse',]
  lam_r1 = diff(c(0,matcht(lamr1$time, lamr1$hazard, tt)))
  lamd0 = lam0[lam0$strata=='death',]
  lam_d0 = diff(c(0,matcht(lamd0$time, lamd0$hazard, tt)))
  lamg0 = lam0[lam0$strata=='gvhd',]
  lam_g0 = diff(c(0,matcht(lamg0$time, lamg0$hazard, tt)))
  lamr0 = lam0[lam0$strata=='relapse',]
  lam_r0 = diff(c(0,matcht(lamr0$time, lamr0$hazard, tt)))
  
  lamd1 = lam1[lam1$strata=='death',]
  lam_d1 = diff(c(0,matcht(lamd1$time, lamd1$hazard, tt)))
  lamg1 = lam1[lam1$strata=='gvhd',]
  lam_g1 = diff(c(0,matcht(lamg1$time, lamg1$hazard, tt)))
  lamr1 = lam1[lam1$strata=='relapse',]
  lam_r1 = diff(c(0,matcht(lamr1$time, lamr1$hazard, tt)))
  lamd0 = lam0[lam0$strata=='death',]
  lam_d0 = diff(c(0,matcht(lamd0$time, lamd0$hazard, tt)))
  lamg0 = lam0[lam0$strata=='gvhd',]
  lam_g0 = diff(c(0,matcht(lamg0$time, lamg0$hazard, tt)))
  lamr0 = lam0[lam0$strata=='relapse',]
  lam_r0 = diff(c(0,matcht(lamr0$time, lamr0$hazard, tt)))
  
  iter = 0; lam = c(lam_d1,lam_g1,lam_r1)
  while(iter<100){
    iter = iter+1; lam0 = lam
    Ee0 = Eeb = 0
    for (b in 1:ngrid){
      e = sqrt(2)*s1*b_list[b]
      Lamd = sapply(1:n, function(i) sum(lam_d1*(tt<=Td[i])*exp(e+Xbd1[i]+delta_gd1*(Tg[i]<tt)+delta_rd1*(Tr[i]<tt))))
      dlamd = sapply(1:n, function(i) sum(lam_d1*(tt==Td[i])*exp(e+Xbd1[i]+delta_gd1*(Tg[i]<tt)+delta_rd1*(Tr[i]<tt))))
      Lamg = sapply(1:n, function(i) sum(lam_g1*(tt<=Tg[i])*exp(e+Xbg1[i]+delta_rg1*(Tr[i]<tt))))
      dlamg = sapply(1:n, function(i) sum(lam_g1*(tt==Tg[i])*exp(e+Xbg1[i]+delta_rg1*(Tr[i]<tt))))
      Lamr = sapply(1:n, function(i) sum(lam_r1*(tt<=Tr[i])*exp(e+Xbr1[i]+delta_gr1*(Tg[i]<tt))))
      dlamr = sapply(1:n, function(i) sum(lam_r1*(tt==Tr[i])*exp(e+Xbr1[i]+delta_gr1*(Tg[i]<tt))))
      dLamd = Dd*log(dlamd); dLamd[dlamd==0] = 0
      dLamg = Dg*log(dlamg); dLamg[dlamg==0] = 0
      dLamr = Dr*log(dlamr); dLamr[dlamr==0] = 0
      lik = exp(-Lamd+dLamd-Lamg+dLamg-Lamr+dLamr)*w_list[b]
      Ee0 = Ee0 + lik
      Eeb = Eeb + lik*exp(e)
    }
    Eeb = Eeb/Ee0
    lam_d1 = sapply(tt, function(t) sum(A*Dd*(Td==t))/sum(A*(Td>=t)*
                                                            Eeb*exp(Xbd1+(Tg<t)*Dg*delta_gd1+(Tr<t)*Dr*delta_rd1)))
    lam_d1[is.nan(lam_d1)] = 0
    lam_g1 = sapply(tt, function(t) sum(A*Dg*(Tg==t))/sum(A*(Tg>=t)*
                                                            Eeb*exp(Xbg1+(Tr<t)*Dr*delta_rg1)))
    lam_g1[is.nan(lam_g1)] = 0
    lam_r1 = sapply(tt, function(t) sum(A*Dr*(Tr==t))/sum(A*(Tr>=t)*
                                                            Eeb*exp(Xbr1+(Tg<t)*Dg*delta_gr1)))
    lam_r1[is.nan(lam_r1)] = 0
    lam = c(lam_d1,lam_g1,lam_r1)
    tol = max(abs(lam-lam0))
    if (tol<1e-8) break
  }
  iter = 0; lam = c(lam_d0,lam_g0,lam_r0)
  while(iter<100){
    iter = iter+1; lam0 = lam
    Ee0 = Eeb = 0
    for (b in 1:ngrid){
      e = sqrt(2)*s0*b_list[b]
      Lamd = sapply(1:n, function(i) sum(lam_d0*(tt<=Td[i])*exp(e+Xbd0[i]+delta_gd0*(Tg[i]<tt)+delta_rd0*(Tr[i]<tt))))
      dlamd = sapply(1:n, function(i) sum(lam_d0*(tt==Td[i])*exp(e+Xbd0[i]+delta_gd0*(Tg[i]<tt)+delta_rd0*(Tr[i]<tt))))
      Lamg = sapply(1:n, function(i) sum(lam_g0*(tt<=Tg[i])*exp(e+Xbg0[i]+delta_rg0*(Tr[i]<tt))))
      dlamg = sapply(1:n, function(i) sum(lam_g0*(tt==Tg[i])*exp(e+Xbg0[i]+delta_rg0*(Tr[i]<tt))))
      Lamr = sapply(1:n, function(i) sum(lam_r0*(tt<=Tr[i])*exp(e+Xbr0[i]+delta_gr0*(Tg[i]<tt))))
      dlamr = sapply(1:n, function(i) sum(lam_r0*(tt==Tr[i])*exp(e+Xbr0[i]+delta_gr0*(Tg[i]<tt))))
      dLamd = Dd*log(dlamd); dLamd[dlamd==0] = 0
      dLamg = Dg*log(dlamg); dLamg[dlamg==0] = 0
      dLamr = Dr*log(dlamr); dLamr[dlamr==0] = 0
      lik = exp(-Lamd+dLamd-Lamg+dLamg-Lamr+dLamr)*w_list[b]
      Ee0 = Ee0 + lik
      Eeb = Eeb + lik*exp(e)
    }
    Eeb = Eeb/Ee0
    lam_d0 = sapply(tt, function(t) sum((1-A)*Dd*(Td==t))/sum((1-A)*(Td>=t)*
                                                                Eeb*exp(Xbd0+(Tg<t)*Dg*delta_gd0+(Tr<t)*Dr*delta_rd0)))
    lam_d0[is.nan(lam_d0)] = 0
    lam_g0 = sapply(tt, function(t) sum((1-A)*Dg*(Tg==t))/sum((1-A)*(Tg>=t)*
                                                                Eeb*exp(Xbg0+(Tr<t)*Dr*delta_rg0)))
    lam_g0[is.nan(lam_g0)] = 0
    lam_r0 = sapply(tt, function(t) sum((1-A)*Dr*(Tr==t))/sum((1-A)*(Tr>=t)*
                                                                Eeb*exp(Xbr0+(Tg<t)*Dg*delta_gr0)))
    lam_r0[is.nan(lam_r0)] = 0
    lam = c(lam_d0,lam_g0,lam_r0)
    tol = max(abs(lam-lam0))
    if (tol<1e-8) break
  }
  
  F_1 = F_0 = F_1o = F_0o = Fd_spe = Fd_int = Fd_cin = 0
  for (b in 1:ngrid){
    e1 = sqrt(2)*s1*b_list[b]
    e0 = sqrt(2)*s0*b_list[b]
    lam_od1 = sapply(1:l, function(t) lam_d1[t]*exp(Xbd1+e1))
    lam_ogd1 = lam_od1 * exp(delta_gd1)
    lam_ord1 = lam_od1 * exp(delta_rd1)
    lam_ogrd1 = lam_orgd1 = lam_od1 * exp(delta_gd1+delta_rd1)
    lam_og1 = sapply(1:l, function(t) lam_g1[t]*exp(Xbg1+e1))
    lam_org1 = lam_og1 *exp(delta_rg1)
    lam_or1 = sapply(1:l, function(t) lam_r1[t]*exp(Xbr1+e1))
    lam_ogr1 = lam_or1 * exp(delta_gr1)
    lam_od0 = sapply(1:l, function(t) lam_d0[t]*exp(Xbd0+e0))
    lam_ogd0 = lam_od0 * exp(delta_gd0)
    lam_ord0 = lam_od0 * exp(delta_rd0)
    lam_ogrd0 = lam_orgd0 = lam_od0 * exp(delta_gd0+delta_rd0)
    lam_og0 = sapply(1:l, function(t) lam_g0[t]*exp(Xbg0+e0))
    lam_org0 = lam_og0 *exp(delta_rg0)
    lam_or0 = sapply(1:l, function(t) lam_r0[t]*exp(Xbr0+e0))
    lam_ogr0 = lam_or0 * exp(delta_gr0)
    
    # potential incidence
    #prev_g = (F_og_1-F_ogr_1-F_ogd_1)/(1-F_od_1-F_or_1-F_ogd_1-F_ogr_1)
    #lam_or1 = lam_ogr1 = prev_g*lam_ogr1+(1-prev_g)*lam_or1
    Lam_o_1 = t(apply(lam_od1+lam_og1+lam_or1, 1, cumsum))
    Lam_og_1 = t(apply(lam_ogd1+lam_ogr1, 1, cumsum))
    Lam_or_1 = t(apply(lam_ord1+lam_org1, 1, cumsum))
    Lam_ogr_1 = t(apply(lam_ogrd1, 1, cumsum))
    Lam_org_1 = t(apply(lam_orgd1, 1, cumsum))
    dF_or_1 = exp(-Lam_o_1)*lam_or1
    dF_og_1 = exp(-Lam_o_1)*lam_og1
    dF_od_1 = exp(-Lam_o_1)*lam_od1
    F_or_1 = t(apply(dF_or_1, 1, cumsum))
    F_og_1 = t(apply(dF_og_1, 1, cumsum))
    F_od_1 = t(apply(dF_od_1, 1, cumsum))
    dF_og__1 = t(apply(dF_og_1*exp(Lam_og_1), 1, cumsum))*exp(-Lam_og_1)
    dF_or__1 = t(apply(dF_or_1*exp(Lam_or_1), 1, cumsum))*exp(-Lam_or_1)
    dF_ogr_1 = dF_og__1*lam_ogr1
    F_ogr_1 = t(apply(dF_ogr_1, 1, cumsum))
    dF_org_1 = dF_or__1*lam_org1
    F_org_1 = t(apply(dF_org_1, 1, cumsum))
    dF_ord_1 = dF_or__1*lam_ord1
    F_ord_1 = t(apply(dF_ord_1, 1, cumsum))
    dF_ogd_1 = dF_og__1*lam_ogd1
    F_ogd_1 = t(apply(dF_ogd_1, 1, cumsum))
    dF_org__1 = t(apply(dF_org_1*exp(Lam_org_1), 1, cumsum))*exp(-Lam_org_1)
    dF_ogr__1 = t(apply(dF_ogr_1*exp(Lam_ogr_1), 1, cumsum))*exp(-Lam_ogr_1)
    dF_orgd_1 = dF_org__1*lam_orgd1
    F_orgd_1 = t(apply(dF_orgd_1, 1, cumsum))
    dF_ogrd_1 = dF_ogr__1*lam_ogrd1
    F_ogrd_1 = t(apply(dF_ogrd_1, 1, cumsum))
    F_1 = F_1 + colMeans(F_od_1+F_ord_1+F_ogd_1+F_orgd_1+F_ogrd_1)*w_list[b]/sqrt(pi)
    
    Lam_o_0 = t(apply(lam_od0+lam_og0+lam_or0, 1, cumsum))
    Lam_og_0 = t(apply(lam_ogd0+lam_ogr0, 1, cumsum))
    Lam_or_0 = t(apply(lam_ord0+lam_org0, 1, cumsum))
    Lam_ogr_0 = t(apply(lam_ogrd0, 1, cumsum))
    Lam_org_0 = t(apply(lam_orgd0, 1, cumsum))
    dF_or_0 = exp(-Lam_o_0)*lam_or0
    dF_og_0 = exp(-Lam_o_0)*lam_og0
    dF_od_0 = exp(-Lam_o_0)*lam_od0
    F_or_0 = t(apply(dF_or_0, 1, cumsum))
    F_og_0 = t(apply(dF_og_0, 1, cumsum))
    F_od_0 = t(apply(dF_od_0, 1, cumsum))
    dF_og__0 = t(apply(dF_og_0*exp(Lam_og_0), 1, cumsum))*exp(-Lam_og_0)
    dF_or__0 = t(apply(dF_or_0*exp(Lam_or_0), 1, cumsum))*exp(-Lam_or_0)
    dF_ogr_0 = dF_og__0*lam_ogr0
    F_ogr_0 = t(apply(dF_ogr_0, 1, cumsum))
    dF_org_0 = dF_or__0*lam_org0
    F_org_0 = t(apply(dF_org_0, 1, cumsum))
    dF_ord_0 = dF_or__0*lam_ord0
    F_ord_0 = t(apply(dF_ord_0, 1, cumsum))
    dF_ogd_0 = dF_og__0*lam_ogd0
    F_ogd_0 = t(apply(dF_ogd_0, 1, cumsum))
    dF_org__0 = t(apply(dF_org_0*exp(Lam_org_0), 1, cumsum))*exp(-Lam_org_0)
    dF_ogr__0 = t(apply(dF_ogr_0*exp(Lam_ogr_0), 1, cumsum))*exp(-Lam_ogr_0)
    dF_orgd_0 = dF_org__0*lam_orgd0
    F_orgd_0 = t(apply(dF_orgd_0, 1, cumsum))
    dF_ogrd_0 = dF_ogr__0*lam_ogrd0
    F_ogrd_0 = t(apply(dF_ogrd_0, 1, cumsum))
    F_0 = F_0 + colMeans(F_od_0+F_ord_0+F_ogd_0+F_orgd_0+F_ogrd_0)*w_list[b]/sqrt(pi)
    
    # counterfactual incidence, separable
    # a_or=a_gr=a_og=a_rg=0,a_od=a_rd=a_gd=1
    lam_og_a = lam_og0
    lam_od_a = lam_od1
    lam_or_a = lam_or0
    lam_ogr_a = lam_ogr0
    lam_ogd_a = lam_ogd1
    lam_org_a = lam_org0
    lam_ord_a = lam_ord1
    lam_ogrd_a = lam_ogrd1
    lam_orgd_a = lam_orgd1
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    Fd_spe = colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)
    
    # counterfactual incidence, conditional intervention
    # a_or=a_gr=0,a_od=a_og=a_rd=a_gd=a_rg=1
    lam_og_a = lam_og1
    lam_od_a = lam_od1
    lam_or_a = lam_or0
    lam_ogr_a = lam_ogr0
    lam_ogd_a = lam_ogd1
    lam_org_a = lam_org1
    lam_ord_a = lam_ord1
    lam_ogrd_a = lam_ogrd1
    lam_orgd_a = lam_orgd1
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    Fd_cin = Fd_cin + colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)*w_list[b]/sqrt(pi)
    
    # overall incidence
    lam_og_a = lam_og1
    lam_od_a = lam_od1
    lam_ogd_a = lam_ogd1
    lam_org_a = lam_org1
    lam_ord_a = lam_ord1
    lam_ogrd_a = lam_ogrd1
    lam_orgd_a = lam_orgd1
    prev_g = (F_og_1-F_ogr_1-F_ogd_1)/(1-F_od_1-F_or_1-F_ogd_1-F_ogr_1)
    lam_or_a = lam_ogr_a = prev_g*lam_ogr1+(1-prev_g)*lam_or1
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    F_1o = F_1o + colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)*w_list[b]/sqrt(pi)
    
    lam_og_a = lam_og0
    lam_od_a = lam_od0
    lam_ogd_a = lam_ogd0
    lam_org_a = lam_org0
    lam_ord_a = lam_ord0
    lam_ogrd_a = lam_ogrd0
    lam_orgd_a = lam_orgd0
    prev_g = (F_og_0-F_ogr_0-F_ogd_0)/(1-F_od_0-F_or_0-F_ogd_0-F_ogr_0)
    lam_or_a = lam_ogr_a = prev_g*lam_ogr0+(1-prev_g)*lam_or0
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    F_0o = F_0o + colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)*w_list[b]/sqrt(pi)
    
    # counterfactual incidence, intervention
    lam_og_a = lam_og1
    lam_od_a = lam_od1
    lam_ogd_a = lam_ogd1
    lam_org_a = lam_org1
    lam_ord_a = lam_ord1
    lam_ogrd_a = lam_ogrd1
    lam_orgd_a = lam_orgd1
    prev_g = (F_og_0-F_ogr_0-F_ogd_0)/(1-F_od_0-F_or_0-F_ogd_0-F_ogr_0)
    lam_or_a = lam_ogr_a = prev_g*lam_ogr0+(1-prev_g)*lam_or0
    Lam_o_a = t(apply(lam_od_a+lam_og_a+lam_or_a, 1, cumsum))
    Lam_og_a = t(apply(lam_ogd_a+lam_ogr_a, 1, cumsum))
    Lam_or_a = t(apply(lam_ord_a+lam_org_a, 1, cumsum))
    Lam_ogr_a = t(apply(lam_ogrd_a, 1, cumsum))
    Lam_org_a = t(apply(lam_orgd_a, 1, cumsum))
    dF_or_a = exp(-Lam_o_a)*lam_or_a
    dF_og_a = exp(-Lam_o_a)*lam_og_a
    dF_od_a = exp(-Lam_o_a)*lam_od_a
    F_or_a = t(apply(dF_or_a, 1, cumsum))
    F_og_a = t(apply(dF_og_a, 1, cumsum))
    F_od_a = t(apply(dF_od_a, 1, cumsum))
    dF_og__a = t(apply(dF_og_a*exp(Lam_og_a), 1, cumsum))*exp(-Lam_og_a)
    dF_or__a = t(apply(dF_or_a*exp(Lam_or_a), 1, cumsum))*exp(-Lam_or_a)
    dF_ogr_a = dF_og__a*lam_ogr_a
    F_ogr_a = t(apply(dF_ogr_a, 1, cumsum))
    dF_org_a = dF_or__a*lam_org_a
    F_org_a = t(apply(dF_org_a, 1, cumsum))
    dF_ord_a = dF_or__a*lam_ord_a
    F_ord_a = t(apply(dF_ord_a, 1, cumsum))
    dF_ogd_a = dF_og__a*lam_ogd_a
    F_ogd_a = t(apply(dF_ogd_a, 1, cumsum))
    dF_org__a = t(apply(dF_org_a*exp(Lam_org_a), 1, cumsum))*exp(-Lam_org_a)
    dF_ogr__a = t(apply(dF_ogr_a*exp(Lam_ogr_a), 1, cumsum))*exp(-Lam_ogr_a)
    dF_orgd_a = dF_org__a*lam_orgd_a
    F_orgd_a = t(apply(dF_orgd_a, 1, cumsum))
    dF_ogrd_a = dF_ogr__a*lam_ogrd_a
    F_ogrd_a = t(apply(dF_ogrd_a, 1, cumsum))
    Fd_int = Fd_int + colMeans(F_od_a + F_ogd_a + F_ord_a + F_ogrd_a + F_orgd_a)*w_list[b]/sqrt(pi)
  }
  
  F_1 = matcht(tt, F_1, xseq)
  F_0 = matcht(tt, F_0, xseq)
  F_1o = matcht(tt, F_1o, xseq)
  F_0o = matcht(tt, F_0o, xseq)
  Fd_int = matcht(tt, Fd_int, xseq)
  Fd_cin = matcht(tt, Fd_cin, xseq)
  F_ide = Fd_int - F_0o
  F_iie = F_1o - Fd_int
  F_cde = Fd_cin - F_0
  F_cie = F_1 - Fd_cin
  res_fr = rbind(F_1, F_0, F_1o, F_0o, Fd_int, Fd_cin, F_ide, F_iie, F_cde, F_cie)
  rownames(res_fr) = c('Potential 1', 'Potential 0', 'Overall 1', 'Overall 0',
                       'Interv', 'Cond Interv',
                       'IDE','IIE','CDE','CIE')
  colnames(res_fr) = xseq
  sigma = c(s1, s0)
  time.end = Sys.time()
  time_fr = difftime(time.end, time.start, unit='secs')
  
  return(list(res_fr=res_fr, time_fr=time_fr, sigma=sigma))
}