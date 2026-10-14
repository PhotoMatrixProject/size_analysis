#setwd("/home/petr_turecek/size4_width")

library(cmdstanr)
library(posterior)
library(loo)
library(nnet)
library(splines)
library(randomcoloR)

# ============================================================
# Data preparation: model absolute width in millimetres with normal likelihood
# ============================================================

d <- read.csv("size_data_2026_10_03.csv")
d$startTime <- as.Date(d$startTime, format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")
d$endTime <- as.Date(d$endTime, format="%Y-%m-%dT%H:%M:%OSZ", tz="UTC")
origin_date <- as.Date("1900-01-01")
d$startTime_numeric <- as.integer(d$startTime-origin_date)+1
d$endTime_numeric <- as.integer(d$endTime-origin_date)+1
d$midTime <- d$startTime_numeric+(d$endTime_numeric-d$startTime_numeric)/2

d_fil <- subset(d,!is.na(startTime) & startTime>=as.Date("1900-01-01") & startTime<=as.Date("1950-12-31"))
d_fil <- subset(d_fil,!is.na(width_mm))
mtime <- mean(d_fil$midTime)
d_fil$midTime <- (d_fil$midTime-mtime)/365
d_fil$CountryLabel <- "Czechia"

K0 <- 5
knots <- seq(min(d_fil$midTime),max(d_fil$midTime),l=K0)
B <- bs(d_fil$midTime, knots=knots[-c(1,K0)], degree=3, intercept=TRUE)

dl <- data.frame(
  IDr=as.integer(as.factor(d_fil$item)),
  width=d_fil$width,
  height=d_fil$height,
  width_mm=d_fil$width_mm,
  page_width_mm=if("page_width_mm" %in% names(d_fil)) d_fil$page_width_mm else NA_real_,
  iss=as.integer(as.factor(d_fil$issue)),
  vol=as.integer(as.factor(d_fil$volume)),
  mag=as.integer(as.factor(d_fil$magazine)),
  country=as.integer(as.factor(d_fil$CountryLabel)),
  time=d_fil$midTime
)

stan_data <- list(N=nrow(dl), y=dl$width_mm, time=dl$time, K=ncol(B), B=unclass(B), N_iss=max(dl$iss), N_vol=max(dl$vol), N_mag=max(dl$mag), iss=dl$iss, vol=dl$vol, mag=dl$mag)

# ============================================================
# Stan model 1: no time, normal on width_mm
# ============================================================

stan_no_time <- '
data {
  int<lower=1> N;
  vector[N] y;

  int<lower=1> N_iss;
  int<lower=1> N_vol;
  int<lower=1> N_mag;

  array[N] int<lower=1,upper=N_iss> iss;
  array[N] int<lower=1,upper=N_vol> vol;
  array[N] int<lower=1,upper=N_mag> mag;
}

parameters {
  real a;

  vector[N_iss] z_iss;
  vector[N_vol] z_vol;
  vector[N_mag] z_mag;

  real<lower=0> sd_iss;
  real<lower=0> sd_vol;
  real<lower=0> sd_mag;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  mu = a + z_iss[iss]*sd_iss + z_vol[vol]*sd_vol + z_mag[mag]*sd_mag;
}

model {
  a ~ normal(100,75);
  z_iss ~ normal(0,1);
  z_vol ~ normal(0,1);
  z_mag ~ normal(0,1);
  sd_iss ~ exponential(0.02);
  sd_vol ~ exponential(0.02);
  sd_mag ~ exponential(0.02);
  sigma ~ exponential(0.02);
  y ~ normal(mu,sigma);
}

generated quantities {
  vector[N] log_lik;
  for(n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n],sigma);
}
'
writeLines(stan_no_time,"normal_width_no_time.stan")

# ============================================================
# Stan model 2: linear time, normal on width_mm
# ============================================================

stan_linear <- '
data {
  int<lower=1> N;
  vector[N] y;
  vector[N] time;

  int<lower=1> N_iss;
  int<lower=1> N_vol;
  int<lower=1> N_mag;

  array[N] int<lower=1,upper=N_iss> iss;
  array[N] int<lower=1,upper=N_vol> vol;
  array[N] int<lower=1,upper=N_mag> mag;
}

parameters {
  real a;
  real b;

  vector[N_iss] z_iss;
  vector[N_vol] z_vol;
  vector[N_mag] z_mag;

  real<lower=0> sd_iss;
  real<lower=0> sd_vol;
  real<lower=0> sd_mag;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  mu = a + b*time + z_iss[iss]*sd_iss + z_vol[vol]*sd_vol + z_mag[mag]*sd_mag;
}

model {
  a ~ normal(100,75);
  b ~ normal(0,5);
  z_iss ~ normal(0,1);
  z_vol ~ normal(0,1);
  z_mag ~ normal(0,1);
  sd_iss ~ exponential(0.02);
  sd_vol ~ exponential(0.02);
  sd_mag ~ exponential(0.02);
  sigma ~ exponential(0.02);
  y ~ normal(mu,sigma);
}

generated quantities {
  vector[N] log_lik;
  for(n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n],sigma);
}
'
writeLines(stan_linear,"normal_width_linear.stan")

# ============================================================
# Stan model 3: spline time, normal on width_mm
# ============================================================

stan_spline <- '
data {
  int<lower=1> N;
  vector[N] y;

  int<lower=1> K;
  matrix[N,K] B;

  int<lower=1> N_iss;
  int<lower=1> N_vol;
  int<lower=1> N_mag;

  array[N] int<lower=1,upper=N_iss> iss;
  array[N] int<lower=1,upper=N_vol> vol;
  array[N] int<lower=1,upper=N_mag> mag;
}

parameters {
  real a;
  vector[K] w;

  vector[N_iss] z_iss;
  vector[N_vol] z_vol;
  vector[N_mag] z_mag;

  real<lower=0> sd_iss;
  real<lower=0> sd_vol;
  real<lower=0> sd_mag;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  mu = a + B*w + z_iss[iss]*sd_iss + z_vol[vol]*sd_vol + z_mag[mag]*sd_mag;
}

model {
  a ~ normal(100,75);
  w ~ normal(0,50);
  z_iss ~ normal(0,1);
  z_vol ~ normal(0,1);
  z_mag ~ normal(0,1);
  sd_iss ~ exponential(0.02);
  sd_vol ~ exponential(0.02);
  sd_mag ~ exponential(0.02);
  sigma ~ exponential(0.02);
  y ~ normal(mu,sigma);
}

generated quantities {
  vector[N] log_lik;
  for(n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n],sigma);
}
'
writeLines(stan_spline,"normal_width_spline.stan")

# ============================================================
# Stan model 4: step time, normal on width_mm
# 4 baselines + 3 estimated cutpoints
# ============================================================

stan_data_step <- stan_data
stan_data_step$t_min <- min(dl$time)
stan_data_step$t_max <- max(dl$time)
stan_data_step$kappa <- 4

stan_step <- '
data {
  int<lower=1> N;
  vector[N] y;
  vector[N] time;
  real t_min;
  real t_max;
  real<lower=0> kappa;

  int<lower=1> N_iss;
  int<lower=1> N_vol;
  int<lower=1> N_mag;

  array[N] int<lower=1,upper=N_iss> iss;
  array[N] int<lower=1,upper=N_vol> vol;
  array[N] int<lower=1,upper=N_mag> mag;
}

parameters {
  vector[4] step_a;
  ordered[3] c_raw;

  vector[N_iss] z_iss;
  vector[N_vol] z_vol;
  vector[N_mag] z_mag;

  real<lower=0> sd_iss;
  real<lower=0> sd_vol;
  real<lower=0> sd_mag;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  vector[3] cut;
  vector[N] s1;
  vector[N] s2;
  vector[N] s3;

  for(j in 1:3) cut[j] = t_min + (t_max-t_min)*inv_logit(c_raw[j]);
  s1 = inv_logit(kappa*(time-cut[1]));
  s2 = inv_logit(kappa*(time-cut[2]));
  s3 = inv_logit(kappa*(time-cut[3]));
  mu = step_a[1]*(1-s1) + step_a[2]*(s1-s2) + step_a[3]*(s2-s3) + step_a[4]*s3 + z_iss[iss]*sd_iss + z_vol[vol]*sd_vol + z_mag[mag]*sd_mag;
}

model {
  step_a ~ normal(100,75);
  c_raw ~ normal(0,1.5);
  z_iss ~ normal(0,1);
  z_vol ~ normal(0,1);
  z_mag ~ normal(0,1);
  sd_iss ~ exponential(0.02);
  sd_vol ~ exponential(0.02);
  sd_mag ~ exponential(0.02);
  sigma ~ exponential(0.02);
  y ~ normal(mu,sigma);
}

generated quantities {
  vector[N] log_lik;
  for(n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n],sigma);
}
'
writeLines(stan_step,"normal_width_step.stan")

# ============================================================
# Compile and fit models
# ============================================================

mod0 <- cmdstan_model("normal_width_no_time.stan")
mod1 <- cmdstan_model("normal_width_linear.stan")
modS <- cmdstan_model("normal_width_spline.stan")
modStep <- cmdstan_model("normal_width_step.stan")

fit0 <- mod0$sample(data=stan_data, chains=4, parallel_chains=4, iter_warmup=1000, iter_sampling=1000, adapt_delta=0.9)
fit1 <- mod1$sample(data=stan_data, chains=4, parallel_chains=4, iter_warmup=1000, iter_sampling=1000, adapt_delta=0.9)
fitS <- modS$sample(data=stan_data, chains=4, parallel_chains=4, iter_warmup=1000, iter_sampling=1000, adapt_delta=0.9)
fitStep <- modStep$sample(data=stan_data_step, chains=4, parallel_chains=4, iter_warmup=1000, iter_sampling=1000, adapt_delta=0.95)

compare_models_waic <- function(models){
  waic_list <- lapply(models,function(model) waic(model$draws("log_lik",format="matrix")))
  waic_values <- sapply(waic_list,function(w) w$estimates["waic","Estimate"])
  waic_se <- sapply(waic_list,function(w) w$estimates["waic","SE"])
  p_waic <- sapply(waic_list,function(w) w$estimates["p_waic","Estimate"])
  p_waic_se <- sapply(waic_list,function(w) w$estimates["p_waic","SE"])
  dWAIC <- waic_values-min(waic_values)
  rel_lik <- exp(-0.5*dWAIC)
  weights <- rel_lik/sum(rel_lik)
  out <- data.frame(Model=names(models),WAIC=waic_values,WAIC_SE=waic_se,pWAIC=p_waic,pWAIC_SE=p_waic_se,dWAIC=dWAIC,Weight=weights)
  out <- out[order(out$WAIC),]
  rownames(out) <- NULL
  out
}

waic_simple <- compare_models_waic(list(no_time=fit0,linear=fit1,spline=fitS,step=fitStep))
print(waic_simple)
write.table(waic_simple,"WAIC_width_normal_models.txt",sep="\t",row.names=FALSE)

# ============================================================
# Plotting setup: absolute-width scale in millimetres
# ============================================================

xdates <- as.Date(c("1900-01-01","1910-01-01","1920-01-01","1930-01-01","1940-01-01","1950-01-01"))
xdati <- as.integer(xdates-origin_date)+1
xdats <- (xdati-mtime)/365
labs <- c("1900","1910","1920","1930","1940","1950")
newx <- seq(min(dl$time),max(dl$time),l=100)
Bnew <- predict(bs(d_fil$midTime, knots=knots[-c(1,K0)], degree=3, intercept=TRUE), newx)

ylim_width <- c(0,241)
y_ticks <- seq(0,240,20)
y_labs <- as.character(y_ticks)
dens_width <- density(dl$width_mm, from=ylim_width[1], to=ylim_width[2], n=1024)

shade_cmdstan <- function(ci,x,col=adjustcolor("black",alpha.f=0.25)){polygon(c(x,rev(x)),c(ci[1,],rev(ci[2,])),col=col,border=NA)}
right_density_panel <- function(){plot(dens_width$y,dens_width$x,type="n",xlab="Density",ylab="",ylim=ylim_width,yaxt="n",xaxs="i",xlim=c(0,max(dens_width$y)*1.10)); axis(2,at=y_ticks,labels=FALSE); polygon(c(0,dens_width$y,0),c(dens_width$x[1],dens_width$x,dens_width$x[length(dens_width$x)]),col="#00808080",border=NA); lines(dens_width$y,dens_width$x,lwd=1); abline(h=y_ticks,lty=2,lwd=0.7)}

set.seed(421)
mcol <- paste0(randomColor(count=length(unique(dl$mag)),luminosity="bright"),"05")
leg.col <- unique(data.frame(mag=as.integer(as.factor(d_fil$magazine)), label=d_fil$magazineLabel))
leg.col <- leg.col[order(leg.col$mag),]

png("tmp.png",width=1,height=1,units="in",res=100,type="cairo")
par(mar=c(0,0,0,0))
plot.new()
tmp <- legend("center",legend=leg.col$label,fill=substr(mcol,1,7),border="black",cex=1.2,title="Magazines",plot=FALSE)
dev.off()
unlink("tmp.png")

w <- tmp$rect$w*2.54+1
h <- tmp$rect$h*2.54+1
png("legend.png",width=w,height=h,units="cm",res=600,type="cairo")
par(mar=c(0,0,0,0))
plot.new()
legend("center",legend=leg.col$label,fill=substr(mcol,1,7),border="black",cex=1.2,title="Magazines")
dev.off()

# ============================================================
# Spline model prediction and figures
# ============================================================

post <- as_draws_df(fitS$draws(c("a","w[1]","w[2]","w[3]","w[4]","w[5]","w[6]","w[7]","sd_mag","sd_vol","sd_iss","sigma")))
W <- as.matrix(post[,grep("^w\\[",names(post))])
pred <- sapply(1:nrow(Bnew),function(i) post$a + W %*% Bnew[i,])
mu <- apply(pred,2,mean)
CI <- apply(pred,2,quantile,probs=c(0.05,0.95))

h <- 12
png("CzechSpline_width_normal1.png",height=h,width=h*1.6,units="cm",res=600,type="cairo")
layout(matrix(c(1,2),nrow=1),widths=c(5,1))
par(mar=c(3.5,5.0,2,1),mgp=c(3.0,0.65,0))
plot(dl$time,dl$width_mm,xlab="Year",ylab="",xaxt="n",yaxt="n",pch=19,col="#00808005",ylim=ylim_width)
axis(1,at=xdats,labels=labs,las=1)
axis(2,at=y_ticks,labels=y_labs,las=2)
mtext("Width in mm",side=2,line=2.7)
title("Czech magazines, normal spline for absolute width",adj=0)
abline(h=y_ticks,lty=2,lwd=0.7)
abline(v=xdats,lty=2,lwd=0.7)
shade_cmdstan(CI,newx)
lines(newx,mu,lwd=1)
par(mar=c(3.5,1,2,0.5),mgp=c(2,0.65,0))
right_density_panel()
dev.off()

png("CzechSpline_width_normal2.png",height=h,width=h*1.6,units="cm",res=600,type="cairo")
layout(matrix(c(1,2),nrow=1),widths=c(5,1))
par(mar=c(3.5,5.0,2,1),mgp=c(3.0,0.65,0))
plot(dl$time,dl$width_mm,xlab="Year",ylab="",xaxt="n",yaxt="n",pch=19,col=mcol[dl$mag],ylim=ylim_width)
axis(1,at=xdats,labels=labs,las=1)
axis(2,at=y_ticks,labels=y_labs,las=2)
mtext("Width in mm",side=2,line=2.7)
title("Czech magazines, normal spline for absolute width",adj=0)
abline(h=y_ticks,lty=2,lwd=0.7)
abline(v=xdats,lty=2,lwd=0.7)
shade_cmdstan(CI,newx)
lines(newx,mu,lwd=1)
par(mar=c(3.5,1,2,0.5),mgp=c(2,0.65,0))
right_density_panel()
dev.off()

# ============================================================
# Variance decomposition figures and numerical tables
# ============================================================

var_tab_from_post <- function(post,model_name,pars=c("sd_mag","sd_vol","sd_iss"),prob=0.9){
  lo <- (1-prob)/2; hi <- 1-lo
  sd_tab <- data.frame(model=model_name,scale="sd_mm",component=pars,mean=sapply(pars,function(p) mean(post[[p]])),median=sapply(pars,function(p) median(post[[p]])),q_low=sapply(pars,function(p) quantile(post[[p]],lo)),q_high=sapply(pars,function(p) quantile(post[[p]],hi)))
  V <- as.matrix(post[,pars])^2
  S <- V/rowSums(V)
  sh_tab <- data.frame(model=model_name,scale="variance_share",component=pars,mean=apply(S,2,mean),median=apply(S,2,median),q_low=apply(S,2,quantile,lo),q_high=apply(S,2,quantile,hi))
  rbind(sd_tab,sh_tab)
}

plot_var_decomp <- function(post,model_name,file,pars=c("sd_mag","sd_vol","sd_iss"),cols=c("#FF0040","#00F080","#2080E0"),sc=0.75){
  h <- 7.2
  x_upper <- max(sapply(pars,function(p) quantile(post[[p]],0.995,na.rm=TRUE)),na.rm=TRUE) * 1.10
  if(!is.finite(x_upper) || x_upper <= 0) x_upper <- 1
  png(file,height=h,width=h*1.7,units="cm",res=600,type="cairo")
  par(mar=c(3.5,5,2,1),mgp=c(2,0.65,0))
  plot(NULL,ylim=c(length(pars),1)+c(1,-1)*0.5,xlim=c(0,x_upper),xlab="Parameter estimate in mm",ylab="",yaxt="n",xaxs="i")
  abline(h=1:length(pars),col="#80808080")
  for(i in 1:length(pars)){
    toplot <- c(post[[pars[i]]])
    dens_i <- density(toplot,from=min(toplot),to=max(toplot))
    polygon(c(dens_i$x,rev(dens_i$x)),i+c(dens_i$y,rev(dens_i$y)*(-1))*sc,col=paste0(cols[i],"80"),border=cols[i])
    lines(PI(toplot,prob=0.9),rep(i,2),lwd=2)
    points(mean(toplot),i,pch=21,bg="white",cex=1.1,col=1,lwd=1.5)
    axis(2,at=i,labels=pars[i],col=cols[i],col.tick=cols[i],las=2)
  }
  abline(v=pretty(c(0,x_upper),n=8),col="#808080",lwd=0.7,lty=2)
  title(paste0("",model_name),adj=0)
  dev.off()
}

post_best <- as_draws_df(fitS$draws(c("sd_mag","sd_vol","sd_iss")))
post_null <- as_draws_df(fit0$draws(c("sd_mag","sd_vol","sd_iss")))
tab_best <- var_tab_from_post(post_best,"spline")
tab_null <- var_tab_from_post(post_null,"no_time")
write.table(tab_best,"Var_Czech_width_normal_spline_table.txt",sep="\t",row.names=FALSE,quote=FALSE)
write.table(tab_null,"Var_Czech_width_normal_no_time_table.txt",sep="\t",row.names=FALSE,quote=FALSE)
plot_var_decomp(post_best,"Czech, Splines","Var_Czech_width_normal_spline.png")
plot_var_decomp(post_null,"Czech, No time effect","Var_Czech_width_normal_no_time.png")

# ============================================================
# Step model prediction, cutpoints, baselines, and figures
# ============================================================

post_step <- as_draws_df(fitStep$draws(c("step_a[1]","step_a[2]","step_a[3]","step_a[4]","c_raw[1]","c_raw[2]","c_raw[3]","sd_mag","sd_vol","sd_iss","sigma")))
A <- as.matrix(post_step[,grep("^step_a\\[",names(post_step))])
Craw <- as.matrix(post_step[,grep("^c_raw\\[",names(post_step))])
Cut <- stan_data_step$t_min + (stan_data_step$t_max-stan_data_step$t_min)*plogis(Craw)

cut_med <- apply(Cut,2,median)
cut_lo <- apply(Cut,2,quantile,0.05)
cut_hi <- apply(Cut,2,quantile,0.95)

pred_step <- sapply(1:length(newx),function(i){
  s1 <- plogis(stan_data_step$kappa*(newx[i]-Cut[,1]))
  s2 <- plogis(stan_data_step$kappa*(newx[i]-Cut[,2]))
  s3 <- plogis(stan_data_step$kappa*(newx[i]-Cut[,3]))
  A[,1]*(1-s1) + A[,2]*(s1-s2) + A[,3]*(s2-s3) + A[,4]*s3
})

mu_step <- apply(pred_step,2,mean)
CI_step <- apply(pred_step,2,quantile,probs=c(0.05,0.95))

cut_years <- origin_date + round(Cut*365 + mtime - 1)
cut_tab <- data.frame(cutpoint=c("cut1","cut2","cut3"), mean=as.Date(apply(cut_years,2,function(x) mean(as.numeric(x))),origin="1970-01-01"), q05=as.Date(apply(cut_years,2,function(x) quantile(as.numeric(x),0.05)),origin="1970-01-01"), q95=as.Date(apply(cut_years,2,function(x) quantile(as.numeric(x),0.95)),origin="1970-01-01"))
print(cut_tab)
write.table(cut_tab,"Step_width_normal_cutpoints_dates.txt",sep="\t",row.names=FALSE,quote=FALSE)

step_a_tab <- data.frame(period=paste0("period_",1:4), mean_width_mm=apply(A,2,mean), q05_width_mm=apply(A,2,quantile,0.05), q95_width_mm=apply(A,2,quantile,0.95))
print(step_a_tab)
write.table(step_a_tab,"Step_width_normal_period_baselines.txt",sep="\t",row.names=FALSE,quote=FALSE)

h <- 12
png("CzechStep_width_normal1.png",height=h,width=h*1.6,units="cm",res=600,type="cairo")
layout(matrix(c(1,2),nrow=1),widths=c(5,1))
par(mar=c(3.5,5.0,2,1),mgp=c(3.0,0.65,0))
plot(dl$time,dl$width_mm,xlab="Year",ylab="",xaxt="n",yaxt="n",pch=19,col="#00808005",ylim=ylim_width)
axis(1,at=xdats,labels=labs,las=1)
axis(2,at=y_ticks,labels=y_labs,las=2)
mtext("Width in mm",side=2,line=2.7)
title("Czech magazines, normal step for absolute width",adj=0)
abline(h=y_ticks,lty=2,lwd=0.7)
abline(v=xdats,lty=2,lwd=0.7)
shade_cmdstan(CI_step,newx)
lines(newx,mu_step,lwd=1)
#for(j in 1:3) rect(cut_lo[j],par("usr")[3],cut_hi[j],par("usr")[4],col="#00000010",border=NA)
abline(v=cut_med,lty=3,lwd=1.2)
par(mar=c(3.5,1,2,0.5),mgp=c(2,0.65,0))
right_density_panel()
dev.off()

png("CzechStep_width_normal2.png",height=h,width=h*1.6,units="cm",res=600,type="cairo")
layout(matrix(c(1,2),nrow=1),widths=c(5,1))
par(mar=c(3.5,5.0,2,1),mgp=c(3.0,0.65,0))
plot(dl$time,dl$width_mm,xlab="Year",ylab="",xaxt="n",yaxt="n",pch=19,col=mcol[dl$mag],ylim=ylim_width)
axis(1,at=xdats,labels=labs,las=1)
axis(2,at=y_ticks,labels=y_labs,las=2)
mtext("Width in mm",side=2,line=2.7)
title("Czech magazines, normal step for absolute width",adj=0)
abline(h=y_ticks,lty=2,lwd=0.7)
abline(v=xdats,lty=2,lwd=0.7)
shade_cmdstan(CI_step,newx)
lines(newx,mu_step,lwd=1)
#for(j in 1:3) rect(cut_lo[j],par("usr")[3],cut_hi[j],par("usr")[4],col="#00000010",border=NA)
abline(v=cut_med,lty=3,lwd=1.2)
par(mar=c(3.5,1,2,0.5),mgp=c(2,0.65,0))
right_density_panel()
dev.off()

post_step_var <- as_draws_df(fitStep$draws(c("sd_mag","sd_vol","sd_iss")))
tab_step_var <- var_tab_from_post(post_step_var,"step")
write.table(tab_step_var,"Var_Czech_width_normal_step_table.txt",sep="\t",row.names=FALSE,quote=FALSE)
plot_var_decomp(post_step_var,"Czech, normal step for absolute width","Var_Czech_width_normal_step.png")

tab_var_all <- rbind(tab_best,tab_null,tab_step_var)
write.table(tab_var_all,"Var_Czech_width_normal_spline_no_time_step_table.txt",sep="\t",row.names=FALSE,quote=FALSE)


# ============================================================
# No-time model prediction
# ============================================================


#Paste variance decomposition images together
# library(magick)
# 
# add_panel_label <- function(file,label,size=48){
#   img <- image_read(file)
#   image_annotate(img,label,size=size,gravity="northeast",location="+35+25",color="black",weight=700)
# }
# 
# imgA <- add_panel_label("Var_Czech_width_normal_spline.png","A")
# imgB <- add_panel_label("Var_Czech_width_normal_no_time.png","B")
# 
# hA <- image_info(imgA)$height
# hB <- image_info(imgB)$height
# if(hA != hB) imgB <- image_scale(imgB,paste0("x",hA))
# 
# combined <- image_append(c(imgA,imgB),stack=FALSE)
# image_write(combined,"Var_Czech_width_normal_spline_vs_no_time.png")



