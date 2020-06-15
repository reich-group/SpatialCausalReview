rm(list=ls())

source("models.R")

# Simulation design

 type  <- 1
 batch <- 1
 nsims <- 12
 iters <- 5000
 burn  <- 1000
 nr    <- 30
 beta  <- 0.5
 phi   <- 0.5
 sige  <- 1
 sigU  <- 2
 sigV  <- 2

 g1 <- function(v,u){u+v}
 g2 <- function(v,u){v+ifelse(u>0,u,0)-0.63}
 g3 <- function(v,u){v + rep(1:30/30,30)*u}

 if(type==1){rhoU  <- 0.99; rhoV  <- 0.99; g<-g1}
 if(type==2){rhoU  <- 0.90; rhoV  <- 0.99; g<-g1}
 if(type==3){rhoU  <- 0.99; rhoV  <- 0.90; g<-g1}
 if(type==4){rhoU  <- 0.90; rhoV  <- 0.90; g<-g1}
 if(type==5){rhoU  <- 0.99; rhoV  <- 0.99; g<-g2}
 if(type==6){rhoU  <- 0.99; rhoV  <- 0.99; g<-g3}


print(paste("Scenario",type,"batch",batch))
filename <- paste0("mod",type,"batch",batch,".txt")


# Set up the spatial covariance matricies

 N     <- nr^2
 s     <- expand.grid(1:nr,1:nr)
 ADJ   <- as.matrix(dist(s))==1
 SU    <- sigU*sigU*solve(diag(rowSums(ADJ))-rhoU*ADJ)
 SV    <- sigV*sigV*solve(diag(rowSums(ADJ))-rhoV*ADJ)
 PU    <- t(chol(SU))
 PV    <- t(chol(SV))

# Start the simulation
 names   <- c("NS","NS+P","S","S+P","S+AIPW","Joint","Joint+cut")
 beta_mn <- matrix(0,nsims,length(names))
 colnames(beta_mn) <- names
 beta_lo <- beta_mn
 beta_hi <- beta_mn
 CPU     <- beta_mn

 for(sim in 1:nsims){     
  set.seed(919*sim*batch)

  print(paste("Starting dataset",sim,"of",nsims))

  # Generate and plot a dataset
  expit <- function(x){1/(1+exp(-x))}
  U     <- as.vector(PU%*%rnorm(N))
  V     <- as.vector(PV%*%rnorm(N))
  eta   <- g(V,phi*U)
  A     <- rbinom(N,1,expit(eta))
  Y     <- beta*A + phi*U + sige*rnorm(N)

  library(fields)
  library(viridis)
  par(mfrow=c(2,2))
  image.plot(1:nr,1:nr,matrix(U,nr,nr),col=viridis(50),main="U",xlab="",ylab="",axes=FALSE)
  image.plot(1:nr,1:nr,matrix(V,nr,nr),col=viridis(50),main="V",xlab="",ylab="",axes=FALSE)
  image.plot(1:nr,1:nr,matrix(A,nr,nr),col=viridis(50),main="A",xlab="",ylab="",axes=FALSE) 
  image.plot(1:nr,1:nr,matrix(Y,nr,nr),col=viridis(50),main="Y",xlab="",ylab="",axes=FALSE)

 # Propensity score regression
  library(splines)
  p_hat  <- spatial_logit(A=A,ADJ=ADJ,iters=iters,burn=burn)$prop_score
  fp_hat <- bs(p_hat,df=5)

  # Response models 
  fit1   <- BLR(Y,A,iters=iters,burn=burn,filename=filename)
  fit2   <- BLR_prop(Y,A,X=fp_hat,iters=iters,burn=burn,filename=filename)
  fit3   <- spatial_gaussian(Y=Y,A=A,ADJ=ADJ,iters=iters,burn=burn,filename=filename)
  fit4   <- spatial_prop(Y=Y,A=A,ADJ=ADJ,X=fp_hat,iters=iters,burn=burn,filename=filename)
  fit5   <- spatial_gaussian(Y=Y,A=A,ADJ=ADJ,prop_score=p_hat,iters=iters,burn=burn,filename=filename)
  fit6   <- joint(Y=Y,A=A,ADJ=ADJ,iters=iters,burn=burn,thin=2,filename=filename,plot=F)
  fit7   <- joint_nofeedback(Y=Y,A=A,ADJ=ADJ,iters=iters,thin=2,burn=burn,filename=filename)

  beta_mn[sim,1] <- mean(fit1$beta)
  beta_mn[sim,2] <- mean(fit2$beta)
  beta_mn[sim,3] <- mean(fit3$beta)
  beta_mn[sim,4] <- mean(fit4$beta)
  beta_mn[sim,5] <- mean(fit5$beta)
  beta_mn[sim,6] <- mean(fit6$beta)
  beta_mn[sim,7] <- mean(fit7$beta)

  beta_lo[sim,1] <- quantile(fit1$beta,0.025)
  beta_lo[sim,2] <- quantile(fit2$beta,0.025)
  beta_lo[sim,3] <- quantile(fit3$beta,0.025)
  beta_lo[sim,4] <- quantile(fit4$beta,0.025)
  beta_lo[sim,5] <- quantile(fit5$beta,0.025)
  beta_lo[sim,6] <- quantile(fit6$beta,0.025)
  beta_lo[sim,7] <- quantile(fit7$beta,0.025)

  beta_hi[sim,1] <- quantile(fit1$beta,0.975)
  beta_hi[sim,2] <- quantile(fit2$beta,0.975)
  beta_hi[sim,3] <- quantile(fit3$beta,0.975)
  beta_hi[sim,4] <- quantile(fit4$beta,0.975)
  beta_hi[sim,5] <- quantile(fit5$beta,0.975)
  beta_hi[sim,6] <- quantile(fit6$beta,0.975)
  beta_hi[sim,7] <- quantile(fit7$beta,0.975)

  CPU[sim,1]     <- fit1$seconds
  CPU[sim,2]     <- fit2$seconds
  CPU[sim,3]     <- fit3$seconds
  CPU[sim,4]     <- fit4$seconds
  CPU[sim,5]     <- fit5$seconds
  CPU[sim,6]     <- fit6$seconds
  CPU[sim,7]     <- fit7$seconds

 }

BIAS <- colMeans(beta_mn-beta)
MSE  <- colMeans((beta_mn-beta)^2)
COV  <- colMeans((beta_lo<=beta) & (beta<=beta_hi))

name <- paste0("Sim",type,"batch",batch,".RData")
save.image(name)

 